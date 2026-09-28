source("environment/music/runtime.R")
options(digits=17)
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==4L, as.character(packageVersion("MuSiC"))=="1.0.0")
manifest <- jsonlite::fromJSON(args[1]); out <- args[2]; triple <- args[3]; mode <- args[4]
stopifnot(mode %in% c("control", "full"))
dir.create(out, recursive=TRUE, showWarnings=FALSE)
started <- proc.time()["elapsed"]
input <- jsonlite::fromJSON(manifest$input_manifest)
types <- input$cell_types; data_dir <- input$data_directory
cells <- read.delim(file.path(data_dir,"cells.tsv"), check.names=FALSE, colClasses="character")
refs <- read.delim(manifest$references_file, check.names=FALSE, colClasses="character")
cases <- read.delim(manifest$cases_file, check.names=FALSE, colClasses="character")
targets <- read.delim(manifest$targets_file, check.names=FALSE, colClasses="character")
all_targets <- read.delim(file.path(data_dir,"targets.tsv"), check.names=FALSE, colClasses="character")
inventory <- read.csv(manifest$inventory_file, check.names=FALSE, colClasses="character")
genes <- readLines(file.path(data_dir,"genes.tsv"))
stopifnot(length(genes)==30172L, !anyDuplicated(genes), length(types)==6L,
          nrow(targets)==112L, all(table(targets$held_out)==8L))
read_matrix <- function(path) {
  h <- gzfile(path,"rb"); on.exit(close(h))
  methods::as(Matrix::readMM(h),"CsparseMatrix")
}
bulk_all <- read_matrix(file.path(data_dir,"targets.mtx.gz"))
dimnames(bulk_all) <- list(genes,all_targets$target_name)
refs <- refs[refs$triple_key==triple,,drop=FALSE]
cases <- cases[cases$triple_key==triple,,drop=FALSE]
stopifnot(nrow(refs)==36L, nrow(cases)==144L, !anyDuplicated(refs$reference_id),
          setequal(refs$level,c("balanced","ratio10")))
needed <- strsplit(refs$reference_donors[1],"|",fixed=TRUE)[[1]]
held <- unique(cases$held_out)
stopifnot(length(needed)==3L, length(held)==4L, !any(held %in% needed),
          all(refs$reference_donors==refs$reference_donors[1]))
targets <- targets[targets$held_out %in% held,,drop=FALSE]
stopifnot(nrow(targets)==32L, !anyDuplicated(targets$target_name),
          all(table(targets$held_out)==8L), all(targets$target_name %in% colnames(bulk_all)))
# Preserve the target order in the frozen selection; do not construct new targets.
bulk_all <- as.matrix(bulk_all[,targets$target_name,drop=FALSE])
counts_by_donor <- lapply(needed,function(d) {
  path <- input$donor_matrices$path[input$donor_matrices$donor==d]
  stopifnot(length(path)==1L)
  x <- read_matrix(path); cm <- cells[cells$donor==d,,drop=FALSE]
  stopifnot(nrow(x)==length(genes), ncol(x)==nrow(cm))
  dimnames(x) <- list(genes,cm$cell_id); x
}); names(counts_by_donor) <- needed
refs <- refs[order(as.integer(refs$block),as.integer(refs$budget),
                   refs$level!="balanced",refs$dominant_donor),,drop=FALSE]
make_sce <- function(ref) {
  ids <- as.integer(strsplit(ref$reference_columns_R,"|",fixed=TRUE)[[1]])
  stopifnot(all(ids>=1L),all(ids<=nrow(cells)))
  cm <- cells[ids,,drop=FALSE]
  quota <- table(cm$donor,cm$cell_type)
  expected <- ifelse(rownames(quota)==ref$dominant_donor,
                     as.integer(ref$major_count),as.integer(ref$minor_count))
  stopifnot(length(ids)==as.integer(ref$budget)*6L,length(unique(ids))==length(ids),
            setequal(unique(cm$donor),needed),all(quota==expected),
            all(colSums(quota)==as.integer(ref$budget)),!any(held %in% cm$donor))
  pieces <- lapply(needed,function(d) counts_by_donor[[d]][,cm$cell_id[cm$donor==d],drop=FALSE])
  x <- do.call(cbind,pieces); x <- x[,cm$cell_id,drop=FALSE]; rownames(cm) <- cm$cell_id
  stopifnot(max(abs(Matrix::colSums(x)-as.numeric(cm$total_counts)))==0)
  SingleCellExperiment(assays=list(counts=x),colData=S4Vectors::DataFrame(cm))
}
truth <- as.matrix(data.frame(lapply(targets[paste0("true_",0:5)],as.numeric)))
colnames(truth) <- types
bulk_genes <- rownames(bulk_all)[rowMeans(bulk_all)>0]
old <- read.csv(file.path("results/music_mechanism/20260917T210142443179Z",triple,"predictions.csv"),check.names=FALSE)
old <- old[old$method=="music_weighted",]
endpoint_diff <- 0
basis_diff <- 0
control_diff <- 0
fit_count <- 0L
deadline <- Sys.time()+3600
fit_components <- function(theta,S,Sigma,support) {
  theta <- theta[support,types,drop=FALSE]
  Sigma <- Sigma[support,types,drop=FALSE]
  S <- S[types]
  D <- sweep(theta,2,S,"*")
  yy <- bulk_all[support,,drop=FALSE]
  yy <- sweep(yy,2,colSums(yy),"/")
  ans <- matrix(NA_real_,nrow(targets),length(types),dimnames=list(targets$target_name,types))
  conv <- character(nrow(targets))
  for(t in seq_len(nrow(targets))) {
    if(Sys.time()>deadline) stop("Fitting time limit reached")
    keep <- yy[,t]>0
    y <- yy[keep,t];names(y) <- support[keep]
    z <- music.iter(Y=y,D=D[keep,,drop=FALSE],S=S,Sigma=Sigma[keep,,drop=FALSE],iter.max=1000,nu=0.0001,eps=0.01,centered=FALSE,normalize=FALSE)
    ans[t,] <- z$p.weight
    conv[t] <- z$converge
    stopifnot(all(is.finite(z$var.p)))
  }
  fit_count <<- fit_count+nrow(targets)
  stopifnot(all(is.finite(ans)),min(ans)>=-1e-10,max(abs(rowSums(ans)-1))<1e-10,
            all(grepl("^Converge at [0-9]+$",conv)))
  list(pred=ans,mae=100*rowMeans(abs(ans-truth)),convergence=conv)
}
get_basis <- function(ref) {
  sce <- make_sce(ref)
  z <- music_basis(sce,non.zero=TRUE,markers=bulk_genes,clusters="cell_type",samples="donor",select.ct=types,cell_size=NULL,ct.cov=FALSE,verbose=FALSE)
  S <- colMeans(z$S,na.rm=TRUE)
  delta <- max(abs(sweep(z$M.theta,2,S,"*")-z$Disgn.mtx))
  basis_diff <<- max(basis_diff,delta)
  stopifnot(delta<=1e-10,all(is.finite(z$Sigma)),all(z$Sigma>=0),all(S>0))
  list(theta=z$M.theta,S=S,Sigma=z$Sigma,support=intersect(rownames(z$Disgn.mtx),bulk_genes))
}
compare_saved <- function(ref,fit) {
  o <- old[old$reference_id==ref$reference_id,]
  o <- o[match(targets$target_name,o$target_name),]
  stopifnot(nrow(o)==nrow(targets),!anyNA(o$target_name))
  p <- as.matrix(o[paste0("pred_",0:5)])
  delta <- max(abs(p-fit$pred));stopifnot(delta<=1e-10)
  control_diff <<- max(control_diff,delta)
}
shapley <- function(scores) {
  ans <- matrix(0,nrow(scores),3,dimnames=list(NULL,c("Theta","S","Sigma")))
  for(j in 0:2) for(mask in 0:7) if(bitwAnd(mask,bitwShiftL(1,j))==0) {
    n <- sum(as.integer(intToBits(mask))[1:3])
    w <- factorial(n)*factorial(2-n)/factorial(3)
    ans[,j+1] <- ans[,j+1]+w*(scores[,bitwOr(mask,bitwShiftL(1,j))+1]-scores[,mask+1])
  }
  stopifnot(max(abs(rowSums(ans)-(scores[,8]-scores[,1])))<=1e-10)
  ans
}
blocks <- if(mode=="control") 0L else 0:2
for(block in blocks) {
  rr <- refs[as.integer(refs$block)==block & as.integer(refs$budget)%in%c(60,300),]
  bases <- lapply(seq_len(nrow(rr)),function(i)get_basis(rr[i,]))
  names(bases) <- rr$reference_id
  common <- Reduce(intersect,lapply(bases,function(x)x$support))
  stopifnot(length(common)>1000)
  endpoint_cache <- list()
  doms <- sort(needed)
  if(mode=="control") doms <- doms[1]
  for(budget in c(60L,300L)) {
    balrow <- rr[as.integer(rr$budget)==budget & rr$level=="balanced",]
    bal <- bases[[balrow$reference_id]]
    for(dom in doms) {
      unrow <- rr[as.integer(rr$budget)==budget & rr$dominant_donor==dom & rr$level=="ratio10",]
      stopifnot(nrow(unrow)==1L)
      un <- bases[[unrow$reference_id]]
      path <- file.path(out,paste0("component_b",block,"_n",budget,"_",dom,".rds"))
      if(file.exists(path)) next
      if(mode=="control") {
        for(refrow in list(balrow,unrow)) {
          q <- bases[[refrow$reference_id]]
          f <- fit_components(q$theta,q$S,q$Sigma,q$support)
          compare_saved(refrow,f)
          standard <- music_prop(bulk.mtx=bulk_all,sc.sce=make_sce(refrow),markers=common,
            clusters="cell_type",samples="donor",select.ct=types,verbose=FALSE)
          same <- fit_components(q$theta,q$S,q$Sigma,common)
          dd <- max(abs(same$pred-standard$Est.prop.weighted[targets$target_name,types]))
          endpoint_diff <- max(endpoint_diff,dd);stopifnot(dd<=1e-10)
        }
      }
      allfits <- vector("list",8)
      for(mask in 0:7) {
        if(mask==0 && !is.null(endpoint_cache[[as.character(budget)]])) {
          allfits[[mask+1]] <- endpoint_cache[[as.character(budget)]]
        } else {
          th <- if(bitwAnd(mask,1L)>0)un$theta else bal$theta
          ss <- if(bitwAnd(mask,2L)>0)un$S else bal$S
          sg <- if(bitwAnd(mask,4L)>0)un$Sigma else bal$Sigma
          allfits[[mask+1]] <- fit_components(th,ss,sg,common)
          if(mask==0) endpoint_cache[[as.character(budget)]] <- allfits[[1]]
        }
      }
      scores <- do.call(cbind,lapply(allfits,function(x)x$mae))
      phi <- shapley(scores)
      original_mae <- function(refrow) {
        o <- old[old$reference_id==refrow$reference_id,]
        100*o$mae[match(targets$target_name,o$target_name)]
      }
      record <- list(triple=triple,block=block,budget=budget,dominant=dom,targets=targets,
        types=types,common_genes=common,original_support=c(balanced=length(bal$support),unequal=length(un$support)),
        fits=allfits,phi=phi,original_balanced_mae=original_mae(balrow),original_unequal_mae=original_mae(unrow))
      saveRDS(record,path)
      cat(basename(path),"saved; elapsed",round(proc.time()["elapsed"]-started,1),"s\n");flush.console()
    }
  }
}
info <- list(mode=mode,triple=triple,completed=TRUE,actual_adapter_fits=fit_count,
             standard_endpoint_control_fits=if(mode=="control")128L else 0L,
             basis_max_difference=basis_diff,saved_prediction_max_difference=control_diff,
             common_endpoint_max_difference=endpoint_diff,elapsed_seconds=unname(proc.time()["elapsed"]-started))
jsonlite::write_json(info,file.path(out,paste0("component_",mode,".json")),auto_unbox=TRUE,pretty=TRUE,digits=17)
capture.output(sessionInfo(),file=file.path(out,"component_session.txt"))
