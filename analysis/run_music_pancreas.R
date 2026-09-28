source("environment/music/runtime.R")
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2L, as.character(packageVersion("MuSiC")) == "1.0.0")
held <- args[1]
mode <- args[2]
stopifnot(mode %in% c("pilot", "full"))
input <- "data/processed/music_pancreas"
out <- file.path("results/music_pancreas", held)
dir.create(out, recursive=TRUE, showWarnings=FALSE)
config <- jsonlite::fromJSON(file.path(input, "design.json"))
stopifnot(held %in% config$donors)
types <- config$cell_types
cells <- read.delim(file.path(input, "cells.tsv"), check.names=FALSE, colClasses="character")
targets <- read.delim(file.path(input, "targets.tsv"), check.names=FALSE)
refs <- read.delim(file.path(input, "references.tsv"), check.names=FALSE, colClasses="character")
genes <- readLines(file.path(input, "genes.tsv"))
read_matrix <- function(path) {
  h <- gzfile(path, "rb")
  on.exit(close(h))
  methods::as(Matrix::readMM(h), "CsparseMatrix")
}
counts <- read_matrix(file.path(input, "counts.mtx.gz"))
dimnames(counts) <- list(genes, cells$cell_id)
bulk <- read_matrix(file.path(input, "targets.mtx.gz"))
dimnames(bulk) <- list(genes, targets$target_name)
targets <- targets[targets$held_out == held, , drop=FALSE]
bulk <- as.matrix(bulk[, targets$target_name, drop=FALSE])
truth <- as.matrix(targets[paste0("true_", types)])
colnames(truth) <- types
stopifnot(nrow(targets) == 60L, max(abs(rowSums(truth) - 1)) < 1e-12)
refs <- refs[refs$held_out == held, , drop=FALSE]
if (mode == "pilot") refs <- refs[refs$block == "0", , drop=FALSE]
.music_trace <- new.env(parent=emptyenv())
.music_trace$items <- list()
invisible(trace("music.iter", where=asNamespace("MuSiC"), print=FALSE, exit=quote({
  z <- returnValue()
  st <- get(".music_trace", envir=.GlobalEnv)
  st$items[[length(st$items) + 1L]] <- list(convergence=z$converge, n_features=nrow(D), variance_finite=all(is.finite(z$var.p)))
})))
for (i in seq_len(nrow(refs))) {
  ref <- refs[i, ]
  path <- file.path(out, paste0(ref$reference_id, ".rds"))
  if (file.exists(path)) next
  ids <- as.integer(strsplit(ref$columns_R, "|", fixed=TRUE)[[1]])
  cm <- cells[ids, , drop=FALSE]
  quota <- table(cm$donor, cm$cell_type)
  expected <- if (ref$dominant == "balanced") rep(as.integer(ref$budget) / 3, 3) else ifelse(rownames(quota) == ref$dominant, as.integer(ref$budget) * 10/12, as.integer(ref$budget) / 12)
  stopifnot(length(ids) == as.integer(ref$budget) * length(types), !anyDuplicated(ids), !(held %in% cm$donor),
            all(quota == expected), all(colSums(quota) == as.integer(ref$budget)))
  x <- counts[, ids, drop=FALSE]
  stopifnot(max(abs(Matrix::colSums(x) - as.numeric(cm$umi_total))) == 0)
  rownames(cm) <- cm$cell_id
  sce <- SingleCellExperiment(assays=list(counts=x), colData=S4Vectors::DataFrame(cm))
  .music_trace$items <- list()
  fit <- suppressMessages(MuSiC::music_prop(bulk.mtx=bulk, sc.sce=sce, markers=NULL, clusters="cell_type", samples="donor", select.ct=types,
               cell_size=NULL, ct.cov=FALSE, verbose=FALSE, iter.max=1000, nu=0.0001, eps=0.01, centered=FALSE, normalize=FALSE))
  estimates <- list(weighted=fit$Est.prop.weighted[targets$target_name, types, drop=FALSE], nnls=fit$Est.prop.allgene[targets$target_name, types, drop=FALSE])
  for (p in estimates) stopifnot(all(is.finite(p)), min(p) >= -1e-10, max(abs(rowSums(p) - 1)) < 1e-10)
  diagnostics <- do.call(rbind, lapply(.music_trace$items, as.data.frame))
  stopifnot(nrow(diagnostics) == nrow(targets))
  rows <- do.call(rbind, lapply(names(estimates), function(method) {
    p <- estimates[[method]]
    tab <- data.frame(reference_id=ref$reference_id, held_out=held, block=as.integer(ref$block), budget=as.integer(ref$budget),
                      dominant=ref$dominant, target_name=targets$target_name, method=method, mae_pp=100 * rowMeans(abs(p-truth)),
                      rmse_pp=100 * sqrt(rowMeans((p-truth)^2)), convergence=diagnostics$convergence,
                      n_features=diagnostics$n_features, variance_finite=diagnostics$variance_finite)
    for (k in seq_along(types)) {
      tab[[paste0("true_", types[k])]] <- truth[, k]
      tab[[paste0("pred_", types[k])]] <- p[, k]
    }
    tab
  }))
  saveRDS(rows, path)
  cat(ref$reference_id, "saved\n")
  flush.console()
}
paths <- file.path(out, paste0(refs$reference_id, ".rds"))
stopifnot(all(file.exists(paths)))
rows <- do.call(rbind, lapply(paths, readRDS))
write.csv(rows, file.path(out, paste0(mode, "_predictions.csv")), row.names=FALSE)
capture.output(sessionInfo(), file=file.path(out, "session_info.txt"))
cat("COMPLETE", held, mode, nrow(rows), "rows\n")
