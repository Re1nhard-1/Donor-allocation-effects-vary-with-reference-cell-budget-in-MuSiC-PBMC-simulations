root <- "results/music_components"
files <- list.files(root,pattern="^component_b.*\\.rds$",recursive=TRUE,full.names=TRUE)
stopifnot(length(files)==252L)
rows <- lapply(files,function(path) {
  x <- readRDS(path)
  do.call(rbind,lapply(0:7,function(mask) {
    p <- x$fits[[mask+1]]$pred
    tm <- x$targets
    z <- data.frame(triple=x$triple,block=x$block,budget=x$budget,dominant=x$dominant,
      target_name=tm$target_name,held_out=tm$held_out,mask=mask,mae=x$fits[[mask+1]]$mae,
      original_balanced=x$original_balanced_mae,original_unequal=x$original_unequal_mae,
      phi_Theta=x$phi[,1],phi_S=x$phi[,2],phi_Sigma=x$phi[,3],
      convergence=x$fits[[mask+1]]$convergence,n_common=length(x$common_genes))
    for(k in 0:5) {
      z[[paste0("pred_",k)]] <- p[,k+1]
      z[[paste0("true_",k)]] <- as.numeric(tm[[paste0("true_",k)]])
    }
    z
  }))
})
z <- do.call(rbind,rows)
stopifnot(nrow(z)==64512L)
for(k in names(z)) if(is.double(z[[k]]))z[[k]] <- sprintf("%.17g",z[[k]])
write.csv(z,file.path(root,"predictions.csv"),row.names=FALSE)
cat(nrow(z),"logical prediction rows exported\n")
