
main = function(data, ncores){
  ngene = nrow(data)

  # fit NB
  fitnb = mclapply(1:ngene, function(i){
    if(i %% 1000 == 0) message(paste("Gene", i))
    fit = try(fitdist(data[i,], "nbinom"), silent = TRUE)
    if(class(fit) == "try-error"){
      mu = NA
      disp = NA
    }else{
      mu = fit$estimate["mu"]
      size = fit$estimate["size"]
      # var.fit = mu + mu^2/size
      disp = 1/size
    }
    return(c(mu, disp))
  }, mc.cores = ncores)
  fitnb = Reduce(rbind, fitnb)
  colnames(fitnb) = c("mu", "dispersion")
  rownames(fitnb) = rownames(data)

  # phitest
  var = apply(data, 1, var)
  mean = rowMeans(data)

  res = as.data.frame(fitnb)
  res$mean = mean
  res$var = var
  ### zero percentage
  d = rowMeans(data == 0)
  res$d = d

  res = res[res$var > res$mean, ]

  qtx = quantile(res$mean, 0.99)
  qty = quantile(res$var, 0.99)
  idx = which(res$mean <= qtx & res$var <= qty)
  y = (res$var[idx] - res$mean[idx])
  x = (res$mean[idx])^2
  lm1 = lm(y ~ x + 0)
  phi = coef(lm1)
  res$phi = phi

  ### zero percentage
  res$dhat.c = (1+phi * res$mean)^(-1/phi)
  res$dhat = (1+res$dispersion * res$mean)^(-1/res$dispersion)

  obj = t.test(res$dhat.c, res$dhat, alternative = "less")
  return(list(pval = obj$p.value, par = res))
}


#' Phitest for analyzing the heterogeneity of single-cell populations
#' @title Appling the Phitest method to a count matrix
#' @param object A matrix of single-cell UMI counts (rows for genes and columns for cells).
#' @param label A character or numeric vector of cluster labels. Length should be the same as
#' cell number and order should match the order in \code{object}.
#' @param ncores Number of cores used for parallel computation.
#' @param min.cell An integer specifying a threshold to filter genes. Genes expressed in
#' fewer than \code{min.cell} cells are filtered out.
#' @return A list of two elements: \code{pval} contains the \emph{P} values, and
#' \code{par} contains the estimated parameters.
#' @export
#' @import parallel
#' @importFrom fitdistrplus fitdist
#' @importFrom stats coef lm quantile t.test
#' @author Wei Vivian Li, \email{vivian.li@rutgers.edu}
phitest = function(object, label, ncores = 1, min.cell = 10){
  data = object[rowSums(object > 0) > min.cell, ]
  ngene = nrow(data)

  if(!(class(label) %in% c("character", "numeric"))) stop("label does not belong to an accepted class!")
  if(length(label) != ncol(data)) stop("Length of label does not match cell number!")

  phires = lapply(unique(label), function(x){
    message(paste("Cluster", x))
    cnt = data[, label == x]
    cnt = cnt[rowSums(cnt > 0) > min.cell, ]

    tp = try(main(cnt, ncores), silent = TRUE)
    if(class(tp) == "try-error") return(NULL)
    return(tp)
  })
  names(phires) = unique(label)
  pval = sapply(phires, function(x) x$pval)
  names(pval) = unique(label)
  par = lapply(phires, function(x) x$par)
  return(list(pval =  pval, par = par))
}

