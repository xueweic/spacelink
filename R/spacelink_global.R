#' Global Spatial Variability Analysis
#'
#' Performs global spatial gene variability analysis across multiple length scales
#' using kernel-based methods.
#'
#' @param normalized_counts Normalized count matrix. Rows and columns indicate genes and spots, respectively.
#' @param spatial_coords Spatial coordinates matrix. Rows indicate spots.
#' @param covariates Spatial features (without intercept). Rows indicate spots.
#' @param lengthscales If NULL, length-scales are calculated using the two-step approach.
#' @param n_lengthscales Number of length-scales (i.e., number of kernels). Default is 5.
#' @param M Controls the minimum length-scale. Minimum length-scale is set to be the minimum distance multiplied by M. Default is 1.
#' @param n_workers Number of workers for parallelization. Default is 1.
#'
#' @return A data frame containing:
#' \describe{
#'   \item{tau.sq}{Nugget variance component}
#'   \item{sigma.sq1, sigma.sq2, ...}{Spatial variance components for each kernel}
#'   \item{raw_ESV}{Effective Spatial Variability score}
#'   \item{pval}{Combined p-value from score tests}
#'   \item{padj}{Benjamini-Hochberg adjusted p-values}
#'   \item{ESV}{ESV adjusted for multiple testing (0 if padj > 0.05)}
#' }
#'
#' @examples
#' # Set up simulated data
#' library(spacelink)
#' set.seed(123)
#' n_spots <- 100
#' n_genes <- 10
#' 
#' # Generate spatial coordinates (2D grid)
#' coords <- expand.grid(
#'   x = seq(0, 10, length.out = sqrt(n_spots)),
#'   y = seq(0, 10, length.out = sqrt(n_spots))
#' )
#' 
#' # Simulate (normalized) expression data with spatial autocorrelation
#' expr_data <- matrix(nrow=n_genes, ncol=n_spots)
#' rownames(expr_data) <- paste0("Gene_", 1:n_genes)
#' D <- as.matrix(dist(coords))
#' K <- exp(-D/2); chol_K <- chol(K)
#' for(i in 1:(n_genes/2)){
#'   expr_data[i,] <- matrix(rnorm(n_spots, 0, 1), nrow=1, ncol=n_spots)
#'   expr_data[i,] <- expr_data[i,] + i*matrix(rnorm(n_spots, 0, 1), 1, n_spots)%*%chol_K
#' }
#' # Simulate (normalized) expression data without spatial autocorrelation
#' for(i in (n_genes/2+1):n_genes){
#'   expr_data[i,] <- matrix(rnorm(n_spots, 0, 1), nrow=1, ncol=n_spots)
#' }
#' 
#' global_results <- spacelink_global(normalized_counts = expr_data, spatial_coords = coords)
#' print(global_results[, c("raw_ESV", "pval", "ESV", "padj")])
#'
#' @importFrom ACAT ACAT
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom stats p.adjust
#' @export
spacelink_global <- function(normalized_counts, spatial_coords, covariates = NULL, lengthscales = NULL, n_lengthscales = 5, M = 1, n_workers = 1){
  Y <- normalized_counts
  X <- covariates
  if(ncol(Y)!=nrow(spatial_coords)){
      stop("The column dimension of normalized_counts and the row dimension of spatial_coords should be the same.")
  }
  if(!is.null(X)){
    if(nrow(spatial_coords)!=nrow(X)){
      stop("The row dimensions of spatial_coords and covariates should be the same.")
    }
  }
  # Create matrix of phi (= 1/lengthscale)
  if(!is.null(lengthscales)){
    if(nrow(lengthscales)!=nrow(Y)){
        stop("The row dimensions of normalized counts and lengthscales should be the same.")
    }
    phi_mat <- 1/lengthscales
    n_lengthscales <- ncol(phi_mat)
  }else{
    out <- select_lengthscales(Y, spatial_coords, n_lengthscales, M = M, is.inverse = TRUE, n_workers)
    phi_mat <- out$phi_mat
    time0 <- out$time
  }
  colnames(phi_mat) <- paste0("phi", 1:n_lengthscales)
  phi_mat <- as.data.frame(phi_mat)

  if(is.null(X)){
    Y <- Y - apply(Y,1,mean)
  }else{
    X <- cbind(rep(1,nrow(X)), X)
    Y <- Y - t(X %*% solve(crossprod(X), crossprod(X,t(Y))))
  }

  out <- bplapply(1:nrow(Y), function(gene_idx) {
    runtime <- system.time({
      y <- as.numeric(Y[gene_idx,])
      phi_seq <- phi_mat[gene_idx,]
      nnls_res <- run_nnls(y, spatial_coords, phi_seq = phi_seq)
      nnls_res <- cbind(nnls_res, phi_seq)

      sigma.sq_vec <- nnls_res[,grep("sigma", colnames(nnls_res))]
      tau.sq <- nnls_res[,grep("tau",colnames(nnls_res))]
      nnls_res$raw_ESV <- calculate_ESV(sigma.sq_vec, tau.sq, spatial_coords, phi_seq = phi_seq)

      test_res <- score_test(y, spatial_coords, phi_seq = phi_seq)
    })

    nnls_res$time <- runtime[["elapsed"]]

    return(list(nnls_res=nnls_res, pval_vec=test_res$pval_vec))
  }, BPPARAM = MulticoreParam(workers = n_workers))

  results <- do.call("rbind", lapply(out, function(x)x$nnls_res))
  if(is.null(lengthscales)){results$time <- results$time + time0}

  pval_mat <- do.call("rbind", lapply(out, function(x)x$pval_vec))
  colnames(pval_mat) <- paste0("pval",1:ncol(pval_mat))
  results <- cbind(results, pval_mat)
  suppressWarnings({results$pval <- apply(pval_mat,1,ACAT)})
  results$padj <- p.adjust(results$pval, method = "BH")

  results$ESV <- results$raw_ESV
  results$ESV[results$padj > 0.05] <- 0

  if(!is.null(rownames(Y))){
    rownames(results) <- make.unique(rownames(Y))
  }

  return(results)
}
