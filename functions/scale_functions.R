# longitudinal standardisation
scale_rm <- function(mat, origin = 1, centre = F) {
  ncols = ncol(mat)
  nrows = nrow(mat)

  origin_rows = which(mat$time == 1)

  if(centre) {
    cent = sapply(3:ncols, function(j) mean(mat[origin_rows,j]))
    mat[,3:ncols] <- sweep(mat[3:ncols], 2L, cent, check.margin=FALSE)
  }
  dist = sapply(3:ncols, function(j) sqrt(sum(mat[,j]**2)/max(1, length(mat[,j]-1))))
  mat_st = sapply(3:ncols, function(j) mat[,j]/dist[j-2])                       # faster than sweep
  # mat_st = sweep(mat, 2, dist, "/")

  mat_st = as.matrix(cbind(mat[,1:2], mat_st), nrow = nrows, ncol = ncols)
  colnames(mat_st) <- colnames(mat)

  return(data.frame(mat_st))
}
