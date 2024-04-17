# longitudinal standardisation
scale_rm <- function(mat, origin = 1, centre = F) {
  ncols = ncol(mat)
  nrows = nrow(mat)
  origin_rows = which(mat$time == origin)
  # scale first columnwise
  cent = sapply(3:ncols, function(j) mean(mat[origin_rows,j]))
  dist   = sapply(3:ncols, function(j) sqrt(sum((mat[,j] - cent[j-2])**2)/max(1, length(mat[,j]-1))))
  mat[,3:ncols] = sapply(3:ncols, function(j) mat[,j]/dist[j-2])                       # faster than sweep
  # centre whole dataset
  if(centre) {
    gen_cent = mean(as.matrix(mat[origin_rows,3:ncols]))
    mat[,3:ncols]  = sapply(3:ncols, function(j) mat[,j] - gen_cent)
  }
  # cent = sapply(3:ncols, function(j) mean(mat[origin_rows,j]))
  # if(centre) {
  #   mat[,3:ncols] <- sapply(3:ncols, function(j) mat[,j] - cent[j-2])
  #   # mat[,3:ncols] <- sweep(mat[3:ncols], 2L, cent, check.margin=FALSE)
  #   dist = sapply(3:ncols, function(j) sqrt(sum(mat[,j]**2)/max(1, length(mat[,j]-1))))
  #
  # } else {
  #   dist = sapply(3:ncols, function(j) sqrt(sum((mat[,j] - cent[j-2])**2)/max(1, length(mat[,j]-1))))
  #
  # }
  # mat_st = sweep(mat[,3:ncols], 2, dist, "/", check.margin = FALSE)
  # mat_st = as.matrix(cbind(mat[,1:2], mat_st), nrow = nrows, ncol = ncols)
  # colnames(mat_st) <- colnames(mat)
  return(data.frame(mat))
}

#
# # longitudinal standardisation
# scale_rm <- function(mat, origin = 1, centre = F) {
#   ncols = ncol(mat)
#   nrows = nrow(mat)
#
#   origin_rows = which(mat$time == origin)
#
#
#   # scale first columnwise
#   cent = sapply(3:ncols, function(j) mean(mat[origin_rows,j]))
#   dist   = sapply(3:ncols, function(j) sqrt(sum((mat[,j] - cent[j-2])**2)/max(1, length(mat[,j]-1))))
#   mat[,3:ncols] = sapply(3:ncols, function(j) mat[,j]/dist[j-2])                       # faster than sweep
#
#   # centre whole dataset
#   if(centre) {
#     gen_cent = mean(as.matrix(mat[origin_rows,3:ncols]))
#     mat[,3:ncols]  = sapply(3:ncols, function(j) mat[,j] - gen_cent)
#   }
#
#
#   # cent = sapply(3:ncols, function(j) mean(mat[origin_rows,j]))
#   # if(centre) {
#   #   mat[,3:ncols] <- sapply(3:ncols, function(j) mat[,j] - cent[j-2])
#   #   # mat[,3:ncols] <- sweep(mat[3:ncols], 2L, cent, check.margin=FALSE)
#   #   dist = sapply(3:ncols, function(j) sqrt(sum(mat[,j]**2)/max(1, length(mat[,j]-1))))
#   #
#   # } else {
#   #   dist = sapply(3:ncols, function(j) sqrt(sum((mat[,j] - cent[j-2])**2)/max(1, length(mat[,j]-1))))
#   #
#   # }
#
#   # mat_st = sweep(mat[,3:ncols], 2, dist, "/", check.margin = FALSE)
#
#   # mat_st = as.matrix(cbind(mat[,1:2], mat_st), nrow = nrows, ncol = ncols)
#   # colnames(mat_st) <- colnames(mat)
#
#   return(data.frame(mat))
# }
