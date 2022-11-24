Rcpp::loadModule("rtreecpp", TRUE)

#' Create an RTree
#'
#' Organizes points in an R-tree.
#'
#' The R-tree is created using the quadratic splitting algorithm,
#' with the maximum number of elements
#' per node set to 16.
#' See \url{https://www.boost.org/doc/libs/1_75_0/libs/geometry/doc/html/geometry/spatial_indexes/introduction.html} for details.
#'
#' @param x A 2-column numeric matrix of point coordinates.
#'
#' @return An RTree S3 object.
#'
#' @export
RTree <- function(x) {

  if (!is.numeric(x)) {
    stop('x must be numeric.')
  }
  if (length(dim(x)) != 2 | dim(x)[2] != 2) {
    stop('x must be a 2-column matrix.')
  }

  rTreeCpp <- new(RTreeCpp, x)
  me <- list(
    rTreeCpp = rTreeCpp
  )

  class(me) <- append(class(me), "RTree")
  return(me)
}

#' Get Points Within Distance
#'
#' For each point \eqn{y_i} in set \code{y}, returns the row-indices
#' of the points indexed in \code{rTree}
#' that are within a given \code{distance} of \eqn{y_i}.
#'
#' @param rTree An \link{RTree} object.
#' @param y A 2-column numeric matrix of point coordinates.
#' @param distance A positive scalar.
#' @return A list with one entry per row of `y`.
#' Each entry is a vector of row indices.
#' @examples
#' # For two sets of points, A and B, find points in B within
#' # a distance of 0.05 to points in A.
#' n <- 10^4
#' A <- cbind(runif(n), runif(n))
#' B <- cbind(runif(n), runif(n))
#'
#' # Make an RTree on B
#' B_rtree <- RTree(B)
#'
#' # Find points of B within 0.05 of each point in A
#' nearby <- withinDistance(B_rtree, A, 0.05)
#'
#' # Each element of nearby is a list of row indices in B
#' print(nearby[[1]])
#' @export
withinDistance <- function(rTree, y, distance) {
  UseMethod("withinDistance", rTree)
}

withinDistance.RTree <- function(rTree, y, distance) {

  if (!inherits(rTree, "RTree")) {
    stop('rTree must be of class RTree.')
  }
  if (!is.numeric(y)) {
    stop('y must be numeric.')
  }
  if (length(dim(y)) != 2 | dim(y)[2] != 2) {
    stop('y must be a 2-column matrix.')
  }
  if (!is.numeric(distance)) {
    stop('distance must be numeric.')
  }
  if (length(distance) != 1) {
    stop('distance must be a scalar.')
  }
  if (distance <= 0) {
    stop('distance must be positive.')
  }

  index.ls <- rTree$rTreeCpp$within_distance_list(y, distance)

  return(index.ls)
}


#' Count Points Within Distance
#'
#' For each point \eqn{y_i} in set \code{y}, returns the number of
#' points indexed in \code{rTree}
#' that are within a given \code{distance} of \eqn{y_i}.
#'
#' @param rTree An \link{RTree} object.
#' @param y A 2-column numeric matrix of point coordinates.
#' @param distance A positive scalar.
#' @return An integer vector of with one entry for each row of `y`.
#' @examples
#' # For two sets of points, A and B, count points in B within
#' # a distance of 0.05 to points in A.
#' n <- 10^4
#' A <- cbind(runif(n), runif(n))
#' B <- cbind(runif(n), runif(n))
#'
#' # Make an RTree on B
#' B_rtree <- RTree(B)
#'
#' # Find points of B within 0.05 of each point in A
#' counts <- countWithinDistance(B_rtree, A, 0.05)
#'
#' # Each element of counts is a count of points in B
#' print(counts[[1]])
#' @export
countWithinDistance <- function(rTree, y, distance) {
  UseMethod("countWithinDistance", rTree)
}

countWithinDistance.RTree <- function(rTree, y, distance) {

  if (!inherits(rTree, "RTree")) {
    stop('rTree must be of class RTree.')
  }
  if (!is.numeric(y)) {
    stop('y must be numeric.')
  }
  if (length(dim(y)) != 2 | dim(y)[2] != 2) {
    stop('y must be a 2-column matrix.')
  }
  if (!is.numeric(distance)) {
    stop('distance must be numeric.')
  }
  if (length(distance) != 1) {
    stop('distance must be a scalar.')
  }
  if (distance <= 0) {
    stop('distance must be positive.')
  }

  index.ls <- rTree$rTreeCpp$count_within_distance_list(y, distance)

  return(index.ls)
}

#' Get Nearest Neighbors
#'
#' For each point \eqn{y_i} in set \code{y}, returns the row-indices of
#' the \code{k} points indexed in \code{rTree}
#' that are closest to \eqn{y_i}.
#'
#' @param rTree An \link{RTree} object.
#' @param y A 2-column numeric matrix of point coordinates.
#' @param k A positive integer.
#' @return A list with one entry per row of `y`.
#' Each entry is a vector of row indices.
#' @examples
#' # For two sets of points, A and B, find the 10 nearest neighbors in B
#' # for each point in A.
#' n <- 10^4
#' A <- cbind(runif(n), runif(n))
#' B <- cbind(runif(n), runif(n))
#'
#' # Make an RTree on B
#' B_rtree <- RTree(B)
#'
#' # Find 10 nearest neighbors in B for each point in A
#' nearest <- knn(B_rtree, A, 10L)
#'
#' # Each element of nearest is a vector of row indices of 10 points in B
#' print(nearest[[1]])
#' @export
knn <- function(rTree, y, k) {
  UseMethod("knn", rTree)
}

knn.RTree <- function(rTree, y, k) {

  if (!inherits(rTree, "RTree")) {
    stop('rTree must be of class RTree.')
  }
  if (!is.numeric(y)) {
    stop('y must be numeric.')
  }
  if (length(dim(y)) != 2 | dim(y)[2] != 2) {
    stop('y must be a 2-column matrix.')
  }
  if (!is.numeric(k)) {
    stop('k must be numeric.')
  }
  if (!is.integer(k)) {
    k <- as.integer(k)
    warning('k was cast to integer, this may lead to unexpected results.')
  }
  if (length(k) != 1) {
    stop('k must be a scalar.')
  }
  if (k <= 0) {
    stop('k must be positive.')
  }

  index.ls <- rTree$rTreeCpp$knn_list(y, k)

  return(index.ls)
}






