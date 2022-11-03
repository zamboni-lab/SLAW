
Rcpp::loadModule("rtreecpp2", TRUE)

#' Create an RTree
#'
#' Organizes points in an R-tree.
#'
#' The R-tree is created using the quadratic splitting algorithm, with the maximum number of elements
#' per node set to 16. See \url{http://www.boost.org/doc/libs/1_63_0/libs/geometry/doc/html/geometry/spatial_indexes/introduction.html} for details.
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

  rTreeCpp <- new(RTreeCpp2, x)
  me <- list(
    rTreeCpp = rTreeCpp
  )

  class(me) <- append(class(me), "RTree")
  return(me)
}

#' Get Points Within Distance
#'
#' For each point \eqn{y_i} in set \code{y}, returns the row-indices of the points indexed in \code{rTree}
#' that are within a given \code{distance} of \eqn{y_i}.
#'
#' @param rTree An \link{RTree} object.
#' @param y A 2-column numeric matrix of point coordinates.
#' @param distance A positive scalar.
#'
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

#' Get Points Within Distance
#'
#' For each point \eqn{y_i} in set \code{y}, returns the row-indices of the points indexed in \code{rTree}
#' that are within a given \code{distance} of \eqn{y_i}.
#'
#' @param rTree An \link{RTree} object.
#' @param y A 2-column numeric matrix of point coordinates.
#' @param distance A positive scalar.
#'
#' @export
withinBox <- function(rTree, y, dx, dy) {
  UseMethod("withinBox", rTree)
}

withinBox.RTree <- function(rTree, y, dx, dy) {

  if (!inherits(rTree, "RTree")) {
    stop('rTree must be of class RTree.')
  }
  if (!is.numeric(y)) {
    stop('y must be numeric.')
  }
  if (length(dim(y)) != 2 | dim(y)[2] != 2) {
    stop('y must be a 2-column matrix.')
  }
  if (!is.numeric(dx)) {
    stop('dx must be numeric.')
  }
  if (length(dx) != 1) {
    stop('dx must be a scalar.')
  }
  if (dx <= 0) {
    stop('dx must be positive.')
  }

  if (!is.numeric(dy)) {
    stop('dy must be numeric.')
  }
  if (length(dy) != 1) {
    stop('dy must be a scalar.')
  }
  if (dy <= 0) {
    stop('dy must be positive.')
  }


  index.ls <- rTree$rTreeCpp$within_box_list(y, dx, dy)

  return(index.ls)
}


#' Get Nearest Neighbors
#'
#' For each point \eqn{y_i} in set \code{y}, returns the row-indices of the \code{k} points indexed in \code{rTree}
#' that are closest to \eqn{y_i}.
#'
#' @param rTree An \link{RTree} object.
#' @param y A 2-column numeric matrix of point coordinates.
#' @param k A positive integer.
#'
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






