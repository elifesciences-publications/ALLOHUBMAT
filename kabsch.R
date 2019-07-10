#===============================================================================
# AlloHubMat::kabsch.R
# Kabsch Algorithm
# (C) 2019 Jamie Macpherson, Jens Kleinjung, Franca Fraternali
#
# This file is part of AlloHubMat.
#
#    AlloHubMat is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Foobar is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Foobar.  If not, see <https://www.gnu.org/licenses/>.
#' 
#' Aligns two sets of points via rotation and translation after
#'   solving the superpositioning problem in an Eigenvalue Ansatz.
#' 
#' Given two sets of points, with one specified as the reference set,
#' the other set will be rotated so that the RMSD between the two is minimized.
#' The format of the matrix is that there should be one row for each of
#' n observations, and the number of columns, d, specifies the dimensionality
#' of the points. The point sets must be of equal size and with the same
#' ordering, i.e. point one of the second matrix is mapped to point one of
#' the reference matrix, point two of the second matrix is mapped to point two 
#' of the reference matrix, and so on.
#'	 
#' @param P n x d matrix of reference points.
#' @param Q n x d matrix of points to align to to \code{pm}
#' @return Matrix \code{qm} rotated and translated so that the ith point 
#'	is aligned to the ith point of \code{pm} in the least-squares sense.
#' @references
#' \url{https://en.wikipedia.org/wiki/Kabsch_algorithm}
#===============================================================================

kabsch = function(Q, P){

	## center objects
	# center Q
	Qc = scale(Q, center=T, scale=F)

	# center P
	Pc = scale(P, center=T, scale=F)

	## determine the cross-covariance matrix
	C = crossprod(Pc, Qc)

	## compute a single value decomposition of the
	## cross-covariance matrix
	C.svd = svd(C)

	## use the sign of the determinant to ensure a right-hand coordinate system
	## (this is only required if the fragments are stereogenic systems, which
	## they aren't)
	d = det(C.svd$u) * det(C.svd$v)

	## determine the optimal rotation matrix
	U = tcrossprod(C.svd$u, C.svd$v)

	## rotate matrix P unto Q
	Prot = Pc %*% U

	## determine the RMSD between rotated P and Q
	# squared differences
	Dsr.sq = (Prot - Qc)**2

	rmsd = sqrt(sum(apply(Dsr.sq, 1, sum)) / dim(Dsr.sq)[1]);

	## spit out the optimal fit RMSD between Q and P
	return(rmsd)
}

#===============================================================================

