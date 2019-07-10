#===============================================================================
# AlloHubMat::MI.R
# Entropy metrics for a pair of character vectors
# Mutual Information, Finite Size Error, Joint Entropy, normalised MI
# See also our paper https://doi.org/10.1096/fj.11-190868.
# (C) 2019 Jens Kleinjung, Jamie Macpherson, Franca Fraternali
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
#===============================================================================
## input data col1 and col2 should be character vectors
## might work with other data types, but that is untested
.mi = function(col1, col2) {
	## nput checks
	stopifnot(length(col1) > 1);
	stopifnot(is.character(col1));
	stopifnot(length(col2) > 1);
	stopifnot(is.character(col2));

	## marginal probabilities
	t1 = table(col1);
	p1 = t1 / sum(t1);
	t1.l = length(t1);

	t2 = table(col2);
	p2 = t2 / sum(t2);
	t2.l = length(t2);

	## vector of letter pair words
	pp.w = paste(c(col1), c(col2), sep = "");
	## vector of pairs ...
	t.pp = table(pp.w);
	t.pp.l = length(t.pp);
	##   ... and their marginal probabilities
	pp = t.pp / sum(t.pp);

	## vector of first letters of pairs ...
	pp.1 = substring(names(pp), 1, 1);
	##   ... and their marginal probabilities
	pp.p1 = p1[pp.1];

	## vector of second letters of pairs ...
	pp.2 = substring(names(pp), 2, 2);
	##   ... and their marginal probabilities
	pp.p2 = p2[pp.2];

	## Mutual Information
	MI = sum(pp * log(pp / (pp.p1 * pp.p2)));

	## Finite Size Error
	FSE = (t.pp.l - t1.l - t2.l + 1) / sum(length(col1), length(col2));

	## Joint Entropy
	JE = -sum(pp * log(pp));

	## normalised Mutual Information
	if (JE > 0) {
		nMI = (MI - FSE) / JE;
	} else {
		nMI = 0;
	}

	return(as.vector(c(MI, FSE, JE, nMI)));
}

#===============================================================================

