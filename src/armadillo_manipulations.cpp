/* Copyright (C) 2014 - 2015  James Balamuta
 *
 * This file is part of GMWM R Methods Package
 *
 * The file uses methods in the r-to-armadillo project and is free software: you can redistribute it and/or modify it
 * under the terms of the MIT License.
 *
 * The r-to-armadillo project is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 */

#include <RcppArmadillo.h>
#include "armadillo_manipulations.h"

//' Sort Matrix by Column
//' 
//' Sorts a given matrix by a specific column while retain the elements in each row.
//' 
//' @param x   A \code{matrix} to sort
//' @param col A \code{int} that indicates the column the matrix should sort by.
//' @details 
//' The functional difference between armadillo's sort() and sort_mat() is straight forward.
//' sort() will sort each column without respect to the rows. 
//' Using sort_matrix will sort only 1 column and retain the other elements to be in the same row.
//' @return The matrix sorted by values in the specified column.
//' @keywords internal
// [[Rcpp::export]]
arma::mat sort_mat(arma::mat x, unsigned int col){
  
  arma::uvec id = arma::sort_index(x.col(col));
  
  for(unsigned int i = 0; i<x.n_cols; i++){
    arma::vec sub = x.col(i);
    x.col(i) = sub.elem(id);
  }
  
  return x;
}

//' @title Reverse Subset Column
//' @description 
//' Subsets the column by going from high indices to low (the reverse of the supported practice)
//' @param x     A \code{matrix} of dimensions M x N
//' @param start A \code{unsigned int} that indicates the starting column.
//' @param end   A \code{unsigned int} that indicates the ending column.
//' @return x A \code{matrix} with matrix rows displayed in reverse order
//' @details Consider a vector x=[[1,2],[3,4]].
//' By setting \code{start=1} and \code{end=0}, the function would output x=[[2,1],[4,1]].
//' Start and end must be valid C++ matrix locations. (e.g. matrix cols start at 0 and not 1)
//' @author JJB
//' @examples
//' x = matrix(c(1,2,3,4), nrow = 2,byrow = TRUE)
//' rev_col_subset(x, 1, 0)
//' @keywords internal
// [[Rcpp::export]]
arma::mat rev_col_subset(arma::mat x, unsigned int start, unsigned int end){
  arma::mat A = arma::mat(x.n_rows, start-end+1);
  for(unsigned int i = 0; i < start-end+1; i++){
    A.col(i) = x.col(start-i);
  }
  return A;
}

//' @title Reverse Subset Row
//' @description Subsets the row by going from high indices to low (the reverse of the supported practice)
//' @param x      A \code{matrix} of dimensions M x N
//' @param start  A \code{unsigned int} that indicates the starting row.
//' @param end    A \code{unsigned int} that indicates the ending row.
//' @return x A \code{matrix} with matrix rows displayed in reversed order
//' @details Consider a vector x=[[1,2],[3,4]], the function would output x=[[3,4],[1,2]].
//' Start and end must be valid C++ matrix locations. (e.g. matrix rows start at 0 and not 1)
//' @author JJB
//' @examples
//' x = matrix(c(1,2,3,4), nrow=2,byrow=TRUE)
//' rev_row_subset(x, 1, 0)
//' @keywords internal
// [[Rcpp::export]]
arma::mat rev_row_subset(arma::mat x, unsigned int start, unsigned int end){
  arma::mat A = arma::mat(start-end+1, x.n_cols);
  for(unsigned int i = 0; i < start-end+1; i++){
    A.row(i) = x.row(start-i);
  }
  return A;
}

//' @title Reverse Armadillo Vector
//' @description Reverses the order of an Armadillo Vector
//' @usage reverse_vec(x)
//' @param x A \code{column vector} of length N
//' @return x A \code{column vector} with its contents reversed.
//' @details Consider a vector x=[1,2,3,4,5], the function would output x=[5,4,3,2,1].
//' @author JJB
//' @examples
//' x = 1:5
//' reverse_vec(x)
//' @keywords internal
// [[Rcpp::export]]
arma::vec reverse_vec(arma::vec x) {
   std::reverse(x.begin(), x.end());
   return x;
}

//' @title Transform an Armadillo field<vec> to a matrix
//' @description Unlists vectors in a field and places them into a matrix
//' @param x A \code{field<vec>}.
//' @return A \code{mat} containing the field elements within a column.
//' @author JJB
//' @examples
//' x=rnorm(100)
//' @keywords internal
// [[Rcpp::export]]
arma::mat field_to_matrix(arma::field<arma::vec> x){
  unsigned int nx = x.n_elem;
  unsigned int row;
  if(nx > 0){
    row = x(0).n_elem;
  }else{
    row = 999999999;
  }
  arma::mat A(row,nx);
  for(unsigned int i =0; i<nx; i++){
    A.col(i) = x(i);
  }
  return A; 
}

//' @title Accumulation of Armadillo field<vec>
//' @description Sums vectors in a field into a single variable.
//' @param x A \code{field<vec>}.
//' @return An \code{mat} containing the field elements within a column.
//' @author JJB
//' @examples
//' x=rnorm(100)
//' @keywords internal
// [[Rcpp::export]]
double sum_field_vec(const arma::field<arma::vec>& x){
  unsigned int nelems = x.n_elem;
  double total_elems = 0;
  
  for(unsigned int i = 0; i < nelems; i++){
    total_elems += sum(x(i));
  }
  
  return total_elems;
}