/*
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2017 Robert McCulloch and Rodney Sparapani
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/GPL-2
 */

#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"
#include "bart.h"
#include "heterbart.h"


using namespace Rcpp;


/*void JointBart(){ //multiple class of bart done


}*/

// [[Rcpp::export]]
void JointBart(const IntegerVector& n, // vector of sample sizes in train
               const size_t& p, //number of variables
               const IntegerVector& np, // vector of sample sizes in test
               const List& x, //x, train, each element p x n_k
               const List& y, //y, train, each element n_k
               const List& xp, //x, test, each element p x np_k,
               const size_t& m, //number of tree
               IntegerVector& nc, //number of cut points same across all the graph
               const size_t& nd, //number of kept draws
               const size_t& burn, //number of burn
               const double& mybeta, //tree prior powe r
               const double& alpha, //tree prior base
               const NumericVector& tau, //sigma prior, vec
               const NumericVector& nu, // sigma prior, vec
               const NumericVector& lambda, // The scale of the prior for variance, vec
               const NumericVector& sigest, // he prior for the error variance, vect
               const List& w, // each element Vector of weights which multiply the standard deviation.
               const bool& dart, //
               const double& theta, //
               const double& omega, //
               IntegerVector& igrp,
               const double& a,
               const double& b,
               const double& rho,
               const bool& aug,
               const size_t& keeptrain,
               const size_t& keeptest,
               const size_t& keeptestme,
               const size_t& keeptreedraws,
               const size_t& printevery,
               const List& iXinfo // each element is a matrix
                 ){

  size_t K = n.size(); // Number of graphs
  Rprintf("K:%d\n", K);

  /*NumericVector xv(as<NumericVector>(x[k])); //k
  Rprintf("xv[0]: %.4f, xv[n[k]-1]: %.4f\n", xv[0], xv[n[k]*p-1]);
  double *ix = &xv[0];*/

  int *numcut = &nc[0];//

  int *grp = &igrp[0];

  //random number generation LH: May need to be modified
  arn gen;

  // start on bart
  std::vector<bart>  mul_bart(K);
  for(size_t k=0; k<K; k++){
    mul_bart[k].setm(m);
    if(iXinfo.size()>0){
      NumericMatrix Xinfo(as<NumericMatrix>(iXinfo[k]));
      xinfo _xi;
      _xi.resize(p);
      for(size_t i=0;i<p;i++) {
        _xi[i].resize(numcut[i]);
        //Rcpp::IntegerVector cutpts(Xinfo[i]);
        for(size_t j=0;j<numcut[i];j++) _xi[i][j]=Xinfo(i, j);
      }

      Rprintf("_xi[0]: %.4f, _xi[]: %.4f\n", _xi[0][0], _xi[p-1][numcut[p-1]-1]);

      mul_bart[k].setxinfo(_xi);
    }

    mul_bart[k].pr();

  }




}


