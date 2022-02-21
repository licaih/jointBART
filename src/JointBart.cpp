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
#include <RcppArmadillo.h>
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
               NumericVector& tau, //sigma prior, vec
               NumericVector& nu, // sigma prior, vec
               NumericVector& lambda, // The scale of the prior for variance, vec
               NumericVector& sigma, // he prior for the error variance, vect
               const List& w, // each element Vector of weights which multiply the standard deviation.
               const bool& dart, //
               const double& theta, //
               const double& omega, //
               IntegerVector& igrp,
               const double& a,
               const double& b,
               const double& rho,
               const bool& aug,
               const List& iXinfo // each element is a matrix
                 ){

  /*
   * initialize varaible
   */
  size_t K = n.size(); // Number of graphs
  Rprintf("K:%d\n", K);
  int *numcut = &nc[0];
  int *grp = &igrp[0];

  NumericVector sdraw(nd+burn);
  arma::cube varprb = arma::zeros<arma::cube>(nd,p,K);
  arma::cube varcnt = arma::zeros<arma::cube>(nd,p,K);

  //random number generation LH: May need to be modified
  //unsigned int n1=111; // additional parameters needed to call from C++
  //unsigned int n2=222;

  arn gen; // ?????
  // sigma
  std::vector<double*> smat(K);

  // start on bart
  std::vector<bart>  mul_bart(K);
  for(size_t k=0; k<K; k++){
    mul_bart[k].setm(m);

    /*
     * set cut off point
     */
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

    /*
     * set other parameters
     */
    NumericVector xv(as<NumericVector>(x[k]));
    //Rprintf("xv[k]: %.4f, xv[n[k]*p-1]: %.4f\n", xv[0], xv[n[k]*p-1]);
    double *ix = &xv[0];

    NumericVector yv(as<NumericVector>(y[k]));
    //Rprintf("yv[k]: %.4f, yv[n[k]-1]: %.4f\n", yv[0], yv[n[k]-1]);
    double *iy = &yv[0];

    mul_bart[k].setprior(alpha,mybeta,tau[k]);
    mul_bart[k].setdata(p,n[k],ix,iy,numcut);
    mul_bart[k].setdart(a,b,rho,aug,dart,theta,omega); // can be removed

    mul_bart[k].pr();

    std::vector<double> wv(as<std::vector<double>>(w[k]));
    mul_bart[k].setw(wv);

    smat[k] = &wv[0];
    for(size_t j=0;j<n[k];j++) smat[k][j] = mul_bart[k].getw(j)*sigma[k]; //
    //Rprintf("smat 0:%f 3:%f \n", smat[k][0], smat[k][3]);
    //Rprintf("wmat 0:%f 3:%f \n", mul_bart[k].getw(0), mul_bart[k].getw(3));

  }

  /*std::stringstream treess;  //string stream to write trees to
  treess.precision(10);*/

  /*
   * MCMC
   */
  Rprintf("\nMCMC\n");

  /*
   * test set probability vector
   * std::vector<double> s = {0.2, 0.01, 0.02,0.03,0.05,0.04,0.02, 0.05,0.09, 0.3};
   std::vector<double> s2 = {0.4, 0.01, 0.02,0.03,0.05,0.04,0.02, 0.05,0.09, 0.2};
   std::vector<double> s3 = {0.2, 0.01, 0.52,0.03,0.05,0.04,0.02, 0.05,0.09, 0.3};

   mul_bart[0].setpv(s);
   mul_bart[1].setpv(s2);
   mul_bart[2].setpv(s3);
   */


  for(size_t iter=0;iter<(nd+burn);iter++) {
    /*
     * update probability of G/adj
     */


    /*
     * update each graph
     */
    for(size_t k=0; k<K; k++){
     mul_bart[k].draw(*smat[k], gen);

      double rss=0.0,restemp=0.0;
      for(size_t j=0; j<n[k]; j++){
        //Rprintf("y %d:%f \n", j, mul_bart[k].gety(j));
        restemp=((mul_bart[k].gety(j))-mul_bart[k].f(j))/(mul_bart[k].getw(j));
        rss += restemp*restemp;
        Rprintf("y j:%d 3:%f \n", j, mul_bart[k].gety(j));

        Rprintf("f j:%d 3:%f \n", j, mul_bart[k].f(j));
        Rprintf("w j:%d 3:%f \n", j, mul_bart[k].getw(j)); // NA values
        Rprintf("rss k:%d 3:%f \n", k, rss);

      }
      sigma[k] = sqrt((nu[k]*lambda[k] + rss)/gen.chi_square(n[k]+nu[k]));
      //Rprintf("chisq k:%d 3:%f \n", k, gen.chi_square(n[k]+nu[k]));
      Rprintf("sigma k:%d 3:%f \n", k, sigma[k]);
      //for(size_t j=0;j<n[k];j++) smat[k][j] = mul_bart[k].getw(j)*sigma[k];
      //Rprintf("smat 0:%f 3:%f \n", smat[k][0], smat[k][3]);
      //Rprintf("wmat 0:%f 3:%f \n", mul_bart[k].getw(0), mul_bart[k].getw(3));
    }

    /*
     * reladeness
     */

    /*
     * edge
     */
  }

  /*
   * end MCMC
   */

  /*
   *   for(size_t k=0; k<K; k++){
   Rprintf("\nEND\n");
   mul_bart[k].pr();
   }
   */


}


