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
List JointBart(const IntegerVector& n, // vector of sample sizes in train
               const size_t& p, //number of variables
               const IntegerVector& np, // vector of sample sizes in test (need to remove)
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
  int *numcut = &nc[0];
  //int *grp = &igrp[0];

  double rss, restemp;
  arma::mat sdraw= arma::zeros<arma::mat>(nd+burn, K);
  arma::cube varprb = arma::zeros<arma::cube>(nd,p,K);
  arma::cube varcnt = arma::zeros<arma::cube>(nd,p,K);
  arma::cube trdraw = arma::zeros<arma::cube>(nd,n[2],K);//need to change

  //random number generation LH: May need to be modified
  //unsigned int n1=111; // additional parameters needed to call from C++
  //unsigned int n2=222;

  arn gen; // ?????
  // sigma
  std::vector<double*> svec(K);

  // start on bart
  std::vector<heterbart>  mul_bart(K);
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

      //Rprintf("_xi[0]: %.4f, _xi[]: %.4f\n", _xi[0][0], _xi[p-1][numcut[p-1]-1]);

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

    //mul_bart[k].pr();

    std::vector<double> wv(as<std::vector<double>>(w[k]));
    mul_bart[k].setw(wv);

    svec[k] = new double[n[k]];
    for(size_t j=0;j<n[k];j++)  svec[k][j] = mul_bart[k].getw(j)*sigma[k];
  }

  std::vector<double> ivarprb (p,0.);
  std::vector<size_t> ivarcnt (p,0);
  /*std::stringstream treess;  //string stream to write trees to
   treess.precision(10);*/

  /*
   * MCMC
   */
  Rprintf("\nMCMC\n");


  for(size_t iter=0;iter<(nd+burn);iter++) {
    /*
     * update probability of G/adj
     */


    /*
     * update each graph
     */
    for(size_t k=0; k<K; k++){
      mul_bart[k].draw(svec[k], gen);
      rss=0.0,restemp=0.0;
      for(size_t j=0; j<n[k]; j++){
        restemp=((mul_bart[k].gety(j))-mul_bart[k].f(j))/(mul_bart[k].getw(j));
        rss += restemp*restemp;

      }
      sigma[k] = sqrt((nu[k]*lambda[k] + rss)/gen.chi_square(n[k]+nu[k]));
      for(size_t j=0;j<n[k];j++)  svec[k][j] = mul_bart[k].getw(j)*sigma[k];
      sdraw(iter, k) = sigma[k];
      //Rprintf("sdraw %f \n", sdraw(iter, k));

      if(iter >= burn){
        for(size_t j=0;j<n[k];j++) trdraw(iter-burn, j, k)+= mul_bart[k].f(j);
        //Rprintf("trdraw %f \n", trdraw(iter-burn, 20, k));
        ivarcnt = mul_bart[k].getnv();
        ivarprb  = mul_bart[k].getpv();

        size_t iter2=(iter-burn);
        for(size_t j=0;j<p;j++){
          varcnt(iter2, j, k)=ivarcnt[j];
          //varcnt(i-burn,j)=ivarcnt[j];
          varprb(iter2, j, k)=ivarprb[j];
          //varprb(i-burn,j)=ivarprb[j];
        }
        //Rprintf("varprb %f \n", varprb(iter2, 5, k));
        //Rprintf("varcnt %f \n", varcnt(iter2, 5, k));
      }

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

  // for(size_t k=0; k<K; k++){
  //   Rprintf("\nEND\n");
  //   mul_bart[k].pr();
  //  }

  Rcpp::List ret;
  ret["sigma"]=sdraw;
  //ret["yhat.train.mean"]=trmean;
  ret["yhat.train"]=trdraw;
  //ret["varcount"]=varcount;
  ret["varcount"]=varcnt;
  ret["varprob"]=varprb;

  return ret;
}


