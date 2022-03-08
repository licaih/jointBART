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

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double mrf_C(const arma::mat& Theta, // KxK matrix of graph similarity
             const double& nui,// p vec of edge-specific prior
             const arma::mat& B // combination of edge inclusion
){
  double mrf_C = 0.0;
  double tmp2, tmp3;

  for(arma::uword i=0;i<B.n_rows;i++){
    tmp3 = accu(B.row(i));
    tmp2 = as_scalar( B.row(i)*Theta*B.row(i).t());
    mrf_C += std::exp(nui*tmp3+tmp2) ;
  }
  return mrf_C;
}

// [[Rcpp::export]]
void up_relatedness(const arma::mat& adj,
                    arma::mat& Theta,
                    const arma::vec& nu,
                    const double& alpha,
                    const double& beta,
                    const double& my_w,
                    const arma::mat& B,
                    arma::mat& accep_gamma,
                    arma::mat& accep_theta,
                    arma::mat& within_model){
  size_t K = Theta.n_rows;
  size_t p = nu.n_elem;

  double alpha_prop = 1.0, beta_prop = 1.0;
  double theta_prop, sum_over_edges, log_ar;
  arma::mat Theta_prop(K,K);

  for(size_t k =0; k<K-1;k++){
    for(size_t kprime =k+1; kprime<K;kprime++){
      /*
       * Between model move
       */
      if(Theta(k, kprime) == 0.0 ){
        theta_prop = R::rgamma(alpha_prop, 1.0/beta_prop);
      }else{
        theta_prop = 0.0;
      }
      Theta_prop = Theta;
      //Rprintf("Theta(k, kprime) %.4f \n",Theta(k, kprime));
      Theta_prop(k, kprime) = theta_prop;
      Theta_prop(kprime, k) = theta_prop;
      //Rprintf("Theta_prop(k, kprime) %.4f \n",Theta_prop(k, kprime));
      //Rprintf("Theta(k, kprime) %.4f\n",Theta(k, kprime));

      // log likelihood
      sum_over_edges = 0.0;
      for(size_t l=0;l<p;l++){
        sum_over_edges += (std::log(mrf_C(Theta, nu(l), B)) +
          2 * (theta_prop - Theta(k, kprime)) * adj(l, k) * adj(l, kprime) -
          log(mrf_C(Theta_prop, nu(l), B)));
      }

      // MH ratio
      log_ar = 0.0;
      if(theta_prop == 0.0){
        log_ar = alpha_prop*(std::log(beta_prop)) - lgamma(alpha_prop) +
          lgamma(alpha) - alpha*std::log(beta) -
          (alpha - alpha_prop)*std::log(Theta(k,kprime)) +
          (beta - beta_prop)*Theta(k, kprime) +
          sum_over_edges + std::log(1-my_w) - std::log(my_w);
        //Rprintf("0 log_ar %d %d :%f\n, Theta(k,kprime): %f\n",
        //        k, kprime, log_ar, Theta(k,kprime));

      }else{
        log_ar = alpha*std::log(beta) - lgamma(alpha) +
          lgamma(alpha_prop) - alpha_prop*std::log(beta_prop) +
          (alpha-alpha_prop)*std::log(theta_prop) -
          (beta-beta_prop)*theta_prop + sum_over_edges +
          std::log(my_w) - std::log(1-my_w);
        //Rprintf("!0 log_ar %d %d :%f\n", k, kprime, log_ar);

      }

      // accept or reject
      if(log_ar > std::log(R::runif(0.0,1.0))){
        Theta(k, kprime) = theta_prop;
        Theta(kprime, k) = theta_prop;
        accep_gamma(k, kprime) += 1;
      }

      /*
       * within model move
       */
      if(Theta(k, kprime) != 0){
        within_model(k, kprime) += 1;
        theta_prop = R::rgamma(alpha_prop, 1.0/beta_prop);
        Theta_prop = Theta;
        Theta_prop(k, kprime) = theta_prop;
        Theta_prop(kprime, k) = theta_prop;

        //log-likelihood
        sum_over_edges = 0.0;
        for(size_t l=0;l<p;l++){
          sum_over_edges += (std::log(mrf_C(Theta, nu(l), B)) +
            2 * (theta_prop - Theta(k, kprime)) * adj(l, k) * adj(l, kprime) -
            std::log(mrf_C(Theta_prop, nu(l), B)));
        }

        // MH ratio
        log_ar = (alpha - alpha_prop)*(std::log(theta_prop) -
          std::log(Theta(k, kprime))) +
          (beta- beta_prop)*(Theta(k, kprime) - theta_prop) + sum_over_edges;
        //Rprintf("within  %d %d :%f\n", k, kprime, log_ar);

        // accept or reject
        if(log_ar > std::log(R::runif(0.0,1.0))){
          Theta(k, kprime) = theta_prop;
          Theta(kprime, k) = theta_prop;
          accep_theta(k, kprime) += 1;
        }
      }
    }
  }
}


// [[Rcpp::export]]
void up_nu(arma::vec& nu,
           const arma::mat& adj,
           const arma::mat& Theta,
           const double& a,
           const double& b,
           const arma::mat& B,
           arma::vec& accep_nu){
  double a_prop = 2.0, b_prop = 4.0;
  double qu, nu_prop, log_ar;
  size_t p = nu.n_elem;

  for(size_t l=0;l<p;l++){
    qu =  R::rgamma(a_prop, 1.0/b_prop);
    nu_prop = std::log(qu) - std::log(1-qu);

    log_ar = (nu_prop - nu(l))*std::log(a - a_prop + sum(adj.row(l))) +
      std::log(mrf_C(Theta, nu(l),B)) +
      (a+b-a_prop-b_prop)*(std::log(1+std::exp(nu(l)))-
      std::log(1+std::exp(nu_prop))) -
      std::log(mrf_C(Theta, nu_prop, B));

    if(log_ar > std::log(R::runif(0.0,1.0))){
      nu(l) = nu_prop;
      accep_nu(l) += 1;
    }

  }
}


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
               const List& iXinfo, // each element is a matrix
               arma::mat& Theta, // graph similarity
               arma::mat& adj,
               arma::vec& graph_nu,
               const arma::mat& B, // combitorial of edge inclusion
               const double& graph_alpha,
               const double& graph_beta,
               const double& my_w,
               const double& graph_a,
               const double& graph_b,
               const bool& Joint){

  /*
   * initialize varaible
   */
  size_t K = n.size(); // Number of graphs
  int *numcut = &nc[0];
  //int *grp = &igrp[0];

  double rss, restemp, adjprobtmp, adjprob, adj_prop, log_ar, diffg, alpha_adj, alpha_adj_prop, sa_prop, sumtmp1;
  //double alpha_ga = 6./p, alpha_gb = 3./p;
  int totalcnt;
  arma::mat sdraw= arma::zeros<arma::mat>(nd+burn, K);
  arma::cube varprb = arma::zeros<arma::cube>(nd,p,K);
  arma::cube varcnt = arma::zeros<arma::cube>(nd,p,K);
  arma::cube trdraw = arma::zeros<arma::cube>(nd,max(n),K);//need to change

  arma::mat accep_gamma  = arma::zeros<arma::mat>(K, K);
  arma::mat accep_theta  = arma::zeros<arma::mat>(K, K);
  arma::mat within_model = arma::zeros<arma::mat>(K, K);
  arma::vec accep_nu     = arma::zeros<arma::vec>(p);
  arma::mat prob         = arma::ones<arma::mat>(p, K);
  prob = prob/p;

  //arma::mat adj          = arma::zeros<arma::mat>(p, K); //
  arma::vec prob_prop  = arma::zeros<arma::vec>(p);
  arma::vec prob1      = arma::zeros<arma::vec>(p);
  arma::cube theta_all = arma::zeros<arma::cube>(nd,K,K);
  arma::cube adj_all   = arma::zeros<arma::cube>(nd,p,K);
  arma::mat nu_all     = arma::zeros<arma::mat>(nd,p);

  //random number generation LH: May need to be modified
  //unsigned int n1=111; // additional parameters needed to call from C++
  //unsigned int n2=222;

  arn gen; // ?????
  // sigma
  std::vector<double*> svec(K);
  std::vector<double>  probvec(p);

  /*
   * set parameters
   */
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

    if(iter == burn && Joint){
      accep_gamma.zeros();
      accep_theta.zeros();
      within_model.zeros();
      accep_nu.zeros();
    }

    /*
     * update each graph
     */
    for(size_t k=0; k<K; k++){
      /*
       * update probability of G/adj
       * Gk and Gk_prime only differs in edge l
       */
      if(Joint && iter > burn/2){

        ivarcnt = mul_bart[k].getnv();
        totalcnt = 0;
        for(size_t j=0;j<p;j++){
          totalcnt += ivarcnt[j];
        }

        for(arma::uword l=0;l<p;l++){
          adjprobtmp = arma::accu(Theta.row(k) % adj.row(l)) -
            Theta(k,k)*adj(l,k);
          adjprobtmp = graph_nu(l) + 2.*adjprobtmp;
          // adjprob    = std::exp(adjprobtmp);
          // adjprob    = adjprob/(1.+adjprob);

          if(k == 1 && l == 10 && iter % 1000 == 0){
            Rprintf("adjprobtmp:%.4f\n", adjprobtmp);
            //Rprintf("curr_prob :%.4f\n", prob(l,k));
            //Rprintf("adjprob:%.4f\n", adjprob);
          }
          adj_prop = 1.-adj(l,k);

          // if(adjprob > R::runif(0.0,1.0)){
          //   adj_prop = 1;
          // }else{
          //   adj_prop = 0;
          // }
          if(k == 1 && l == 10 && iter % 1000 == 0) Rprintf("adj_propose:%.4f\n", adj_prop);
          if(adj(l,k) != adj_prop){

            diffg = adj_prop - adj(l,k);
            prob_prop    = adj.col(k);
            prob_prop(l) = adj_prop;
            prob_prop    = (prob_prop + 1.)*0.5;
            prob1        = (adj.col(k) + 1.)*0.5;
            //sa_prop    = adj_prop*alpha_ga + (1.-adj_prop)*alpha_gb + ivarcnt[l];
            //prob_prop(l) =  R::rgamma(sa_prop, 1.0);
            //if(l == 10 && iter % 1000 == 0) Rprintf("prob_prop(l):%.4e\n", prob_prop(l));

            //MH

            sumtmp1 = adjprobtmp +
              lgamma(arma::accu(prob_prop)) - lgamma(arma::accu(prob1)) +
              lgamma(prob1(l)) -
              lgamma(prob_prop(l))+
              lgamma(arma::accu(prob_prop)+totalcnt) -
              lgamma(arma::accu(prob1)+totalcnt)  +
              lgamma(prob_prop(l)+ivarcnt[l]) -
              lgamma(prob1(l)+ivarcnt[l]) ;
            log_ar  = sumtmp1;
            if(k == 1 && l == 10 && iter % 1000 == 0)
              Rprintf("log_ar:%.4f--, %.4f, %.4f, cnt:%d, total:%d \n",
                      log_ar, lgamma(arma::accu(prob_prop)),
                      lgamma(arma::accu(prob1)), ivarcnt[l], totalcnt);

            // log_ar = totalcnt*(std::log(sumtmp1) -
            //   std::log(sumtmp1-prob(l,k)+prob_prop(l))) +
            //   diffg*(-lgamma(alpha_ga) + lgamma(alpha_gb) -
            //   lgamma(alpha_gb+ ivarcnt[l]) + lgamma(alpha_ga+ ivarcnt[l]));
            // if(l == 10 && iter % 1000 == 0) Rprintf("log_ar:%.4f, %4f, %4f\n",
            //    diffg*adjprobtmp,  totalcnt*(std::log(sumtmp1) -
            //      std::log(sumtmp1-prob(l,k)+prob_prop(l))),
            //      diffg*(-lgamma(alpha_ga) + lgamma(alpha_gb) -
            //        lgamma(alpha_gb+ ivarcnt[l]) + lgamma(alpha_ga+ ivarcnt[l])));

            if(log_ar > std::log(R::runif(0.0,1.0))){
              adj(l,k)  = adj_prop;
              if(k == 1 && l == 10 && iter % 1000 == 0) Rprintf("--------accept!!!\n");
            }


            // else{
            //   /*
            //    * update probability
            //    */
            //   sa_prop      =  adj(l,k)*alpha_ga + (1.- adj(l,k))*alpha_gb + ivarcnt[l];
            //   prob_prop(l) =  R::rgamma(sa_prop, 1.0);
            //
            //   // MH
            //   sumtmp1 = 0.0;
            //   sumtmp1 = arma::accu(prob.col(k));
            //   log_ar = 0.0;
            //   log_ar = totalcnt*(std::log(sumtmp1) -
            //     std::log(sumtmp1-prob(l,k)+prob_prop(l))) -
            //     ivarcnt[l]*(std::log(prob_prop(l)) - std::log(prob(l,k))) ;
            //
            //   if(log_ar > std::log(R::runif(0.0,1.0))){
            //     prob(l,k)  = prob_prop(l);
            //     if(l == 10 && iter % 1000 == 0) Rprintf("accept new prob:%.4f\n", prob(l,k));
            //
            //   }
            // }

          }

        }

        for(size_t j=0;j<p;j++) probvec[j] = (adj(j, k) + 1.)*0.5 + (double)ivarcnt[j];
        if(k == 1 &&iter % 1000 == 0) Rprintf("--------------------probvec:%.4e!!!\n", probvec[10]);
        probvec    = gen.log_dirichlet(probvec);
        if(k == 1 &&iter % 1000 == 0) Rprintf("--------------------probvec:%.4e!!!\n", probvec[10]);
        for(size_t j=0;j<p;j++) probvec[j] = std::exp(probvec[j]);
        if(k == 1 && iter % 1000 == 0) Rprintf("--------------------prob11:%.4e!!!\n", probvec[10]);
        //prob1   = prob.col(k);
        //probvec = arma::conv_to<std::vector<double>>::from(prob1);
        // update prob vector
        mul_bart[k].setpv(probvec);
      }



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

      ivarcnt = mul_bart[k].getnv();
      ivarprb  = mul_bart[k].getpv();

      if(iter >= burn){
        for(size_t j=0;j<n[k];j++) trdraw(iter-burn, j, k)+= mul_bart[k].f(j);
        //Rprintf("trdraw %f \n", trdraw(iter-burn, 20, k));

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
    if(Joint && iter > burn/2) up_relatedness(adj, Theta, graph_nu, graph_alpha, graph_beta, my_w,
       B, accep_gamma, accep_theta, within_model);
    /*
     * edge-specific
     */
    if(Joint && iter > burn/2) up_nu(graph_nu, adj, Theta, graph_a, graph_b, B, accep_nu);


    //record nu, adj, Theta
    if(iter >= burn && Joint){
      size_t iter2=(iter-burn);

      theta_all.row(iter2) = Theta;
      adj_all.row(iter2)   = adj;
      nu_all.row(iter2)    = graph_nu.t();
    }

  }

  /*
   * end MCMC
   */
  Rprintf("\nEND\n");
  // for(size_t k=0; k<K; k++){
  //   Rprintf("\nEND\n");
  //   mul_bart[k].pr();
  // }

  Rcpp::List ret;
  ret["sigma"]=sdraw;
  //ret["yhat.train.mean"]=trmean;
  ret["yhat.train"]=trdraw;
  //ret["varcount"]=varcount;
  ret["varcount"]=varcnt;
  ret["varprob"]=varprb;
  if(Joint){
    ret["theta_all"]=theta_all;
    ret["adj_all"]=adj_all;
    ret["nu_all"]=nu_all;
    ret["accep_gamma"] = accep_gamma;
    ret["accep_nu"] = accep_nu;
    ret["accep_theta"] = accep_theta;
  }

  return ret;
}


