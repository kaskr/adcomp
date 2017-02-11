/*##########################################################################
#' GMRF with TMB
#' ==========================================================================
#'
#' by Mark R Payne
#' DTU-Aqua, Charlottenlund, Denmark
#' mpa@aqua.dtu.dk
#'
#' Tue Aug  5 09:15:38 2014
#'
#' TMB code to solves the classic FELSPINE problem that is used as an example in the mgcv
#' package (particularly in relation to soap bubble smoothing) using a GMRF. 
#
#  This work is subject to a Creative Commons "Attribution" "ShareALike" License.
#  You are largely free to do what you like with it, so long as you "attribute" 
#  me for my contribution. See the fine print at the end for exact details.
#
#  To do:
#  
#  Notes:
#
########################################################################## */
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{ 
  using namespace density;
  
  /* DATA SECTION */ 
  DATA_ARRAY(grid);     /*Grid structure info: dim(grid)=c(2,ngrid) */
  DATA_VECTOR(obs);     //Observation data
  DATA_INTEGER(n_obs);  //Nunber of observations
  DATA_IVECTOR(obs_idx);   /*Grid index corresponding to observation*/
  
  /* FIXED EFFECTS / PARAMETERS */
  PARAMETER(logdelta);    /*delta>0 !*/
  PARAMETER(logscale);      // scaling of the GMRF
  PARAMETER(logSigma);   // sd deviation of the observations
  
  /* RANDOM EFFECTS  */
  PARAMETER_VECTOR(eta); /* Random effects on the GMRF, length(eta)=ngrid */

  /* OTHER */
  Type jnll=Type(0);                        //Objective function

  /* Objective function*/
  /* --------------------*/
  //Setup GMRF density object 
  GMRF_t<Type> nldens=GMRF(grid,exp(logdelta));
  
  //JNLL contribution for the grid (including the scaling)
  jnll=Type(0);
  jnll += SCALE(nldens,exp(logscale))(eta);

  //JNLL for the observations
  for(int i=0; i<n_obs; i++){
    jnll -= dnorm(obs[i], eta(obs_idx[i]-1), exp(logSigma),true);
  }
  
  return jnll;
}

/*
#' This work by Mark R Payne is licensed under a  Creative Commons
#' Attribution-NonCommercial-ShareAlike 3.0 Unported License. 
#' For details, see http://creativecommons.org/licenses/by-nc-sa/3.0/deed.en_US
#' Basically, this means that you are free to "share" and "remix" for 
#' non-commerical purposes as you see fit, so long as you "attribute" me for my
#' contribution. Derivatives can be distributed under the same or 
#' similar license.
#'
#' This work comes with ABSOLUTELY NO WARRANTY or support.
#'
#' This work should also be considered as BEER-WARE. For details, see
#' http://en.wikipedia.org/wiki/Beerware
#' 
*/

