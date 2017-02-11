/*##########################################################################
#' GMRF-AR1 spacetime example
#' ==========================================================================
#'
#' by Mark R Payne
#' DTU-Aqua, Charlottenlund, Denmark
#' mpa@aqua.dtu.dk
#'
#' Thu Aug  7 08:58:26 2014
#'
#' Solves a simple space-time example using TMB. The spatial
#' dimension is modelled as a GMRF, while the temporal is modelled as an
#' AR1 process.
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

/* Parameter transform */
template <class Type>
Type f(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}

/* Objective function */
template<class Type>
Type objective_function<Type>::operator() ()
{ 
  using namespace density;
  
  /* DATA SECTION */ 
  DATA_ARRAY(spacegrid);     /*Grid structure info: dim(spacegrid)=c(2,ngrid) */
  DATA_VECTOR(obs);     //Observation data
  DATA_IVECTOR(spacetime_idx);   /*Spacetitme grid index corresponding to observation*/
    
  /* FIXED EFFECTS / PARAMETERS */
  PARAMETER(phi_trans);    /*delta>0 !*/
  PARAMETER(logdelta);    /*delta>0 !*/
  PARAMETER(logscale);    // scaling of the spacettime field
  PARAMETER(logsigma);   // sd deviation of the observations
  
  /* RANDOM EFFECTS  */
  PARAMETER_ARRAY(eta); // Random effects on the spacetime field 

  /* OTHER */
  Type jnll=Type(0);                        //Objective function

  /* Objective function*/
  /* --------------------*/
  //Setup GMRF density object
  Type delta=exp(logdelta);
  GMRF_t<Type> spat_nldens=GMRF(spacegrid,delta);
  
  //JNLL contribution for the spacetime (including the scaling amd AR1)
  Type phi=f(phi_trans);
  jnll += SCALE(AR1(phi,spat_nldens),exp(logscale))(eta);

  //JNLL for the observations
  for(int i=0; i<obs.size(); i++){
    jnll -= dnorm(obs[i], eta[spacetime_idx[i]-1], exp(logsigma),true);
  }
  
  //Finish up
  ADREPORT(delta);
  
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

