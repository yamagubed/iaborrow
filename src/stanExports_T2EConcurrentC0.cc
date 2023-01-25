// Generated by rstantools.  Do not edit by hand.

#include <Rcpp.h>
using namespace Rcpp ;
#include "stanExports_T2EConcurrentC0.h"

RCPP_MODULE(stan_fit4T2EConcurrentC0_mod) {


    class_<rstan::stan_fit<stan_model, boost::random::ecuyer1988> >("rstantools_model_T2EConcurrentC0")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<stan_model, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<stan_model, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<stan_model, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<stan_model, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<stan_model, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<stan_model, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<stan_model, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<stan_model, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<stan_model, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<stan_model, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<stan_model, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<stan_model, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<stan_model, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<stan_model, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<stan_model, boost::random::ecuyer1988> ::constrained_param_names)
    .method("standalone_gqs", &rstan::stan_fit<stan_model, boost::random::ecuyer1988> ::standalone_gqs)
    ;
}
