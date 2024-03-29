# Generated by rstantools.  Do not edit by hand.

# names of stan models
stanmodels <- c("BinConcurrent", "BinExternal", "BinHybrid", "ContConcurrent", "ContExternal", "ContHybrid", "T2EConcurrent", "T2EConcurrentC0", "T2EExternal", "T2EHybrid")

# load each stan module
Rcpp::loadModule("stan_fit4BinConcurrent_mod", what = TRUE)
Rcpp::loadModule("stan_fit4BinExternal_mod", what = TRUE)
Rcpp::loadModule("stan_fit4BinHybrid_mod", what = TRUE)
Rcpp::loadModule("stan_fit4ContConcurrent_mod", what = TRUE)
Rcpp::loadModule("stan_fit4ContExternal_mod", what = TRUE)
Rcpp::loadModule("stan_fit4ContHybrid_mod", what = TRUE)
Rcpp::loadModule("stan_fit4T2EConcurrent_mod", what = TRUE)
Rcpp::loadModule("stan_fit4T2EConcurrentC0_mod", what = TRUE)
Rcpp::loadModule("stan_fit4T2EExternal_mod", what = TRUE)
Rcpp::loadModule("stan_fit4T2EHybrid_mod", what = TRUE)

# instantiate each stanmodel object
stanmodels <- sapply(stanmodels, function(model_name) {
  # create C++ code for stan model
  stan_file <- if(dir.exists("stan")) "stan" else file.path("inst", "stan")
  stan_file <- file.path(stan_file, paste0(model_name, ".stan"))
  stanfit <- rstan::stanc_builder(stan_file,
                                  allow_undefined = TRUE,
                                  obfuscate_model_name = FALSE)
  stanfit$model_cpp <- list(model_cppname = stanfit$model_name,
                            model_cppcode = stanfit$cppcode)
  # create stanmodel object
  methods::new(Class = "stanmodel",
               model_name = stanfit$model_name,
               model_code = stanfit$model_code,
               model_cpp = stanfit$model_cpp,
               mk_cppmodule = function(x) get(paste0("rstantools_model_", model_name)))
})
