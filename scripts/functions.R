
# bootstrapping CIs for GLMMs -------------------------------------------------------

# from https://github.com/rikudoukarthik/covid-ebirding/blob/main/scripts/functions.R

boot_conf_GLMM = function(model, 
                          new_data, # separately specify dataframe with vars for model
                          new_data_string, # string for clusterExport()
                          re_form = NA,
                          nsim = 1000,
                          pred_type = "link")
{
  
  require(tidyverse)
  require(glue)
  require(lme4)
  require(VGAM)
  require(parallel) # to parallelise bootstrap step
  
  pred_fun <- function(model) {
    predict(model, newdata = new_data, type = pred_type, re.form = re_form, 
            allow.new.levels = TRUE)
    # not specifying type = "response" because will later transform prediction along with SE
  }
  
  par_cores <- max(1, floor(detectCores()/2))
  par_cluster <- makeCluster(rep("localhost", par_cores), outfile = "log.txt")
  clusterEvalQ(par_cluster, library(lme4))
  clusterExport(par_cluster, varlist = new_data_string)
  
  print(glue("Using {par_cores} cores."))
  
  pred_bootMer <- bootMer(model, nsim = nsim, FUN = pred_fun,
                          parallel = "snow", 
                          use.u = FALSE, type = "parametric", 
                          ncpus = par_cores, cl = par_cluster)
  
  stopCluster(par_cluster)
  
  return(pred_bootMer$t)
  
}

