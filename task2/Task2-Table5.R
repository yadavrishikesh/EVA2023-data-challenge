setwd(this.path::here())
library(evgam)
load("../outputs/threshold_levels_lists.Rdata")
load("../outputs/models_and_data_imputed.Rdata")
source("Task2-functions.R")


mods_id <- c(1,18, 42:46)
thresh_id <- 1:4 
ncomplete <- nrow(na.omit(Amaurot))
probs_new <- revdbayes::rDir(n = 1, alpha = rep(1, ncomplete))
newdata <-  na.omit(Amaurot) |> dplyr::select(!Y)
if(!"outputs/C2_results.Rdata" %in% list.files("..", recursive = TRUE)){
results_thr_GP_models <- matrix(nrow = length(thresh_id), ncol = length(mods_id))
for(i in seq_along(thresh_id)){
   for(j in seq_along(mods_id)){
         ### loop over the GP models
         loss <-
            try(fun_evgam_gpd(
               fit_thr_ald = thr.levels.lists[[i]],
               ### output of the fitted thresholds models i
               model.gpd = model_gpd.list[[mods_id[j]]],
               ### the formula for the GPD models j
               data_train = Amaurot.imp.final,
               newdata = newdata,
               bayesboot = FALSE
            ))
         if(!inherits(loss, "try-error")){
            results_thr_GP_models[i,j] <- loss
         }
         print(paste0("Threshold ", i, ", model ", j, ": ", Sys.time()))
   }
}
 save(results_thr_GP_models, file = "../outputs/C2_results.Rdata")
} else{
   load("../outputs/C2_results.Rdata")
}
retlev <- matrix(unlist(unlist(results_thr_GP_models)), nrow = length(mods_id))
tab <- round(retlev, 1)
rownames(tab) <- 1:nrow(tab)
cat(
   knitr::kable(tab,
                booktabs = TRUE,
                format = "latex",
                col.names = paste0("$q_{", c(0.95, 0.96, 0.97, 0.98), "}$"),
                linesep = "",
                row.names = TRUE,
                escape = FALSE),
   file = "../tables/Table5.tex")

