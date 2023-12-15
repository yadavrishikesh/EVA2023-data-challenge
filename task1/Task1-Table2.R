devtools::install_github("nicolagnecco/erf")
library(erf)
library(quantregForest)
library(evgam)

setwd(this.path::here())
load("../outputs/Amaurot_imputed.Rdata")

Y <- Amaurot.imp.final$Y
V1 <- Amaurot.imp.final$V1
V2 <- Amaurot.imp.final$V2
V3 <- Amaurot.imp.final$V3
V4 <- Amaurot.imp.final$V4
WindDirection <- Amaurot.imp.final$WindDirection
WindSpeed <- Amaurot.imp.final$WindSpeed
Atmosphere <- Amaurot.imp.final$Atmosphere
Season <- ifelse(Amaurot.imp.final$Season == "S2", 1, 0)
CP <- ifelse(Amaurot.imp.final$CP == "1", 1, 0)

data <- data.frame(Y = Y,
                   V1 = V1,
                   V2 = V2,
                   V3 = V3,
                   V4 = V4,
                   Season = Season,
                   WindDirection = WindDirection,
                   WindSpeed = WindSpeed,
                   Atmosphere = Atmosphere,
                   CP = CP)

#-------------------------------------------------------------------------------

prob.choices <- seq(0.95, 0.99, 0.01)

X <- cbind(V1, V2, V3, V4, WindDirection, WindSpeed, Atmosphere, Season, CP)

# quantregForest

quantregforest.object <- quantregForest(X, Y, nthreads = 4, keep.inbag = TRUE)

thr.levels.quantregforest <- predict(quantregforest.object, what = prob.choices)

props.quantregforest <- as.vector(round(100 * colMeans(Y > thr.levels.quantregforest), 2))

# ERF

erf.object <- erf(X, Y)

thr.levels.erf <- predict(erf.object, newdata = X, quantiles = prob.choices)

props.erf <- as.vector(round(100 * colMeans(Y > thr.levels.erf), 2))

# evgam

thr.levels.evgam <- sapply(prob.choices, function(zeta){
  model_thr_ald <- list(Y ~ Season + CP + s(V1, bs = "cs") + 
                          s(V2, bs = "cs") + s(V3, bs = "cs") + s(V4, bs = "cs") +
                          s(WindDirection, bs = "cs") + s(WindSpeed, bs = "cs") + 
                          s(Atmosphere, bs = "cs"),
                        ~ Season + CP + s(V1, bs = "cs") + 
                          s(V2, bs = "cs") + s(V3, bs = "cs") + s(V4, bs = "cs") +
                          s(WindDirection, bs = "cs") + s(WindSpeed, bs = "cs") + 
                          s(Atmosphere, bs = "cs"))
  
  fit_thr_ald <- evgam(model_thr_ald, data, family = "ald", ald.args = list(tau = zeta))
  fit_thr_ald$location$fitted})

props.evgam <- as.vector(round(100 * colMeans(Y > thr.levels.evgam), 2))

save(thr.levels.quantregforest, props.quantregforest,
     thr.levels.erf, props.erf,
     thr.levels.evgam, props.evgam,
     file = "../outputs/thrmodel_comparison.Rdata")

#-------------------------------------------------------------------------------

props.table <- rbind(props.quantregforest, props.erf, props.evgam)

props.table <- apply(props.table, 1:2, function(x){sprintf(x, fmt = "%.2f")})
rownames(props.table) <- c("\\texttt{quantregForest}", "\\texttt{erf}", "\\texttt{evgam}")
colnames(props.table) <- c("$q_{0.95}$","$q_{0.96}$","$q_{0.97}$","$q_{0.98}$","$q_{0.99}$")
cat(knitr::kable(props.table, 
                 format = "latex", 
                 escape = FALSE, 
                 booktabs = TRUE),
    file = "../tables/Table2.tex")
