## Load libraries
setwd(this.path::here())
library(mev)
data <- read.csv("../data/Coputopia.csv", header = TRUE)
d12 <- mev::taildep(data = with(data, cbind(Y1, Y2)[Season == "S1",]),
             u = seq(0.8, 0.995, by = 0.01),
             method = list(eta = "betacop",
                           chi = "betacop"),
             confint = "lrt",
             plot = FALSE)
d13 <- mev::taildep(data = cbind(data$Y1, data$Y3),
                    u = seq(0.8, 0.995, by = 0.01),
                    method = list(eta = "betacop",
                                  chi = "betacop"),
                    confint = "lrt",
                    plot = FALSE)
d23 <- mev::taildep(data = cbind(data$Y2, data$Y3),
                    u = seq(0.8, 0.995, by = 0.01),
                    method = list(eta = "betacop",
                                  chi = "betacop"),
                    confint = "lrt",
                    plot = FALSE)

pdf("../figures/Figure6.pdf", width = 8, height = 4)
par(mfrow = c(1,2), mar = c(4,4,1,1), bty = "l")

matplot(x = d12$u,
        y = d12$eta,
        lty = c(1,2,2),
        col = 1,
        ylim = c(0.5,1),
        yaxs = "i",
        type = "l",
        ylab = expression(eta),
        xlab = "quantile level")
matplot(x = d13$u,
        y = d13$eta,
        lty = c(1,2,2),
        col = "gray90",
        add = TRUE,
        type = "l")
matplot(x = d23$u,
        y = d23$eta,
        lty = c(1,2,2),
        add = TRUE,
        col = "gray50",
        type = "l")
matplot(x = d12$u,
        y = d12$chi,
        lty = c(1,2,2),
        col = 1,
        ylim = c(0,0.5),
        yaxs = "i",
        type = "l",
        ylab = expression(chi),
        xlab = "quantile level")
matplot(x = d13$u,
        y = d13$chi,
        lty = c(1,2,2),
        col = "gray90",
        add = TRUE,
        type = "l")
matplot(x = d23$u,
        y = d23$chi,
        lty = c(1,2,2),
        add = TRUE,
        col = "gray50",
        type = "l")
dev.off()
