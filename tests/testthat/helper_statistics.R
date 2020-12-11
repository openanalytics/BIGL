library(BIGL)

data("directAntivirals", package = "BIGL")
data <- subset(directAntivirals, experiment == 1)[, c("effect", "d1", "d2")]

transforms <- list(
  "BiolT" = function(y, args) with(args, N0*exp(y*time.hours)),
  "InvBiolT" = function(T, args) with(args, 1/time.hours*log(T/N0)),
  "PowerT" = function(y, args) with(args, log(y)),
  "InvPowerT" = function(T, args) with(args, exp(T)),
  "compositeArgs" = list("N0" = 1,
                         "time.hours" = 72)
)

fit <- fitMarginals(data, transforms = transforms, method = "nlslm")
rs <- fitSurface(data, fit, transforms = transforms,
                 B.CP = 2, B.B = NULL, parallel = FALSE,
                 statistic = "both")
R <- rs$maxR$Ymean$R
reps <- with(rs$offAxisTable, tapply(effect, d1d2, length))
n1 <- rs$meanR$n1
