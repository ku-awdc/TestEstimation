## Sample size simulations

library("tidyverse")
library("pbapply")
source("notebooks/tests_sim.R")

#Rcpp::sourceCpp("notebooks/hw_2t.cpp")

tribble(~Population, ~Prevalence, ~SampleProportion,
  "1", 0.05, 1,
  "2", 0.1, 1,
  "3", 0.25, 1,

) ->
  populations


tribble(~Test, ~Parameter, ~Alpha, ~Beta, ~Estimate,
  "G1", "Se", NA_real_, NA_real_, 0.999,
  "G1", "Sp", NA_real_, NA_real_, 0.999,
  "A1", "Se", NA_real_, NA_real_, 0.8,
  "A1", "Sp", NA_real_, NA_real_, 0.99,
  "A2", "Se", NA_real_, NA_real_, 0.5,
  "A2", "Sp", NA_real_, NA_real_, 0.99,
  "B1", "Se", NA_real_, NA_real_, 0.9,
  "B1", "Sp", NA_real_, NA_real_, 0.95,
  "B2", "Se", NA_real_, NA_real_, 0.8,
  "B2", "Sp", NA_real_, NA_real_, 0.9,
  "C1", "Se", NA_real_, NA_real_, 0.9,
  "C1", "Sp", NA_real_, NA_real_, 0.95,
  "C2", "Se", NA_real_, NA_real_, 0.8,
  "C2", "Sp", NA_real_, NA_real_, 0.9,
  "D1", "Se", NA_real_, NA_real_, 0.8,
  "D1", "Sp", NA_real_, NA_real_, 0.9,
  "D2", "Se", NA_real_, NA_real_, 0.8,
  "D2", "Sp", NA_real_, NA_real_, 0.9,
) ->
  tests

tribble(~"TestA", ~"TestB", ~"Parameter", ~"Correlation",
  "D1", "A1", "Se", 0.75,
  "D1", "A1", "Sp", 0.75,
  "D2", "A1", "Se", 0.25,
  "D2", "A1", "Sp", 0.25,
  "D1", "A2", "Se", 0.75,
  "D1", "A2", "Sp", 0.75,
  "D2", "A2", "Se", 0.25,
  "D2", "A2", "Sp", 0.25,
) ->
  correlations

source("notebooks/tests_sim.R")
ts <- TestSim$new(populations, tests |> filter(Test %in% c("A1","C1")), correlations |> filter(TestA %in% c("A1","C1"), TestB %in% c("A1","C1")))
tally <- (ts$simulate_data(1000, FALSE))

tally <- tally |> simplify2array()
dim(tally) <- c(4,3)
tally



library("TMB")
try(setwd("notebooks/TMB"))

try(dyn.unload(dynlib("hw_2t_cor")))
compile("hw_2t_cor.cpp")
dyn.load(dynlib("hw_2t_cor"))



try(dyn.unload(dynlib("hw_2t")))
compile("hw_2t.cpp")
dyn.load(dynlib("hw_2t"))

Params <- list(selg=rep(1,2L), splg=rep(1,2L), prevlg=rep(0,3L))
Data <- list(INF=Inf, N=1L, P=3L, tally=tally, se_alpha=rep(2,2), se_beta=rep(1,2), sp_alpha=rep(2,2), sp_beta=rep(1,2))
obj <- MakeADFun(data=Data, parameters=Params, DLL="hw_2t")

library("optimx")
t1 <- system.time(p1 <- plogis(do.call("optim", obj)$par))
t2 <- system.time(p2 <- plogis(nlminb(unlist(Params), obj$fn, obj$gr)$par))
obj$method <- "CG"
t3 <- system.time(p3 <- plogis(do.call("optim", obj)$par))
t4 <- system.time(p4 <- plogis(optimr(unlist(Params), obj$fn, obj$gr, obj$he, method="nlm")$par))
t5 <- system.time(p5 <- plogis(optimr(unlist(Params), obj$fn, obj$gr, obj$he, method="bobyqa")$par))
t6 <- system.time(p6 <- plogis(optimr(unlist(Params), obj$fn, obj$gr, obj$he, method="slsqp")$par))
t7 <- system.time(p7 <- plogis(optimr(unlist(Params), obj$fn, obj$gr, obj$he, method="tnewt")$par))
t8 <- system.time(p8 <- plogis(optimr(unlist(Params), obj$fn, obj$gr, obj$he, method="CG")$par))
t1;p1
t2;p2
t3;p3
t4;p4
t5;p5
t6;p6
t7;p7
t8;p8
# checkallsolvers()

## Conclusions:
# CG, bobyqa, slsqp, tnewt, nlminb and optimr(nlm) agree (but then the last 2 are the same method??)
# BFGS is often slighly different?
# SANN and Nelder-Mead don't work very well (not shown)
# TODO: proper study to get time differences!



par; obj$fn(par)
obj$fn(obj$par)
p1 <- plogis(do.call("optim", obj)$par)
p2 <- plogis(nlminb(unlist(Params), obj$fn, obj$gr)$par)
#p1 <- c(p1, se=1-(p1[1:2]*p1[3:4]), sp=p1[1:2]*(p1[3:4]+1))
#p2 <- c(p2, se=1-(p2[1:2]*p2[3:4]), sp=p2[1:2]*(p2[3:4]+1))
p1; p2
obj$fn(p1[1:7])
obj$fn(p2[1:7])

plogis(nlminb(qlogis(p1), obj$fn, obj$gr, obj$he)$par)
plogis(nlminb(qlogis(p2), obj$fn, obj$gr, obj$he)$par)
p2

obj$par <- qlogis(p1); plogis(do.call("optim", obj)$par)
obj$par <- qlogis(p2); plogis(do.call("optim", obj)$par)
p1
p2

obj$fn(p1)
obj$fn(p2)


obj$method <- "SANN"
(p3 <- plogis(do.call("optim", obj)$par))
obj$fn(p1[1:7])
obj$fn(p2[1:7])
obj$fn(p3[1:7])



yj <- rbeta(1000,1,1)
roc <- rbeta(1000,1,1)
plot(1-(yj*roc), 1-(yj*(1-roc)))
plot.ecdf(yj*roc)
plot.ecdf(yj*(1-roc))

tp <- unlist(Params)
tp[3] <- 1.5
optim(tp, obj$fn, obj$gr, method="Nelder-Mead")$par
optim(tp, obj$fn, obj$gr, method="BFGS")$par
nlminb(tp, obj$fn, obj$gr)
do.call("optim", obj)$par

optim(obj$par, obj$fn, obj$gr, method="L-BFGS-B", lower=0.001, upper=0.999, control=list(trace=TRUE))
obj$method <- "L-BFGS-B"
obj$lower <- -Inf
obj$upper <- Inf
do.call("optim", obj)$par

par; obj$fn(par)

#obj$hessian <- TRUE
par <- do.call("optim", obj)$par
obj$fn(par)
tests |> filter(Test %in% c("A1","C1"))
populations

setwd("....//")

obj$fn(obj$par)

