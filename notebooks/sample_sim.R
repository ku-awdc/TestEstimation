## Sample size simulations

library("tidyverse")
library("pbapply")
source("notebooks/tests_sim.R")

#Rcpp::sourceCpp("notebooks/hw_2t.cpp")


library("TMB")
try(setwd("notebooks"))
try(dyn.unload(dynlib("hw_2t")))
compile("hw_2t.cpp")
dyn.load(dynlib("hw_2t"))
setwd("../")

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

Params <- list(se=rep(0.75,2L), sp=rep(0.75,2L), covse=0.0, covsp=0.0, prev=rep(0.5,3L))
Data <- list(N=1L, P=3L, tally=tally)
obj <- MakeADFun(data=Data, parameters=Params, DLL="hw_2t")
obj$hessian <- TRUE
do.call("optim", obj)

obj$fn(obj$par)

rm(ts); rm(TestSim); gc()
source("notebooks/tests_sim.R")
ts <- TestSim$new(populations, tests, correlations)
tally <- (ts$simulate_data(1000))
ts$method_A()
ts$method_B()
ts$method_C("m")$Pars
#ts$method_C("p")
ts$method_D()$Pars

source("notebooks/tests_sim.R")
ts <- TestSim$new(populations, tests |> filter(Test!="B"),  correlations)
tally <- (ts$simulate_data(1000, FALSE))
ts$method_A()
ts$method_B()
ts$method_C("m")$Pars
ts$method_C("p")$Pars
ts$method_C("n")$Pars
ts$method_D()$Pars

ts <- TestSim$new(populations, tests, correlations)


source("notebooks/tests_sim.R")

using <- c("A1","B1","D1")
tt <- tests |> filter(Test %in% using)
tc <- correlations |> filter(TestA %in% using, TestB %in% using)
ntests <- nrow(tt)/2L
ts <- TestSim$new(populations, tt, tc)
ts$simulate_data(100L, FALSE)
iterations <- 250L

expand_grid(TotalN = c(100L,250L,500L), Iteration = seq_len(iterations)) |>
  rowwise() |>
  group_split() |>
  pblapply(function(x){
    ts$simulate_data(x$TotalN, FALSE)
    x |>
      bind_cols(
        bind_rows(
          ts$method_A(),
          ts$method_B(),
          if(ntests==2) ts$method_C("m"),
          if(FALSE && ntests==2) ts$method_C("p"),
          ts$method_C("n"),
          if(FALSE && ntests==2) ts$method_D()
        )
      )
  }, cl=6L) |>
  bind_rows() ->
  res

res |>
  pivot_longer(Se:Sp, names_to="Parameter", values_to="Estimate") |>
  ggplot(aes(x=factor(TotalN), y=Estimate, col=Method, group=str_c(TotalN,Method))) +
  geom_hline(aes(yintercept=Estimate), tt |> filter(Test%in%c("C1","C2","D1","D2"))) +
  geom_boxplot() +
  facet_wrap(~Parameter)
