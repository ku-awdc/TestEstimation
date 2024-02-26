## Ideas for sample size stuff

library("TeachingDemos")
library("coda")
library("runjags")
library("rjags")
library("tidyverse")
library("pbapply")

runjags.options(silent.jags=TRUE, silent.runjags=TRUE, inits.warning=FALSE, predraw.plots=FALSE)
theme_set(theme_light())
set.seed(2024-02-23)

## With gold standard:

tibble(
  RefSe = 1,
  RefSp = 1,
)

## With known reference standard:

tribble(~N, ~Prev,
  100, 0.2,
#  100, 0.3,
  ) ->
  populations

tribble(~TestNumber, ~Parameter, ~Alpha, ~Beta,
  1, "Se", 36+1, 4+1,
  1, "Sp", 99+1, 1+1,
) |>
  mutate(Estimate = (Alpha-1)/(Alpha+Beta-2), TestName="Reference") ->
  ref
ref

tribble(~TestNumber, ~Parameter, ~Alpha, ~Beta,
  1, "Se", 20+1, 20+1,
  1, "Sp", 95+1, 5+1,
) |>
  mutate(Estimate = (Alpha-1)/(Alpha+Beta-2), TestName="Reference") ->
  ref
ref


## TODO:
# 1. Incorporate ppp to get CI for true prevalence
# 2. Implement CI for se/sp by bootstrapping from distn. of input ref se/sp and estimated true prevalence
# 3. Speed up special case of 1 population
# 4. Implement correlated reference and comparator
# 5. Implement 2x ref tests
# 6. Re-estimate prevalence from comparator test using the estimated se and sp,
#    and check against original estimate - this gives us some estimate of correlation

source("notebooks/tests_R6.R")
test_eval <- TestEval$new(ref, 500L)
(tally <- test_eval$simulate_data(0.2, c(0.8,0.9)))[1,,]
(tally <- test_eval$simulate_data(0.2, c(0.8,0.9), reassign_prob=0.75))[1,,]
test_eval$estimate_prev()$par
test_eval$estimate_comp()$par
test_eval$re_estimate_prev()$par

test_eval$get_estimates()

move1 <- rbinom(dim(tally)[1],tally[,2,1],0.5)
move2 <- rbinom(dim(tally)[1],tally[,1,2],0.5)
tally[,1,1] <- tally[,1,1] +move1
tally[,2,1] <- tally[,2,1] -move1
tally[,1,2] <- tally[,1,2] +move2
tally[,2,2] <- tally[,2,2] -move2
test_eval$set_data(tally)

test_eval$estimate_prev()$par
test_eval$estimate_comp()$par

## Use ppp to get CI for true prevalence?
ppp_prev <- c(
  apply(test_eval$calculate_ppp(use_comp = FALSE) * apply(tally,c(1,3),sum),1,sum),
  apply((1-test_eval$calculate_ppp(use_comp = FALSE)) * apply(tally,c(1,3),sum),1,sum)
) |> as.matrix()
stopifnot(all.equal(apply(ppp_prev,2,sum), (apply(tally,c(1,3),sum) |> apply(1,sum))))
ppp_prev[1,] / apply(ppp_prev,2,sum)
qbeta(c(0.025,0.975), ppp_prev[1,1]+1, ppp_prev[2,1]+1)

test_eval$calculate_ppp()
test_eval$get_estimates()


source("notebooks/tests_R6.R")
test_eval <- TestEval$new(ref, c(500L, 250L))
(tally <- test_eval$simulate_data(c(0.1,0.2), c(0.8,0.9), reassign_prob = 0.95))
tally[1,,]
tally[2,,]
test_eval$estimate_prev()$par
test_eval$estimate_comp()$par

test_eval$re_estimate_prev()$par
test_eval$estimate_hw()$par


move1 <- rbinom(dim(tally)[1],tally[,2,1],0.5)
move2 <- rbinom(dim(tally)[1],tally[,1,2],0.5)
tally[,1,1] <- tally[,1,1] +move1
tally[,2,1] <- tally[,2,1] -move1
tally[,1,2] <- tally[,1,2] +move2
tally[,2,2] <- tally[,2,2] -move2
test_eval$set_data(tally)

test_eval$estimate_prev()$par
test_eval$estimate_comp()$par

test_eval$estimate_prev()$par
test_eval$estimate_comp()$par

## Use ppp to get CI for true prevalence?
ppp_prev <- bind_rows(
  apply(test_eval$calculate_ppp(use_comp = FALSE) * apply(tally,c(1,3),sum),1,sum),
  apply((1-test_eval$calculate_ppp(use_comp = FALSE)) * apply(tally,c(1,3),sum),1,sum)
) |> as.matrix()
stopifnot(apply(ppp_prev,2,sum) == (apply(tally,c(1,3),sum) |> apply(1,sum)))
ppp_prev[1,] / apply(ppp_prev,2,sum)
qbeta(c(0.025,0.975), ppp_prev[1,1]+1, ppp_prev[2,1]+1)

test_eval$estimate_comp(include_prev=TRUE)$par
test_eval$estimate_hw()$par
test_eval$calculate_ppp()
test_eval$get_estimates()




source("notebooks/tests_R6.R")
test_eval <- TestEval$new(ref, c(500L, 250L))
tibble(Iteration = 1:1000) |>
  rowwise() |>
  group_split() |>
  pblapply(function(x){
    test_eval$simulate_data(c(0.1,0.2), c(0.8,0.9), reassign_prob=0.75)
    #test_eval$simulate_data(c(0.1), c(0.8,0.9))
    test_eval$estimate_prev()$par
    test_eval$estimate_comp()$par
    x |> bind_cols(test_eval$get_estimates(beta_n=0))
  }, cl=6L) |>
  bind_rows() ->
  res

ggplot(res, aes(x=Parameter, y=Estimate, col=Method)) +
  geom_hline(aes(yintercept=Estimate), tibble(Parameter=c("Se","Sp"),Estimate=c(0.8,0.9))) +
  geom_boxplot()

ggplot(res2, aes(x=Parameter, y=Estimate, col=Method)) +
  geom_hline(aes(yintercept=Estimate), tibble(Parameter=c("Se","Sp"),Estimate=c(0.8,0.9))) +
  geom_boxplot()




stop("OLDER BELOW HERE, but some of it is still useful - JAGS code etc")


tribble(
  ~TestNumber, ~Se, ~Sp,
  1, 0.8, 0.9,
) ->
  comp


## Using point estimates:

ref |>
  mutate(Estimate = (Alpha-1)/(Alpha+Beta-2), Type="Reference") |>
  select(Type, TestNumber, Parameter, Estimate) |>
  pivot_wider(names_from="Parameter", values_from="Estimate") |>
  bind_rows(
    comp |> mutate(Type = "Comparator")
  ) |>
  mutate(TestIndex = str_c("Test_", row_number())) ->
  tests

se <- tests$Se
sp <- tests$Sp
p <- 1

make_probs_fun <- function(n_tests, n_pops){
  stopifnot(n_pops==1L, n_tests==2L)

  compiler::cmpfun(
    function(prev, test1, test2){
      se <- c(test1[1], test2[1])
      sp <- c(test1[2], test2[2])
      probs <- array(dim=c(1,2,2,2), dimnames = list("Pop1",c("Se","Sp"),str_c("Test1",c("-","+")),str_c("Test2",c("-","+"))))
      probs["Pop1","Se","Test1-","Test2-"] <- prev[p]*(1-se[1])*(1-se[2])
      probs["Pop1","Se","Test1+","Test2-"] <- prev[p]*se[1]*(1-se[2])
      probs["Pop1","Se","Test1-","Test2+"] <- prev[p]*(1-se[1])*se[2]
      probs["Pop1","Se","Test1+","Test2+"] <- prev[p]*se[1]*se[2]
      probs["Pop1","Sp","Test1-","Test2-"] <- (1-prev[p])*sp[1]*sp[2]
      probs["Pop1","Sp","Test1+","Test2-"] <- (1-prev[p])*(1-sp[1])*sp[2]
      probs["Pop1","Sp","Test1-","Test2+"] <- (1-prev[p])*sp[1]*(1-sp[2])
      probs["Pop1","Sp","Test1+","Test2+"] <- (1-prev[p])*(1-sp[1])*(1-sp[2])
      probs
    }
  )
}

make_ppp_fun <- function(n_tests, n_pops){
  stopifnot(n_pops==1L, n_tests==2L)

  compiler::cmpfun(
    function(prev, test1, test2){
      se <- c(test1[1], test2[1])
      sp <- c(test1[2], test2[2])
      probs <- array(dim=c(1,2,2,2), dimnames = list("Pop1",c("Se","Sp"),str_c("Test1",c("-","+")),str_c("Test2",c("-","+"))))
      probs["Pop1","Se","Test1-","Test2-"] <- prev[p]*(1-se[1])*(1-se[2])
      probs["Pop1","Se","Test1+","Test2-"] <- prev[p]*se[1]*(1-se[2])
      probs["Pop1","Se","Test1-","Test2+"] <- prev[p]*(1-se[1])*se[2]
      probs["Pop1","Se","Test1+","Test2+"] <- prev[p]*se[1]*se[2]
      probs["Pop1","Sp","Test1-","Test2-"] <- (1-prev[p])*sp[1]*sp[2]
      probs["Pop1","Sp","Test1+","Test2-"] <- (1-prev[p])*(1-sp[1])*sp[2]
      probs["Pop1","Sp","Test1-","Test2+"] <- (1-prev[p])*sp[1]*(1-sp[2])
      probs["Pop1","Sp","Test1+","Test2+"] <- (1-prev[p])*(1-sp[1])*(1-sp[2])
      probs
    }
  )
}

make_prob <- make_probs_fun(2,1)
make_prob(0.2, test1, test2)

tally <- rmultinom(1, size=populations$N[p], prob=make_prob(populations$Prev[p], test1, test2) |> apply(c(3,4), sum))

optimise(function(x) dmultinom(tally, populations$N[p], make_prob(x, test1, test2) |> apply(c(3,4), sum), log=TRUE), c(0,1), maximum=TRUE)

reftally <- tally[c(1,2)] + tally[c(3,4)]
opt <- optimise(function(x) dmultinom(reftally, populations$N[p], make_prob(x, test1, test2) |> apply(c(3), sum), log=TRUE), c(0,1), maximum=TRUE)
(estprev <- opt$maximum)


make_ppp <-

probs <- make_prob(estprev, test1, test2) |> apply(c(2,3), sum)
ppp <- probs["Se",] / (probs["Se",] + probs["Sp",])


optim(c(se=0.9, sp=0.9), function(x) -dmultinom(tally, populations$N[p], make_prob(estprev, test1, x) |> apply(c(3,4), sum), log=TRUE), method="L-BFGS-B", lower=1e-6, upper=1-1e-6)
optim(c(se=0.9, sp=0.9), function(x) -dmultinom(tally, populations$N[p], make_prob(estprev, test1, x) |> apply(c(3,4), sum), log=TRUE))



qbeta(0.5, (tally[4]*ppp[2] + tally[3]*ppp[1])+1, (tally[2]*ppp[2] + tally[1]*ppp[1])+1)
tally[4]/(tally[4]+tally[3])

qbeta(0.5, (1 - tally[4]*ppp[2] + tally[3]*ppp[1])+1, (tally[2]*ppp[2] + tally[1]*ppp[1])+1)
tally[4]/(tally[4]+tally[3])

reframe(
  Se = rbeta(beta_n, sum(PPP[Result==1])+1, sum(PPP[Result==0])+1),
  Sp = rbeta(beta_n, sum(1-PPP[Result==0])+1, sum(1-PPP[Result==1])+1)
)



## Not identifiable as df < params:
ctally <- tally[c(1,3)] + tally[c(2,4)]
optim(c(se=0.9, sp=0.9), function(x) -dmultinom(ctally, populations$N[p], make_prob(estprev, test1, x) |> apply(c(4), sum), log=TRUE), method="L-BFGS-B", lower=0.00001, upper=0.99999)
optim(c(se=0.9, sp=0.9), function(x) -dmultinom(ctally, populations$N[p], make_prob(estprev, test1, x) |> apply(c(4), sum), log=TRUE))



models <- list(
  readLines("notebooks/models/tests_1B.txt") |> str_c(collapse="\n"),
  readLines("notebooks/models/tests_2.txt") |> str_c(collapse="\n")
)



## Shows bias, and CI are not the same as when incorporating uncertainty in PPP:
cifun <- function(populations, ref, comp, include_comparator = FALSE, beta_n = 1000L){

  stopifnot(
    is.data.frame(populations),
    nrow(populations) >= 1L,
    c("N","Prev") %in% names(populations),
    populations$N >= 1L,
    populations$N %% 1 == 0,
    populations$Prev >= 0, populations$Prev <= 1
  )

  stopifnot(is_tibble(ref), c("TestNumber","Parameter","Alpha","Beta") %in% names(ref), ref[["Parameter"]] %in% c("Se","Sp"))
  stopifnot(is_tibble(comp), c("TestNumber", "Se","Sp") %in% names(comp))
  stopifnot(length(N)==1, N>0, N%%1==0L)
  stopifnot(length(prevalence)==1, prevalence>=0, prevalence<=1)

  ref <- ref |> arrange(TestNumber, Parameter=="Sp")
  populations <- populations |> mutate(Population = fct(str_c("Pop", row_number())))


  ## Simulate data:
  ref |>
    mutate(Estimate = (Alpha-1)/(Alpha+Beta-2), Type="Reference") |>
    select(Type, TestNumber, Parameter, Estimate) |>
    pivot_wider(names_from="Parameter", values_from="Estimate") |>
    bind_rows(
      comp |> mutate(Type = "Comparator")
    ) |>
    mutate(TestIndex = str_c("Test_", row_number())) ->
    tests

  stopifnot(
    is.data.frame(tests),
    nrow(tests) >= 2L,
    c("TestIndex", "Se","Sp") %in% names(tests),
    tests$Se >= 0, tests$Se <= 1,
    tests$Sp >= 0, tests$Sp <= 1
  )

  populations |>
    group_by(Population, N, Prev) |>
    reframe(Individual = str_c(Population, "_", seq_len(N))) |>
    ungroup() |>
    mutate(Status = rbinom(n(), 1L, prevalence)) |>
    expand_grid(
      tests |> select(TestIndex, Se, Sp)
    )|>
    mutate(Result = rbinom(n(), 1L, Status*Se + (1-Status)*(1-Sp))) |>
    select(Population, Individual, Status, TestIndex, Result) ->
    simdata

  simdata |>
    semi_join(
      tests |> {\(x) if(include_comparator) identity(x) else filter(x, Type=="Reference")}(),
      by="TestIndex"
    ) ->
    longrefdata

  longrefdata |>
    pivot_wider(names_from="TestIndex", values_from="Result") ->
    refdata

  (if(include_comparator){
    ref |>
      mutate(Estimate = (Alpha-1)/(Alpha+Beta-2)) |>
      bind_rows(
        comp |>
          pivot_longer(Se:Sp, names_to="Parameter", values_to="Estimate") |>
          mutate(TestNumber = TestNumber + max(ref$TestNumber), Alpha=2, Beta=1)
      )
  }else{
    ref |> mutate(Estimate = (Alpha-1)/(Alpha+Beta-2))
  }) |>
    arrange(Parameter, TestNumber) ->
    tpriors

  Alpha <- tpriors[["Alpha"]] |> matrix(ncol=2)
  Beta <- tpriors[["Beta"]] |> matrix(ncol=2)
  inits <- list(
    prev=populations$Prev,
    se=tpriors |> filter(Parameter=="Se") |> pull(Estimate),
    sp=tpriors |> filter(Parameter=="Sp") |> pull(Estimate)
  )

  ntests <- ncol(refdata)-3L
  stopifnot(nrow(refdata)==sum(populations$N))

  if(ntests == 1L){

    populations |>
      select(Population) |>
      expand_grid(Combo=c("0","1")) |>
      arrange(Population, Combo) ->
      combos

    refdata |>
      group_by(Population) |>
      summarise(Positive = sum(Test_1), .groups="drop") |>
      arrange(Population) |>
      pull(Positive) ->
      Positive

    mm <- run.jags(models[[1]], data=list(Positive=Positive, N=populations$N, P=nrow(populations), Alpha=Alpha, Beta=Beta), inits=inits, n.chain=1, adapt=1000, burnin=1000, sample=0)

  }else if(ntests == 2L){

    populations |>
      select(Population) |>
      arrange(Population) |>
      expand_grid(
        Combo = expand_grid(Test1=c("0","1"),Test2=c("0","1"))[,2:1] |> apply(1,str_c,collapse="")
      ) ->
      combos

    refdata |>
      count(Population, Test_1, Test_2) |>
      complete(Population, Test_1, Test_2, fill = list(n = 0L)) |>
      arrange(Population, Test_2, Test_1) ->
      tabdata

    stopifnot(nrow(tabdata)==4L*nrow(populations), sum(tabdata$n) == sum(populations$N))
    Tally <- tabdata |> pull(n) |> matrix(ncol=nrow(populations))

    mm <- run.jags(models[[2]], data=list(Tally=Tally, N=populations$N, P=nrow(populations), Alpha=Alpha, Beta=Beta), inits=inits, n.chain=1, adapt=1000, burnin=1000, sample=0)

  }else{
    stop("Unhandled number of reference tests")
  }

  res <- mm |> as.jags() |> jags.samples("ppp", 10000, type="mean", force.list=TRUE, progress.bar="none")

  longrefdata |>
    group_by(Population, Individual) |>
    arrange(Individual, TestIndex) |>
    summarise(Combo = str_c(Result, collapse=""), .groups="drop") |>
    full_join(
      combos |> mutate(PPP = as.numeric(res$mean$ppp)),
      by=c("Combo","Population")
    ) |>
    full_join(
      simdata,
      by=c("Individual","Population")
    ) |>
    inner_join(
      tests |> filter(Type=="Comparator"),
      by="TestIndex"
    ) |>
    group_by(TestIndex, TestNumber) |>
    reframe(
      Se = rbeta(beta_n, sum(PPP[Result==1])+1, sum(PPP[Result==0])+1),
      Sp = rbeta(beta_n, sum(1-PPP[Result==0])+1, sum(1-PPP[Result==1])+1)
    ) |>
    pivot_longer(Se:Sp, names_to="Parameter", values_to="Estimate") |>
    group_by(TestIndex, TestNumber, Parameter) |>
    summarise(Mean = mean(Estimate), Median = median(Estimate), Lower95 = coda::HPDinterval(coda::as.mcmc(Estimate))[1], Upper95 = coda::HPDinterval(coda::as.mcmc(Estimate))[2], .groups="drop") |>
    inner_join(
      tests |> pivot_longer(Se:Sp, names_to="Parameter", values_to="True"),
      by=c("TestIndex","TestNumber","Parameter")
    )




  tibble(
    Status = c(0,1),

  )

  ppp <-



  |>
    count(Population, Test1, Test2) |>
    complete(Population, Test1, Test2, fill = list(n = 0L)) |>
    arrange(Population, Test2, Test1) ->
    data
  stopifnot(nrow(data)==8L, sum(data$n) == sum(populations$N))
  return(data |> pull(n) |> matrix(ncol=2))



}



## Simulate data:
stopifnot(
  is.data.frame(populations),
  nrow(populations) >= 2L,
  c("N","Prev") %in% names(populations),
  populations$N >= 1L,
  populations$N %% 1 == 0,
  populations$Prev >= 0, populations$Prev <= 1
)
stopifnot(
  is.data.frame(tests),
  nrow(tests) >= 2L,
  c("Se","Sp") %in% names(tests),
  tests$Se >= 0, tests$Se <= 1,
  tests$Sp >= 0, tests$Sp <= 1
)

populations |>
  mutate(Population = fct(str_c("Pop", row_number()))) |>
  group_by(Population, N, Prev) |>
  reframe(Individual = str_c(Population, "_", seq_len(N))) |>
  ungroup() |>
  mutate(Status = rbinom(n(), 1L, Prev)) |>
  expand_grid(
    tests |>
      mutate(TestNum = fct(str_c("Test", row_number())))
  ) |>
  mutate(Result = rbinom(n(), 1L, Status*Se + (1-Status)*(1-Sp))) |>
  select(Population, Individual, TestNum, Result) |>
  pivot_wider(names_from="TestNum", values_from=Result) |>
  count(Population, Test1, Test2) |>
  complete(Population, Test1, Test2, fill = list(n = 0L)) |>
  arrange(Population, Test2, Test1) ->
  data
stopifnot(nrow(data)==8L, sum(data$n) == sum(populations$N))
return(data |> pull(n) |> matrix(ncol=2))

