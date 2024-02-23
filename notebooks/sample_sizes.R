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

tribble(~TestNumber, ~Parameter, ~Alpha, ~Beta,
  1, "Se", 20, 2,
  1, "Sp", 40, 1,
) ->
  ref


tribble(
  ~TestNumber, ~Se, ~Sp,
  1, 0.8, 0.9,
) ->
  comp

models <- list(
  readLines("notebooks/models/tests_1.txt") |> str_c(collapse="\n")
)


cifun <- function(N, prevalence, ref, comp, include_comparator = FALSE, beta_n = 1000L){

  stopifnot(is_tibble(ref), c("TestNumber","Parameter","Alpha","Beta") %in% names(ref), ref[["Parameter"]] %in% c("Se","Sp"))
  stopifnot(is_tibble(comp), c("TestNumber", "Se","Sp") %in% names(comp))
  stopifnot(length(N)==1, N>0, N%%1==0L)
  stopifnot(length(prevalence)==1, prevalence>=0, prevalence<=1)

  ref <- ref |> arrange(TestNumber, Parameter=="Sp")

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

  tibble(Individual = seq_len(N)) |>
    mutate(Status = rbinom(n(), 1L, prevalence)) |>
    expand_grid(
      tests |> select(TestIndex, Se, Sp)
    )|>
    mutate(Result = rbinom(n(), 1L, Status*Se + (1-Status)*(1-Sp))) |>
    select(Individual, Status, TestIndex, Result) ->
    simdata

  simdata |>
    semi_join(
      tests |> filter(Type=="Reference"),
      by="TestIndex"
    ) ->
    refdata

  if(include_comparator){

  }


  if(nrow(ref)/2L==1L){

    combos <- tibble(Combo=c("0","1"))

    stopifnot(nrow(refdata)==N)
    Positive <- sum(refdata[["Result"]])
    Alpha <- ref[["Alpha"]] |> array(dim=c(1,2))
    Beta <- ref[["Beta"]] |> array(dim=c(1,2))

    mm <- run.jags(models[[1]], data=list(Positive=Positive, N=N, Alpha=Alpha, Beta=Beta), adapt=1000, burnin=1000, sample=0)

  }else if(nrow(ref)/2L==2L){

  }else{
    stop("Unhandled number of reference tests")
  }

  res <- mm |> as.jags() |> jags.samples("ppp", 5000, type="mean", force.list=TRUE, progress.bar="none")

  refdata |>
    group_by(Individual) |>
    arrange(Individual, TestIndex) |>
    summarise(Combo = str_c(Result, collapse=""), .groups="drop") |>
    full_join(
      combos |> mutate(PPP = apply(res$mean$ppp,1,mean)),
      by="Combo"
    ) |>
    full_join(
      simdata,
      by="Individual"
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
    summarise(Mean = mean(Estimate), Median = median(Estimate), Lower95 = coda::HPDinterval(coda::as.mcmc(Estimate))[1], Upper95 = coda::HPDinterval(coda::as.mcmc(Estimate))[2], .groups="drop")




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

