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
  100, 0.3,
  ) ->
  populations

tribble(~TestNumber, ~Parameter, ~Alpha, ~Beta,
  1, "Se", 50, 50,
  1, "Sp", 100, 1,
) ->
  ref

tribble(
  ~TestNumber, ~Se, ~Sp,
  1, 0.8, 0.9,
) ->
  comp

models <- list(
  readLines("notebooks/models/tests_1B.txt") |> str_c(collapse="\n"),
  readLines("notebooks/models/tests_2.txt") |> str_c(collapse="\n")
)


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

