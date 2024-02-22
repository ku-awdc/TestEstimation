## Ideas for sample size stuff

library("tidyverse")
library("TeachingDemos")
library("coda")
library("runjags")

## With gold standard:

tibble(
  RefSe = 1,
  RefSp = 1,

)

## With known reference standard:

tibble(
  Hyperprior = c("Alpha","Beta"),
  Se = c(20,2),
  Sp = c(40,1)
) ->
  test_ref


tibble(
  Se = 0.8,
  Sp = 0.98
) ->
  test_comp

