library("tidyverse")
library("R6")
library("compiler")
library("coda")

library("optimx")
# install.packages(checkallsolvers())


# This is a helper class specifically for the simulation study

# Inputs to c'tor are:
# Scenario number 1-6
# Population prevalences and sample proportions
# Test performances - either 2 or 3 tests (maybe with Alpha/Beta for ref test(s), which is unused for now) and correlation strength vs test #1 (expressed as max of possible, where currently only 1 can be non-zero)

# Iterate method simulates data, fits 4 methods, returns estimates
# Simple simulation only i.e. Dendukuri method
# Start with optim but look at TMB/optimr for speed

cor_to_cov <- function(test, cor){
  stopifnot(length(test)==2L, length(cor)==1L)
  if(cor == 0){
    0.0
  }else if(cor < 0){
    mincor <- (test[1L]-1)*(1-test[2L])
    -cor * mincor
  }else{
    maxcor <- pmin(test[1L],test[2L])-test[1L]*test[2L]
    cor * maxcor
  }
}
find_iss <- function(se, sp, cor_p, cor_n){
  # TODO: find limits for cor_n and cor_p and impose at start of function
  # Note: cor_n can't be equal to cor_p but that should be OK as se/sp limits it anyway
  stopifnot(cor_n!=cor_p)

  isp <- (cor_n*(-se) + cor_n - sp*cor_p) / (cor_n-cor_p)
  ise <- (-cor_n*se + se - sp*cor_p + sp + cor_p -1) / (cor_p-cor_n)

  stopifnot(isp>=0, isp<=1, ise>=0, isp<=1)

  #re_est_se <- cor_p*ise + (1-cor_p)*(1-isp)
  #re_est_sp <- (1-cor_n)*isp + cor_n*(1-ise)
  c(ise=ise, isp=isp) #, re_est_se=re_est_se, re_est_sp=re_est_sp)
}
#find_iss(0.8, 0.9, 0.9, 0.1)
find_cors <- function(se, sp, ise, isp){
  # TODO: find limits for ise and isp and impose at start of function
  # Note: ise+isp must be >0
  stopifnot((ise+isp) >0, (se+sp) >0)

  cor_p <- (isp + se -1) / (isp + ise -1)
  cor_n <- (isp-sp) / (isp+ise-1)

  stopifnot(cor_p>=0, cor_p<=1, cor_n>=0, cor_n<=1)

  #re_est_se <- cor_p*ise + (1-cor_p)*(1-isp)
  #re_est_sp <- (1-cor_n)*isp + cor_n*(1-ise)
  c(cor_p=cor_p, cor_n=cor_n)#, re_est_se=re_est_se, re_est_sp=re_est_sp)
}
#find_cors(0.8, 0.5, 0.8, 0.9)
find_ss <- function(ise, isp, cor_p, cor_n){
  se <- cor_p*ise + (1-cor_p)*(1-isp)
  sp <- (1-cor_n)*isp + cor_n*(1-ise)
  c(se=se, sp=sp)
}
#find_ss(1, 0.9, 0.9, 0.9)
create_hw_function <- function(n_tests, cov=FALSE){

  str_c("Test",seq_len(n_tests)) |>
    as.list() |>
    set_names() |>
    lapply(\(x) c(0,1)) |>
    do.call("expand_grid", args=_) |>
    (\(x) mutate(x, Combo=apply(x,1,str_c,collapse="")))() |>
    select(Combo, everything()) ->
    status

  status |>
    mutate(Comp = Test1) |>
    pivot_longer(!c(Combo,Comp), names_to="Test", values_to="Status") |>
    mutate(Index = str_replace(Test, "Test", "")) |>
    mutate(SeText = if_else(Status==1L, str_c("se[",Index,"]"), str_c("(1-se[",Index,"])"))) |>
    mutate(SpText = if_else(Status==1L, str_c("(1-sp[",Index,"])"), str_c("sp[",Index,"]"))) |>
    mutate(CovSeText = case_when(
      !cov ~ "",
      Test=="Test1" ~ "",
      Status==Comp ~ str_c("+dse[", Index, "]"),
      Status!=Comp ~ str_c("-dse[", Index, "]"),
    )) |>
    mutate(CovSpText = case_when(
      !cov ~ "",
      Test=="Test1" ~ "",
      Status==Comp ~ str_c("+dsp[", Index, "]"),
      Status!=Comp ~ str_c("-dsp[", Index, "]"),
    )) |>
    group_by(Combo) |>
    mutate(SeText = str_c(str_c(SeText, collapse="*"), " ", str_c(CovSeText, collapse=""))) |>
    mutate(SpText = str_c(str_c(SpText, collapse="*"), " ", str_c(CovSpText, collapse=""))) |>
    ungroup() |>
    select(-CovSeText, -CovSpText, -Index, -Comp) |>
    pivot_wider(names_from=Test, values_from=Status) |>
    # Re-order to match the Tally order:
    arrange(!!!rlang::parse_exprs(str_c("Test",rev(seq_len(n_tests))))) ->
    test_txt

  if(FALSE){
  str_c("T",seq_len(n_tests)) |>
    as.list() |>
    lapply(\(x) str_c(x, "_", c(0,1))) |>
    c(list(str_c("P",seq_len(n_pops)), c("Se","Sp"))) ->
    dn
  }

  fun_txt <- c(
    if(cov) "function(se, sp, dse, dsp){" else "function(se, sp){",
    #str_c("stopifnot(length(prev)==", n_pops, "L, length(se)==", n_tests, "L, length(sp)==", n_tests, "L)"),
    #str_c("stopifnot(length(corse)==", n_tests, "L, length(corsp)==", n_tests, "L)"),
    if(cov) "if(!all(dse>-1, dse<1, dsp>-1, dsp<1, se>=0, se<=1, sp>=0, sp<=1)) return(FALSE)",
    if(!cov) "if(!all(se>=0, se<=1, sp>=0, sp<=1)) return(FALSE)",
    str_c("se_probs=sp_probs <- numeric(", nrow(test_txt), "L)"),
    str_c("se_probs[", seq_len(nrow(test_txt)), "] <- ", test_txt$SeText),
    str_c("sp_probs[", seq_len(nrow(test_txt)), "] <- ", test_txt$SpText),
    if(cov && n_tests>2L) "if(!all(se_probs >= 0 & se_probs <= 1 & sp_probs >= 0 & sp_probs <= 1)) return(FALSE)",
    "return(list(se_probs,sp_probs))",
    #str_c("prev <- rep(prev, each=", nrow(test_txt), "L)"),
    #str_c("probs <- array(c(prev*se_probs, (1-prev)*sp_probs), dim=c(", str_c(rep("2",n_tests),collapse=","), ",",n_pops, ",2))"),
    # str_c("attr(probs, 'dimnames') <- ", str_c(capture.output(dput(dn)), collapse="")),
    "return(probs)",
    "}"
  )
  #str_c(fun_txt, collapse="\n") |> cat()
  str_c(fun_txt, collapse="\n") |> parse(text=_) |> eval() |> cmpfun()
}



TestSim <- R6Class(
  "TestSim",
  public = list(

    initialize = function(populations, tests, correlations=tribble(~"TestA", ~"TestB", ~"Parameter", ~"Correlation")){

      stopifnot(
        is.data.frame(populations), c("Population","Prevalence","SampleProportion") %in% names(populations),
        populations[["SampleProportion"]] >= 0, populations[["SampleProportion"]] <= 1
      )
      stopifnot(
        is.data.frame(tests),
        c("Test","Parameter","Alpha","Beta") %in% names(tests),
        tests[["Parameter"]] %in% c("Se","Sp")
      )
      stopifnot(
        is.data.frame(tests),
        c("Test","Parameter","Alpha","Beta") %in% names(tests),
        tests[["Parameter"]] %in% c("Se","Sp")
      )

      ## Normalise:
      populations[["SampleProportion"]] <- populations[["SampleProportion"]] / sum(populations[["SampleProportion"]])

      if(!"Estimate" %in% names(tests)) tests[["Estimate"]] <- NA_real_
      tests |>
        mutate(Estimate = case_when(
          is.na(Estimate) ~ (Alpha-1)/(Alpha+Beta-2),
          .default = Estimate
        )) ->
        tests
      stopifnot(tests[["Estimate"]] >= 0, tests[["Estimate"]] <= 1)

      if(!"Type" %in% names(tests)){
        tests |>
          mutate(Type = case_when(
            Test %in% c("C1","C2","D1","D2") ~ "Comparator",
            .default = "Reference"
          )) ->
          tests
      }
      comptest <- tests |> filter(Type=="Comparator") |> pull(Test)
      stopifnot(sum(tests[["Type"]]=="Comparator") == 2L, length(comptest)==2L, comptest[1L]==comptest[2L])
      stopifnot(!is.na(tests$Test), tests |> count(Test) |> pull(n) == 2L)

      ## Re-arrange so the comparator test is always first:
      tests |>
        mutate(Test = as.character(Test)) |>
        arrange(desc(Type=="Comparator"), Test, Parameter) ->
        tests
      testnames <- tests |> distinct(Test) |> pull()

      stopifnot(
        tests |> count(Test, Parameter) |> pull(n) == 1L,
        tests |> count(Test) |> pull(n) == 2L
      )

      n_tests <- (nrow(tests)/2L)
      stopifnot(n_tests %in% 2L:3L)
      n_pops <- nrow(populations)

      se <- tests |> filter(Parameter=="Se") |> pull(Estimate)
      sp <- tests |> filter(Parameter=="Sp") |> pull(Estimate)
      names(se)=names(sp) <- testnames

      ## Get covariances from correlations
      stopifnot(correlations$TestA == testnames[1L], correlations$TestB %in% testnames[-1L], correlations$Parameter %in% c("Se","Sp"))

      getcov <- function(tt, parameter){
        if(tt==testnames[1L]) return(0.0)
        cor <- correlations |> filter(Parameter==parameter, TestB==tt) |> pull(Correlation)
        if(length(cor)==0L) return(0.0)
        cor_to_cov(if(parameter=="Se") se[c(testnames[1L], tt)] else sp[c(testnames[1L], tt)], cor) |> as.numeric()
      }
      covse <- sapply(testnames, getcov, parameter="Se")
      covsp <- sapply(testnames, getcov, parameter="Sp")

      private$parameters <- list(
        n_tests=n_tests,
        n_combos=2^n_tests,
        se=se,
        covse=covse,
        sp=sp,
        covsp=covsp,
        n_pops=n_pops,
        prev=populations$Prevalence,
        tests=tests,
        populations=populations
      )

      ## Generate the probability functions (separate for 2 test, 3 test, 3 test with cor):
      # TODO: use C++, and maybe hard-code ref test performance for efficiency???
      private$prob_funs <- list(
        partial = create_hw_function(n_tests-1L, cov=FALSE),
        nocov = create_hw_function(n_tests, cov=FALSE),
        full = create_hw_function(n_tests, cov=TRUE)
      )

      ## Then do a data sim
      private$data <- list(
        N = 0L
      )

      invisible(self)
    },

    simulate_data = function(total_n, sample_pars=TRUE){

      N <- round(private$parameters$populations$SampleProportion * total_n)
      if(sum(N)!=total_n) N[length(N)] <- total_n - sum(N[-length(N)])
      stopifnot(N>0, sum(N)==total_n)

      if(sample_pars){
        stop("Implement me")
        ## Resample se and sp of ref tests, as well as covse and covsp
      }else{
        probs <- private$prob_funs$full(private$parameters$se, private$parameters$sp, private$parameters$covse, private$parameters$covsp)
      }

      stopifnot(!isFALSE(probs))

      prev <- private$parameters$prev
      seq_along(prev) |>
        lapply(function(p){
          rmultinom(1, N[p], prev[p]*probs[[1]] + (1-prev[p])*probs[[2]])
        }) ->
        tally_full

      indexes <- seq(2L,private$parameters$n_combos,by=2L)
      tally_full |>
        lapply(function(x) x[indexes-1L] + x[indexes]) ->
        tally_partial
      tally_full |>
        vapply(function(x) x[indexes-1L], numeric(private$parameters$n_combos/2L)) ->
        tally_comp_neg
      tally_full |>
        vapply(function(x) x[indexes], numeric(private$parameters$n_combos/2L)) ->
        tally_comp_pos

      names(tally_full)=names(tally_partial) <- str_c("Pop",seq_along(prev))
      dimnames(tally_comp_neg)[[2]]=dimnames(tally_comp_pos)[[2]] <- str_c("Pop",seq_along(prev))

      private$data$N <- N
      private$data$tally_full <- tally_full
      private$data$tally_partial <- tally_partial
      private$data$tally_comp <- list(neg=tally_comp_neg, pos=tally_comp_pos)
      private$est_partial <- list()

      invisible(tally_full)
    },

    estimate_partial = function(){
      if(length(private$est_partial)>0L) return(private$est_partial)
      N <- private$data$N
      probs <- private$prob_funs$partial(private$parameters$se[-1L], private$parameters$sp[-1L])
      sapply(seq_along(N), function(pp){
        optimise(
          function(x) -dmultinom(private$data$tally_partial[[pp]], N[pp], x*probs[[1]] + (1-x)*probs[[2]], log=TRUE),
          interval=c(0,1)
        )$minimum
      }) -> prevs
      private$est_partial <- list(prevs=prevs, probs=probs)
      private$est_partial
    },

    method_A = function(){

      tally_comp <- private$data$tally_comp |> lapply(function(x) apply(x,1,sum))

      ## If we have >2 tests then the gold standard is any positive:
      if(private$parameters$n_tests>2L){
        tally_comp |>
          lapply(function(x) c(x[1], sum(x)-x[1])) ->
          tally_comp
      }

      estse <- tally_comp$pos[2] / (tally_comp$pos[2] + tally_comp$neg[2])
      estsp <- tally_comp$neg[1] / (tally_comp$pos[1] + tally_comp$neg[1])

      tibble(Method = "A", Se = estse, Sp = estsp)
    },

    method_B = function(opt_method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B")){

      opt_method <- match.arg(opt_method)
      lower <- -Inf
      upper <- Inf
      if(opt_method=="L-BFGS-B"){
        lower <- 1e-6
        upper <- 1.0 - 1e-6
      }

      ## First get PPP based on ref tests:
      estpart <- self$estimate_partial()
      prevs <- estpart$prevs
      probs <- estpart$probs
      vapply(seq_along(prevs), function(pp){
        se_ppp <- probs[[1]]*prevs[pp]
        sp_ppp <- probs[[2]]*(1-prevs[pp])
        ppp <- se_ppp / (se_ppp + sp_ppp)
        rep(ppp, each=2L)
        }, numeric(private$parameters$n_combos)
      ) -> pppmat

      ## Then optimise for se and sp of comparator using the trick I sent to Hayley:
      pars <- c(
        se=0.75,
        sp=0.75
      )
      optim(pars, function(x){
        if(any(x<0) || any(x>1) || sum(x)<0) return(Inf)
        seq_len(private$parameters$n_pops) |>
          sapply(function(pp){
            probs <- c(1-x[1L],x[1L])*pppmat[,pp] + c(x[2L],1-x[2L])*(1-pppmat[,pp])
            -dmultinom(private$data$tally_full[[pp]], private$data$N[pp], probs, log=TRUE)
          }) |> sum()
      })$par |>
        as.numeric() ->
        pars

      return(tibble(Method = "B", Se = pars[1L], Sp=pars[2L]))


      ## Equivalent:
      estpart <- self$estimate_partial()
      prevs <- estpart$prevs

      pars <- c(
        se=private$parameters$se[1L],
        se=private$parameters$sp[1L]
      )
      zs <- rep(0, private$parameters$n_tests)
      optim(pars, function(x){
        se <- private$parameters$se
        se[1] <- x[1L]
        sp <- private$parameters$sp
        sp[1] <- x[2L]
        self$ll_hw(prevs, se, sp, zs, zs)
      }, method=opt_method, lower=lower, upper=upper)$par ->
        pars

      return(tibble(Method = "B", Se = pars[1L], Sp=pars[2L]))

    },

    method_C_TMB = function(method=c("HW","HWpriors"), opt_method="nlm"){

      st <- Sys.time()
      stopifnot(private$parameters$n_test==2L)

      method <- match.arg(method)

      tally <- private$data$tally_full |> do.call("bind_cols", args=_) |> as.matrix()
      tally[] <- as.double(tally)
      Params <- list(selg=rep(1,2), splg=rep(1,2), prevlg=rep(0,private$parameters$n_test))
      Data <- list(P=private$parameters$n_test, tally=tally, se_alpha=rep(2,2), se_beta=rep(1,2), sp_alpha=rep(2,2), sp_beta=rep(1,2))

      if(method=="HWpriors"){
        Params$selg[2] <- qlogis(private$parameters$se[2])
        Params$splg[2] <- qlogis(private$parameters$sp[2])
        sepr <- abs(1000*(c(0,1)-private$parameters$se[2]))
        sppr <- abs(1000*(c(0,1)-private$parameters$sp[2]))
        Data$se_alpha <- c(1,sepr[1])
        Data$se_beta[2] <- sepr[2]
        Data$sp_alpha <- c(1,sppr[1])
        Data$sp_beta[2] <- sppr[2]
      }
      obj <- MakeADFun(data=Data, parameters=Params, DLL="hw_2t")
      pars <- plogis(optimr(unlist(Params), obj$fn, obj$gr, obj$he, method=opt_method)$par)

      ests <- pars[c(1L,3L)]
      dt <- as.numeric(Sys.time() - st, units="secs")

      tibble(Method = method, Se = ests[1L], Sp = ests[2L], Time = dt, Pars = list(pars))

    },

    method_C_cheat = function(correlation=c("model","post-hoc","none"), opt_method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B")){

      correlation <- match.arg(correlation)

      opt_method <- match.arg(opt_method)
      lower <- -Inf
      upper <- Inf
      if(opt_method=="L-BFGS-B"){
        lower <- 1e-6
        upper <- 1.0 - 1e-6
      }

      estpart <- self$estimate_partial()
      prevs <- estpart$prevs

      ## Need to replace with TMB!!
      stopifnot(private$parameters$n_tests==2L)

      if(correlation=="model"){

        method <- "C1"

        ## Fit a HW model with correlation, but ref test performance is fixed:
        pars <- c(
          se=0.75,
          sp=0.75,
          covsehp1=rep(0.5,private$parameters$n_test-1L),
          covsphp1=rep(0.5,private$parameters$n_test-1L)
        )
        indexes <- list(
          se=1L,
          sp=2L,
          covsehp1=2L+seq_len(private$parameters$n_tests-1L),
          covsphp1=2L+private$parameters$n_tests-1L+seq_len(private$parameters$n_tests-1L)
        )
        stopifnot(all(unlist(indexes) %in% seq_along(unlist(pars))), all(table(unlist(indexes))==1L))

        optim(pars, fn=function(x){
          se <- private$parameters$se
          se[1L] <- x[indexes$se]
          sp <- private$parameters$sp
          sp[1L] <- x[indexes$sp]
          self$ll_hw(prevs, se, sp, c(0, 2*(x[indexes$covsehp1]-0.5)), c(0, 2*(x[indexes$covsphp1]-0.5)))
        }, method=opt_method, lower=lower, upper=upper)$par ->
          pars

        pars[c(indexes$covsehp1, indexes$covsphp1)] <- 2*(pars[c(indexes$covsehp1, indexes$covsphp1)]-0.5)
        ests <- pars[c(indexes$se[1L], indexes$sp[1L])]

      }else if(correlation=="post-hoc"){
        ## Fit a HW model without correlation:

        method <- "C2"

        stopifnot(private$parameters$n_tests==2L)
        ## Maybe relax this by adding a negative correlation between ref tests - would have to swap test order (for 3 tests max)

        pars <- c(
          se=rep(0.75, private$parameters$n_tests),
          se=rep(0.75, private$parameters$n_tests)
        )
        indexes <- list(
          se=seq_len(private$parameters$n_tests),
          sp=private$parameters$n_tests+seq_len(private$parameters$n_tests)
        )
        #stopifnot(all(unlist(indexes) %in% seq_along(unlist(pars))), all(table(unlist(indexes))==1L))
        zs <- rep(0, private$parameters$n_tests)
        optim(pars, function(x){
          self$ll_hw(prevs, x[indexes$se], x[indexes$sp], zs, zs)
        }, method=opt_method, lower=lower, upper=upper)$par ->
          pars

        truese <- private$parameters$se
        truesp <- private$parameters$sp
        estse <- pars[indexes$se]
        estsp <- pars[indexes$sp]

        cors <- find_cors(truese[2], truesp[2], max(truese[2],estse[2]), max(truesp[2],estsp[2]))
        ests <- find_ss(estse[1], estsp[1], cors[1], cors[2])
        names(ests) <- c("se","sp")

      }else if(correlation=="none"){

        ## Fit a HW model without correlation and with fixed parameters:

        method <- "C3"

        pars <- c(
          se=0.75,
          se=0.75
        )
        indexes <- list(
          se=1L,
          sp=2L
        )
        #stopifnot(all(unlist(indexes) %in% seq_along(unlist(pars))), all(table(unlist(indexes))==1L))
        zs <- rep(0, private$parameters$n_tests)
        optim(pars, function(x){
          se <- private$parameters$se
          se[1] <- x[indexes$se]
          sp <- private$parameters$sp
          sp[1] <- x[indexes$sp]
          self$ll_hw(prevs, se, sp, zs, zs)
        }, method=opt_method, lower=lower, upper=upper)$par ->
          pars

        ests <- pars[c(indexes$se, indexes$sp)]
        names(ests) <- c("se","sp")

      }else{
        stop("Unrecognised correlation type")
      }

      tibble(Method = method, Se = ests[1L], Sp = ests[2L], Pars = list(pars))

    },

    method_C = function(correlation=c("model","post-hoc","none"), opt_method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B")){

      correlation <- match.arg(correlation)

      opt_method <- match.arg(opt_method)
      lower <- -Inf
      upper <- Inf
      if(opt_method=="L-BFGS-B"){
        lower <- 1e-6
        upper <- 1.0 - 1e-6
      }

      estpart <- self$estimate_partial()
      prevs <- estpart$prevs

      ## Need to replace with TMB!!
      #stopifnot(private$parameters$n_tests==2L)

      if(correlation=="model"){

        method <- "C1"

        ## Fit a HW model with correlation, but ref test performance is fixed:
        pars <- c(
          prev=prevs,
          se=0.75,
          sp=0.75,
          covsehp1=rep(0.5,private$parameters$n_test-1L),
          covsphp1=rep(0.5,private$parameters$n_test-1L)
        )
        indexes <- list(
          prevs=seq_len(private$parameters$n_pops),
          se=private$parameters$n_pops+1L,
          sp=private$parameters$n_pops+2L,
          covsehp1=private$parameters$n_pops+2L+seq_len(private$parameters$n_tests-1L),
          covsphp1=private$parameters$n_pops+2L+private$parameters$n_tests-1L+seq_len(private$parameters$n_tests-1L)
        )
        stopifnot(all(unlist(indexes) %in% seq_along(unlist(pars))), all(table(unlist(indexes))==1L))

        optim(pars, fn=function(x){
          se <- private$parameters$se
          se[1L] <- x[indexes$se]
          sp <- private$parameters$sp
          sp[1L] <- x[indexes$sp]
          self$ll_hw(x[indexes$prevs], se, sp, c(0, 2*(x[indexes$covsehp1]-0.5)), c(0, 2*(x[indexes$covsphp1]-0.5)))
        }, method=opt_method, lower=lower, upper=upper)$par ->
          pars

        pars[c(indexes$covsehp1, indexes$covsphp1)] <- 2*(pars[c(indexes$covsehp1, indexes$covsphp1)]-0.5)
        ests <- pars[c(indexes$se[1L], indexes$sp[1L])]

      }else if(correlation=="post-hoc"){
        ## Fit a HW model without correlation:

        method <- "C2"

        stopifnot(private$parameters$n_tests==2L)
        ## Maybe relax this by adding a negative correlation between ref tests - would have to swap test order (for 3 tests max)

        pars <- c(
          prev=prevs,
          se=rep(0.75, private$parameters$n_tests),
          se=rep(0.75, private$parameters$n_tests)
        )
        indexes <- list(
          prevs=seq_len(private$parameters$n_pops),
          se=private$parameters$n_pops+seq_len(private$parameters$n_tests),
          sp=private$parameters$n_pops+private$parameters$n_tests+seq_len(private$parameters$n_tests)
        )
        #stopifnot(all(unlist(indexes) %in% seq_along(unlist(pars))), all(table(unlist(indexes))==1L))
        zs <- rep(0, private$parameters$n_tests)
        optim(pars, function(x){
          self$ll_hw(x[indexes$prevs], x[indexes$se], x[indexes$sp], zs, zs)
        }, method=opt_method, lower=lower, upper=upper)$par ->
          pars

        truese <- private$parameters$se
        truesp <- private$parameters$sp
        estse <- pars[indexes$se]
        estsp <- pars[indexes$sp]

        cors <- find_cors(truese[2], truesp[2], max(truese[2],estse[2]), max(truesp[2],estsp[2]))
        ests <- find_ss(estse[1], estsp[1], cors[1], cors[2])
        names(ests) <- c("se","sp")

      }else if(correlation=="none"){

        ## Fit a HW model without correlation and with fixed parameters:

        method <- "C3"

        pars <- c(
          prev=prevs,
          se=0.75,
          se=0.75
        )
        indexes <- list(
          prevs=seq_len(private$parameters$n_pops),
          se=private$parameters$n_pops+1L,
          sp=private$parameters$n_pops+2L
        )
        #stopifnot(all(unlist(indexes) %in% seq_along(unlist(pars))), all(table(unlist(indexes))==1L))
        zs <- rep(0, private$parameters$n_tests)
        optim(pars, function(x){
          se <- private$parameters$se
          se[1] <- x[indexes$se]
          sp <- private$parameters$sp
          sp[1] <- x[indexes$sp]
          self$ll_hw(x[indexes$prevs], se, sp, zs, zs)
        }, method=opt_method, lower=lower, upper=upper)$par ->
          pars

        ests <- pars[c(indexes$se, indexes$sp)]
        names(ests) <- c("se","sp")

      }else{
        stop("Unrecognised correlation type")
      }

      tibble(Method = method, Se = ests[1L], Sp = ests[2L], Pars = list(pars))

    },

    method_D = function(opt_method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B")){

      opt_method <- match.arg(opt_method)
      lower <- -Inf
      upper <- Inf
      if(opt_method=="L-BFGS-B"){
        lower <- 1e-6
        upper <- 1.0 - 1e-6
      }

      if(private$parameters$n_tests==2L){

        ## Fit a HW model without correlation, as correlation is unidentifiable anyway:
        pars <- c(
          prev=rep(0.5, private$parameters$n_pop),
          se=rep(0.75, private$parameters$n_test),
          sp=rep(0.75, private$parameters$n_test)
        )
        indexes <- list(
          prev = seq_len(private$parameters$n_pops),
          se=private$parameters$n_pops+seq_len(private$parameters$n_tests),
          sp=private$parameters$n_pops+private$parameters$n_tests+seq_len(private$parameters$n_tests)
        )
        #stopifnot(all(unlist(indexes) %in% seq_along(unlist(pars))), all(table(unlist(indexes))==1L))
        zs <- rep(0, private$parameters$n_tests)
        optim(pars, function(x){
          self$ll_hw(x[indexes$prev], x[indexes$se], x[indexes$sp], zs, zs)
        }, method=opt_method, lower=lower, upper=upper)$par ->
          pars

        pars[c(indexes$covsehp1, indexes$covsphp1)] <- 2*(pars[c(indexes$covsehp1, indexes$covsphp1)]-0.5)

      }else{

        ## Fit a HW model with correlation:
        pars <- c(
          prev=rep(0.5, private$parameters$n_pop),
          se=rep(0.75, private$parameters$n_test),
          sp=rep(0.75, private$parameters$n_test),
          covsehp1=rep(0.5,private$parameters$n_test-1L),
          covsphp1=rep(0.5,private$parameters$n_test-1L)
        )
        indexes <- list(
          prev = seq_len(private$parameters$n_pops),
          se=private$parameters$n_pops+seq_len(private$parameters$n_tests),
          sp=private$parameters$n_pops+private$parameters$n_tests+seq_len(private$parameters$n_tests),
          covsehp1=private$parameters$n_pops+private$parameters$n_tests*2L+seq_len(private$parameters$n_tests-1L),
          covsphp1=private$parameters$n_pops+private$parameters$n_tests*2L+private$parameters$n_tests-1L+seq_len(private$parameters$n_tests-1L)
        )
        #stopifnot(all(unlist(indexes) %in% seq_along(unlist(pars))), all(table(unlist(indexes))==1L))
        zs <- rep(0, private$parameters$n_tests)
        optim(pars, function(x){
          self$ll_hw(x[indexes$prev], x[indexes$se], x[indexes$sp], c(0, 2*(x[indexes$covsehp1]-0.5)), c(0, 2*(x[indexes$covsphp1]-0.5)))
          #self$ll_hw(x[indexes$prev], x[indexes$se], x[indexes$sp], zs, zs)
        }, method=opt_method, lower=lower, upper=upper)$par ->
          pars

        pars[c(indexes$covsehp1, indexes$covsphp1)] <- 2*(pars[c(indexes$covsehp1, indexes$covsphp1)]-0.5)
      }

      tibble(Method = "D", Se = pars[indexes$se[1]], Sp = pars[indexes$sp[1]], Pars = list(pars))

    },


    ll_hw = function(prev, se, sp, covse, covsp){
      if(any(prev<0) || any(prev>1) || any(se<0) || any(se>1) || any(sp<0) || any(sp>1) ||
          any(covse<(-1)) || any(covse>1) || any(covsp<(-1)) || any(covsp>1) || any((se+sp)<1)) return(Inf)
      probs <- private$prob_funs$full(se, sp, covse, covsp)
      if(any(probs[[1]]<0) || any(probs[[1]]>1) || abs(sum(probs[[1]])-1)>sqrt(.Machine$double.eps)) return(Inf)
      if(any(probs[[2]]<0) || any(probs[[2]]>1) || abs(sum(probs[[2]])-1)>sqrt(.Machine$double.eps)) return(Inf)
      sapply(seq_along(prev), function(pp){
        -dmultinom(private$data$tally_full[[pp]], private$data$N[pp], prev[pp]*probs[[1]] + (1-prev[pp])*probs[[2]], log=TRUE)
      }) |> sum() -> ll
      ll
      #lp <- sum(dbeta(se, 2, 1, log=TRUE) + dbeta(sp, 2, 1, log=TRUE))
      #ll+lp
    },



    method_broken_ppp = function(){

      estpart <- self$estimate_partial()
      prevs <- estpart$prevs
      probs <- estpart$probs

      ## Then calculate PPPs:
      vapply(seq_along(prevs), function(pp){
        se_ppp <- probs[[1]]*prevs[pp] / (probs[[1]]*prevs[pp] + (1-probs[[1]])*prevs[pp])
        sp_ppp <- probs[[2]]*(1-prevs[pp]) / (probs[[2]]*(1-prevs[pp]) + (1-probs[[2]])*(1-prevs[pp]))
        ppp <- se_ppp / (se_ppp + sp_ppp)
      }, numeric(private$parameters$n_combos/2L)
      ) -> pppmat

      ppv_p <- (se*prevs) / (se*prevs + (1-sp)*(1-prevs))
      ppv_n <- 1 - ((sp*(1-prevs)) / (sp*(1-prevs) + (1-se)*prevs))

      ppv_p * private$data$tally_comp$pos[2,]
      ppv_n * private$data$tally_comp$neg[1,]

      se <- sum(ppv_p * private$data$tally_comp$pos[2,]) / sum(ppv_p * private$data$tally_comp$pos[2,] + ppv_n * private$data$tally_comp$neg[1,])
      sp <- 1 - (sum((1-ppv_p) * private$data$tally_comp$pos[1,]) / sum((1-ppv_p) * private$data$tally_comp$pos[1,] + (1-ppv_n) * private$data$tally_comp$neg[2,]))

      return(tibble(Method = "B", Se = se, Sp=sp))


      (ppv_p * t(private$data$tally_comp$pos)) / ((ppv_p * t(private$data$tally_comp$pos)) + (ppv_n * t(private$data$tally_comp$neg)))
      ppv_n * t(private$data$tally_comp$neg)

      ppv_p
      ppv_n

      ppp_pos <- sum(private$data$tally_comp$pos * pppmat)
      ppp_neg <- sum(private$data$tally_comp$neg * pppmat)
      ppp_pos_1m <- sum(private$data$tally_comp$pos * (1-pppmat))
      ppp_neg_1m <- sum(private$data$tally_comp$neg * (1-pppmat))

      se <- (ppp_pos+1) / (ppp_pos+ppp_neg+1)
      sp <- (ppp_neg_1m+1) / (ppp_neg_1m+ppp_pos_1m+1)




      tibble(Method = "B", Se = se, Sp=sp)

    },

    method_C_old = function(correlation=c("model","post-hoc","none","sep_prev"), opt_method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B")){

      correlation <- match.arg(correlation)

      opt_method <- match.arg(opt_method)
      lower <- -Inf
      upper <- Inf
      if(opt_method=="L-BFGS-B"){
        lower <- 1e-6
        upper <- 1.0 - 1e-6
      }

      if(correlation=="model"){

        method <- "C1"

        ## Fit a HW model with correlation, but ref test performance is fixed:
        pars <- c(
          prev=rep(0.5,private$parameters$n_pops),
          se=0.75,
          sp=0.75,
          covsehp1=rep(0.5, private$parameters$n_tests-1L),
          covsphp1=rep(0.5, private$parameters$n_tests-1L)
        )
        indexes <- list(
          prev = seq_len(private$parameters$n_pops),
          se=private$parameters$n_pops+1L,
          sp=private$parameters$n_pops+2L,
          covsehp1=private$parameters$n_pops+2L+seq_len(private$parameters$n_tests-1L),
          covsphp1=private$parameters$n_pops+2L+private$parameters$n_tests-1L+seq_len(private$parameters$n_tests-1L)
        )
        #stopifnot(all(unlist(indexes) %in% seq_along(unlist(pars))), all(table(unlist(indexes))==1L))
        optim(pars, function(x){
          se <- private$parameters$se
          se[1L] <- x[indexes$se]
          sp <- private$parameters$sp
          sp[1L] <- x[indexes$sp]
          self$ll_hw(x[indexes$prev], se, sp, c(0, 2*(x[indexes$covsehp1]-0.5)), c(0, 2*(x[indexes$covsphp1]-0.5)))
        }, method=opt_method, lower=lower, upper=upper)$par ->
          pars

        pars[c(indexes$covsehp1, indexes$covsphp1)] <- 2*(pars[c(indexes$covsehp1, indexes$covsphp1)]-0.5)
        ests <- pars[c("se","sp")]


      }else if(correlation=="post-hoc"){
        ## Fit a HW model without correlation:

        method <- "C2"

        stopifnot(private$parameters$n_tests==2L)
        ## Maybe relax this by adding a negative correlation between ref tests - would have to swap test order (for 3 tests max)

        ## First estimate prevalence:
        N <- private$data$N
        probs <- private$prob_funs$partial(private$parameters$se[-1L], private$parameters$sp[-1L])
        if(any(probs[[1]]<0) || any(probs[[1]]>1) || abs(sum(probs[[1]])-1)>sqrt(.Machine$double.eps)) return(Inf)
        if(any(probs[[2]]<0) || any(probs[[2]]>1) || abs(sum(probs[[2]])-1)>sqrt(.Machine$double.eps)) return(Inf)
        sapply(seq_along(N), function(pp){
          optimise(
            function(x) -dmultinom(private$data$tally_partial[[pp]], N[pp], x*probs[[1]] + (1-x)*probs[[2]], log=TRUE),
            interval=c(0,1)
          )$minimum
        }) ->
          prevs

        pars <- c(
          se=rep(0.75, private$parameters$n_tests),
          se=rep(0.75, private$parameters$n_tests)
        )
        indexes <- list(
          se=seq_len(private$parameters$n_tests),
          sp=private$parameters$n_tests+seq_len(private$parameters$n_tests)
        )
        #stopifnot(all(unlist(indexes) %in% seq_along(unlist(pars))), all(table(unlist(indexes))==1L))
        zs <- rep(0, private$parameters$n_tests)
        optim(pars, function(x){
          self$ll_hw(prevs, x[indexes$se], x[indexes$sp], zs, zs)
        }, method=opt_method, lower=lower, upper=upper)$par ->
          pars

        truese <- private$parameters$se
        truesp <- private$parameters$sp
        estse <- pars[indexes$se]
        estsp <- pars[indexes$sp]

        cors <- find_cors(truese[2], truesp[2], max(truese[2],estse[2]), max(truesp[2],estsp[2]))
        ests <- find_ss(estse[1], estsp[1], cors[1], cors[2])
        names(ests) <- c("se","sp")

      }else if(correlation=="none"){

        ## Fit a HW model without correlation and with fixed parameters:

        method <- "C3"

        pars <- c(
          prev=rep(0.5,private$parameters$n_pops),
          se=0.75,
          se=0.75
        )
        indexes <- list(
          prev = seq_len(private$parameters$n_pops),
          se=private$parameters$n_pops+1L,
          sp=private$parameters$n_pops+2L
        )
        #stopifnot(all(unlist(indexes) %in% seq_along(unlist(pars))), all(table(unlist(indexes))==1L))
        zs <- rep(0, private$parameters$n_tests)
        optim(pars, function(x){
          se <- private$parameters$se
          se[1] <- x[indexes$se]
          sp <- private$parameters$sp
          sp[1] <- x[indexes$sp]
          self$ll_hw(x[indexes$prev], se, sp, zs, zs)
        }, method=opt_method, lower=lower, upper=upper)$par ->
          pars

        ests <- pars[c(indexes$se, indexes$sp)]
        names(ests) <- c("se","sp")

      }else if(correlation=="sep_prev"){

        ## First estimate prevalence:
        N <- private$data$N
        probs <- private$prob_funs$partial(private$parameters$se[-1L], private$parameters$sp[-1L])
        if(any(probs[[1]]<0) || any(probs[[1]]>1) || abs(sum(probs[[1]])-1)>sqrt(.Machine$double.eps)) return(Inf)
        if(any(probs[[2]]<0) || any(probs[[2]]>1) || abs(sum(probs[[2]])-1)>sqrt(.Machine$double.eps)) return(Inf)
        sapply(seq_along(N), function(pp){
          optimise(
            function(x) -dmultinom(private$data$tally_partial[[pp]], N[pp], x*probs[[1]] + (1-x)*probs[[2]], log=TRUE),
            interval=c(0,1)
          )$minimum
        }) ->
          prevs

        ## Then fit a HW model without correlation and with fixed parameters:

        method <- "C4"

        pars <- c(
          se=0.75,
          se=0.75
        )
        indexes <- list(
          se=1L,
          sp=2L
        )
        #stopifnot(all(unlist(indexes) %in% seq_along(unlist(pars))), all(table(unlist(indexes))==1L))
        zs <- rep(0, private$parameters$n_tests)
        optim(pars, function(x){
          se <- private$parameters$se
          se[1] <- x[indexes$se]
          sp <- private$parameters$sp
          sp[1] <- x[indexes$sp]
          self$ll_hw(prevs, se, sp, zs, zs)
        }, method=opt_method, lower=lower, upper=upper)$par ->
          pars

        ests <- pars[c(indexes$se, indexes$sp)]
        names(ests) <- c("se","sp")
      }else{
        stop("Unrecognised correlation type")
      }

      tibble(Method = method, Se = ests["se"], Sp = ests["sp"], Pars = list(pars))

    },


    set_data = function(tally){
      stopifnot(tally>=0, dim(tally)==c(rep(2L,length=length(private$dims$ref)+1L), length(private$sample_size)), apply(tally,length(dim(tally)),sum)==private$sample_size)
      stopifnot(tally>=0)
      private$tally <- tally
    },

    estimate_cor = function(with_cor=TRUE, method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B")){

      stopifnot(length(private$tally)==length(private$sample_size) * 2^(length(private$dims$ref)+1))

      if(length(private$dims$ref)==1){
        rv <- self$estimate_cor_2t(with_cor, method)
      }else if(length(private$dims$ref)==2){
        rv <- self$estimate_cor_3t(with_cor, method)
      }else{
        stop("Number of tests (>3) not supported")
      }
      rv
    },

    estimate_cor_2t = function(with_cor=TRUE, method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B")){
      stopifnot(length(private$dims$ref)==1, length(private$sample_size)==2)  # hard-coded indexes below!

      method <- match.arg(method)
      lower <- -Inf
      upper <- Inf
      if(method=="L-BFGS-B"){
        lower <- 1e-6
        upper <- 1.0 - 1e-6
      }

      nref <- length(private$dims$ref)
      pars <- c(se=0.75, sp=0.75, rep(0.5, nref*2), rep(0.5,length(private$sample_size)))
      names(pars) <- c("se","sp", str_c("corse",seq_len(nref)+1),str_c("corsp",seq_len(nref)+1), str_c("prev",seq_along(private$sample_size)))
      if(!with_cor) pars <- pars[-(3:4)]
      est_comp <- optim(pars, function(x){
        if(any(x<=0) || any(x>=1) || sum(x[1:2])<=1) return(Inf)
        if(with_cor){
          probs <- private$make_probs$sim(x[5:6], x[1:2], x[3]-0.5, x[4]-0.5)
        }else{
          probs <- private$make_probs$sim(x[3:4], x[1:2], 0, 0)
        }
        if(!attr(probs, "is_valid")) return(Inf)
        seq_along(private$sample_size) |>
          as.list() |>
          sapply(function(p){
            -dmultinom(private$tally[,,p], private$sample_size[p], probs[,p], log=TRUE)
          }) |>
          sum()
      }, method=method, lower=lower, upper=upper)
      if(with_cor) est_comp$par[3:4] <- est_comp$par[3:4]-0.5
      est_comp$par
    },

    estimate_cor_3t = function(with_cor=TRUE, method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B")){
      stopifnot(length(private$dims$ref)==2, length(private$sample_size)==1)  # hard-coded indexes below!

      method <- match.arg(method)
      lower <- -Inf
      upper <- Inf
      if(method=="L-BFGS-B"){
        lower <- 1e-6
        upper <- 1.0 - 1e-6
      }

      nref <- length(private$dims$ref)
      pars <- c(se=0.75, sp=0.75, rep(0.5, nref*2), rep(0.5,length(private$sample_size)))
      names(pars) <- c("se","sp", str_c("corse",seq_len(nref)+1),str_c("corsp",seq_len(nref)+1), str_c("prev",seq_along(private$sample_size)))
      if(!with_cor) pars <- pars[-(3:6)]
      est_comp <- optim(pars, function(x){
        if(any(x<=0) || any(x>=1) || sum(x[1:2])<=1) return(Inf)
        if(with_cor){
          probs <- private$make_probs$sim(x[7], x[1:2], c(0,x[3:4]-0.5), c(0,x[5:6]-0.5))
        }else{
          probs <- private$make_probs$sim(x[3], x[1:2], 0, 0)
        }
        if(!attr(probs, "is_valid")) return(Inf)
        seq_along(private$sample_size) |>
          as.list() |>
          sapply(function(p){
            -dmultinom(private$tally[], private$sample_size[p], probs[], log=TRUE)
          }) |>
          sum()
      }, method=method, lower=lower, upper=upper)
      if(with_cor) est_comp$par[3:6] <- est_comp$par[3:6]-0.5
      est_comp$par
    },

    estimate_prev = function(){
      stopifnot(length(private$tally)==length(private$sample_size) * 2^(length(private$ref_dims)+1))
      reftally <- apply(private$tally,c(1L,private$ref_dims-1L),sum)
      stopifnot(length(reftally)==length(private$sample_size) * 2^length(private$ref_dims))

      if(length(private$sample_size)==1L){
        private$est_prev <- optimise(function(x){
          if(x<=0 || x>=1) return(Inf)
          probs <- private$make_probs(x,c(0.5,0.5)) |> apply(c(1,private$ref_dims), sum)
          seq_along(private$sample_size) |>
            as.list() |>
            sapply(function(p){
              -dmultinom(reftally[p,], private$sample_size[p], probs[p,], log=TRUE)
            }) |>
            sum()
        }, c(0,1))
        private$est_prev$par <- private$est_prev$minimum
      }else if(length(private$sample_size)==2L){
        private$est_prev <- optim(c(prev1=0.5, prev2=0.5), function(x){
          if(any(x<=0) || any(x>=1)) return(Inf)
          probs <- private$make_probs(x,c(0.5,0.5)) |> apply(c(1,private$ref_dims), sum)
          seq_along(private$sample_size) |>
            as.list() |>
            sapply(function(p){
              -dmultinom(reftally[p,], private$sample_size[p], probs[p,], log=TRUE)
            }) |>
            sum()
        })
      }

      invisible(private$est_prev)
    },

    estimate_comp = function(include_prev=FALSE, method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B")){

      stopifnot(!is.null(private$est_prev$par))
      method <- match.arg(method)
      lower <- -Inf
      upper <- Inf
      if(method=="L-BFGS-B"){
        lower <- 1e-6
        upper <- 1.0 - 1e-6
      }

      if(include_prev){

        pars <- c(se=0.75, sp=0.75, rep(0.5,length(private$sample_size)))
        names(pars) <- c("se","sp",str_c("prev",1:length(private$sample_size)))
        est_comp <- optim(pars, function(x){
          if(any(x<=0) || any(x>=1) || sum(x[1:2])<=1) return(Inf)
          probs <- private$make_probs(x[-(1:2)], x[1:2]) |> apply(c(1L,3L,private$ref_dims), sum)
          seq_along(private$sample_size) |>
            as.list() |>
            sapply(function(p){
              -dmultinom(private$tally[p,,], private$sample_size[p], probs[p,,], log=TRUE)
            }) |>
            sum()
        }, method=method, lower=lower, upper=upper)

      }else{
        estprev <- private$est_prev$par

        est_comp <- optim(c(se=0.75, sp=0.75), function(x){
          if(any(x<=0) || any(x>=1) || sum(x)<=1) return(Inf)
          probs <- private$make_probs(estprev, x) |> apply(c(1L,3L,private$ref_dims), sum)
          seq_along(private$sample_size) |>
            as.list() |>
            sapply(function(p){
              -dmultinom(private$tally[p,,], private$sample_size[p], probs[p,,], log=TRUE)
            }) |>
            sum()
        }, method=method, lower=lower, upper=upper)

      }

      private$est_comp <- est_comp
      invisible(est_comp)
    },

    estimate_comp_ind = function(method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B")){

      stopifnot(!is.null(private$est_prev$par))
      method <- match.arg(method)
      lower <- -Inf
      upper <- Inf
      if(method=="L-BFGS-B"){
        lower <- 1e-6
        upper <- 1.0 - 1e-6
      }

      estprev <- private$est_prev$par
      comptally <- private$tally |> apply(c(1L,2L), sum)

      est_comp <- optim(c(se=0.75, sp=0.75), function(x){
        if(any(x<=0) || any(x>=1) || sum(x)<=1) return(Inf)
        probs <- private$make_probs(estprev, x) |> apply(c(1L,3L), sum)
        seq_along(private$sample_size) |>
          as.list() |>
          sapply(function(p){
            -dmultinom(comptally[p,], private$sample_size[p], probs[p,], log=TRUE)
          }) |>
          sum()
      }, method=method, lower=lower, upper=upper)

      private$est_comp <- est_comp
      invisible(est_comp)
    },

    re_estimate_prev = function(){
      stopifnot(length(private$tally)==length(private$sample_size) * 2^(length(private$ref_dims)+1))
      comptally <- apply(private$tally,c(1L,2L),sum)
      stopifnot(length(comptally)==length(private$sample_size) * 2^length(private$ref_dims))
      browser()

      if(length(private$sample_size)==1L){
        est_prev <- optimise(function(x){
          if(x<=0 || x>=1) return(Inf)
          probs <- private$make_probs(x,private$est_comp$par) |> apply(c(1,3), sum)
          seq_along(length(private$sample_size)) |>
            as.list() |>
            sapply(function(p){
              -dmultinom(comptally[p,], private$sample_size[p], probs[p,], log=TRUE)
            }) |>
            sum()
        }, c(0,1))
        est_prev$par <- private$est_prev$minimum
      }else if(length(private$sample_size)==2L){
        est_prev <- optim(c(prev1=0.5, prev2=0.5), function(x){
          if(any(x<=0) || any(x>=1)) return(Inf)
          probs <- private$make_probs(x,private$est_comp$par) |> apply(c(1,3), sum)
          seq_along(private$sample_size) |>
            as.list() |>
            sapply(function(p){
              -dmultinom(comptally[p,], private$sample_size[p], probs[p,], log=TRUE)
            }) |>
            sum()
        })
      }

      invisible(est_prev)
    },

    estimate_hw = function(method=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B")){

      stopifnot(length(private$sample_size)>1L)
      method <- match.arg(method)
      lower <- -Inf
      upper <- Inf
      if(method=="L-BFGS-B"){
        lower <- 1e-6
        upper <- 1.0 - 1e-6
      }

      pars <- c(se1=0.75, sp1=0.75, se2=0.75, sp2=0.75, rep(0.5,length(private$sample_size)))
      names(pars) <- c("se1","sp1", "se2","sp2", str_c("prev",1:length(private$sample_size)))
      est_comp <- optim(pars, function(x){
        if(any(x<=0) || any(x>=1) || sum(x[1:2])<=1) return(Inf)
        probs <- private$make_probs(x[-(1:4)], x[1:2], refse=x[3], refsp=x[4]) |> apply(c(1L,3L,private$ref_dims), sum)
        seq_along(private$sample_size) |>
          as.list() |>
          sapply(function(p){
            -dmultinom(private$tally[p,,], private$sample_size[p], probs[p,,], log=TRUE)
          }) |>
          sum()
      }, method=method, lower=lower, upper=upper)

      invisible(est_comp)
    },

    calculate_ppp = function(use_comp=TRUE){

      stopifnot(!identical(private$est_comp, numeric(0)))
      stopifnot(!identical(private$est_comp, numeric(0)))

      ppp <- private$make_probs(private$est_prev$par, private$est_comp$par[1:2])

      if(use_comp){
        ppp <- apply(ppp, c(1L,2L,3L,private$ref_dims), sum)
        ppp <- ppp[,"Se",,] / (ppp[,"Se",,]+ppp[,"Sp",,])
      }else{
        ppp <- apply(ppp, c(1L,2L,private$ref_dims), sum)
        ppp <- ppp[,"Se",] / (ppp[,"Se",]+ppp[,"Sp",])
      }

      private$ppp <- ppp
      ppp
    },

    get_estimates = function(beta_n = 1e4, boot_n=0L){

      stopifnot(boot_n==0L)

      ppp <- self$calculate_ppp(use_comp=TRUE)
      dim(ppp) <- dim(private$tally)
      se <- (ppp*private$tally) |> apply(2,sum) |> rev()
      sp <- ((1-ppp)*private$tally) |> apply(2,sum)

      if(boot_n > 0 && beta_n > 0){
        stop("Need to bootstrap, otherwise the below can be replaced with TeachingDemos::hpd anyway")
        seci <- rbeta(beta_n, se[1]+1, se[2]+1) |> as.mcmc() |> HPDinterval() |> as.numeric()
        spci <- rbeta(beta_n, sp[1]+1, sp[2]+1) |> as.mcmc() |> HPDinterval() |> as.numeric()
      }else{
        seci=spci <- rep(NA_real_, 2L)
      }

      gs <- c(
        private$tally[4]/apply(private$tally,2,sum)[2],
        private$tally[1]/apply(private$tally,2,sum)[1]
      )

      pars <- self$estimate_comp(include_prev=FALSE)

      tribble(~Method, ~Parameter, ~Estimate, ~LowerCI, ~UpperCI,
        "GoldStandard", "Se", gs[1], NA_real_, NA_real_,
        "GoldStandard", "Sp", gs[2], NA_real_, NA_real_,
        "MLE-SepPrev", "Se", pars$par["se"], NA_real_, NA_real_,
        "MLE-SepPrev", "Sp", pars$par["sp"], NA_real_, NA_real_,
        "Mode-WithComp", "Se", se[1]/sum(se), seci[1], seci[2],
        "Mode-WithComp", "Sp", sp[1]/sum(sp), spci[1], spci[2],
      ) -> rv


      ppp <- self$calculate_ppp(use_comp=FALSE)
      dim(ppp) <- dim(private$tally)[-2L]
      se <- (ppp*private$tally[,-2L,]) |> apply(2,sum) |> rev()
      sp <- ((1-ppp)*private$tally[,-2L,]) |> apply(2,sum)

      if(boot_n > 0 && beta_n > 0){
        stop("Need to bootstrap")
        seci <- rbeta(beta_n, se[1]+1, se[2]+1) |> as.mcmc() |> HPDinterval() |> as.numeric()
        spci <- rbeta(beta_n, sp[1]+1, sp[2]+1) |> as.mcmc() |> HPDinterval() |> as.numeric()
      }else{
        seci=spci <- rep(NA_real_, 2L)
      }

      pars <- self$estimate_comp(include_prev=TRUE)

      if(length(private$sample_size)>1L){
        est_comp <- self$estimate_hw()
        hwe <- est_comp$par[1:2]
      }else{
        hwe <- rep(NA_real_, 2L)
      }

      rv |>
        bind_rows(
          tribble(~Method, ~Parameter, ~Estimate, ~LowerCI, ~UpperCI,
            "Mode-NoComp", "Se", se[1]/sum(se), seci[1], seci[2],
            "Mode-NoComp", "Sp", sp[1]/sum(sp), spci[1], spci[2],
            "MLE-CombPrev", "Se", pars$par["se"], NA_real_, NA_real_,
            "MLE-CombPrev", "Sp", pars$par["sp"], NA_real_, NA_real_,
            "MLE-HW", "Se", hwe[1], NA_real_, NA_real_,
            "MLE-HW", "Sp", hwe[2], NA_real_, NA_real_,
          )
        )

    }

  ),

  private = list(
    parameters = list(),
    prob_funs = list(),
    data = list(),
    est_partial = list(),


    reference_tests = tibble(),
    sample_size = integer(),
    dims = integer(),
    dimnames = list(),
    tally = array(),
    est_comp = double(),
    ppp = double(),
    make_probs = list()
  )
)

