library("tidyverse")
library("R6")
library("compiler")
library("coda")

TestEval <- R6Class(
  "TestEval",
  public = list(

    initialize = function(reference_tests, sample_size){
      stopifnot(c("Parameter","Estimate","TestName") %in% names(reference_tests), c("Se","Sp") %in% reference_tests$Parameter)

      n_tests <- (nrow(reference_tests)/2L) + 1L
      stopifnot(n_tests %in% 2L:3L)
      n_pops <- length(sample_size)

      private$reference_tests <- reference_tests
      private$sample_size <- sample_size
      dims <- list(comp=1L, ref=seq_len(n_tests-1L)+1L, pop=n_tests+1L, type=n_tests+2L)
      private$dims <- dims

      ref_se <- reference_tests |> filter(Parameter=="Se") |> arrange(TestName) |> pull("Estimate")
      ref_sp <- reference_tests |> filter(Parameter=="Sp") |> arrange(TestName) |> pull("Estimate")

      cov_to_cor <- function(test, cov){
        mincor <- (test[1]-1)*(1-test[])
        maxcor <- pmin(test[1],test[])-test[1]*test[]
        cor <- cov
        cor[cov < 0] <- -cov[cov < 0] / mincor[cov < 0]
        cor[cov > 0] <- cov[cov > 0] / maxcor[cov > 0]
        cor[1] <- 0
        cor
      }
      cor_to_cov <- function(test, cor){
        mincor <- (test[1]-1)*(1-test[])
        maxcor <- pmin(test[1],test[])-test[1]*test[]
        cov <- cor
        cov[cor < 0] <- -cor[cor < 0] * mincor[cov < 0]
        cov[cor > 0] <- cor[cor > 0] * maxcor[cov > 0]
        cov[1] <- 0
        cov
      }

      str_c("Test",seq_len(n_tests)) |>
        as.list() |>
        set_names() |>
        lapply(\(x) c(0,1)) |>
        do.call("expand_grid", args=_) |>
        (\(x) mutate(x, Combo=apply(x,1,str_c,collapse="")))() |>
        select(Combo, everything()) ->
        status

      status

      status |>
        mutate(Comp = Test1) |>
        pivot_longer(!c(Combo,Comp), names_to="Test", values_to="Status") |>
        mutate(Index = str_replace(Test, "Test", "")) |>
        mutate(SeText = if_else(Status==1L, str_c("se[",Index,"]"), str_c("(1-se[",Index,"])"))) |>
        mutate(SpText = if_else(Status==1L, str_c("(1-sp[",Index,"])"), str_c("sp[",Index,"]"))) |>
        mutate(CovSeText = case_when(
          Test=="Test1" ~ "",
          Status==Comp ~ str_c("+dse[", Index, "]"),
          Status!=Comp ~ str_c("-dse[", Index, "]"),
        )) |>
        mutate(CovSpText = case_when(
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

      str_c("T",seq_len(n_tests)) |>
        as.list() |>
        lapply(\(x) str_c(x, "_", c(0,1))) |>
        c(list(str_c("P",seq_len(n_pops)), c("Se","Sp"))) ->
        dn
      private$dimnames <- dn
      names(private$dimnames) <- names(private$dims)

      fun_txt <- c(
        "function(prev, se, sp, corse, corsp){",
        str_c("stopifnot(length(prev)==", n_pops, "L, length(se)==", n_tests, "L, length(sp)==", n_tests, "L)"),
        str_c("stopifnot(length(corse)==", n_tests, "L, length(corsp)==", n_tests, "L)"),
        "is_valid <- all(corse>-1, corse<1, corsp>-1, corsp<1, se>0, se<1, sp>0, sp<1, prev>0, prev<1)",
        str_c("dse <- corse#cor_to_cov(se, corse)"),
        str_c("dsp <- corsp#cor_to_cov(sp, corsp)"),
        str_c("se_probs=sp_probs <- numeric(", nrow(test_txt), "L)"),
        str_c("se_probs[", seq_len(nrow(test_txt)), "] <- ", test_txt$SeText),
        str_c("sp_probs[", seq_len(nrow(test_txt)), "] <- ", test_txt$SpText),
        "is_valid <- is_valid && all(se_probs >= 0 & se_probs <= 1 & sp_probs >= 0 & sp_probs <= 1)",
        str_c("prev <- rep(prev, each=", nrow(test_txt), "L)"),
        str_c("probs <- array(c(prev[]*se_probs, (1-prev[])*sp_probs), dim=c(", str_c(rep("2",n_tests),collapse=","), ",",n_pops, ",2))"),
        "attr(probs, 'is_valid') <- is_valid",
        str_c("attr(probs, 'dimnames') <- ", str_c(capture.output(dput(dn)), collapse="")),
        "return(probs)",
        "}"
      )
      str_c(fun_txt, collapse="\n") |> parse(text=_) |> eval() |> cmpfun() -> test_fun

      private$make_probs$full <- function(prev, comparator, corse=0, corsp=0, refse=ref_se, refsp=ref_sp){
        se <- c(comparator[1], refse)
        sp <- c(comparator[2], refsp)
        if(length(corse)==1L) corse <- rep(corse, n_tests)
        if(length(corsp)==1L) corsp <- rep(corsp, n_tests)
        probs <- test_fun(prev, se, sp, corse, corsp)
        probs
      }

      private$make_probs$sim <- function(prev, comparator, corse=0, corsp=0){
        probs <- private$make_probs$full(prev, comparator, corse=corse, corsp=corsp)
        rv <- probs |> apply(private$dims[c("comp","ref","pop")] |> unlist(), sum)
        dim(rv) <- c(2^(length(private$dims$ref)+1L), length(prev))
        attr(rv, "is_valid") <- attr(probs, "is_valid", exact=TRUE)
        rv
      }
    },

    simulate_data = function(prev, comparator, corse=0, corsp=0){
      stopifnot(length(prev)==length(private$sample_size))
      probs <- private$make_probs$sim(prev,comparator,corse=corse,corsp=corsp)
      if(!attr(probs, "is_valid")) stop("Invalid probabilities generated")

      tally <- array(dim=c(2^(length(private$dims$ref)+1L), length(prev)))
      for(p in seq_len(length(prev))){
        tally[,p] <- rmultinom(1L, size=private$sample_size[p], prob=probs[,p])
      }
      dim(tally) <- c(rep(2L,length=length(private$dims$ref)+1L), length(prev))
      dimnames(tally) <- private$dimnames[c("comp","ref","pop")] |> set_names("")

      self$set_data(tally)
      tally
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
    reference_tests = tibble(),
    sample_size = integer(),
    dims = integer(),
    dimnames = list(),
    tally = array(),
    est_prev = double(),
    est_comp = double(),
    ppp = double(),
    make_probs = list()
  )
)

