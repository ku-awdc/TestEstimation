library("R6")
library("compiler")
library("coda")

TestEval <- R6Class(
  "TestEval",
  public = list(

    initialize = function(reference_tests, sample_size){
      stopifnot(nrow(reference_tests)==2L, c("Parameter","Estimate","TestName") %in% names(reference_tests), c("Se","Sp") %in% reference_tests$Parameter)

      private$reference_tests <- reference_tests
      private$sample_size <- sample_size
      private$ref_dims <- 4L

      ref_se <- reference_tests |> filter(Parameter=="Se") |> arrange(TestName) |> pull("Estimate")
      ref_sp <- reference_tests |> filter(Parameter=="Sp") |> arrange(TestName) |> pull("Estimate")

      private$make_probs <- function(prev, comparator, refse=ref_se, refsp=ref_sp){
        se <- c(comparator[1], refse)
        sp <- c(comparator[2], refsp)
        probs <- array(dim=c(length(sample_size),2,2,2), dimnames = list(str_c("Pop",seq_along(sample_size)),c("Se","Sp"),str_c("Comp",c("-","+")),str_c("Ref1",c("-","+"))))
        probs[,"Se","Comp-","Ref1-"] <- prev[]*(1-se[1])*(1-se[2])
        probs[,"Se","Comp+","Ref1-"] <- prev[]*se[1]*(1-se[2])
        probs[,"Se","Comp-","Ref1+"] <- prev[]*(1-se[1])*se[2]
        probs[,"Se","Comp+","Ref1+"] <- prev[]*se[1]*se[2]
        probs[,"Sp","Comp-","Ref1-"] <- (1-prev[])*sp[1]*sp[2]
        probs[,"Sp","Comp+","Ref1-"] <- (1-prev[])*(1-sp[1])*sp[2]
        probs[,"Sp","Comp-","Ref1+"] <- (1-prev[])*sp[1]*(1-sp[2])
        probs[,"Sp","Comp+","Ref1+"] <- (1-prev[])*(1-sp[1])*(1-sp[2])
        probs
      }
    },

    simulate_data = function(prev, comparator, reassign_prob=0){
      stopifnot(length(prev)==length(private$sample_size))
      probs <- private$make_probs(prev,comparator) |> apply(c(1L,3L, private$ref_dims), sum)
      tally <- array(dim=c(length(prev), 2^(length(private$ref_dims)+1L)))
      for(p in seq_len(length(prev))){
        tally[p,] <- rmultinom(1L, size=private$sample_size[p], prob=probs[p,,])
      }
      dim(tally) <- c(length(prev), rep(2L,length=length(private$ref_dims)+1L))

      if(reassign_prob>0){
        move1 <- rbinom(dim(tally)[1],tally[,2,1],reassign_prob)
        move2 <- rbinom(dim(tally)[1],tally[,1,2],reassign_prob)
        tally[,1,1] <- tally[,1,1] +move1
        tally[,2,1] <- tally[,2,1] -move1
        tally[,1,2] <- tally[,1,2] -move2
        tally[,2,2] <- tally[,2,2] +move2
      }

      self$set_data(tally)
      tally
    },

    set_data = function(tally){
      stopifnot(dim(tally)==c(length(private$sample_size),rep(2L,length=length(private$ref_dims)+1L)), apply(tally,1,sum)==private$sample_size)
      stopifnot(tally>=0)
      private$tally <- tally
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
    ref_dims = integer(),
    tally = array(),
    est_prev = double(),
    est_comp = double(),
    ppp = double(),
    make_probs = NULL
  )
)

