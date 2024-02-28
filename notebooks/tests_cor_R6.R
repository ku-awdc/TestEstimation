library("tidyverse")
library("vctrs")
library("R6")
library("compiler")
library("coda")

CorTestEval <- R6Class(
  "CorTestEval",
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

      ise <- c(0.9,0.99)
      isp <- c(0.99,0.99)
      cor_p <- 0.8
      cor_n <- 0.1
      trueprev <- 0.25

      find_ss(ise[1],isp[1],cor_p,cor_n)
      find_ss(ise[2],isp[2],cor_p,cor_n)

      test_fun <- function(prev, se, sp, cor_p, cor_n){
        iss1 <- find_iss(se[1], sp[1], cor_p, cor_n)
        iss2 <- find_iss(se[2], sp[2], cor_p, cor_n)
        ise <- c(iss1[1], iss2[1])
        isp <- c(iss1[2], iss2[2])
        if(any(ise>1) || any(ise<0) || any(isp>1) || any(isp<0)) stop('Invalid parameters')
        prev <- cor_p*prev + cor_n*(1-prev)
        prob <- numeric(4)
        prob[1] <- prev * ise[1]*ise[2] + (1-prev) *(1-isp[1])*(1-isp[2])
        prob[2] <- prev * (1-ise[1])*ise[2] + (1-prev) *isp[1]*(1-isp[2])
        prob[3] <- prev * (ise[1])*(1-ise[2]) + (1-prev)*(1-isp[1])*(isp[2])
        prob[4] <- prev * (1-ise[1])*(1-ise[2]) + (1-prev)*(isp[1])*(isp[2])
        dim(prob) <- c(2,2)
        prob
      }

      ise <- matrix(NA_real_, ncol=3, nrow=3)
      ise[1,2] <- 0.8
      ise[1,3] <- 0.9
      ise[2,2] <- 0.9
      ise[3,3] <- 0.95
      isp <- 1-0.25*(1-ise)

      tibble(Ind=seq_len(N[1])) |>
        mutate(Status=rbinom(n(),1,prev[1])) |>
        mutate(C12 = Status*cor_p[2] + (1-Status)*cor_n[2]) |>
        mutate(C13 = Status*cor_p[3] + (1-Status)*cor_n[3]) |>
        mutate(Test1a = rbinom(n(), 1, ise[1,2]*C12 + (1-isp[1,2])*(1-C12))) |>
        mutate(Test1b = rbinom(n(), 1, ise[1,3]*C13 + (1-isp[1,3])*(1-C13))) |>
        mutate(Test1 = pmax(Test1a, Test1b)) |>  ## Is this right?!? Maybe
        mutate(Test2 = rbinom(n(), 1, ise[2,2]*C12 + (1-isp[2,2])*(1-C12))) |>
        mutate(Test3 = rbinom(n(), 1, ise[3,3]*C13 + (1-isp[3,3])*(1-C13)))

      # p(+++ | C12+, C13+) = p(Test1 | C12+, C13+) x p(Test2 | C12+, C13+) p(Test3 | C12+, C13+)
      # p(Test1 | C12+, C13+)
      #       = p(Test1 | C12+) + p(Test1 | C13+) - p(Test1 | C12+) x p(Test1 | C13+)
      # p(Test2 | C12+, C13+) = p(Test2 | C12+) + p(Test2 | C13+) - p(Test2 | C12+) x p(Test2 | C13+)
      # p(Test2 | C13+) = p(Test2)
      # p(Test3 | C12+, C13+) = p(Test3 | C13+) ...



      test_fun_3t <- function(prev, se, sp, cor_p, cor_n){

        iss12 <- find_iss(se[1], sp[1], cor_p[2], cor_n[2])
        iss13 <- find_iss(se[1], sp[1], cor_p[3], cor_n[3])
        iss2 <- find_iss(se[2], sp[2], cor_p[2], cor_n[2])
        iss3 <- find_iss(se[3], sp[3], cor_p[3], cor_n[3])
        ise12 <- c(iss12[1], iss2[1])
        ise13 <- c(iss13[1], iss3[1])
        isp12 <- c(iss12[2], iss2[2])
        isp13 <- c(iss13[2], iss3[2])

        if(any(ise12>1) || any(ise12<0) || any(isp12>1) || any(isp12<0)) stop('Invalid parameters')
        if(any(ise13>1) || any(ise13<0) || any(isp13>1) || any(isp13<0)) stop('Invalid parameters')

        prev12 <- cor_p[2]*prev + cor_n[2]*(1-prev)
        prev13 <- cor_p[3]*prev + cor_n[3]*(1-prev)

        prob <- numeric(8)

        cp <- function(a,b) a+b-a*b

        prev *
          (
            cor_p[2]*ise12[1]*ise12[2] * cor_p[3]*ise13[1]*ise13[2] +
            (1-cor_p[2])*(1-isp12[1])*(1-isp12[2]) * cor_p[3]*ise13[1]*ise13[2] +
            (1-cor_p[2])*(1-isp12[1])*(1-isp12[2]) * cor_p[3]*ise13[1]*ise13[2] +

        cp(prev12*ise12[1]*ise12[2], prev13*ise13[1]*ise13[2])

        prob[1] <- prev12 * ise12[1]*ise12[2]* +
                  (1-prev12) *(1-isp[1])*(1-isp[2])
        prob[2] <- prev * (1-ise[1])*ise[2] + (1-prev) *isp[1]*(1-isp[2])
        prob[3] <- prev * (ise[1])*(1-ise[2]) + (1-prev)*(1-isp[1])*(isp[2])
        prob[4] <- prev * (1-ise[1])*(1-ise[2]) + (1-prev)*(isp[1])*(isp[2])



        dim(prob) <- c(2,2)
        prob
      }

      se <- c(0.7,0.5)
      sp <- c(0.9,0.95)
      prev <- c(0.1,0.5)
      N <- c(500L,500L)
      cor_p <- max(se)
      cor_n <- 1-max(sp)
      test_fun(prev[2], se, sp, cor_p, cor_n)
      test_fun(prev[2], se, sp, 1, 0)

      t1 <- find_ss(0.9,0.99,0.8,0.01)
      t2 <- find_ss(0.6,0.99,0.8,0.01)
      t1; t2

      test_fun(prev[2], c(t1[1],t2[1]), c(t1[2],t2[2]), 0.8, 0.01)
      test_fun(prev[2], c(t1[1],t2[1]), c(t1[2],t2[2]), 1, 0)


      test_fun_old(prev[1], se, sp, 1, 0)
      test_fun_old(prev[1], se, sp, 0.7, 0.05)

      #cors <- c(1,0)
      #cors <- c(0.8,0.05)
      #cors <- c(0.9,0.05)

      cors <- c(0.5,0.04)
      test_fun(prev[1],se,sp,cors[1],cors[2]); test_fun(prev[1],se,sp,1,0);
      test_fun(prev[2],se,sp,cors[1],cors[2]); test_fun(prev[2],se,sp,1,0);



      test_fun(0.2, se, sp, cor_p, cor_n)

      sum(prob)

      prev <- trueprev
      #prob <- numeric(4)
      prob[1] <- cor_p*prev * ise[1]*ise[2] + cor_n*(1-prev) * ise[1]*ise[2] +
                cor_p*prev * (1-isp[1])*(1-isp[2]) + cor_n*(1-prev) * (1-isp[1])*(1-isp[2])
      #prob[2] <- prev * (1-ise[1])*ise[2] + (1-prev) *isp[1]*(1-isp[2])
      #prob[3] <- prev * (ise[1])*(1-ise[2]) + (1-prev)*(1-isp[1])*(isp[2])
      #prob[4] <- prev * (1-ise[1])*(1-ise[2]) + (1-prev)*(isp[1])*(isp[2])
      dim(prob) <- c(2,2)
      prob
      sum(prob)


      fun_txt <- vec_c(
        "function(prev, se, sp, cor_p, cor_n){",
        "iss1 <- find_iss(se[1], sp[1], cor_p, cor_n)",
        "iss2 <- find_iss(se[2], sp[2], cor_p, cor_n)",
        "ise <- c(iss1[1], iss2[1])",
        "isp <- c(iss1[2], iss2[2])",
        "if(any(ise>1) || any(ise<0) || any(isp>1) || any(isp<0)) stop('Invalid parameters')",
        "prob <- numeric(4)",
        "prob[1] <- prev * (cor_p*ise[1]*ise[2] + (1-cor_p)*(1-isp[1])*(1-isp[2])) +
        (1-prev) * (cor_n*(1-ise[1])*(1-ise[2]) + (1-cor_n)*isp[1]*isp[2])",
        "prob[2] <- prev * (cor_p*(1-ise[1])*ise[2] + (1-cor_p)*isp[1]*(1-isp[2])) +
        (1-prev) * (cor_n*ise[1]*(1-ise[2]) + (1-cor_n)*(1-isp[1])*isp[2])",
        "prob[3] <- prev * (cor_p*(ise[1])*(1-ise[2]) + (1-cor_p)*(1-isp[1])*(isp[2])) +
        (1-prev) * (cor_n*(1-ise[1])*(ise[2]) + (1-cor_n)*(isp[1])*(1-isp[2]))",
        "prob[4] <- prev * (cor_p*(1-ise[1])*(1-ise[2]) + (1-cor_p)*(isp[1])*(isp[2])) +
        (1-prev) * (cor_n*ise[1]*ise[2] + (1-cor_n)*(1-isp[1])*(1-isp[2]))",
        "dim(prob) <- c(2,2)",
        "prob",
        "}"
      )
      fun_txt |> cat(sep="\n")
      str_c(fun_txt, collapse="\n") |> parse(text=_) |> eval() |> cmpfun() -> test_fun_old

      hw_txt <- vec_c(
        "function(prev, ise, isp){",
        "prob <- numeric(4)",
        "prob[1] <- prev * ise[1]*ise[2] + (1-prev) *(1-isp[1])*(1-isp[2])",
        "prob[2] <- prev * (1-ise[1])*ise[2] + (1-prev) *isp[1]*(1-isp[2])",
        "prob[3] <- prev * (ise[1])*(1-ise[2]) + (1-prev)*(1-isp[1])*(isp[2])",
        "prob[4] <- prev * (1-ise[1])*(1-ise[2]) + (1-prev)*(isp[1])*(isp[2])",
        "dim(prob) <- c(2,2)",
        "prob",
        "}"
      )
      hw_txt |> cat(sep="\n")
      str_c(hw_txt, collapse="\n") |> parse(text=_) |> eval() |> cmpfun() -> hw_fun

      probs <- test_fun(0.2,c(0.5,0.8),c(0.9,0.9),0.8,0.1)
      test_fun(0.2,c(0.5,0.8),c(0.9,0.9),1,0)

      tt <- function(i){

        se <- c(0.7,0.5)
        sp <- c(0.9,0.95)
        prev <- c(0.1,0.5)
        N <- c(500L,500L)

        cors <- c(1,0)
        cors <- c(0.8,0.05)
        #cors <- c(0.9,0.05)

        #cors <- c(0.85,0.04)
        test_fun(prev[1],se,sp,cors[1],cors[2]); test_fun(prev[1],se,sp,1,0);
        test_fun(prev[2],se,sp,cors[1],cors[2]); test_fun(prev[2],se,sp,1,0);

      data <- matrix(c(
        rmultinom(1L, N[1], test_fun(prev[1],se,sp,cors[1],cors[2])),
        rmultinom(1L, N[2], test_fun(prev[2],se,sp,cors[1],cors[2]))
      ), ncol=2)

      pars <- c(prev1=0.5, prev2=0.5, se1=0.75, sp1=0.75, cor_p=1, cor_n=0)
      optim(pars, function(pars){
        if(any(pars<0) || any(pars>1) || ((pars["se1"]+pars["sp1"]) < 1)) return(Inf)
#        if(pars["cor_p"] < pars["se1"] || pars["cor_p"] < se[2]) return(Inf)
#        if((1-pars["cor_n"]) < pars["sp1"] || (1-pars["cor_n"]) < sp[2]) return(Inf)

        probs <- try(test_fun(pars["prev1"], c(pars["se1"], se[2]), c(pars["sp1"], sp[2]), pars["cor_p"], pars["cor_n"]), silent=TRUE)
        if(inherits(probs, "try-error")) return(Inf)
        if(any(probs<0) || any(probs>1)) return(Inf)
        ll1 <- -dmultinom(data[,1], N[1], probs, log=TRUE) |> sum()

        probs <- try(test_fun(pars["prev2"], c(pars["se1"], se[2]), c(pars["sp1"], sp[2]), pars["cor_p"], pars["cor_n"]), silent=TRUE)
        if(any(probs<0) || any(probs>1)) return(Inf)
        if(inherits(probs, "try-error")) return(Inf)
        ll2 <- -dmultinom(data[,2], N[2], probs, log=TRUE) |> sum()

        ll1+ll2
      })$par |> as.list() |> as_tibble() |> mutate(mod="fancy", iteration=i) ->
        res1

      pars <- c(prev1=0.5, prev2=0.5, se1=0.75, sp1=0.75, se2=0.75, sp2=0.75)
      optim(pars, function(pars){
        if(any(pars<0) || any(pars>1) || ((pars["se1"]+pars["sp1"]) < 1)) return(Inf)
        if((pars["se2"]+pars["sp2"]) < 1) return(Inf)

        probs <- hw_fun(pars["prev1"], c(pars["se1"], pars["se2"]), c(pars["sp1"], pars["sp2"]))
        ll1 <- -dmultinom(data[,1], N[1], probs, log=TRUE) |> sum()

        probs <- hw_fun(pars["prev2"], c(pars["se1"], pars["se2"]), c(pars["sp1"], pars["sp2"]))
        ll2 <- -dmultinom(data[,2], N[2], probs, log=TRUE) |> sum()

        ll1+ll2
      })$par |> as.list() |> as_tibble() |> mutate(mod="hw", iteration=i) ->
        res2

      cors <- find_cors(se[2], sp[2], max(se[2],res2$se2), max(sp[2],res2$sp2))
      ss <- find_ss(res2$se1, res2$sp1, cors[1], cors[2])
      res3 <- res2
      res3$se1 <- ss[1]
      res3$sp1 <- ss[2]
      res3$mod <- "posthoc"

      bind_rows(res1,res2,res3)

      }

      out <- pbapply::pblapply(1:500, tt, cl=6L) |> bind_rows()
      summary(out |> filter(mod=="hw"))
      summary(out |> filter(mod=="fancy"))
      summary(out |> filter(mod=="posthoc"))





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

