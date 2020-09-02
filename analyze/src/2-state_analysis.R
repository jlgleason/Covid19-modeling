#!/usr/bin/env Rscript --vanilla
# set expandtab ts=4 sw=4 ai fileencoding=utf-8
#
# Author: KL JJ
# Maintainer(s): PB, MG 
# License: (c) HRDAG 2020, GPL v2 or newer
#
# Covid19-undercount/analyze/src/2-state_analysis.R
# -----------------------------------------------------------
#

if (!require('pacman')){install.packages('pacman')}
pacman::p_load(argparse, readr, tidyr, dplyr, stringr, 
               ggplot2, lubridate, glue, pracma)

stopifnot(str_ends(getwd(), "Covid19-modeling/analyze"))

set.seed(19481210)

#Restaurant: 
# TX 21, CA 15, FL 20, GA 24, OH 15, NY 16, PA 16, AZ 20, MI 16, IL 16
#School:
# TX 21, CA 13, FL 17, GA 18, OH 17, NY 18, PA 16, AZ 16, MI 16, IL 17
# Convergence Errors:
# NY, MI

##-----GLOBAL vars
##-----set initial conditions to get reasonable starting point for MCMC
States.sd <- data.frame(
  State=c("TX", "CA", "FL", "GA", "OH", "NY", "PA", "AZ", "MI", "IL"),
  date.sd=c("20200321", "20200315", "20200320", "20200324", "20200317",
              "20200318", "20200316", "20200320", "20200316", "20200317"),
  date.sd.prison= c("20200421", "20200415", "20200420", "20200424", "20200417",
                    "20200418", "20200416", "20200420", "20200416", "20200417"),
  T0=rep(30, 10),
  R0=c(2.5, 2.3, 2.5, 2.5, 2.5, 2.8, 2.5, 2.5, 2.5, 2.5),
  phi=rep(0.4, 10),
  phi.prison=rep(0.6, 10)
)

# State <- "CA"
smooth.deaths <- TRUE

states <- c("AZ", "IL")
ps <- c(0.005, 0.01)

for (State in states) {
  for (p in ps) {
    
    T0 <- States.sd$T0[States.sd$State == State]
    R0.init <- States.sd$R0[States.sd$State == State]
    p_str <- gsub('\\.', '_', p)
    
    print(paste(State, p_str))
    
    getargs <- function() {
      # sets file paths in one place, could be used in command line call
      # works in RStudio, too
      parser <- argparse::ArgumentParser()
      parser$add_argument("--funs", default="src/functions.R")
      parser$add_argument("--pop", default="staff") # 'general', incarcerated', 'staff'
      parser$add_argument("--theta", default="output_general/theta.rds")
      parser$parse_args()
    }
    args <- getargs()
    
    args[["mcmc_output"]] <- glue("output_{args$pop}/mcmc_sim_{tolower(State)}_{p_str}.RData")
    args[["posterior_hists"]] <- glue("output_{args$pop}/posterior_hists_{tolower(State)}_{p_str}.png")
    args[["posterior_trace"]] <- glue("output_{args$pop}/posterior_trace_{tolower(State)}_{p_str}.png")
    
    source(args$funs)
     
    if (args$pop == 'general') {
      datas <- load.data.general(State = State)
      date.sd <- States.sd$date.sd[States.sd$State == State]
      phi <- States.sd$phi[States.sd$State == State]
    } else if (args$pop == 'incarcerated') {
      datas <- load.data.incarcerated(State = State)
      date.sd <- States.sd$date.sd.prison[States.sd$State == State]
      phi <- States.sd$phi.prison[States.sd$State == State]
    } else if (args$pop == 'staff') {
      datas <- load.data.jail.staff(State = State)
      date.sd <- States.sd$date.sd.prison[States.sd$State == State]
      phi <- States.sd$phi.prison[States.sd$State == State]
    } else {
      stop("'pop' argument must be one of 'general', 'incarcerated', or 'staff'")
    }
    deaths <- datas$deaths
    pops <- datas$pops
    N <- datas$N
    
    n.to.add <- as.numeric(min(deaths$date) - ymd("20200101"))
    Ds <- c(rep(0, n.to.add), deaths$death_increase)
    
    T1 <- as.numeric(n.to.add + ymd(date.sd) - min(deaths$date))
    I <- 1 / N
    gammar <- 1 / 6
    zeta <- c(mean=6.4, sd=1.5, a=3.4, b=9.4) # limits of uniform prior on gammar^(-1)
    beta <- gammar * R0.init
    xi <- c(mean=2.5, sd=1.5, a=1, b=4) # prior on beta | gammar  is Uniform(xi[1]*gammar,xi[2]*gammar) on this interval
    zeta.T0 <- c(1, 90)
    zeta.phi <- c(0.01, 0.99)
    
    theta <- readRDS(args$theta)
    T <- length(Ds)
    t.max <- length(Ds)
    i.max <- length(theta)
    max.time <- length(theta) - 1
    
    n.thin <- 1000
    
    sim1 <- simulate.sir.timebreak(beta, gammar, 1 / N, N, T0, T, T1, phi)
    Nus <- sim1$Nus
    Ds.sim <- calculate.death.means(Nus, theta, p)
    
    plot(Ds.sim[1:T], Ds, xlim=c(0, max(Ds.sim[1:T], Ds)))
    abline(0, 1) #looks pretty good!
    
    c0 <- c(0.002^2, 0.002^2, 0.75^2, 0.002^2)
    c1 <- (2.38 / 4) * 0.5 
    nmc <- 100000
    burn <- 50000
    k.adapt <- 20000
    thinned <- round(seq(burn, nmc, length.out=n.thin))
    
    adaptive <- TRUE
    
    # if you aren't sure you have the tuning parameters right, then set to TRUE to see trace plots
    # the plotting slows things down, so if you are confident it will mix well, then don't make the plots as it runs
    plotting <- FALSE
    disp.int <- 5000

    mcmc.out <- run.mcmc.state(nmc, theta, beta, gammar, T, I, Ds, Nus, xi, zeta,
                               disp.int, k.adapt=k.adapt, c0=c0, c1=c1, T0=T0,
                               zeta.T0=zeta.T0, phi, T1, zeta.phi, plotting)
    
    save(file=args$mcmc_output, mcmc.out)
    load(file=args$mcmc_output)
    
    df <- data.frame(beta=mcmc.out$BETA,
                     gammar=mcmc.out$GAMMAR,
                     T0=mcmc.out$T0S,
                     phi=mcmc.out$PHI,
                     SIt=mcmc.out$SIS[dim(mcmc.out$SIS)[1],1,])
    df <- df %>% mutate(R0=beta / gammar, RT1=(R0 * phi) * SIt)
    df$iter <- seq(nmc)
    
    df <- df %>% pivot_longer(-iter)
    
    png(args$posterior_hists, width=600, height=400)
    print({ggplot(df %>% filter(iter >= burn, name %in% c("R0", "RT1", "T0", "phi")),
                  aes(x=value)) +
        geom_histogram() + 
        facet_wrap(~name, scales='free') +
        ggtitle(pops$State[pops$abbreviation == State])})
    dev.off()
    
    png(args$posterior_trace, width=600, height=400)
    print({ggplot(df %>% filter(name %in% c("beta", "gammar", "R0")),
                  aes(x=iter, y=value)) +
        geom_point(size=0.2) + 
        facet_wrap(~name, scales='free')})
    dev.off()

  }
}

# done.
