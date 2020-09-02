#!/usr/bin/env Rscript --vanilla
# set expandtab ts=4 sw=4 ai fileencoding=utf-8
#
# Author: JJ
# Maintainer(s): MG, PB, JJ, KL
# License: (c) HRDAG 2020, GPL v2 or newer
#
# Covid19-undercount/analyze/src/4-make_state_tables.R
# -----------------------------------------------------------
#

if (!require('pacman')){install.packages('pacman')}
pacman::p_load(xtable, lubridate, dplyr, ggplot2,
               argparse, glue, stringr, readr, tidyr, zoo)

stopifnot(str_ends(getwd(), "Covid19-modeling/analyze"))

set.seed(19481210)

# states <- c("OH", "AZ", "IL", "CA", "FL", "GA")
# state_group <- "1"
states <- c("TX", "PA")
state_group <- "2"
groups <- c("general", "incarcerated", "staff")

ngroups <- length(groups)
nstate <- length(states)
states <- sort(states)

parser <- ArgumentParser()
parser$add_argument("--p", default="0_01") # "0_01"
args <- parser$parse_args()

args[["state_RT1"]] <- glue("output_figures/state{state_group}_RT1_{args$p}.tex")
args[["state_R0"]] <- glue("output_figures/state{state_group}_R0_{args$p}.tex")
args[["state_undercount_graph"]] <- glue("output_figures/state{state_group}_undercount_daily_{args$p}.png")
args[["state_undercount_table_stub"]] <- glue("output_figures/state{state_group}_undercount_{args$p}")

# general population data
state.dat <- read_csv("input/us-states-20200821.csv")
abbrev <- read_csv("input/state_abbrev.csv")
pops <- read_csv("input/pop_data.csv")
pops <- pops %>% select(State,Pop2018)
state.dat <- state.dat %>%
  inner_join(abbrev, by="state") %>%
  inner_join(pops, by=c("state"="State"))
state.dat.general <- state.dat %>%
  filter(is.element(abbreviation, states)) %>%
  rename(pop = Pop2018) %>%
  select(date, abbreviation, pop, cases)
state.dat.general$group <- "general"

# incarcerated data
state.dat <- read_csv("input/covid_prison_cases.csv")
pops <- read_csv("input/prison_populations.csv")
state.dat$date <- mdy(as.character(state.dat$as_of_date))
pops <- pops %>% select(abbreviation,march_pop,april_pop,june_pop)
pops$pop <- round(
  rowMeans(pops[,c('march_pop', 'april_pop', 'june_pop')]),
  0
)
state.dat <- state.dat %>%
  filter(total_prisoner_cases > 0) %>%
  inner_join(pops, by=c("abbreviation"="abbreviation"))
state.dat.incarcerated <- state.dat %>%
  filter(is.element(abbreviation, states)) %>%
  rename(cases = total_prisoner_cases) %>%
  select(date, abbreviation, pop, cases)
state.dat.incarcerated$group <- "incarcerated"

# jail staff data
state.dat <- read_csv("input/covid_prison_cases.csv")
pops <- read_csv("input/staff_populations.csv")
state.dat$date <- mdy(as.character(state.dat$as_of_date))
pops <- pops %>% select(abbreviation,april_pop)
state.dat <- state.dat %>%
  filter(total_staff_cases > 0) %>%
  inner_join(pops, by=c("abbreviation"="abbreviation"))
state.dat.staff <- state.dat %>%
  filter(is.element(abbreviation, states)) %>%
  rename(cases = total_staff_cases, pop = april_pop) %>%
  select(date, abbreviation, pop, cases)
state.dat.staff$group <- "staff"

state.dat <- rbind(state.dat.general, state.dat.incarcerated, state.dat.staff)

state.dat <- state.dat %>%
  arrange(abbreviation, date)

RT1.mat <- matrix('', ngroups, nstate)
rownames(RT1.mat) <- groups
colnames(RT1.mat) <- states

R0.mat <- matrix('', ngroups, nstate)
rownames(R0.mat) <- groups
colnames(R0.mat) <- states

ctr <- 0
for (state in states) {
  for (group in groups) {
    ctr <- ctr + 1
    print(paste(state, group))
    if (exists('mcmc.out')) {rm(mcmc.out)}
    load(glue("output_{group}/mcmc_sim_{tolower(state)}_{args$p}.RData"))
    T <- length(mcmc.out$SIS[,1,1])
    RT1 <- (mcmc.out$BETA 
            * mcmc.out$PHIS 
            * mcmc.out$SIS[dim(mcmc.out$SIS)[1],1,]
            / (mcmc.out$GAMMAR))
    R0 <-  mcmc.out$BETA / mcmc.out$GAMMAR
    burn <- floor(length(mcmc.out$BETA) / 2)
    nmc <- length(mcmc.out$BETA)
    
    RT1.mat[which(group == groups), which(state == states)] <- paste(sprintf("%.2f", round(mean(RT1[burn:nmc]),2)), " (",
                                                             sprintf("%.2f", round(quantile(RT1[burn:nmc], 0.025), 2)), ",",
                                                             sprintf("%.2f", round(quantile(RT1[burn:nmc], 0.975), 2)), ")", sep='')
    
    R0.mat[which(group == groups), which(state == states)] <- paste(sprintf("%.2f", round(mean(R0[burn:nmc]), 2)), " (",
                                                            sprintf("%.2f", round(quantile(R0[burn:nmc], 0.025), 2)), ",",
                                                            sprintf("%.2f", round(quantile(R0[burn:nmc], 0.975), 2)), ")", sep='')
    
    df.tmp <- data.frame(date=date(ymd("20200101") + seq(from=0, to=T - 1)),
                         abbreviation=state,
                         group=group,
                         S.mean=apply(mcmc.out$SIS[,1,burn:nmc], 1, mean),
                         S.025=apply(mcmc.out$SIS[,1,burn:nmc], 1, quantile, probs=0.025),
                         S.975=apply(mcmc.out$SIS[,1,burn:nmc], 1, quantile, probs=0.975))
    state.dat.tmp <- state.dat %>%
      inner_join(df.tmp, by=c("abbreviation"="abbreviation", "date"="date", "group" = "group"))
    if (ctr==1) {
      all.state.dat <- state.dat.tmp
    } else {
      all.state.dat <- rbind(all.state.dat, state.dat.tmp)
    }

  }
}

all.state.dat <- all.state.dat %>%
  mutate(S.mean=S.mean * pop, S.025=S.025 * pop, S.975=S.975 * pop) %>%
  mutate(ucount.mean=(pop - S.mean) / cases,
         ucount.025=(pop - S.975) / cases,
         ucount.975=(pop - S.025) / cases,
         ucount.mean.lead5=(lag(pop, 5) - lag(S.mean, 5)) / cases,
         ucount.mean.lead10=(lag(pop, 10) - lag(S.mean, 10)) / cases)

print(xtable(RT1.mat,
             caption="Estimated posterior mean of $R_{T_1}$, with 95 percent posterior credible interval shown in parentheses 
             \\label{tab:RT1}"),
      file=args$state_RT1,
      size='footnotesize')

print(xtable(R0.mat,
             caption="Estimated posterior mean of $R_{0}$, with 95 percent posterior credible interval shown in parentheses 
             \\label{tab:R0}"),
      file=args$state_R0,
      size='footnotesize')

ulim <- all.state.dat %>%
  filter(date >= ymd("20200325")) %>% 
  summarise(ulim=max(ucount.mean))
plt <- ggplot(all.state.dat %>% filter(date >= ymd("20200325"), abbreviation %in% states),
              aes(x=date, y=ucount.mean, col=factor(group))) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~abbreviation, scales='free') + ylim(0, ulim$ulim)

png(args$state_undercount_graph, width=900, height=600)
print(plt)
dev.off()

ucount_ests_start <- all.state.dat %>% 
  filter(date >= ymd("20200325")) %>% 
  filter(date <= ymd("20200410")) %>% 
  select(abbreviation, group, ucount.mean, ucount.025, ucount.975) %>%
  group_by(abbreviation, group) %>%
  summarise_if(is.numeric, mean)
ucount_mns <- ucount_ests_start %>%
  select(abbreviation, group, ucount.mean) %>% 
  pivot_wider(names_from=abbreviation, values_from=ucount.mean) %>%
  select(-group)
ucount_025 <- ucount_ests_start %>%
  select(abbreviation, group, ucount.025) %>%
  pivot_wider(names_from=abbreviation, values_from=ucount.025) %>%
  select(-group)
ucount_975 <- ucount_ests_start %>%
  select(abbreviation, group, ucount.975) %>%
  pivot_wider(names_from=abbreviation, values_from=ucount.975) %>%
  select(-group)

ucount.mat <- matrix('', nstate, ngroups)
rownames(ucount.mat) <- states
colnames(ucount.mat) <- groups

for (group in groups) {
  for (state in states) {
    row <- which(group == groups)
    col <- which(state == states)
    ucount.mat[col, row] <- paste(sprintf("%.2f", round(ucount_mns[row,col], 2)), " (",
                                  sprintf("%.2f", round(ucount_025[row,col], 2)), ",",
                                  sprintf("%.2f", round(ucount_975[row,col], 2)), ")", sep='')
  }
}

print(xtable(ucount.mat,
             caption="Estimated posterior mean of undercount, with 95 percent posterior credible interval shown in parentheses 
             \\label{tab:undercount}"),
      file=glue("{args$state_undercount_table_stub}_start.tex"),
      size='footnotesize')

ucount_ests_end <- all.state.dat %>% 
  filter(date >= ymd("20200801")) %>% 
  filter(date <= ymd("20200820")) %>% 
  select(abbreviation, group, ucount.mean, ucount.025, ucount.975) %>%
  group_by(abbreviation, group) %>%
  summarise_if(is.numeric, mean)
ucount_mns <- ucount_ests_end %>%
  select(abbreviation, group, ucount.mean) %>% 
  pivot_wider(names_from=abbreviation, values_from=ucount.mean) %>%
  select(-group)
ucount_025 <- ucount_ests_end %>%
  select(abbreviation, group, ucount.025) %>%
  pivot_wider(names_from=abbreviation, values_from=ucount.025) %>%
  select(-group)
ucount_975 <- ucount_ests_end %>%
  select(abbreviation, group, ucount.975) %>%
  pivot_wider(names_from=abbreviation, values_from=ucount.975) %>%
  select(-group)

ucount.mat <- matrix('', nstate, ngroups)
rownames(ucount.mat) <- states
colnames(ucount.mat) <- groups

for (group in groups) {
  for (state in states) {
    row <- which(group == groups)
    col <- which(state == states)
    ucount.mat[col, row] <- paste(sprintf("%.2f", round(ucount_mns[row,col], 2)), " (",
                                  sprintf("%.2f", round(ucount_025[row,col], 2)), ",",
                                  sprintf("%.2f", round(ucount_975[row,col], 2)), ")", sep='')
  }
}

print(xtable(ucount.mat,
             caption="Estimated posterior mean of undercount, with 95 percent posterior credible interval shown in parentheses 
             \\label{tab:undercount}"),
      file=glue("{args$state_undercount_table_stub}_end.tex"),
      size='footnotesize')

# done.
