#!/usr/bin/env Rscript --vanilla
# set expandtab ts=4 sw=4 ai fileencoding=utf-8
#
# Author: JG
#
# Covid19-undercount/analyze/src/prison-population-sort.R
# -----------------------------------------------------------
#

# US States with highest sum of incarcerated population (mean from March/April/June)
# and staff population (from April):
# TX (170459)
# CA (159690)
# FL (113527)
# GA (61690) 
# OH (60160)
# NY (59944) 
# PA (59268) 
# AZ (50255) 
# MI (49183) 
# IL (46296) 

if (!require('pacman')){install.packages('pacman')}
pacman::p_load(readr, dplyr)

set.seed(19481210)

incarcerated_pop <- read_csv("input/prison_populations.csv") 
staff_pop <- read_csv("input/staff_populations.csv") 

incarcerated_pop$mean_pop <- rowMeans(
        incarcerated_pop[,c('march_pop', 'april_pop', 'june_pop')]
)
incarcerated_pop <- incarcerated_pop %>% select(abbreviation, mean_pop)

pop <- incarcerated_pop %>% 
        inner_join(staff_pop, by = 'abbreviation') %>%
        rename(incarcerated_pop = mean_pop, staff_pop = april_pop)
        
pop$total_pop <- pop$incarcerated_pop + pop$staff_pop
pop <- pop[order(-pop$total_pop),]
print(head(pop, 11))

# cases <- read_csv("input/covid_prison_cases.csv") 
# cases <- cases %>%
#         group_by(abbreviation) %>%
#         mutate(max_staff_deaths = max(total_staff_deaths, na.rm = TRUE)) %>%
#         mutate(max_incarcerated_deaths = max(total_prisoner_deaths, na.rm = TRUE)) %>%
#         select(abbreviation, max_staff_deaths, max_incarcerated_deaths) %>%
#         distinct()
# prev <- pop %>% inner_join(cases, by = 'abbreviation')
# prev$staff_prev <- prev$max_staff_deaths / prev$staff_pop
# prev$incarcerated_prev <- prev$max_incarcerated_deaths / prev$incarcerated_pop
# prev$prev_ratio <- prev$incarcerated_prev / prev$staff_prev 
# prev <- prev %>% select(abbreviation, total_pop, staff_prev, incarcerated_prev, prev_ratio)
# prev <- prev[order(-prev$total_pop),]
# print(head(prev, 11))
