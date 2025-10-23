#-------------------------------------------------------------#
# MSc Health Data Science 2023                                #
# Dataset: Mtb Within-host genetic & tx data                  #
# Frequency of resistance acquisition                         #
#-------------------------------------------------------------#

# FC: code edited on 22/10/2025 to calculate DR acquisition rates in different subsets (pan-S, MDR and XDR populations)

# install.packages("tidyverse", type = "source")
# install.packages("freqtables", type = "source")

library(tidyverse)
library(freqtables)

# setwd("~/Desktop/LSHTM/Summer Project/analysis")
setwd("/Users/francesccoll/ibv/projects/7.mtb_wh_paper/manuscript/tables_and_figures/acquisition_rates_DR/")

# 1. load data----
# data <- readRDS('datasets/tx_fixedmut_tbprof.Rds')
data <- readRDS('tx_fixedmut_tbprof.Rds')

# write.table(data, file = "tx_fixedmut_tbprof.csv", sep = "\t", col.names = T, row.names = F, quote = F)

# 2. create mini dataframe for each drug/gene profile----
# amikacin
am_freq <- data %>% 
  select(strain_id, Am, Rv0529, Rv2416c, Rv3106, MTB000019, Rv3805c, Rv3862c, Rv3197A) %>%
  mutate(num_am_mutations = rowSums(.[3:9]),
         num_strains_am_mutations = case_when(Rv0529 == 1 | Rv2416c == 1 | Rv3106 == 1 | MTB000019 == 1 | Rv3805c == 1 | Rv3862c == 1 | Rv3197A == 1 ~ 1,
                                              T ~ 0)) %>% 
  filter(Am==1) %>% 
  freq_table(num_strains_am_mutations)

# bedaquiline
bdq_freq <- data %>% 
  select(strain_id, Bdq, Rv1305, Rv2535c, Rv1979c, Rv0678, Rv0676c, Rv0677c) %>%
  mutate(num_bdq_mutations = rowSums(.[3:8]),
         num_strains_bdq_mutations = case_when(Rv1305 == 1 | Rv2535c == 1 | Rv1979c == 1 | Rv0678 == 1 | Rv0676c == 1 | Rv0677c == 1 ~ 1,
                                              T ~ 0)) %>% 
  filter(Bdq==1) %>% 
  freq_table(num_strains_bdq_mutations)

# capreomycin
cm_freq <- data %>% 
  select(strain_id, Cm, Rv0529, Rv3106, MTB000019, Rv0668, Rv3805c, Rv3862c) %>%
  mutate(num_cm_mutations = rowSums(.[3:8]),
         num_strains_cm_mutations = case_when(Rv0529 == 1 | Rv3106 == 1 | MTB000019 == 1 | Rv0668 == 1 | Rv3805c == 1 | Rv3862c == 1 ~ 1,
                                               T ~ 0)) %>% 
  filter(Cm==1) %>% 
  freq_table(num_strains_cm_mutations)

# clofazimine
cfz_freq <- data %>% 
  select(strain_id, Cfz, Rv2535c, Rv1979c, Rv0678, Rv0676c, Rv0677c) %>%
  mutate(num_cfz_mutations = rowSums(.[3:7]),
         num_strains_cfz_mutations = case_when(Rv2535c == 1 | Rv1979c == 1 | Rv0678 == 1 | Rv0676c == 1 | Rv0677c == 1 ~ 1,
                                              T ~ 0)) %>% 
  filter(Cfz==1) %>% 
  freq_table(num_strains_cfz_mutations)

# cycloserine
cs_freq <- data %>% 
  select(strain_id, Cs, Rv2780, Rv3423c) %>%
  mutate(num_cs_mutations = rowSums(.[3:4]),
         num_strains_cs_mutations = case_when(Rv2780 == 1 | Rv3423c == 1 ~ 1,
                                               T ~ 0)) %>% 
  filter(Cs==1) %>% 
  freq_table(num_strains_cs_mutations)

# delamanid
# dld_freq <- data %>% 
#   select(strain_id, Dld, Rv3547, Rv3261, Rv3262, Rv1173, Rv2983, Rv0407) %>%
#   mutate(num_dld_mutations = rowSums(.[3:8]),
#          num_strains_dld_mutations = case_when(Rv3547 == 1 | Rv3261 == 1 | Rv3262 == 1 | Rv1173 == 1 | Rv2983 == 1 | Rv0407 == 1 ~ 1,
#                                               T ~ 0)) %>% 
#   filter(Dld==1) %>% 
#   freq_table(num_strains_dld_mutations)

# ethambutol
eth_freq <- data %>% 
  select(strain_id, E, Rv3795, Rv3793, Rv3794, Rv3792, Rv1267c, Rv3806c) %>%
  mutate(num_e_mutations = rowSums(.[3:8]),
         num_strains_e_mutations = case_when(Rv3795 == 1 | Rv3793 == 1 | Rv3794 == 1 | Rv3792 == 1 | Rv1267c == 1 | Rv3806c == 1 ~ 1,
                                               T ~ 0)) %>% 
  filter(E==1) %>% 
  freq_table(num_strains_e_mutations)

# ethionamide
eto_freq <- data %>% 
  select(strain_id, Eto, Rv3083, Rv1484, Rv0486, Rv1854c, Rv3854c, Rv3855) %>%
  mutate(num_eto_mutations = rowSums(.[3:8]),
         num_strains_eto_mutations = case_when(Rv3083 == 1 | Rv1484 == 1 | Rv0486 == 1 | Rv1854c == 1 | Rv3854c == 1 | Rv3855 == 1 ~ 1,
                                             T ~ 0)) %>% 
  filter(Eto==1) %>% 
  freq_table(num_strains_eto_mutations)

# isoniazid
inh_freq <- data %>% 
  select(strain_id, H, Rv1908c, Rv2428, Rv1484, Rv0486, Rv1854c, Rv1258c, Rv2752c) %>%
  mutate(num_h_mutations = rowSums(.[3:9]),
         num_strains_h_mutations = case_when(Rv1908c == 1 | Rv2428 == 1 | Rv1484 == 1 | Rv0486 == 1 | Rv1854c == 1 | Rv1258c == 1 | Rv2752c == 1 ~ 1,
                                               T ~ 0)) %>% 
  filter(H==1) %>% 
  freq_table(num_strains_h_mutations)

# kanamycin
km_freq <- data %>% 
  select(strain_id, Km, MTB000019, Rv2461c, Rv3197A) %>%
  mutate(num_km_mutations = rowSums(.[3:5]),
         num_strains_km_mutations = case_when(MTB000019 == 1 | Rv2461c == 1 | Rv3197A == 1 ~ 1,
                                             T ~ 0)) %>% 
  filter(Km==1) %>% 
  freq_table(num_strains_km_mutations)

# fluoroquinolones
fq_freq <- data %>% 
  select(strain_id, Fq, Rv0006, Rv0005) %>%
  mutate(num_fq_mutations = rowSums(.[3:4]),
         num_strains_fq_mutations = case_when(Rv0006 == 1 | Rv0005 == 1 ~ 1,
                                              T ~ 0)) %>% 
  filter(Fq==1) %>% 
  freq_table(num_strains_fq_mutations)

# linezolid
lzd_freq <- data %>% 
  select(strain_id, Lzd, Rv0701, MTB000020) %>%
  mutate(num_lzd_mutations = rowSums(.[3:4]),
         num_strains_lzd_mutations = case_when(Rv0701 == 1 | MTB000020 == 1 ~ 1,
                                               T ~ 0)) %>% 
  filter(Lzd==1) %>% 
  freq_table(num_strains_lzd_mutations)

# moxifloxacin
# mfx_freq <- data %>% 
#   select(strain_id, Mfx, Rv0006, Rv0005) %>%
#   mutate(num_mfx_mutations = rowSums(.[3:4]),
#          num_strains_mfx_mutations = case_when(Rv0006 == 1 | Rv0005 == 1 ~ 1,
#                                                T ~ 0)) %>% 
#   filter(Mfx==1) %>% 
#   freq_table(num_strains_mfx_mutations)

# pyrazinamide 
pyr_freq <- data %>% 
  select(strain_id, Z, Rv3596c, Rv3601c, Rv2043c, Rv1918c, Rv1258c, Rv3236c) %>%
  mutate(num_z_mutations = rowSums(.[3:8]),
         num_strains_z_mutations = case_when(Rv3596c == 1 | Rv3601c == 1 | Rv2043c == 1 | Rv1918c == 1 | Rv1258c == 1 | Rv3236c == 1  ~ 1,
                                               T ~ 0)) %>%
  filter(Z==1) %>% 
  freq_table(num_strains_z_mutations)

# rifampicin
rif_freq <- data %>% 
  select(strain_id, R, Rv0667, Rv0957, Rv3457c, Rv2752c) %>%
  mutate(num_r_mutations = rowSums(.[3:6]),
         num_strains_r_mutations = case_when(Rv0667 == 1 | Rv0957 == 1 | Rv3457c == 1 | Rv2752c == 1  ~ 1,
                                             T ~ 0)) %>% 
  filter(R==1) %>% 
  freq_table(num_strains_r_mutations)

# streptomycin
str_freq <- data %>% 
  select(strain_id, S, Rv3919c, Rv0682, MTB000019, Rv1258c, Rv3862c, Rv3197A) %>%
  mutate(num_s_mutations = rowSums(.[3:6]),
         num_strains_s_mutations = case_when(Rv3919c == 1 | Rv0682 == 1 | MTB000019 == 1 | Rv1258c == 1  | Rv3862c == 1 | Rv3197A == 1 ~ 1,
                                             T ~ 0)) %>% 
  filter(S==1) %>% 
  freq_table(num_strains_s_mutations)


all_freq <- rbind(am_freq,
                  bdq_freq,
                  cm_freq,
                  cfz_freq,
                  cs_freq,
                  eth_freq,
                  eto_freq,
                  fq_freq,
                  inh_freq,
                  km_freq,
                  lzd_freq,
                  #mfx_freq,
                  pyr_freq,
                  rif_freq,
                  str_freq)

# write.csv(all_freq, 'output/treatment_fixed_mutation/frequency_table.csv')
write.csv(all_freq, 'frequency_table.csv')

# table(data$resistance_type)
# MDR  mono-resistant pan-susceptible       resistant             XDR 
# 158              35             174              54             129 

sort(table(data$geography_recode))
# Azerbaijan       Norway        Spain       Sweden       Russia       Uganda      Romania        China       Ghana        Malawi      Georgia South Africa      Moldova      Belarus 
#         1            1            1            1            4            5            6           14           20           31           75          106          114          171 

# 3. create frequency table stratified by baseline resistance (pan-susceptible vs not)
data <- readRDS('tx_fixedmut_tbprof.Rds')
# data <- data %>% 
#   mutate(pansusceptible = case_when(resistance_type == 'pan-susceptible' ~ 1, T ~ 0))
data <- data %>% 
  mutate(pansusceptible = case_when(resistance_type == 'pan-susceptible' ~ 1, T ~ 0))
data <- data %>% 
  mutate(mdr = case_when(resistance_type == 'MDR' ~ 1, T ~ 0))
data <- data %>% 
  mutate(xdr = case_when(resistance_type == 'XDR' ~ 1, T ~ 0))


# amikacin
am_freq_strat <- data %>% 
  select(strain_id, pansusceptible, Am, Rv0529, Rv2416c, Rv3106, MTB000019, Rv3805c, Rv3862c, Rv3197A) %>%
  mutate(num_am_mutations = rowSums(.[3:9]),
         num_strains_am_mutations = case_when(Rv0529 == 1 | Rv2416c == 1 | Rv3106 == 1 | MTB000019 == 1 | Rv3805c == 1 | Rv3862c == 1 | Rv3197A == 1 ~ 1,
                                              T ~ 0)) %>% 
  filter(Am==1) %>% 
  freq_table(pansusceptible, num_strains_am_mutations) 

# bedaquiline
bdq_freq_strat <- data %>% 
  select(strain_id, pansusceptible, Bdq, Rv1305, Rv2535c, Rv1979c, Rv0678, Rv0676c, Rv0677c) %>%
  mutate(num_bdq_mutations = rowSums(.[3:8]),
         num_strains_bdq_mutations = case_when(Rv1305 == 1 | Rv2535c == 1 | Rv1979c == 1 | Rv0678 == 1 | Rv0676c == 1 | Rv0677c == 1 ~ 1,
                                               T ~ 0)) %>% 
  filter(Bdq==1) %>% 
  freq_table(pansusceptible, num_strains_bdq_mutations)

# capreomycin
cm_freq_strat <- data %>% 
  select(strain_id, pansusceptible, Cm, Rv0529, Rv3106, MTB000019, Rv0668, Rv3805c, Rv3862c) %>%
  mutate(num_cm_mutations = rowSums(.[3:8]),
         num_strains_cm_mutations = case_when(Rv0529 == 1 | Rv3106 == 1 | MTB000019 == 1 | Rv0668 == 1 | Rv3805c == 1 | Rv3862c == 1 ~ 1,
                                              T ~ 0)) %>% 
  filter(Cm==1) %>% 
  freq_table(pansusceptible, num_strains_cm_mutations)

# clofazimine
cfz_freq_strat <- data %>% 
  select(strain_id, pansusceptible, Cfz, Rv2535c, Rv1979c, Rv0678, Rv0676c, Rv0677c) %>%
  mutate(num_cfz_mutations = rowSums(.[3:7]),
         num_strains_cfz_mutations = case_when(Rv2535c == 1 | Rv1979c == 1 | Rv0678 == 1 | Rv0676c == 1 | Rv0677c == 1 ~ 1,
                                               T ~ 0)) %>% 
  filter(Cfz==1) %>% 
  freq_table(pansusceptible, num_strains_cfz_mutations)

# cycloserine
cs_freq_strat <- data %>% 
  select(strain_id, pansusceptible, Cs, Rv2780, Rv3423c) %>%
  mutate(num_cs_mutations = rowSums(.[3:4]),
         num_strains_cs_mutations = case_when(Rv2780 == 1 | Rv3423c == 1 ~ 1,
                                              T ~ 0)) %>% 
  filter(Cs==1) %>% 
  freq_table(pansusceptible, num_strains_cs_mutations)

# delamanid
# dld_freq <- data %>% 
#   select(strain_id, Dld, Rv3547, Rv3261, Rv3262, Rv1173, Rv2983, Rv0407) %>%
#   mutate(num_dld_mutations = rowSums(.[3:8]),
#          num_strains_dld_mutations = case_when(Rv3547 == 1 | Rv3261 == 1 | Rv3262 == 1 | Rv1173 == 1 | Rv2983 == 1 | Rv0407 == 1 ~ 1,
#                                               T ~ 0)) %>% 
#   filter(Dld==1) %>% 
#   freq_table(num_strains_dld_mutations)

# ethambutol
eth_freq_strat <- data %>% 
  select(strain_id, pansusceptible, E, Rv3795, Rv3793, Rv3794, Rv3792, Rv1267c, Rv3806c) %>%
  mutate(num_e_mutations = rowSums(.[3:8]),
         num_strains_e_mutations = case_when(Rv3795 == 1 | Rv3793 == 1 | Rv3794 == 1 | Rv3792 == 1 | Rv1267c == 1 | Rv3806c == 1 ~ 1,
                                             T ~ 0)) %>% 
  filter(E==1) %>% 
  freq_table(pansusceptible, num_strains_e_mutations)

# ethionamide
eto_freq_strat <- data %>% 
  select(strain_id, pansusceptible, Eto, Rv3083, Rv1484, Rv0486, Rv1854c, Rv3854c, Rv3855) %>%
  mutate(num_eto_mutations = rowSums(.[3:8]),
         num_strains_eto_mutations = case_when(Rv3083 == 1 | Rv1484 == 1 | Rv0486 == 1 | Rv1854c == 1 | Rv3854c == 1 | Rv3855 == 1 ~ 1,
                                               T ~ 0)) %>% 
  filter(Eto==1) %>% 
  freq_table(pansusceptible, num_strains_eto_mutations)

# isoniazid
inh_freq_strat <- data %>% 
  select(strain_id, pansusceptible, H, Rv1908c, Rv2428, Rv1484, Rv0486, Rv1854c, Rv1258c, Rv2752c) %>%
  mutate(num_h_mutations = rowSums(.[3:9]),
         num_strains_h_mutations = case_when(Rv1908c == 1 | Rv2428 == 1 | Rv1484 == 1 | Rv0486 == 1 | Rv1854c == 1 | Rv1258c == 1 | Rv2752c == 1 ~ 1,
                                             T ~ 0)) %>% 
  filter(H==1) %>% 
  freq_table(pansusceptible, num_strains_h_mutations)

# kanamycin
km_freq_strat <- data %>% 
  select(strain_id, pansusceptible, Km, MTB000019, Rv2461c, Rv3197A) %>%
  mutate(num_km_mutations = rowSums(.[3:5]),
         num_strains_km_mutations = case_when(MTB000019 == 1 | Rv2461c == 1 | Rv3197A == 1 ~ 1,
                                              T ~ 0)) %>% 
  filter(Km==1) %>% 
  freq_table(pansusceptible, num_strains_km_mutations)

# fluoroquinolones
fq_freq_strat <- data %>% 
  select(strain_id, pansusceptible, Fq, Rv0006, Rv0005) %>%
  mutate(num_fq_mutations = rowSums(.[3:4]),
         num_strains_fq_mutations = case_when(Rv0006 == 1 | Rv0005 == 1 ~ 1,
                                              T ~ 0)) %>% 
  filter(Fq==1) %>% 
  freq_table(pansusceptible, num_strains_fq_mutations)

# linezolid
lzd_freq_strat <- data %>% 
  select(strain_id, pansusceptible, Lzd, Rv0701, MTB000020) %>%
  mutate(num_lzd_mutations = rowSums(.[3:4]),
         num_strains_lzd_mutations = case_when(Rv0701 == 1 | MTB000020 == 1 ~ 1,
                                               T ~ 0)) %>% 
  filter(Lzd==1) %>% 
  freq_table(pansusceptible, num_strains_lzd_mutations)

# pyrazinamide 
pyr_freq_strat <- data %>% 
  select(strain_id, pansusceptible, Z, Rv3596c, Rv3601c, Rv2043c, Rv1918c, Rv1258c, Rv3236c) %>%
  mutate(num_z_mutations = rowSums(.[3:8]),
         num_strains_z_mutations = case_when(Rv3596c == 1 | Rv3601c == 1 | Rv2043c == 1 | Rv1918c == 1 | Rv1258c == 1 | Rv3236c == 1  ~ 1,
                                             T ~ 0)) %>%
  filter(Z==1) %>% 
  freq_table(pansusceptible, num_strains_z_mutations)

# rifampicin
rif_freq_strat <- data %>% 
  select(strain_id, pansusceptible, R, Rv0667, Rv0957, Rv3457c, Rv2752c) %>%
  mutate(num_r_mutations = rowSums(.[3:6]),
         num_strains_r_mutations = case_when(Rv0667 == 1 | Rv0957 == 1 | Rv3457c == 1 | Rv2752c == 1  ~ 1,
                                             T ~ 0)) %>% 
  filter(R==1) %>% 
  freq_table(pansusceptible, num_strains_r_mutations)

# streptomycin
str_freq_strat <- data %>% 
  select(strain_id, pansusceptible, S, Rv3919c, Rv0682, MTB000019, Rv1258c, Rv3862c, Rv3197A) %>%
  mutate(num_s_mutations = rowSums(.[3:6]),
         num_strains_s_mutations = case_when(Rv3919c == 1 | Rv0682 == 1 | MTB000019 == 1 | Rv1258c == 1  | Rv3862c == 1 | Rv3197A == 1 ~ 1,
                                             T ~ 0)) %>% 
  filter(S==1) %>% 
  freq_table(pansusceptible, num_strains_s_mutations)

all_freq_strat <- rbind(am_freq_strat,
                        bdq_freq_strat,
                        cm_freq_strat,
                        cfz_freq_strat,
                        cs_freq_strat,
                        eth_freq_strat,
                        eto_freq_strat,
                        fq_freq_strat,
                        inh_freq_strat,
                        km_freq_strat,
                        lzd_freq_strat,
                        # mfx_freq_strat,
                        pyr_freq_strat,
                        rif_freq_strat,
                        str_freq_strat)

write.csv(all_freq_strat, 'frequency_table_strat.panS.csv')



# amikacin
am_freq_strat <- data %>% 
  select(strain_id, mdr, Am, Rv0529, Rv2416c, Rv3106, MTB000019, Rv3805c, Rv3862c, Rv3197A) %>%
  mutate(num_am_mutations = rowSums(.[3:9]),
         num_strains_am_mutations = case_when(Rv0529 == 1 | Rv2416c == 1 | Rv3106 == 1 | MTB000019 == 1 | Rv3805c == 1 | Rv3862c == 1 | Rv3197A == 1 ~ 1,
                                              T ~ 0)) %>% 
  filter(Am==1) %>% 
  freq_table(mdr, num_strains_am_mutations) 

# bedaquiline
bdq_freq_strat <- data %>% 
  select(strain_id, mdr, Bdq, Rv1305, Rv2535c, Rv1979c, Rv0678, Rv0676c, Rv0677c) %>%
  mutate(num_bdq_mutations = rowSums(.[3:8]),
         num_strains_bdq_mutations = case_when(Rv1305 == 1 | Rv2535c == 1 | Rv1979c == 1 | Rv0678 == 1 | Rv0676c == 1 | Rv0677c == 1 ~ 1,
                                               T ~ 0)) %>% 
  filter(Bdq==1) %>% 
  freq_table(mdr, num_strains_bdq_mutations)

# capreomycin
cm_freq_strat <- data %>% 
  select(strain_id, mdr, Cm, Rv0529, Rv3106, MTB000019, Rv0668, Rv3805c, Rv3862c) %>%
  mutate(num_cm_mutations = rowSums(.[3:8]),
         num_strains_cm_mutations = case_when(Rv0529 == 1 | Rv3106 == 1 | MTB000019 == 1 | Rv0668 == 1 | Rv3805c == 1 | Rv3862c == 1 ~ 1,
                                              T ~ 0)) %>% 
  filter(Cm==1) %>% 
  freq_table(mdr, num_strains_cm_mutations)

# clofazimine
cfz_freq_strat <- data %>% 
  select(strain_id, mdr, Cfz, Rv2535c, Rv1979c, Rv0678, Rv0676c, Rv0677c) %>%
  mutate(num_cfz_mutations = rowSums(.[3:7]),
         num_strains_cfz_mutations = case_when(Rv2535c == 1 | Rv1979c == 1 | Rv0678 == 1 | Rv0676c == 1 | Rv0677c == 1 ~ 1,
                                               T ~ 0)) %>% 
  filter(Cfz==1) %>% 
  freq_table(mdr, num_strains_cfz_mutations)

# cycloserine
cs_freq_strat <- data %>% 
  select(strain_id, mdr, Cs, Rv2780, Rv3423c) %>%
  mutate(num_cs_mutations = rowSums(.[3:4]),
         num_strains_cs_mutations = case_when(Rv2780 == 1 | Rv3423c == 1 ~ 1,
                                              T ~ 0)) %>% 
  filter(Cs==1) %>% 
  freq_table(mdr, num_strains_cs_mutations)

# delamanid
# dld_freq <- data %>% 
#   select(strain_id, Dld, Rv3547, Rv3261, Rv3262, Rv1173, Rv2983, Rv0407) %>%
#   mutate(num_dld_mutations = rowSums(.[3:8]),
#          num_strains_dld_mutations = case_when(Rv3547 == 1 | Rv3261 == 1 | Rv3262 == 1 | Rv1173 == 1 | Rv2983 == 1 | Rv0407 == 1 ~ 1,
#                                               T ~ 0)) %>% 
#   filter(Dld==1) %>% 
#   freq_table(num_strains_dld_mutations)

# ethambutol
eth_freq_strat <- data %>% 
  select(strain_id, mdr, E, Rv3795, Rv3793, Rv3794, Rv3792, Rv1267c, Rv3806c) %>%
  mutate(num_e_mutations = rowSums(.[3:8]),
         num_strains_e_mutations = case_when(Rv3795 == 1 | Rv3793 == 1 | Rv3794 == 1 | Rv3792 == 1 | Rv1267c == 1 | Rv3806c == 1 ~ 1,
                                             T ~ 0)) %>% 
  filter(E==1) %>% 
  freq_table(mdr, num_strains_e_mutations)

# ethionamide
eto_freq_strat <- data %>% 
  select(strain_id, mdr, Eto, Rv3083, Rv1484, Rv0486, Rv1854c, Rv3854c, Rv3855) %>%
  mutate(num_eto_mutations = rowSums(.[3:8]),
         num_strains_eto_mutations = case_when(Rv3083 == 1 | Rv1484 == 1 | Rv0486 == 1 | Rv1854c == 1 | Rv3854c == 1 | Rv3855 == 1 ~ 1,
                                               T ~ 0)) %>% 
  filter(Eto==1) %>% 
  freq_table(mdr, num_strains_eto_mutations)

# isoniazid
inh_freq_strat <- data %>% 
  select(strain_id, mdr, H, Rv1908c, Rv2428, Rv1484, Rv0486, Rv1854c, Rv1258c, Rv2752c) %>%
  mutate(num_h_mutations = rowSums(.[3:9]),
         num_strains_h_mutations = case_when(Rv1908c == 1 | Rv2428 == 1 | Rv1484 == 1 | Rv0486 == 1 | Rv1854c == 1 | Rv1258c == 1 | Rv2752c == 1 ~ 1,
                                             T ~ 0)) %>% 
  filter(H==1) %>% 
  freq_table(mdr, num_strains_h_mutations)

# kanamycin
km_freq_strat <- data %>% 
  select(strain_id, mdr, Km, MTB000019, Rv2461c, Rv3197A) %>%
  mutate(num_km_mutations = rowSums(.[3:5]),
         num_strains_km_mutations = case_when(MTB000019 == 1 | Rv2461c == 1 | Rv3197A == 1 ~ 1,
                                              T ~ 0)) %>% 
  filter(Km==1) %>% 
  freq_table(mdr, num_strains_km_mutations)

# fluoroquinolones
fq_freq_strat <- data %>% 
  select(strain_id, mdr, Fq, Rv0006, Rv0005) %>%
  mutate(num_fq_mutations = rowSums(.[3:4]),
         num_strains_fq_mutations = case_when(Rv0006 == 1 | Rv0005 == 1 ~ 1,
                                              T ~ 0)) %>% 
  filter(Fq==1) %>% 
  freq_table(mdr, num_strains_fq_mutations)

# linezolid
lzd_freq_strat <- data %>% 
  select(strain_id, mdr, Lzd, Rv0701, MTB000020) %>%
  mutate(num_lzd_mutations = rowSums(.[3:4]),
         num_strains_lzd_mutations = case_when(Rv0701 == 1 | MTB000020 == 1 ~ 1,
                                               T ~ 0)) %>% 
  filter(Lzd==1) %>% 
  freq_table(mdr, num_strains_lzd_mutations)

# moxifloxacin
# mfx_freq_strat <- data %>% 
#   select(strain_id, mdr, Mfx, Rv0006, Rv0005) %>%
#   mutate(num_mfx_mutations = rowSums(.[3:4]),
#          num_strains_mfx_mutations = case_when(Rv0006 == 1 | Rv0005 == 1 ~ 1,
#                                                T ~ 0)) %>% 
#   filter(Mfx==1) %>% 
#   freq_table(mdr, num_strains_mfx_mutations)

# pyrazinamide 
pyr_freq_strat <- data %>% 
  select(strain_id, mdr, Z, Rv3596c, Rv3601c, Rv2043c, Rv1918c, Rv1258c, Rv3236c) %>%
  mutate(num_z_mutations = rowSums(.[3:8]),
         num_strains_z_mutations = case_when(Rv3596c == 1 | Rv3601c == 1 | Rv2043c == 1 | Rv1918c == 1 | Rv1258c == 1 | Rv3236c == 1  ~ 1,
                                             T ~ 0)) %>%
  filter(Z==1) %>% 
  freq_table(mdr, num_strains_z_mutations)

# rifampicin
rif_freq_strat <- data %>% 
  select(strain_id, mdr, R, Rv0667, Rv0957, Rv3457c, Rv2752c) %>%
  mutate(num_r_mutations = rowSums(.[3:6]),
         num_strains_r_mutations = case_when(Rv0667 == 1 | Rv0957 == 1 | Rv3457c == 1 | Rv2752c == 1  ~ 1,
                                             T ~ 0)) %>% 
  filter(R==1) %>% 
  freq_table(mdr, num_strains_r_mutations)

# streptomycin
str_freq_strat <- data %>% 
  select(strain_id, mdr, S, Rv3919c, Rv0682, MTB000019, Rv1258c, Rv3862c, Rv3197A) %>%
  mutate(num_s_mutations = rowSums(.[3:6]),
         num_strains_s_mutations = case_when(Rv3919c == 1 | Rv0682 == 1 | MTB000019 == 1 | Rv1258c == 1  | Rv3862c == 1 | Rv3197A == 1 ~ 1,
                                             T ~ 0)) %>% 
  filter(S==1) %>% 
  freq_table(mdr, num_strains_s_mutations)

all_freq_strat <- rbind(am_freq_strat,
                        bdq_freq_strat,
                        cm_freq_strat,
                        cfz_freq_strat,
                        cs_freq_strat,
                        eth_freq_strat,
                        eto_freq_strat,
                        fq_freq_strat,
                        inh_freq_strat,
                        km_freq_strat,
                        lzd_freq_strat,
                        # mfx_freq_strat,
                        pyr_freq_strat,
                        rif_freq_strat,
                        str_freq_strat)


write.csv(all_freq_strat, 'frequency_table_strat.MDR.csv')



# amikacin
am_freq_strat <- data %>% 
  select(strain_id, xdr, Am, Rv0529, Rv2416c, Rv3106, MTB000019, Rv3805c, Rv3862c, Rv3197A) %>%
  mutate(num_am_mutations = rowSums(.[3:9]),
         num_strains_am_mutations = case_when(Rv0529 == 1 | Rv2416c == 1 | Rv3106 == 1 | MTB000019 == 1 | Rv3805c == 1 | Rv3862c == 1 | Rv3197A == 1 ~ 1,
                                              T ~ 0)) %>% 
  filter(Am==1) %>% 
  freq_table(xdr, num_strains_am_mutations) 

# bedaquiline
bdq_freq_strat <- data %>% 
  select(strain_id, xdr, Bdq, Rv1305, Rv2535c, Rv1979c, Rv0678, Rv0676c, Rv0677c) %>%
  mutate(num_bdq_mutations = rowSums(.[3:8]),
         num_strains_bdq_mutations = case_when(Rv1305 == 1 | Rv2535c == 1 | Rv1979c == 1 | Rv0678 == 1 | Rv0676c == 1 | Rv0677c == 1 ~ 1,
                                               T ~ 0)) %>% 
  filter(Bdq==1) %>% 
  freq_table(xdr, num_strains_bdq_mutations)

# capreomycin
cm_freq_strat <- data %>% 
  select(strain_id, xdr, Cm, Rv0529, Rv3106, MTB000019, Rv0668, Rv3805c, Rv3862c) %>%
  mutate(num_cm_mutations = rowSums(.[3:8]),
         num_strains_cm_mutations = case_when(Rv0529 == 1 | Rv3106 == 1 | MTB000019 == 1 | Rv0668 == 1 | Rv3805c == 1 | Rv3862c == 1 ~ 1,
                                              T ~ 0)) %>% 
  filter(Cm==1) %>% 
  freq_table(xdr, num_strains_cm_mutations)

# clofazimine
cfz_freq_strat <- data %>% 
  select(strain_id, xdr, Cfz, Rv2535c, Rv1979c, Rv0678, Rv0676c, Rv0677c) %>%
  mutate(num_cfz_mutations = rowSums(.[3:7]),
         num_strains_cfz_mutations = case_when(Rv2535c == 1 | Rv1979c == 1 | Rv0678 == 1 | Rv0676c == 1 | Rv0677c == 1 ~ 1,
                                               T ~ 0)) %>% 
  filter(Cfz==1) %>% 
  freq_table(xdr, num_strains_cfz_mutations)

# cycloserine
cs_freq_strat <- data %>% 
  select(strain_id, xdr, Cs, Rv2780, Rv3423c) %>%
  mutate(num_cs_mutations = rowSums(.[3:4]),
         num_strains_cs_mutations = case_when(Rv2780 == 1 | Rv3423c == 1 ~ 1,
                                              T ~ 0)) %>% 
  filter(Cs==1) %>% 
  freq_table(xdr, num_strains_cs_mutations)

# delamanid
# dld_freq <- data %>% 
#   select(strain_id, Dld, Rv3547, Rv3261, Rv3262, Rv1173, Rv2983, Rv0407) %>%
#   mutate(num_dld_mutations = rowSums(.[3:8]),
#          num_strains_dld_mutations = case_when(Rv3547 == 1 | Rv3261 == 1 | Rv3262 == 1 | Rv1173 == 1 | Rv2983 == 1 | Rv0407 == 1 ~ 1,
#                                               T ~ 0)) %>% 
#   filter(Dld==1) %>% 
#   freq_table(num_strains_dld_mutations)

# ethambutol
eth_freq_strat <- data %>% 
  select(strain_id, xdr, E, Rv3795, Rv3793, Rv3794, Rv3792, Rv1267c, Rv3806c) %>%
  mutate(num_e_mutations = rowSums(.[3:8]),
         num_strains_e_mutations = case_when(Rv3795 == 1 | Rv3793 == 1 | Rv3794 == 1 | Rv3792 == 1 | Rv1267c == 1 | Rv3806c == 1 ~ 1,
                                             T ~ 0)) %>% 
  filter(E==1) %>% 
  freq_table(xdr, num_strains_e_mutations)

# ethionamide
eto_freq_strat <- data %>% 
  select(strain_id, xdr, Eto, Rv3083, Rv1484, Rv0486, Rv1854c, Rv3854c, Rv3855) %>%
  mutate(num_eto_mutations = rowSums(.[3:8]),
         num_strains_eto_mutations = case_when(Rv3083 == 1 | Rv1484 == 1 | Rv0486 == 1 | Rv1854c == 1 | Rv3854c == 1 | Rv3855 == 1 ~ 1,
                                               T ~ 0)) %>% 
  filter(Eto==1) %>% 
  freq_table(xdr, num_strains_eto_mutations)

# isoniazid
inh_freq_strat <- data %>% 
  select(strain_id, xdr, H, Rv1908c, Rv2428, Rv1484, Rv0486, Rv1854c, Rv1258c, Rv2752c) %>%
  mutate(num_h_mutations = rowSums(.[3:9]),
         num_strains_h_mutations = case_when(Rv1908c == 1 | Rv2428 == 1 | Rv1484 == 1 | Rv0486 == 1 | Rv1854c == 1 | Rv1258c == 1 | Rv2752c == 1 ~ 1,
                                             T ~ 0)) %>% 
  filter(H==1) %>% 
  freq_table(xdr, num_strains_h_mutations)

# kanamycin
km_freq_strat <- data %>% 
  select(strain_id, xdr, Km, MTB000019, Rv2461c, Rv3197A) %>%
  mutate(num_km_mutations = rowSums(.[3:5]),
         num_strains_km_mutations = case_when(MTB000019 == 1 | Rv2461c == 1 | Rv3197A == 1 ~ 1,
                                              T ~ 0)) %>% 
  filter(Km==1) %>% 
  freq_table(xdr, num_strains_km_mutations)

# fluoroquinolones
fq_freq_strat <- data %>% 
  select(strain_id, xdr, Fq, Rv0006, Rv0005) %>%
  mutate(num_fq_mutations = rowSums(.[3:4]),
         num_strains_fq_mutations = case_when(Rv0006 == 1 | Rv0005 == 1 ~ 1,
                                              T ~ 0)) %>% 
  filter(Fq==1) %>% 
  freq_table(xdr, num_strains_fq_mutations)

# linezolid
lzd_freq_strat <- data %>% 
  select(strain_id, xdr, Lzd, Rv0701, MTB000020) %>%
  mutate(num_lzd_mutations = rowSums(.[3:4]),
         num_strains_lzd_mutations = case_when(Rv0701 == 1 | MTB000020 == 1 ~ 1,
                                               T ~ 0)) %>% 
  filter(Lzd==1) %>% 
  freq_table(xdr, num_strains_lzd_mutations)

# moxifloxacin
# mfx_freq_strat <- data %>% 
#   select(strain_id, xdr, Mfx, Rv0006, Rv0005) %>%
#   mutate(num_mfx_mutations = rowSums(.[3:4]),
#          num_strains_mfx_mutations = case_when(Rv0006 == 1 | Rv0005 == 1 ~ 1,
#                                                T ~ 0)) %>% 
#   filter(Mfx==1) %>% 
#   freq_table(xdr, num_strains_mfx_mutations)

# pyrazinamide 
pyr_freq_strat <- data %>% 
  select(strain_id, xdr, Z, Rv3596c, Rv3601c, Rv2043c, Rv1918c, Rv1258c, Rv3236c) %>%
  mutate(num_z_mutations = rowSums(.[3:8]),
         num_strains_z_mutations = case_when(Rv3596c == 1 | Rv3601c == 1 | Rv2043c == 1 | Rv1918c == 1 | Rv1258c == 1 | Rv3236c == 1  ~ 1,
                                             T ~ 0)) %>%
  filter(Z==1) %>% 
  freq_table(xdr, num_strains_z_mutations)

# rifampicin
rif_freq_strat <- data %>% 
  select(strain_id, xdr, R, Rv0667, Rv0957, Rv3457c, Rv2752c) %>%
  mutate(num_r_mutations = rowSums(.[3:6]),
         num_strains_r_mutations = case_when(Rv0667 == 1 | Rv0957 == 1 | Rv3457c == 1 | Rv2752c == 1  ~ 1,
                                             T ~ 0)) %>% 
  filter(R==1) %>% 
  freq_table(xdr, num_strains_r_mutations)

# streptomycin
str_freq_strat <- data %>% 
  select(strain_id, xdr, S, Rv3919c, Rv0682, MTB000019, Rv1258c, Rv3862c, Rv3197A) %>%
  mutate(num_s_mutations = rowSums(.[3:6]),
         num_strains_s_mutations = case_when(Rv3919c == 1 | Rv0682 == 1 | MTB000019 == 1 | Rv1258c == 1  | Rv3862c == 1 | Rv3197A == 1 ~ 1,
                                             T ~ 0)) %>% 
  filter(S==1) %>% 
  freq_table(xdr, num_strains_s_mutations)

all_freq_strat <- rbind(am_freq_strat,
                        bdq_freq_strat,
                        cm_freq_strat,
                        cfz_freq_strat,
                        cs_freq_strat,
                        eth_freq_strat,
                        eto_freq_strat,
                        fq_freq_strat,
                        inh_freq_strat,
                        km_freq_strat,
                        lzd_freq_strat,
                        # mfx_freq_strat,
                        pyr_freq_strat,
                        rif_freq_strat,
                        str_freq_strat)


write.csv(all_freq_strat, 'frequency_table_strat.XDR.csv')


