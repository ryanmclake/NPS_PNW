set.seed(42)

ZOOP_all <- env_zoop_data %>%
  group_by(site_code, event_year, park_code)%>%
  summarize_all(funs(mean), na.rm = F)%>%
  arrange(event_year)%>%
  ungroup(.)%>%
  select(event_year, site_code, park_code, COPE, CLAD, MICRO, RAP, Ca, lake_temp, Chlorophyll, fish)%>%
  mutate(site_code = case_when(
    site_code == "LH15" ~ 1,
    site_code == "Allen" ~ 2,
    site_code == "LP19" ~ 3,
    site_code == "Deadwood" ~ 4,
    site_code == "Blue" ~ 5,
    site_code == "Blum" ~ 6,
    site_code == "Silent" ~ 7,
    site_code == "EasyRidge" ~ 8,
    site_code == "East" ~ 9,
    site_code == "Bowan" ~ 10,
    site_code == "Triplet" ~ 11,
    site_code == "Gladys" ~ 12,
    site_code == "Ferry" ~ 13,
    site_code == "Heather" ~ 14,
    site_code == "Crazy" ~ 15,
    site_code == "Milk" ~ 16,
    site_code == "LaCrosse" ~ 17,
    site_code == "Sunup" ~ 18,
    site_code == "Connie" ~ 19))%>%
  melt(., id.vars = c("event_year","site_code","park_code","Ca","lake_temp","Chlorophyll","fish"))

ZOOP_fishless <- ZOOP_all %>% filter(fish == 0)
ZOOP_fish <- ZOOP_all %>% filter(fish == 1)


# GLOBAL ZOOPLANKTON MODEL
#
zoop_taxa <- c("COPE", "CLAD", "MICRO", "RAP")
for(s in 1:length(zoop_taxa)){

  ZOOP_jags <- ZOOP_all %>% filter(variable == zoop_taxa[s])
  sites <- length(levels(as.factor(ZOOP_jags$site_code)))

  jagsdata <- with(ZOOP_jags, list(value = value, lake_temp = lake_temp, Ca = Ca, site_code = site_code,
                              Chlorophyll = Chlorophyll, fish = fish,  N = length(ZOOP_jags), sites = sites))

  init_values <- function(){
    list(alpha = rnorm(4), beta = rnorm(4), omega = rnorm(4), phi = rnorm(4), sigma = runif(1), sigma_a = runif(1))
  }

  params <- c("alpha", "beta", "omega", "phi", "sigma", "sigma_a")

  fit_lm_global <- jags(data = jagsdata, inits = init_values,
                      parameters.to.save = params, model.file = global_lme_jags,
                      n.chains = 3, n.iter = 200000, n.burnin = 2000, n.thin = 10,
                      DIC = F)

  fit_lm_global <- as.mcmc(fit_lm_global)
  gelman.diag(fit_lm_global)

  global_parameters <- fit_lm_global %>%
    spread_draws(sigma, sigma_a, `alpha[1]`, `alpha[2]`, `alpha[3]`, `alpha[4]`,
                 `beta[1]`, `beta[2]`, `beta[3]`, `beta[4]`,
                 `omega[1]`, `omega[2]`, `omega[3]`, `omega[4]`,
                 `phi[1]`, `phi[2]`, `phi[3]`, `phi[4]`) %>%
    rename(ensemble = .iteration)

  saveRDS(global_parameters, paste0("./output/global_parameters",zoop_taxa[s],".rds"))
}

# FISHLESS ZOOPLANKTON MODEL
#
zoop_taxa <- c("COPE", "CLAD", "MICRO", "RAP")
for(s in 1:length(zoop_taxa)){

  ZOOP_jags <- ZOOP_fishless %>% filter(variable == zoop_taxa[s])
  sites <- length(levels(as.factor(ZOOP_jags$site_code)))

  jagsdata <- with(ZOOP_jags, list(value = value, lake_temp = lake_temp, Ca = Ca, site_code = site_code,
                                   Chlorophyll = Chlorophyll,  N = length(ZOOP_jags), sites = sites))

  init_values <- function(){
    list(alpha = rnorm(1), beta = rnorm(1), omega = rnorm(1), phi = rnorm(1), sigma = runif(1), sigma_a = runif(1))
  }

  params <- c("alpha", "beta", "omega", "phi", "sigma", "sigma_a")

  fit_lm_fishless <- jags(data = jagsdata, inits = init_values,
                        parameters.to.save = params, model.file = fishless_lme_jags,
                        n.chains = 3, n.iter = 200000, n.burnin = 2000, n.thin = 10,
                        DIC = F)

  fit_lm_global <- as.mcmc(fit_lm_global)
  gelman.diag(fit_lm_global)
}


