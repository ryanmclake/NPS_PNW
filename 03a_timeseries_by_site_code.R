macro_join_ts <- macro_join %>%
  group_by(site_code)%>%
  mutate(Amphipoda_lag = lag(Amphipoda),
         Ntoda_lag = lag(Ntoda),
         Ephemeroptera_lag = lag(Ephemeroptera),
         Nmorpha_lag = lag(Nmorpha),
         Veneroida_lag = lag(Veneroida),
         Diptera_lag = lag(Diptera),
         Megaloptera_lag = lag(Megaloptera),
         Trichoptera_lag = lag(Trichoptera),
         Basommatophora_lag = lag(Basommatophora),
         Hemiptera_lag = lag(Hemiptera),
         Platy_lag = lag(Platy),
         Coleoptera_lag = lag(Coleoptera),
         Hirudinida_lag = lag(Hirudinida),
         Ntopda_lag = lag(Ntopda),
         Plecoptera_lag = lag(Plecoptera),
         Acari_lag = lag(Acari),
         Diplostraca_lag = lag(Diplostraca),
         Isopoda_lag = lag(Isopoda),
         Odonata_lag = lag(Odonata),
         Poridera_lag = lag(Poridera))%>%
  mutate_at(vars(-site_code, -park_code, -event_year),
            funs(imputeTS::na_interpolation(., option = "spline")))

sites <- c(unique(macro_join_ts$site_code))

taxa <- c("Amphipoda",
          "Ntoda",
          "Ephemeroptera",
          "Nmorpha",
          "Veneroida",
          "Diptera",
          "Megaloptera",
          "Trichoptera",
          "Basommatophora",
          "Hemiptera",
          "Platy",
          "Coleoptera",
          "Hirudinida",
          "Ntopda",
          "Plecoptera",
          "Acari",
          "Diplostraca",
          "Isopoda",
          "Odonata",
          "Poridera")


drivers <- data.frame(matrix(NA, nrow = length(sites),ncol = length(taxa)))

for(h in 1:length(taxa)){
for(s in 1:length(sites)){

  ts <- macro_join_ts %>%
    filter(site_code == sites[s])%>%
    select(taxa[h],paste0(taxa[h],"_lag"), delta_temp,`Total P`,`Total N`, Ca, Chlorophyll)%>%
    rename(TP = `Total P`,
           TN = `Total N`,
           macro = taxa[h],
           macro_lag = paste0(taxa[h],"_lag"))

  if(sum(ts$macro != 0)){

  all <- glm(macro ~ macro_lag+
               delta_temp+Chlorophyll+
               TN+TP+Ca,
             na.action = "na.fail",
             data = ts)

  all_models <- dredge(all, extra = c("R^2","adjR^2"),
                       fixed = "macro_lag", rank = "AICc")

  best_model <- as.data.frame(subset(all_models, delta <= 2))

  best_fits <- best_model[,colSums(is.na(best_model))<nrow(best_model)] %>%
    summarise_all(funs(mean))

  best_fits <- best_fits %>% select(-`R^2`, -`adjR^2`, -df, -logLik, -AICc, -delta, -weight)%>%
    cbind(., sites[s])%>%
    reshape2::melt(., id.vars = "sites[s]")%>%
    select(variable)

  drivers[s,h] <- as.data.frame(cbind(paste(best_fits$variable,collapse=", ")))
}
}
}

drivers <- drivers %>% rename(Amphipoda = `X1`,
                              Ntoda = `X2`,
                              Ephemeroptera = `X3`,
                              Nmorpha = `X4`,
                              Veneroida = `X5`,
                              Diptera = `X6`,
                              Megaloptera = `X7`,
                              Trichoptera = `X8`,
                              Basommatophora = `X9`,
                              Hemiptera = `X10`,
                              Platy = `X11`,
                              Coleoptera = `X12`,
                              Hirudinida = `X13`,
                              Ntopda = `X14`,
                              Plecoptera = `X15`,
                              Acari = `X16`,
                              Diplostraca = `X17`,
                              Isopoda = `X18`,
                              Odonata = `X19`,
                              Poridera = `X20`)%>%
  mutate(site_code = as.character(sites))

map_drivers <- left_join(drivers, mapping_zoop, by = "site_code")



cca_df_ts <- cca_df %>%
  group_by(site_code)%>%
  mutate(RAP_lag = lag(RAP),
         MICRO_lag = lag(MICRO),
         CLAD_lag = lag(CLAD),
         COPE_lag = lag(COPE))%>%
  mutate_at(vars(-site_code, -park_code, -event_year),
            funs(imputeTS::na_interpolation(., option = "spline")))

sites <- c(unique(cca_df_ts$site_code))

taxa <- c("CLAD",
          "COPE",
          "MICRO",
          "RAP")

drivers <- data.frame(matrix(NA, nrow = length(sites),ncol = length(taxa)))

for(h in 1:length(taxa)){
  for(s in 1:length(sites)){

    ts <- cca_df_ts %>%
      filter(site_code == sites[s])%>%
      select(taxa[h],paste0(taxa[h],"_lag"), delta_temp,`Total P`,`Total N`, Ca, Chlorophyll)%>%
      rename(TP = `Total P`,
             TN = `Total N`,
             macro = taxa[h],
             macro_lag = paste0(taxa[h],"_lag"))

    if(sum(ts$macro != 0)){

      all <- glm(macro ~ macro_lag+
                   delta_temp+Chlorophyll+
                   TN+TP+Ca,
                 na.action = "na.fail",
                 data = ts)

      all_models <- dredge(all, extra = c("R^2","adjR^2"),
                           fixed = "macro_lag", rank = "AICc")

      best_model <- as.data.frame(subset(all_models, delta <= 2))

      best_fits <- best_model[,colSums(is.na(best_model))<nrow(best_model)] %>%
        summarise_all(funs(mean))

      best_fits <- best_fits %>% select(-`R^2`, -`adjR^2`, -df, -logLik, -AICc, -delta, -weight)%>%
        cbind(., sites[s])%>%
        reshape2::melt(., id.vars = "sites[s]")%>%
        select(variable)

      drivers[s,h] <- as.data.frame(cbind(paste(best_fits$variable,collapse=", ")))
    }
  }
}

drivers <- drivers %>% rename(CLAD = `X1`,
                              COPE = `X2`,
                              MICRO = `X3`,
                              RAP = `X4`)%>%
  mutate(site_code = as.character(sites))

zoop_map_drivers <- left_join(drivers, mapping_zoop, by = "site_code")



