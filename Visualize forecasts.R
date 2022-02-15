rmse = function(m, o){
  sqrt(mean((m - o)^2))
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DATA PATH NEEDS TO BE CHANGED to ./data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

data_path <- "./output"


site_prediction_jags <- list.files(data_path, pattern = "MAIN_MCMC_JAGS_PREDICTION")%>%
  map(~ readRDS(file.path(data_path, .))) %>%
  data.table::rbindlist(fill = T)%>%
  na.omit(.)%>%
  left_join(., env_zoop_data, by = c("site_code","event_year", "park_code", "fish"))

evals <- site_prediction_jags %>%
  group_by(site_code, taxa)%>%
  summarize(RMSE = rmse(mean, value),
            var = mean(var))%>%
  left_join(., env_zoop_data, by = "site_code")

ggplot(evals, aes(taxa, var))+
  geom_boxplot(aes(fill=park_code))

ggplot(evals, aes(var, RMSE, group = taxa))+
  geom_point(aes(color = park_code))+
  facet_wrap(~taxa)

ggplot(site_prediction_jags, aes(mean, value, group = taxa))+
  geom_point(aes(color = site_code))+
  facet_wrap(~taxa)

CLAD_TS <- site_prediction_jags%>%filter(taxa == "CLAD")%>%
           ggplot(.)+
             geom_line(aes(x = event_year, y = mean, group = park_code, color = as.character(fish)), lwd = 1, color = "midnightblue")+
             geom_ribbon(aes(x = event_year, ymin = lower_95, ymax = upper_95, group = park_code), alpha = 0.2, fill = "midnightblue") +
             geom_point(aes(x = event_year, y = value, group = park_code), cex=4, color = "black", pch = 21, bg = "grey70")+
             geom_point(aes(x = event_year, y = mean, fill = as.character(fish)), cex = 4, color = "black", pch = 21)+
             labs(y = expression(paste("log_10(Abundance)")), x = "", title = "CLAD predictions")+
             theme_bw()+
             theme(text = element_text(size=15, color = "black"),
                   axis.text = element_text(size = 15, color = "black"),
                   axis.text.x = element_text(angle = 45, hjust = 1))+
             facet_wrap(~site_code)

site_prediction_jags %>% filter(taxa == "COPE")%>%
  ggplot(.)+
  geom_line(aes(x = event_year, y = mean, group = park_code, color = as.character(fish)), lwd = 1, color = "midnightblue")+
  geom_ribbon(aes(x = event_year, ymin = lower_95, ymax = upper_95, group = park_code), alpha = 0.2, fill = "midnightblue") +
  geom_point(aes(x = event_year, y = value, group = park_code), cex=4, color = "black", pch = 21, bg = "grey70")+
  geom_point(aes(x = event_year, y = mean, fill = as.character(fish)), cex = 4, color = "black", pch = 21)+
  labs(y = expression(paste("log_10(Abundance)")), x = "", title = "COPE predictions")+
  theme_bw()+
  theme(text = element_text(size=15, color = "black"),
        axis.text = element_text(size = 15, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~site_code)

site_prediction_jags %>% filter(taxa == "RAP")%>%
  ggplot(.)+
  geom_line(aes(x = event_year, y = mean, group = park_code, color = as.character(fish)), lwd = 1, color = "midnightblue")+
  geom_ribbon(aes(x = event_year, ymin = lower_95, ymax = upper_95, group = park_code), alpha = 0.2, fill = "midnightblue") +
  geom_point(aes(x = event_year, y = value, group = park_code), cex=4, color = "black", pch = 21, bg = "grey70")+
  geom_point(aes(x = event_year, y = mean, fill = as.character(fish)), cex = 4, color = "black", pch = 21)+
  labs(y = expression(paste("log_10(Abundance)")), x = "", title = "RAP predictions")+
  theme_bw()+
  theme(text = element_text(size=15, color = "black"),
        axis.text = element_text(size = 15, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~site_code)

site_prediction_jags %>% filter(taxa == "MICRO")%>%
  ggplot(.)+
  geom_line(aes(x = event_year, y = mean, group = park_code, color = as.character(fish)), lwd = 1, color = "midnightblue")+
  geom_ribbon(aes(x = event_year, ymin = lower_95, ymax = upper_95, group = park_code), alpha = 0.2, fill = "midnightblue") +
  geom_point(aes(x = event_year, y = value, group = park_code), cex=4, color = "black", pch = 21, bg = "grey70")+
  geom_point(aes(x = event_year, y = mean, fill = as.character(fish)), cex = 4, color = "black", pch = 21)+
  labs(y = expression(paste("log_10(Abundance)")), x = "", title = "MICRO predictions")+
  theme_bw()+
  theme(text = element_text(size=15, color = "black"),
        axis.text = element_text(size = 15, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~site_code)


##### UNCERTAINTY PARTITIONING DATA ####

TEMP <- list.files(data_path, pattern = "TEMPERATURE_")%>%
  map(~ readRDS(file.path(data_path, .))) %>%
  data.table::rbindlist(fill = T)

trap_all_SWE <- list.files(data_path, pattern = "SWE_UNC")%>%
  map(~ readRDS(file.path(data_path, .))) %>%
  data.table::rbindlist(fill = T)

trap_all_SWE <- cbind(trap_all_SWE, forecast_type, forecast_number)

trap_all_Parameter <- list.files(data_path, pattern = "Parameter_UNC")%>%
  map(~ readRDS(file.path(data_path, .))) %>%
  data.table::rbindlist(fill = T)

trap_all_Parameter <- cbind(trap_all_Parameter, forecast_type, forecast_number)

trap_all_PROCESS <- list.files(data_path, pattern = "PROCESS_UNC")%>%
  map(~ readRDS(file.path(data_path, .))) %>%
  data.table::rbindlist(fill = T)

trap_all_PROCESS <- cbind(trap_all_PROCESS, forecast_type, forecast_number)

partition <- cbind(trap_all_IC[,c(1,8,9,7)], trap_all_PROCESS[,7], trap_all_SWE[,7], trap_all_Parameter[,7])
names(partition) <- c("event_year", "forecast_type", "forecast_number", "var_ic", "var_pro", "var_dri", "var_par")

all_partitioned <- partition %>%
  mutate(sum_var = var_ic+var_pro+var_dri+var_par)%>%
  mutate(`Initial condition` = var_ic/sum_var)%>%
  mutate(`Driver data` = var_dri/sum_var)%>%
  mutate(`Model parameter` = var_par/sum_var)%>%
  mutate(`Model process` = var_pro/sum_var)%>%
  mutate(sum_check = `Model process`+`Model parameter`+`Driver data`+`Initial condition`)%>%
  select(event_year,forecast_type, forecast_number, `Initial condition`, `Driver data`, `Model parameter`, `Model process`, sum_check)%>%
  ungroup(.)%>%
  filter(forecast_type == "f")

all_partitioned_melt <- all_partitioned%>%
  select(-sum_check)%>%
  melt(., id = c("event_year","forecast_type","forecast_number"))

c <- all_partitioned_melt%>%
  group_by(variable)%>%
  ggplot(., aes(x = event_year, y = value, group=variable, fill = variable))+
  geom_area(aes(fill = variable))+
  scale_fill_manual(values = c("#E69F00", "#D55E00", "#CC79A7", "#56B4E9"))+
  theme_bw()+
  labs(title = "A: Partitioned uncertainty of each year forecast")+
  ylab("Proportion of total variance")+
  xlab("")+
  theme(axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        title = element_text(size = 15), legend.position = "top",
        legend.text = element_text(size = 10, color = "black"))
c

trap_all_SWE_parm <- list.files(data_path, pattern = "SWE_PARM")%>%
  map(~ readRDS(file.path(data_path, .))) %>%
  data.table::rbindlist(fill = T)

trap_all_SWE_parm <- cbind(trap_all_SWE_parm, forecast_type, forecast_number)


trap_all_AR_parm <- list.files(data_path, pattern = "AR_PARM")%>%
  map(~ readRDS(file.path(data_path, .))) %>%
  data.table::rbindlist(fill = T)

trap_all_AR_parm <- cbind(trap_all_AR_parm, forecast_type, forecast_number)

partition <- cbind(trap_all_SWE_parm[,c(1,8,9,7)], trap_all_AR_parm[,7])
names(partition) <- c("event_year", "forecast_type", "forecast_number", "var_swe", "var_ar")

all_partitioned <- partition %>%
  mutate(sum_var = var_swe+var_ar)%>%
  mutate(`SWE Sensitivity` = var_swe/sum_var)%>%
  mutate(`AR Sensitivity` = var_ar/sum_var)%>%
  mutate(sum_check = `SWE Sensitivity`+`AR Sensitivity`)%>%
  select(event_year,forecast_type, forecast_number, `SWE Sensitivity`, `AR Sensitivity`, sum_check)%>%
  ungroup(.)%>%
  filter(forecast_type == "f")

all_partitioned_melt <- all_partitioned%>%
  select(-sum_check)%>%
  melt(., id = c("event_year","forecast_type","forecast_number"))


d <- all_partitioned_melt%>%
  group_by(variable)%>%
  ggplot(., aes(x = event_year, y = value, group=variable, fill = variable))+
  geom_area(aes(fill = variable))+
  scale_fill_manual(values = c("#E69F00", "#D55E00", "#CC79A7", "#56B4E9"))+
  theme_bw()+
  labs(title = "A: Partitioned uncertainty of each year forecast")+
  ylab("Proportion of total variance")+
  xlab("")+
  theme(axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        title = element_text(size = 15), legend.position = "top",
        legend.text = element_text(size = 10, color = "black"))
d
