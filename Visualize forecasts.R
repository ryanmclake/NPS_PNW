rmse = function(m, o){
  sqrt(mean((m - o)^2))
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DATA PATH NEEDS TO BE CHANGED to ./data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

data_path <- "./data"

forecast_type <- c("o", "f","o", "f","o", "f","o", "f","o", "f","o", "f",
                   "o", "f","o", "f","o", "f","o", "f","o", "f")
forecast_number <- c("1", "1","2", "2","3", "3","4", "4","5", "5","6", "6",
                     "7", "7","8", "8","9", "9","10", "10","11", "11")
clad_noca_forecasts <- list.files(data_path, pattern = "CLAD_NOCA_forecast_wDA_")%>%
  map(~ readRDS(file.path(data_path, .))) %>%
  data.table::rbindlist(fill = T)

clad_noca_forecasts <- cbind(clad_noca_forecasts, forecast_type, forecast_number)

CLAD_NOCA_observed <- CLAD_NOCA %>%
  select(site_code, event_year, CLAD)%>%
  group_by(event_year)%>%
  summarize(mean_CLAD = mean(CLAD),
            sd_CLAD = sd(CLAD))

CLAD_NOCA_all <- CLAD_NOCA %>%
  select(site_code, event_year, CLAD)

CLAD_NOCA_forecasts_wDA <- clad_noca_forecasts %>%
  ggplot(., aes(x = event_year, y = mean, group = forecast_number)) +
  geom_point(data = CLAD_NOCA_all, aes(x = event_year, y = CLAD, fill = "midnightblue", color = "black"), inherit.aes = F, cex = 1, alpha = 1)+
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), alpha = 0.2, fill = "midnightblue") +
  geom_line(color = "purple4", size = 1, alpha = 0.6)+
  geom_point(data = CLAD_NOCA_observed, aes(x = event_year, y = mean_CLAD), inherit.aes = FALSE, pch = 21, color = "black", fill = "red", cex = 3) +
  theme_bw()+
  labs(title = "A: NOCA CLAD FORECASTS")+
  ylab(expression(paste("log10(CLAD Abundance + 1)")))+
  xlab("")+
  theme(axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        title = element_text(size = 15),legend.position = "none",
        legend.text = element_text(size = 16, color = "black"))


##### UNCERTAINTY PARTITIONING DATA ####

trap_all_IC <- list.files(data_path, pattern = "IC_UNC")%>%
  map(~ readRDS(file.path(data_path, .))) %>%
  data.table::rbindlist(fill = T)

trap_all_IC <- cbind(trap_all_IC, forecast_type, forecast_number)

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
