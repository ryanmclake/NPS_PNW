# Analyze Time Series of Temperature in the National Park Lakes
se <- function(x) sd(x) / sqrt(length(x))
normalize <- function(x){(x-min(x))/(max(x)-min(x))}
standard_error <- function(x) {sd(x) / sqrt(length(x))}
percent_change <- function(x) {((x - lead(x))/(x))*100}
zscore <- function(x) {(x - mean(x))/sd(x)}
slope <- function(x, y){
  mean_x <- mean(x)
  mean_y <- mean(y)
  nom <- sum((x - mean_x)*(y-mean_y))
  denom <- sum((x - mean_x)^2)
  m <- nom / denom
  return(m)
}

path = "./data/temp_time_series/"

lakes <- c("Connie", "Crazy", "Ferry", "Gladys", "Heather", "LaCrosse", "Milk", "Sunup",
           "Blue", "Allen", "Bowan", "Deadwood", "EasyRidge", "LH15", "LowerBlum", "LowerEast",
           "LowerSilent", "LP19", "Pallisades", "UpperTriplet")
lakes

for(s in 1:length(lakes)){

  Temp_CV <- list.files(path, pattern = paste0("_",lakes[s],"_")) %>%
  map(~read_csv(paste0(path,.), col_types = cols(), col_names = FALSE, id = "file_name", skip = 4))%>%
  data.table::rbindlist(fill = T)%>%
  rename(date_time = X1,
         temp_C = X2)%>%
  select(-X3,-X4,-X5)%>%
  mutate(file_name = case_when(
    file_name == paste0("./data/temp_time_series/AirTemp_",lakes[s],"_POR.csv") ~ "AirTemp",
    file_name == paste0("./data/temp_time_series/BottomTemp_",lakes[s],"_POR.csv") ~ "BottomTemp",
    file_name == paste0("./data/temp_time_series/MidTemp_",lakes[s],"_POR.csv") ~ "MidTemp",
    file_name == paste0("./data/temp_time_series/SurfaceTemp_",lakes[s],"_POR.csv") ~ "SurfTemp",
    TRUE ~ NA_character_))%>%
  mutate(date_time = lubridate::mdy_hms(date_time))%>%
  filter(file_name!="AirTemp")%>%
  mutate(date = lubridate::date(date_time))%>%
  group_by(date)%>%
  summarize(daily_temp_mean = mean(temp_C, na.rm = T),
            daily_temp_sd = sd(temp_C, na.rm = T),
            daily_temp_CV = abs(daily_temp_sd/daily_temp_mean))%>%
  reshape2::melt(., id = "date")%>%
  filter(variable == "daily_temp_CV")%>%
  mutate(month = lubridate::month(date),
         year = lubridate::year(date))%>%
  ggplot(., aes(date, value))+
  geom_line()+
  theme_bw()

  Temp_CV

  ggsave(path = ".", filename = paste0("./figures/Temp_CV_",lakes[s],".jpg"),
         width = 18, height = 12, device = "jpg", dpi = 400)

  }

for(s in 1:length(lakes)){

  Temp_mean <- list.files(path, pattern = paste0("_",lakes[s],"_")) %>%
    map(~read_csv(paste0(path,.), col_types = cols(), col_names = FALSE, id = "file_name", skip = 4))%>%
    data.table::rbindlist(fill = T)%>%
    rename(date_time = X1,
           temp_C = X2)%>%
    select(-X3,-X4,-X5)%>%
    mutate(file_name = case_when(
      file_name == paste0("./data/temp_time_series/AirTemp_",lakes[s],"_POR.csv") ~ "AirTemp",
      file_name == paste0("./data/temp_time_series/BottomTemp_",lakes[s],"_POR.csv") ~ "BottomTemp",
      file_name == paste0("./data/temp_time_series/MidTemp_",lakes[s],"_POR.csv") ~ "MidTemp",
      file_name == paste0("./data/temp_time_series/SurfaceTemp_",lakes[s],"_POR.csv") ~ "SurfTemp",
      TRUE ~ NA_character_))%>%
    mutate(date_time = lubridate::mdy_hms(date_time))%>%
    filter(file_name!="AirTemp")%>%
    mutate(date = lubridate::date(date_time))%>%
    group_by(date)%>%
    summarize(daily_temp_mean = mean(temp_C, na.rm = T),
              daily_temp_sd = sd(temp_C, na.rm = T),
              daily_temp_CV = abs(daily_temp_sd/daily_temp_mean))%>%
    reshape2::melt(., id = "date")%>%
    filter(variable == "daily_temp_mean")%>%
    mutate(month = lubridate::month(date),
           year = lubridate::year(date))

  write_csv(Temp_mean, paste0("./data/Temp_mean_",lakes[s],".csv"))

    temp_mean_plot <- ggplot(Temp_mean, aes(date, value))+
    geom_line()+
    theme_bw()

    temp_mean_plot

  ggsave(path = ".", filename = paste0("./figures/Temp_mean_",lakes[s],".jpg"),
         width = 18, height = 12, device = "jpg", dpi = 400)


}


for(s in 1:length(lakes)){

  AirTemp_mean <- list.files(path, pattern = paste0("_",lakes[s],"_")) %>%
    map(~read_csv(paste0(path,.), col_types = cols(), col_names = FALSE, id = "file_name", skip = 4))%>%
    data.table::rbindlist(fill = T)%>%
    rename(date_time = X1,
           temp_C = X2)%>%
    select(-X3,-X4,-X5)%>%
    mutate(file_name = case_when(
      file_name == paste0("./data/temp_time_series/AirTemp_",lakes[s],"_POR.csv") ~ "AirTemp",
      file_name == paste0("./data/temp_time_series/BottomTemp_",lakes[s],"_POR.csv") ~ "BottomTemp",
      file_name == paste0("./data/temp_time_series/MidTemp_",lakes[s],"_POR.csv") ~ "MidTemp",
      file_name == paste0("./data/temp_time_series/SurfaceTemp_",lakes[s],"_POR.csv") ~ "SurfTemp",
      TRUE ~ NA_character_))%>%
    mutate(date_time = lubridate::mdy_hms(date_time))%>%
    filter(file_name=="AirTemp")%>%
    mutate(date = lubridate::date(date_time))%>%
    group_by(date)%>%
    summarize(daily_temp_mean = mean(temp_C, na.rm = T),
              daily_temp_sd = sd(temp_C, na.rm = T),
              daily_temp_CV = abs(daily_temp_sd/daily_temp_mean))%>%
    reshape2::melt(., id = "date")%>%
    filter(variable == "daily_temp_mean")%>%
    mutate(month = lubridate::month(date),
           year = lubridate::year(date))%>%
    ggplot(., aes(date, value))+
    geom_line()+
    theme_bw()

  AirTemp_mean

  ggsave(path = ".", filename = paste0("./figures/AirTemp_mean_",lakes[s],".jpg"),
         width = 18, height = 12, device = "jpg", dpi = 400)

}



for(s in 1:length(lakes)){

  Temp_mean <- list.files(path, pattern = paste0("_",lakes[s],"_")) %>%
    map(~read_csv(paste0(path,.), col_types = cols(), col_names = FALSE, id = "file_name", skip = 4))%>%
    data.table::rbindlist(fill = T)%>%
    rename(date_time = X1,
           temp_C = X2)%>%
    select(-X3,-X4,-X5)%>%
    mutate(file_name = case_when(
      file_name == paste0("./data/temp_time_series/AirTemp_",lakes[s],"_POR.csv") ~ "AirTemp",
      file_name == paste0("./data/temp_time_series/BottomTemp_",lakes[s],"_POR.csv") ~ "BottomTemp",
      file_name == paste0("./data/temp_time_series/MidTemp_",lakes[s],"_POR.csv") ~ "MidTemp",
      file_name == paste0("./data/temp_time_series/SurfaceTemp_",lakes[s],"_POR.csv") ~ "SurfTemp",
      TRUE ~ NA_character_))%>%
    mutate(date_time = lubridate::mdy_hms(date_time))%>%
    dplyr::filter(file_name %in% c("BottomTemp", "MidTemp", "SurfTemp"))%>%
    mutate(date = lubridate::date(date_time))%>%
    group_by(date)%>%
    summarize(daily_temp_mean = mean(temp_C, na.rm = T),
              daily_temp_sd = sd(temp_C, na.rm = T),
              daily_temp_CV = abs(daily_temp_sd/daily_temp_mean))%>%
    reshape2::melt(., id = "date")%>%
    filter(variable == "daily_temp_mean")%>%
    mutate(month = lubridate::month(date),
           year = lubridate::year(date))%>%
    ggplot(., aes(date, value))+
    geom_line()+
    theme_bw()

  Temp_mean

  ggsave(path = ".", filename = paste0("./figures/Temp_mean_",lakes[s],".jpg"),
         width = 18, height = 12, device = "jpg", dpi = 400)

}




BottomTemp_Lake <- list.files(path, pattern = paste0("BottomTemp_")) %>%
    map(~read_csv(paste0(path,.), col_types = cols(), col_names = FALSE, id = "file_name", skip = 4))%>%
    data.table::rbindlist(fill = T)%>%
    rename(date_time = X1,
           temp_C = X2)%>%
    select(-X3,-X4,-X5)%>%
    mutate(file_name = case_when(
      file_name == paste0("./data/temp_time_series/BottomTemp_Connie_POR.csv") ~ "Connie",
      file_name == paste0("./data/temp_time_series/BottomTemp_Crazy_POR.csv") ~ "Crazy",
      file_name == paste0("./data/temp_time_series/BottomTemp_Ferry_POR.csv") ~ "Ferry",
      file_name == paste0("./data/temp_time_series/BottomTemp_Gladys_POR.csv") ~ "Gladys",
      file_name == paste0("./data/temp_time_series/BottomTemp_Heather_POR.csv") ~ "Heather",
      file_name == paste0("./data/temp_time_series/BottomTemp_LaCrosse_POR.csv") ~ "LaCrosse",
      file_name == paste0("./data/temp_time_series/BottomTemp_Milk_POR.csv") ~ "Milk",
      file_name == paste0("./data/temp_time_series/BottomTemp_Sunup_POR.csv") ~ "Sunup",
      file_name == paste0("./data/temp_time_series/BottomTemp_Allen_POR.csv") ~ "Allen",
      file_name == paste0("./data/temp_time_series/BottomTemp_Blue_POR.csv") ~ "Blue",
      file_name == paste0("./data/temp_time_series/BottomTemp_Bowan_POR.csv") ~ "Bowan",
      file_name == paste0("./data/temp_time_series/BottomTemp_Deadwood_POR.csv") ~ "Deadwood",
      file_name == paste0("./data/temp_time_series/BottomTemp_EasyRidge_POR.csv") ~ "EasyRidge",
      file_name == paste0("./data/temp_time_series/BottomTemp_LH15_POR.csv") ~ "LH15",
      file_name == paste0("./data/temp_time_series/BottomTemp_LowerBlum_POR.csv") ~ "LowerBlum",
      file_name == paste0("./data/temp_time_series/BottomTemp_LowerEast_POR.csv") ~ "LowerEast",
      file_name == paste0("./data/temp_time_series/BottomTemp_LowerSilent_POR.csv") ~ "LowerSilent",
      file_name == paste0("./data/temp_time_series/BottomTemp_LP19_POR.csv") ~ "LP19",
      file_name == paste0("./data/temp_time_series/BottomTemp_Pallisades_POR.csv") ~ "Pallisades",
      file_name == paste0("./data/temp_time_series/BottomTemp_UpperTriplet_POR.csv") ~ "UpperTriplet",
      TRUE ~ NA_character_))%>%
    mutate(date_time = lubridate::mdy_hms(date_time))%>%
    mutate(date = lubridate::date(date_time))%>%
    group_by(date, file_name)%>%
    summarize(bottom_temp_mean = mean(temp_C, na.rm = T),
              bottom_temp_se = se(temp_C))%>%
  arrange(file_name)%>%
  mutate(year = lubridate::year(date))%>%
  group_by(file_name, year)%>%
  summarize(max_temp = max(bottom_temp_mean))%>%
  summarize

ggplot(BottomTemp_Lake, aes(year, max_temp))+
  geom_point(size = 2)+
  geom_smooth(method = "lm")+
  theme_bw()+
  facet_wrap(~file_name, scales = "free_y")

AirTemp_Lake <- list.files(path, pattern = paste0("AirTemp_")) %>%
  map(~read_csv(paste0(path,.), col_types = cols(), col_names = FALSE, id = "file_name", skip = 4))%>%
  data.table::rbindlist(fill = T)%>%
  rename(date_time = X1,
         temp_C = X2)%>%
  select(-X3,-X4,-X5)%>%
  mutate(file_name = case_when(
    file_name == paste0("./data/temp_time_series/AirTemp_Connie_POR.csv") ~ "Connie",
    file_name == paste0("./data/temp_time_series/AirTemp_Crazy_POR.csv") ~ "Crazy",
    file_name == paste0("./data/temp_time_series/AirTemp_Ferry_POR.csv") ~ "Ferry",
    file_name == paste0("./data/temp_time_series/AirTemp_Gladys_POR.csv") ~ "Gladys",
    file_name == paste0("./data/temp_time_series/AirTemp_Heather_POR.csv") ~ "Heather",
    file_name == paste0("./data/temp_time_series/AirTemp_LaCrosse_POR.csv") ~ "LaCrosse",
    file_name == paste0("./data/temp_time_series/AirTemp_Milk_POR.csv") ~ "Milk",
    file_name == paste0("./data/temp_time_series/AirTemp_Sunup_POR.csv") ~ "Sunup",
    file_name == paste0("./data/temp_time_series/AirTemp_Allen_POR.csv") ~ "Allen",
    file_name == paste0("./data/temp_time_series/AirTemp_Blue_POR.csv") ~ "Blue",
    file_name == paste0("./data/temp_time_series/AirTemp_Bowan_POR.csv") ~ "Bowan",
    file_name == paste0("./data/temp_time_series/AirTemp_Deadwood_POR.csv") ~ "Deadwood",
    file_name == paste0("./data/temp_time_series/AirTemp_EasyRidge_POR.csv") ~ "EasyRidge",
    file_name == paste0("./data/temp_time_series/AirTemp_LH15_POR.csv") ~ "LH15",
    file_name == paste0("./data/temp_time_series/AirTemp_LowerBlum_POR.csv") ~ "LowerBlum",
    file_name == paste0("./data/temp_time_series/AirTemp_LowerEast_POR.csv") ~ "LowerEast",
    file_name == paste0("./data/temp_time_series/AirTemp_LowerSilent_POR.csv") ~ "LowerSilent",
    file_name == paste0("./data/temp_time_series/AirTemp_LP19_POR.csv") ~ "LP19",
    file_name == paste0("./data/temp_time_series/AirTemp_Pallisades_POR.csv") ~ "Pallisades",
    file_name == paste0("./data/temp_time_series/AirTemp_UpperTriplet_POR.csv") ~ "UpperTriplet",
    TRUE ~ NA_character_))%>%
  mutate(date_time = lubridate::mdy_hms(date_time))%>%
  mutate(date = lubridate::date(date_time))%>%
  group_by(date, file_name)%>%
  summarize(air_temp_mean = mean(temp_C, na.rm = T),
            air_temp_se = se(temp_C))%>%
  arrange(file_name)%>%
  mutate(year = lubridate::year(date))%>%
  group_by(file_name, year)%>%
  summarize(max_temp = max(air_temp_mean))

ggplot(AirTemp_Lake, aes(year, max_temp))+
  geom_point(size = 2)+
  geom_smooth(method = "lm")+
  theme_bw()+
  facet_wrap(~file_name, scales = "free_y")

data_join <- left_join(BottomTemp_Lake, AirTemp_Lake, by = c("year", "file_name"))

