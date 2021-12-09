#########################
# ZOOPLANKTON COMPILATION
#########################

df <- read_excel("./data/NCCN_Zoops_Combined_For_FCA_May_2019.xlsx")
df[1487, "Taxa"] <- "ROT"
df[1487, "Taxonid"] <- 5
zoop_no_rotifers <- df %>%
  filter(Taxa != "ROT")%>%
  filter(Taxa != "INS")

#Partition Rotifers into functional feeding groups
rot_wFFG <- df %>% filter(Taxa == "ROT")
rot_wFFG$GenSp <- sub(" .*", "", rot_wFFG$GenSp)

rot_wFFG <- rot_wFFG %>%
  mutate(trophi = case_when(
    GenSp == "Conochilius" ~ "malleoramate",
    GenSp == "Kellicottia" ~ "malleate",
    GenSp == "Keratella" ~ "malleate",
    GenSp == "Lecane" ~ "malleate",
    GenSp == "Monostyla" ~ "malleate",
    GenSp == "Polyarthra" ~ "virgate",
    GenSp == "Trichotria" ~ "malleate",
    GenSp == "Synchaeta" ~ "virgate",
    GenSp == "Collotheca" ~ "uncinate",
    GenSp == "Notholca" ~ "malleate",
    GenSp == "Brachionus" ~ "malleate",
    GenSp == "Gastropus" ~ "virgate",
    GenSp == "Philodina" ~ "malleate",
    GenSp == "Ascomorpha" ~ "virgate",
    GenSp == "Filinia" ~ "malleoramate",
    GenSp == "Trichocerca" ~ "virgate",
    GenSp == "Ploesoma" ~ "virgate",
    GenSp == "Conochiloides" ~ "malleoramate",
    TRUE ~ NA_character_),
    feed_type = case_when(
      trophi == "uncinate" ~ "RAP",
      trophi == "virgate" ~ "RAP",
      trophi == "malleate" ~ "MICRO",
      trophi == "malleoramate" ~ "MICRO",
      TRUE ~ NA_character_))%>%
  select(`Park&Year&Lake`,Park, Year, Lake, Code, Density, GenSp, feed_type, Taxonid)%>%
  rename(Taxa = feed_type)

match_sites <- readRDS(file = "./data/name_site_matches.rds")

zoop <- bind_rows(zoop_no_rotifers, rot_wFFG)%>%
  left_join(x = ., y = match_sites, by = c("Lake" = "old_name"))%>%
  clean_names()%>%
  group_by(park, site_code, lake, year, taxa) %>%
  summarise(density = sum(density))%>%
  ungroup() %>%
  mutate(park_site_year = paste0(park, "_", lake, "_", year)) %>%
  select(park_site_year, taxa, density) %>%
  data.frame() %>%
  matrify()%>%
  select(-HAR)%>%
  rownames_to_column(.) %>%
  separate(col = rowname, into = c("park_code", "site_code", "event_year"), sep = "_")%>%
  mutate(event_year = as.character(event_year))%>%
  mutate(RAP = log10(RAP+1),
         MICRO = log10(MICRO+1),
         CLAD = log10(CLAD+1),
         COPE = log10(COPE+1))

########################################
# PHYSICAL AND CHEMICAL DATA COMPILATION
########################################

env_dat_all <- readRDS(file = "./data/bigjoin.rds")
env_dat_yr <- env_dat_all %>%
  filter(variable %in% c("Ca", "Chlorophyll", "DOC","Total N","Total P",
                         "BotTemp","MidTemp","SurfTemp", "ice_out_doy",
                         "BRK_ever","CCT_ever","RBT_ever","TSS_ever","WCT_ever")) %>%
  unique() %>%
  select(park_code, Lake, site_code, event_year, variable, value, solar_jas, Elevation_m, Depth_max, lon, lat)%>%
  spread(key = variable, value = value)%>%
  mutate(event_year = as.character(event_year))%>%
  select(-site_code)%>%
  rename(site_code = Lake)

env_dat_yr$site_code <- gsub("Lake ", "", env_dat_yr$site_code)
env_dat_yr$site_code <- gsub(" Lake", "", env_dat_yr$site_code)
env_dat_yr$site_code <- gsub("Upper ", "", env_dat_yr$site_code)
env_dat_yr$site_code <- gsub("Lower ", "", env_dat_yr$site_code)
env_dat_yr$site_code <- gsub("Easy Ridge", "EasyRidge", env_dat_yr$site_code)
env_dat_yr$site_code <- gsub("Lower ", "", env_dat_yr$site_code)
env_dat_yr$site_code <- gsub("La Crosse", "LaCrosse", env_dat_yr$site_code)

########################
# JOIN DATA FOR ANALYSIS
########################

env_zoop_data <- env_dat_yr %>% left_join(., zoop, by = c("park_code","site_code","event_year")) %>%
  filter(site_code != "La Crosse" & site_code != "Palisades" & site_code != "Hoh")%>%
  group_by(site_code, park_code)%>%
  mutate(fish = ifelse((CCT_ever+TSS_ever+WCT_ever+RBT_ever+BRK_ever)>0,1,0))%>%
  select(-CCT_ever, -TSS_ever, -WCT_ever, -RBT_ever, -BRK_ever)%>%
  mutate(lake_temp = (SurfTemp+BotTemp+MidTemp)/3,
         stability = SurfTemp - BotTemp)%>%
  select(-SurfTemp,-BotTemp,-MidTemp)%>%
  mutate_at(vars(-site_code, -park_code, -event_year),
            funs(imputeTS::na_interpolation(., option = "spline")))%>%
  mutate_at(vars(Ca, Chlorophyll, DOC, ice_out_doy,`Total N`,
                 `Total P`,lake_temp,stability),
            funs((log10(.+1))))%>%
  mutate_at(vars(-site_code, -park_code, -event_year),
            funs(imputeTS::na_interpolation(., option = "spline")))%>%
  ungroup(.)

viz_dat <- vis_dat(env_zoop_data)
viz_dat
ggsave("./figures/MISSING_DATA_CHECK.jpg", width = 10, height = 8, units = "in")

mapping_zoop <- env_zoop_data %>% select(park_code, site_code, lon, lat, Elevation_m, Depth_max, fish, solar_jas)
