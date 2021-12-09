#Random Forest Script
set.seed(71)

directory <- here::here()
config <- yaml::read_yaml(file.path(paste0(directory,"/config/zoop_RF_config.yml")))

zoop_data <- read_csv(config$file_path$zoop_data)%>%
  mutate(event_year = as.character(event_year))
  # mutate(RAP = log10(RAP+1),
  #        MICRO = log10(MICRO+1),
  #        CLAD = log10(CLAD+1),
  #        COPE = log10(COPE+1))

env_zoop_data <- read_csv(config$file_path$environment_data)%>%
  mutate(event_year = as.character(event_year))%>%
  left_join(., zoop_data, by = c("park_code","site_code","event_year"))

write_csv(env_zoop_data, "ZoopData.csv")

park = config$attributes$park_code
taxa = config$attributes$taxa


RF_data <- env_zoop_data%>%
  filter(site_code != "La Crosse" & site_code != "Palisades" & site_code != "Hoh")%>%
  group_by(site_code, park_code)%>%
  #summarize_all(funs(mean))%>%
  mutate(fish = ifelse((CCT_ever+TSS_ever+WCT_ever+RBT_ever+BRK_ever)>0,1,0))%>%
  select(-CCT_ever, -TSS_ever, -WCT_ever, -RBT_ever, -BRK_ever)%>%
  mutate(lake_temp = (SurfTemp+BotTemp+MidTemp)/3,
         stability = SurfTemp - BotTemp)%>%
  select(-SurfTemp,-BotTemp,-MidTemp)%>%
  mutate_at(vars(-site_code, -park_code, -event_year),
            funs(imputeTS::na_interpolation(., option = "linear")))%>%
  ungroup(.)

com = RF_data[,13:16]
m_com = as.matrix(com)

ord = metaMDS(m_com, distance = "bray", trace = F, autotransform = T)
nmds_scores <- cbind(as_tibble(scores(ord)),
                     RF_data[1:2])

NMDS = data.frame(MDS1 = ord$points[,1], MDS2 = ord$points[,2], lake = RF_data[,2], park = RF_data[,1])
simp <- simper(m_com, RF_data$park_code)
summary(simp)
vec.sp<-envfit(ord$points, RF_data[,-c(1,2,3,13,14,15,16)], strata = RF_data$site_code, perm=1000)
vec.sp.df<-as.data.frame(vec.sp$vectors$arrows*sqrt(vec.sp$vectors$r))
vec.sp.df$species<-rownames(vec.sp.df)

zoop_cap <- capscale(formula = select(RF_data, CLAD:RAP) ~
                       stability + Chlorophyll + ice_out_doy + Elevation_m +
                       Ca + solar_jas + DOC + lake_temp +`Total N` +
                       `Total P` + Depth_max,
                     data = RF_data, distance = "bray")

zoop_null <- capscale(formula = select(RF_data, CLAD:RAP) ~ 1,
                      data = RF_data, distance = "bray")

mods <- ordistep(zoop_null, scope = formula(zoop_cap), trace = 0)

mods$anova

vif.cca(zoop_cap)



custom_site_points <- bind_cols(data.frame(scores(zoop_cap)$sites),
                                RF_data)
library(ggvegan)

all_data_model <- autoplot(zoop_cap, layers = c("biplot", "species")) +
  geom_point(data = custom_site_points,
             aes(x = CAP1, y = CAP2, color = as.character(site_code), alpha = 0.5), size = 2)




plot_test <- ggplot(data = NMDS, aes(MDS1, MDS2)) +
  geom_point(data = nmds_scores, aes(NMDS1, NMDS2, color = site_code ,inherit_aes=F), size = 3)+
  #stat_ellipse(geom = "polygon", aes(group = park_code, color = park_code, fill = park_code), alpha = 0.1) +
  geom_segment(data=vec.sp.df,aes(x=0,xend=MDS1,y=0,yend=MDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="black",inherit_aes=F) +
  geom_text(data=vec.sp.df,aes(x=MDS1,y=MDS2,label=species),size=3)+
  xlim(c(-1,1.1))+
  ylim(c(-0.6,0.6))+
  theme_bw()

plot_test <- ggplot(data = NMDS, aes(MDS1, MDS2)) +
  geom_point(data = nmds_scores, aes(NMDS1, NMDS2, color = park_code ,inherit_aes=F), size = 3)+
  #stat_ellipse(geom = "polygon", aes(group = park_code, color = park_code, fill = park_code), alpha = 0.1) +
  geom_segment(data=vec.sp.df,aes(x=0,xend=MDS1,y=0,yend=MDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="black",inherit_aes=F) +
  geom_text(data=vec.sp.df,aes(x=MDS1,y=MDS2,label=species),size=3)+
  xlim(c(-1,1.1))+
  ylim(c(-0.6,0.6))+
  theme_bw()

#stress 0.128
