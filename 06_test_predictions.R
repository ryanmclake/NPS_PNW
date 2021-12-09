#Prediction 1 Warmer temperatures in the water column will favor cladocerans and rotifers

clad <- ggplot(env_zoop_data, aes(lake_temp, CLAD))+
  geom_point(aes(color = park_code))+
  geom_smooth(method = "lm")+
  labs(title = "slope coefficient = 6.14")+theme_classic()+
  theme(legend.position = c(0.2,0.8))
m1 <- glm(CLAD~lake_temp, data = env_zoop_data)
summary(m1)

rap <- ggplot(env_zoop_data, aes(lake_temp, RAP))+
  geom_point(aes(color = park_code))+
  geom_smooth(method = "lm")+
  labs(title = "slope coefficient = 2.27")+
  theme(legend.position = "none")+theme_classic()
m2 <- glm(RAP~lake_temp, data = env_zoop_data)
summary(m2)

micro <- ggplot(env_zoop_data, aes(lake_temp, MICRO))+
  geom_point(aes(color = park_code))+
  geom_smooth(method = "lm")+
  labs(title = "slope coefficient = 6.18")+theme_classic()+
  theme(legend.position = "none")
m3 <- glm(MICRO~lake_temp, data = env_zoop_data)
summary(m3)

cope <- ggplot(env_zoop_data, aes(lake_temp, COPE))+
  geom_point(aes(color = park_code))+
  geom_smooth(method = "lm")+
  labs(title = "slope coefficient = 1.83")+theme_classic()+
  theme(legend.position = "none")
m4 <- glm(COPE~lake_temp, data = env_zoop_data)
summary(m4)

zoops_temp <- (clad+cope)/(micro+rap)
zoops_temp
ggsave("./figures/Zoop_temp_relations.jpg", width = 15, height = 15, units = "in")


# Fish abundance will favor rotifers and copepods
cope_fish <- t.test(COPE ~ as.character(fish), data = env_zoop_data)
cope_fish

rap_fish <- t.test(RAP ~ as.character(fish), data = env_zoop_data)
rap_fish

micro_fish <- t.test(MICRO ~ as.character(fish), data = env_zoop_data)
micro_fish

clad_fish <- t.test(CLAD ~ as.character(fish), data = env_zoop_data)
clad_fish


library(ggpubr)

cope_fish <- ggplot(env_zoop_data, aes(as.character(fish), COPE, group = as.character(fish), color = park_code, size = Depth_max))+
  geom_boxplot(fill = "grey90")+
  geom_jitter(width = 0.05)+
  stat_compare_means(label = "p.signif", method = "t.test",label.x = 1.5, size = 10, paired = F)+
  xlab("Fish? 0=Y 1=N")+
  theme_classic2()+
  theme(legend.position = "none")


clad_fish <- ggplot(env_zoop_data, aes(as.character(fish), CLAD, group = as.character(fish), color = park_code, size = Depth_max))+
  geom_boxplot(fill = "grey90")+
  geom_jitter(width = 0.05)+
  stat_compare_means(label = "p.signif", method = "t.test",label.x = 1.5, size = 10, paired = F)+
  xlab("Fish? 0=Y 1=N")+
  theme_classic2()+
  theme(legend.position = "none")

rap_fish <- ggplot(env_zoop_data, aes(as.character(fish), RAP, group = as.character(fish), color = park_code, size = Depth_max))+
  geom_boxplot(fill = "grey90")+
  geom_jitter(width = 0.05)+
  stat_compare_means(label = "p.signif", method = "t.test",label.x = 1.5, size = 10, paired = F)+
  xlab("Fish? 0=Y 1=N")+
  theme_classic2()+
  theme(legend.position = "none")

micro_fish <- ggplot(env_zoop_data, aes(as.character(fish), MICRO, group = as.character(fish), color = park_code, size = Depth_max))+
  geom_boxplot(fill = "grey90")+
  geom_jitter(width = 0.05)+
  stat_compare_means(label = "p.signif", method = "t.test",label.x = 1.5, size = 10, paired = F)+
  xlab("Fish? 0=Y 1=N")+
  theme_classic2()+
  theme(legend.position = "none")

zoops_fish <- (clad_fish+cope_fish)/(micro_fish+rap_fish)
zoops_fish
ggsave("./figures/Zoop_fish_relations.jpg", width = 15, height = 15, units = "in")


# Dominance of microphagous rotifers coincided with a decrease of cladocerans and
# raptorial rotifers will coincide with an increase in cladocerans. Following Obertregger 2011

env_zoop_data %>% group_by(park_code, event_year)%>%
  summarize_all(funs(mean))%>%

  ggplot(., aes(as.numeric(event_year), CLAD, color = park_code))+
  geom_line()
