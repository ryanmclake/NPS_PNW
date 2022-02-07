set.seed(71)


### Perform an RDA on the Zooplankton
zoop_cap <- capscale(formula = select(env_zoop_data, CLAD:RAP) ~
                       Chlorophyll + Elevation_m +
                       Ca + solar_jas + DOC + lake_temp +`Total N` +
                       `Total P` + fish, data = env_zoop_data, distance = "bray")

zoop_null <- capscale(formula = select(env_zoop_data, CLAD:RAP) ~ 1,
                      data = env_zoop_data, distance = "bray")

mods <- ordiR2step(zoop_null, scope = formula(zoop_cap), trace = 0, direction = c("forward"),
                 permutations = how(nperm = 999), steps = 100)
mods$anova
model_site_points <- bind_cols(data.frame(scores(zoop_cap)$sites),env_zoop_data)

all_data_model <- autoplot(zoop_cap, layers = c("biplot", "species")) +
  geom_point(data = model_site_points,
             aes(x = CAP1, y = CAP2, color = as.character(park_code), alpha = 0.5), size = 4)+
  theme_classic()

all_data_model
ggsave("./figures/CAPSCALE_model_output.jpg", width = 10, height = 10, units = "in")



### Perform an RDA on the Macro invertebrates
macro_cap <- capscale(formula = select(macro_join, Acari:Veneroida) ~
                        Chlorophyll + Elevation_m +
                        Ca + solar_jas + DOC + lake_temp +`Total N` +
                        `Total P` + fish, data = macro_join, distance = "bray")

macro_null <- capscale(formula = select(macro_join, Acari:Veneroida) ~ 1,
                      data = macro_join, distance = "bray")

mods_macro <- ordiR2step(macro_null, scope = formula(macro_cap), trace = 0, direction = c("forward"),
                   permutations = how(nperm = 999), steps = 100)
mods_macro$anova
model_site_points <- bind_cols(data.frame(scores(macro_cap)$sites),macro_join)

all_data_model <- autoplot(macro_cap, layers = c("biplot", "species")) +
  geom_point(data = model_site_points,
             aes(x = CAP1, y = CAP2, color = as.character(park_code), alpha = 0.5), size = 2.5)+
  theme_classic()+
  labs(title = "Aquatic Macroinvertebrate RDA by Taxon")

all_data_model
ggsave("./figures/CAPSCALE_model_output_macros.jpg", width = 10, height = 10, units = "in")




##############################################################
# RUN CCA ANALYSIS OF DRIVER SIGNIFICATNLY DIFFERENT FROM NULL
##############################################################

Y <- env_zoop_data %>% select(CLAD:RAP)
X <- env_zoop_data %>% select(lake_temp,solar_jas,stability,Elevation_m,ice_out_doy,Depth_max)

cc_results <- cancor(X,Y) # that is it

str(cc_results)
cc_results$xcoef
cc_results$ycoef
cc_results$cor

CC1_X <- as.matrix(X) %*% cc_results$xcoef[, 1]
CC2_X <- as.matrix(X) %*% cc_results$xcoef[, 2]
CC3_X <- as.matrix(X) %*% cc_results$xcoef[, 3]
CC4_X <- as.matrix(X) %*% cc_results$xcoef[, 4]
CC5_X <- as.matrix(X) %*% cc_results$xcoef[, 5]
CC6_X <- as.matrix(X) %*% cc_results$xcoef[, 6]

CC1_Y <- as.matrix(Y) %*% cc_results$ycoef[, 1]
CC2_Y <- as.matrix(Y) %*% cc_results$ycoef[, 2]
CC3_Y <- as.matrix(Y) %*% cc_results$ycoef[, 3]
CC4_Y <- as.matrix(Y) %*% cc_results$ycoef[, 4]

cca_df <- env_zoop_data %>%
  mutate(CC1_X=CC1_X,
         CC1_Y=CC1_Y,
         CC2_X=CC2_X,
         CC2_Y=CC2_Y,
         CC3_X=CC3_X,
         CC3_Y=CC3_Y,
         CC4_X=CC4_X,
         CC4_Y=CC4_Y)

CC1X_fig <- cca_df %>%
  ggplot(aes(x=park_code,y=CC1_X, color=park_code))+
  geom_boxplot(width=0.5)+
  geom_jitter(width=0.15)+
  theme(legend.position="none")+
  theme_classic()

CC1X_fig
ggsave("./figures/CC1X_boxplot_parks.jpg", width = 10, height = 10, units = "in")

CC1Y_fig <-cca_df %>%
  ggplot(aes(x=park_code,y=CC1_Y, color=park_code))+
  geom_boxplot(width=0.5)+
  geom_jitter(width=0.15)+
  theme(legend.position="none")+
  theme_classic()

CC1Y_fig
ggsave("./figures/CC1Y_boxplot_parks.jpg", width = 10, height = 10, units = "in")

CC1X_CC1Y_park <- cca_df %>%
  ggplot(aes(x=CC1_X,y=CC1_Y, color=park_code))+
  geom_point()+
  theme_classic()

CC1X_CC1Y_park
ggsave("./figures/CC1X_vs_CC1Y_parks.jpg", width = 10, height = 10, units = "in")

CC1X_CC1Y_lake <- cca_df %>%
  ggplot(aes(x=CC1_X,y=CC1_Y, color=site_code))+
  geom_point()+
  theme_classic()

CC1X_CC1Y_lake
ggsave("./figures/CC1X_vs_CC1Y_lakes.jpg", width = 10, height = 10, units = "in")
