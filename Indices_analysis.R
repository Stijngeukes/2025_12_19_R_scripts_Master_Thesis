#-----------------LIBRARIES_AND_WD----------------

# Set working directory.
setwd('C:/Users/stijn/Documents/Drive/Earth Sciences/Thesis')

# Load libraries
library(reshape2)
library(suncalc)
library(tidyverse)
library(readxl)
library(stringr)
library(dplyr)
library(emmeans)   
library(writexl)
library(glmmTMB)
library(scales)
library(purrr)
library(patchwork)

#Load outputs
files <- list.files("Indices_outputs", pattern = "^acoustic.*\\.csv$", full.names = TRUE)
extrafiles <- list.files("Indices_outputs/acoustic_indices_bugg", pattern = "^acoustic.*\\.csv$", full.names = TRUE)

#Create list of indices and useable columns
indices = c('BI', 'ACI', 'AEI', 'NDSI', 'ADI', 'H', 'CENT', 'TQ')
column_keep = c('bugg', 'File', 'Date', 'BI', 'ACI', 'AEI', 'NDSI', 'ADI', 'H', 'CENT', 'TQ')

habitat_colours = c('#33A02C', '#FF7F00', '#1F78B4')
#Load deployment info
deployment = read_xlsx('Deployment_locations.xlsx')

#------------COMBINE_FILES------------------

#merge old and new files into a complete df

cutoff_date <- ymd("2025-08-31")

original_files_df <- files %>%
  lapply(function(f) read.csv(f, header = TRUE, sep = c(",", '.'))) %>%
  bind_rows()

original_files_df <- original_files_df %>%
  mutate(Date_only = as_date(Date))

original_part <- original_files_df %>%
  filter(Date_only < cutoff_date)

extra_files_df <- extrafiles %>%
  lapply(function(f) read.csv(f, header = TRUE, sep = c(",", '.'))) %>%
  bind_rows()

complete_df = bind_rows(original_part, extra_files_df)

#Compute acoustic entropy index from spectral and temporal entropy
complete_df$H = complete_df$Hf * complete_df$Ht

#remove redundant columns
filtered_df = complete_df[,column_keep]

# Couple device ID to values
filtered_df <- filtered_df %>%
  mutate(Device_ID = str_sub(bugg, -8, -1))

#Merge to deployment data
merged_df <- filtered_df %>%
  left_join(deployment, by = "Device_ID")

#Bin for 15 minutes
merged_df$Date <- ymd_hms(as.character(merged_df$Date))
merged_df$Time15 <- floor_date(merged_df$Date, "15 minutes")
merged_df$Time15_ofday <- format(merged_df$Time15, "%H:%M:%S")
merged_df$day = as.Date(merged_df$Date)

#-------------SEASONAL_PERIODS------------------
merged_df <- merged_df %>%
  # Extract the day of the year. 
  mutate(DayOfYear = lubridate::yday(Date)) %>%
  
  # March 21st is Day 80, June 21st is Day 172, September 21st is Day 264
  mutate(season = case_when(
    DayOfYear >= 80 & DayOfYear < 172 ~ "spring", 
    DayOfYear >= 172 & DayOfYear <= 264 ~ "summer", 
    TRUE ~ 'other'
  ))

#----------------Diel_part--------------------------
#Compute sunrise, sunset and nautical twilight times for every timepoint and location
sun_times <- merged_df %>%
  distinct(Device_ID, Lat, Long, day) %>%
  rowwise() %>%
  mutate(
    times = list(getSunlightTimes(
      date = day,
      lat = Lat,
      lon = Long,
      keep = c("nauticalDawn", "nauticalDusk", "sunrise", "sunset")
    ))
  ) %>%
  unnest(cols = c(times))

#join the files
merged_df <- merged_df %>%
  left_join(sun_times, by = c("Device_ID", "Lat", "Long", "day"))

#Compute diel part based on suntimes
merged_df <- merged_df %>%
  mutate(
    diel = case_when(
      Time15 >= nauticalDusk | Time15 < nauticalDawn ~ "night",
      Time15 >= nauticalDawn & Time15 < sunrise ~ "dawn",
      Time15 >= sunrise & Time15 < sunset ~ "day",
      Time15 >= sunset | Time15 < nauticalDusk ~ "dusk",
      TRUE ~ NA_character_
    )
  )

#Clean up df
merged_df = merged_df %>%
  select(-DayOfYear, -lat, -lon, -nauticalDawn, -nauticalDusk, -sunrise, -sunset)

merged_df = merged_df %>%
  mutate(File = str_remove(File, "\\.mp3$"))
#----------------Standardised_indices------------------------
standard_df <- merged_df
N <- nrow(merged_df)  # sample size

# Linearly transform each index to (0, 1) using the formula from Smithson and Verkuilen 2006:
# y' = (y - a) / (b - a)
# y"  = [y'(N - 1) + 1/2] / N
for (col in indices) {
  
  # For NDSI, first adjust range from -1-1 to 0-1 before scaling
  if (col == "NDSI") {
    temp <- (merged_df[[col]] + 1) / 2
  } else {
    temp <- merged_df[[col]]
  }
  
  # Get min and max of the column
  a <- min(temp, na.rm = TRUE)
  b <- max(temp, na.rm = TRUE)
  
  # Compute y' and y
  y_prime <- (temp - a) / (b - a)
  y_scaled <- ((y_prime * (N - 1)) + 0.5) / N
  
  # Assign scaled values
  standard_df[[col]] <- y_scaled
}

#----------------stratified random sample--------------------------------------
files_to_exclude = read.csv("Subsampling/subsampled_audio_files_complete.csv", sep = ';')


subsample_df_extra <- merged_df %>%
  anti_join(files_to_exclude, by = "File") %>%
  group_by(Device_ID, season, diel) %>%
  slice_sample(n = 3, replace = FALSE) %>%
  ungroup()  %>%
  select(File, Device_ID, season, diel)

# Save to CSV
write.csv(subsample_df_extra, "subsampled_audio_files_extra.csv")


#----------------AIC comparison between model complexities---------------------
AIChab2 = list()
AIChabnest = list()
AIChab = list()
AICnest = list()
AICint = list()
AICnonint = list()
AICseas = list()
AICdiel = list()

for (i in indices){
  #formulas for each model
  formulahab2 = as.formula(paste(i, '~diel*season + diel * Habitat + season * Habitat + (1|Device_ID)'))
  formulahabnest = as.formula(paste(i, '~diel * season + Habitat + (1|Cluster/Device_ID)'))
  formulahab = as.formula(paste(i, '~diel * season + Habitat + (1|Device_ID)'))
  formulanest = as.formula(paste(i, '~diel * season + (1|Cluster/Device_ID)'))
  formulaint = as.formula(paste(i, '~diel * season + (1|Device_ID)'))
  formulanonint  = as.formula(paste(i, '~diel + season + (1|Device_ID)'))
  formulaseas  = as.formula(paste(i, '~season + (1|Device_ID)'))
  formuladiel  = as.formula(paste(i, '~diel + (1|Device_ID)'))
  
  hab2mod = glmmTMB(formulahab2,
                    data = standard_df,
                    family = beta_family(link = 'logit'),
                    control = glmmTMBControl(optArgs = list(iter.max = 1000, eval.max = 1000))  )
  habnestmod = glmmTMB(formulahabnest,
                    data = standard_df,
                    family = beta_family(link = 'logit'),
                    control = glmmTMBControl(optArgs = list(iter.max = 1000, eval.max = 1000))  )
  habmod = glmmTMB(formulahab,
                    data = standard_df,
                    family = beta_family(link = 'logit'),
                   control = glmmTMBControl(optArgs = list(iter.max = 1000, eval.max = 1000))  )
  nestmod = glmmTMB(formulanest,
                    data = standard_df,
                    family = beta_family(link = 'logit'),
                    control = glmmTMBControl(optArgs = list(iter.max = 1000, eval.max = 1000)) ) 
  intmod = glmmTMB(formulaint,
                   data = standard_df,
                   family = beta_family(link = 'logit'),
                   control = glmmTMBControl(optArgs = list(iter.max = 1000, eval.max = 1000)) ) 
  nonintmod = glmmTMB(formulanonint,
                      data = standard_df,
                      family = beta_family(link = 'logit'),
                      control = glmmTMBControl(optArgs = list(iter.max = 1000, eval.max = 1000)))  
  dimod = glmmTMB(formuladiel,
                  data = standard_df,
                  family = beta_family(link = 'logit'),
                  control = glmmTMBControl(optArgs = list(iter.max = 1000, eval.max = 1000)))  
  semod = glmmTMB(formulaseas,
                  data = standard_df,
                  family = beta_family(link = 'logit'),
                  control = glmmTMBControl(optArgs = list(iter.max = 1000, eval.max = 1000)))  
  #Put AIC difference in list.
  AIChab2 [[i]] = AIC(hab2mod)
  AIChabnest[[i]] = AIC(habnestmod)
  AIChab[[i]] = AIC(habmod)
  AICnest[[i]] = AIC(nestmod)
  AICint[[i]] = AIC(intmod)
  AICnonint[[i]] = AIC(nonintmod)
  AICseas[[i]] = AIC(semod)
  AICdiel[[i]] = AIC(dimod)
  
  #If AIC is lower, model fits better.
  #This means that if the values in the list are negative, the habnestmod performs best. 
  #If positive, the subtracted model performs best.
  
  message('Finished:', i)
}

AICcombined = list(
  eight = AIChab2,
  seven = AIChabnest,
  six = AIChab,
  five = AICnest,
  four = AICint,
  three = AICnonint,
  two = AICseas,
  one = AICdiel
)

AICdf = as.data.frame(do.call(cbind, AICcombined))
cols_to_convert <- c('one', 'two', 'three', 
                    'four', 'five', 'six', 'seven', 'eight')
AICdf[cols_to_convert] <- lapply(AICdf[cols_to_convert], as.numeric)

write_xlsx(AICdf, path = 'analysis/AIC_comparison.xlsx')

beepr::beep(sound = 3)


#-----------run the best fitted model-------------

list_glmer = list()
list_emmeans = list()
list_pairs = list()

list_emshabitatdiel = list()
list_emshabitatseason = list()
list_emsdielseason = list()

list_pairshabitatdiel = list()
list_pairshabitatseason = list()
list_pairsdielseason = list()

for (i in indices){
  
  #Formula of the model
  formulaint = as.formula(paste(i, '~diel * season + Habitat * season + Habitat * diel + (1|Device_ID)'))
  
  #Run the model
  intmod = glmmTMB(formulaint,
                   data = standard_df,
                   family = beta_family('logit'))
  
  #Run emmeans and pairs on all interactions
  sum_intmod = summary(intmod)
  ems_intmod = emmeans(intmod, specs = ~ diel * season + Habitat * season + Habitat * diel)
  pairs_intmod = pairs(ems_intmod, adjust = 'tukey')
  
  #Run emmeans and pairs on specific interactions
  emshabitatdiel = emmeans(intmod, specs = ~ Habitat * diel)
  emshabitatseason = emmeans(intmod, specs = ~ Habitat * season)
  emsdielseason = emmeans(intmod, specs = ~ diel * season)

  pairshabitatdiel = pairs(emshabitatdiel, adjust = 'tukey')
  pairshabitatseason = pairs(emshabitatseason, adjust = 'tukey')
  pairsdielseason = pairs(emsdielseason, adjust = 'tukey')
  
  #Save all results in dataframes
  sumdf = as.data.frame(sum_intmod$coefficients$cond)
  emsdf = as.data.frame(ems_intmod)
  pairsdf = as.data.frame(pairs_intmod)
  
  emshabitatdieldf = as.data.frame(emshabitatdiel)
  emshabitatseasondf = as.data.frame(emshabitatseason)
  emsdielseasondf = as.data.frame(emsdielseason)
  
  pairshabitatdieldf = as.data.frame(pairshabitatdiel)
  pairshabitatseasondf = as.data.frame(pairshabitatseason)
  pairsdielseasondf = as.data.frame(pairsdielseason)

  sumdf$index = i
  emsdf$index = i
  pairsdf$index = i
  
  emshabitatdieldf$index = i
  emshabitatseasondf$index = i
  emsdielseasondf$index = i

  pairshabitatdieldf$index = i
  pairshabitatseasondf$index = i
  pairsdielseasondf$index = i
  
  list_glmer[[i]] = sumdf
  list_emmeans[[i]] = emsdf
  list_pairs[[i]] = pairsdf
  
  list_emshabitatdiel[[i]] = emshabitatdieldf
  list_emshabitatseason[[i]] = emshabitatseasondf
  list_emsdielseason[[i]] = emsdielseasondf
  
  list_pairshabitatdiel[[i]] = pairshabitatdieldf
  list_pairshabitatseason[[i]] = pairshabitatseasondf
  list_pairsdielseason[[i]] = pairsdielseasondf

  message("Finished:", i)
}

df_glmer = bind_rows(list_glmer)
df_glmer$term <- rownames(df_glmer)
df_emmeans = bind_rows(list_emmeans)
df_pairs = bind_rows(list_pairs)

df_emmeanshabitatdiel = bind_rows(list_emshabitatdiel)
df_emmeanshabitatseason = bind_rows(list_emshabitatseason)
df_emmeansdielseason = bind_rows(list_emsdielseason)

df_pairshabitatdiel = bind_rows(list_pairshabitatdiel)
df_pairshabitatseason = bind_rows(list_pairshabitatseason)
df_pairsdielseason = bind_rows(list_pairsdielseason)

write_xlsx(as.data.frame(df_glmer), path = 'analysis/model/glmmTMB_summary.xlsx')
write_xlsx(as.data.frame(df_emmeans), path = 'analysis/model/emmeans_all_interactions.xlsx')
write_xlsx(as.data.frame(df_pairs), path = 'analysis/model/pairwise_tukey.xlsx')

write_xlsx(as.data.frame(df_emmeanshabitatdiel), path = 'analysis/model/emmeans_habitat_diel.xlsx')
write_xlsx(as.data.frame(df_emmeanshabitatseason), path = 'analysis/model/emmeans_habitat_season.xlsx')
write_xlsx(as.data.frame(df_emmeansdielseason), path = 'analysis/model/emmeans_diel_season.xlsx')

write_xlsx(as.data.frame(df_pairshabitatdiel), path = 'analysis/model/pairs_habitat_diel.xlsx')
write_xlsx(as.data.frame(df_pairshabitatseason), path = 'analysis/model/pairs_habitat_season.xlsx')
write_xlsx(as.data.frame(df_pairsdielseason), path = 'analysis/model/pairs_diel_season.xlsx')
beepr::beep(sound = 3)

#---------------CLEVELAND_DOTPLOTS--------------------------

#Create df in long format
long_df <- standard_df %>%
  select(Device_ID, Cluster, all_of(indices)) %>%
  pivot_longer(
    cols = all_of(indices),
    names_to = "index",
    values_to = "value"
  )


for (idx in indices) {
  df_idx <- long_df %>% filter(index == idx)
  
  cleve_plot <- ggplot(df_idx, aes(x = value, y = Device_ID)) +
    geom_point(size = 1.5, alpha = 0.6) +
    theme_light(base_size = 14) +
    labs(y = "Device", x = "Value")
  
  ggsave(filename = paste0("plots/Cleveland_outputs/", idx, "_cleveland_dotplot.png"),
         plot = cleve_plot, width = 12, height = 8)
}

#-------------HISTOGRAMS------------------------------

for (idx in indices) {
  df_idx <- long_df %>% filter(index == idx)
  
  hist_plot <- ggplot(df_idx, aes(x = value)) +
    geom_histogram(bins = 30, fill = "grey", color = "white", alpha = 0.7) +
    facet_wrap(~Device_ID, scales = "free_y") +
    theme_light(base_size = 14) +
    labs(x = "Value", y = "Count")
  
  ggsave(filename = paste0("plots/Histograms_outputs/", idx, "_histograms.png"),
         plot = hist_plot, width = 12, height = 8)
}

#-----------------CORRELATION_MATRIX------------------------

num_df <- merged_df[,indices]
cor_matrix <- cor(num_df, use = "pairwise.complete.obs", method = "pearson")
cor_melt <- melt(cor_matrix)

cor_all = ggplot(cor_melt, aes(x=Var1, y=Var2, fill = value)) +
  geom_tile(color = 'white') +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="Correlation") +
  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3) +
  xlab('')+
  ylab('')+
  theme_light() +
  theme(panel.border = element_blank(),
        panel.grid = element_blank()) +
  coord_fixed()
                 
cor_all
ggsave(filename = paste0('plots/correlation_matrix_all_data.png'), width = 10, height = 10)



#-----------Means per index-------------

# Summarise: mean per 15-min bin across all days
summary_df <- merged_df %>%
  pivot_longer(cols = all_of(indices), names_to = "Index", values_to = "Value") %>%
  group_by(Index, Time15_ofday, Habitat, season) %>%
  summarise(
    median_val = median(Value, na.rm = TRUE),
    mean_val = mean(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Time15 = as.POSIXct(Time15_ofday, format = "%H:%M:%S"))

summary_df_diel <- merged_df %>%
  pivot_longer(cols = all_of(indices), names_to = "Index", values_to = "Value") %>%
  group_by(Index, Time15_ofday) %>%
  summarise(
    median_val = median(Value, na.rm = TRUE),
    mean_val = mean(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Time15 = as.POSIXct(Time15_ofday, format = "%H:%M:%S"))

summary_df_diel_season <- merged_df %>%
  pivot_longer(cols = all_of(indices), names_to = "Index", values_to = "Value") %>%
  group_by(Index, Time15_ofday, season) %>%
  summarise(
    median_val = median(Value, na.rm = TRUE),
    mean_val = mean(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Time15 = as.POSIXct(Time15_ofday, format = "%H:%M:%S"))

summary_df_diel_habitat_season <- merged_df %>%
  pivot_longer(cols = all_of(indices), names_to = "Index", values_to = "Value") %>%
  group_by(Index, Time15_ofday, Habitat, season) %>%
  summarise(
    median_val = median(Value, na.rm = TRUE),
    mean_val = mean(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Time15 = as.POSIXct(Time15_ofday, format = "%H:%M:%S"))


combined_diel = ggplot(summary_df_diel, aes(x = Time15, y = mean_val))+
  geom_line(linewidth = 0.7) +
  facet_wrap(~Index, scales = 'free_y', ncol = 4)+
  geom_point(size = 0.4)+
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "6 hour") +
  labs(x = "Time of day", y = "Mean value")+
  theme_light(base_size = 20)+
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black",  size = 16, face = 'bold.italic'),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

ggsave(combined_diel, filename = 'plots/means_all_indices_diurnal.png', width = 18, height = 10)

combined_habitat = ggplot(summary_df_diel_habitat, aes(x = Time15, y = mean_val, color = Habitat))+
  geom_line(linewidth = 1)+
  facet_wrap(~Index, scales = 'free_y')+
  labs(x = "Time of day", y = "Mean index value", color = 'Habitat')+
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "4 hour") +
  geom_point(size = 0.4)+
  theme_light(base_size = 20)+
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black",  size = 16, face = 'bold.italic'),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

ggsave(combined_habitat, filename = 'plots/means_all_indices_habitat.png', width = 18, height = 10)

combined_diel_season = ggplot(summary_df_diel_season, aes(x = Time15, y = mean_val, color = season))+
  geom_line(linewidth = 0.7)+
  facet_wrap(~Index, scales = 'free_y', ncol = 4)+
  labs(x = "Time of day", y = "Mean value", color = 'Season')+
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "6 hour") +
  geom_point(size = 0.4)+
  theme_light(base_size = 20)+
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black",  size = 16, face = 'bold.italic'),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

ggsave(combined_diel_season, filename = 'plots/means_all_indices_diurnal_seasonal.png', width = 12, height = 8)

combined_season_habitat = ggplot(summary_df_diel_habitat_season, aes(x = Time15, y = mean_val, color = Habitat))+
  geom_line(linewidth = 1)+
  facet_wrap(Index~season, scales = 'free_y', ncol = 2)+
  scale_color_manual(values = habitat_colours)+
  labs(x = "Time of day", y = "Mean index value")+
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "6 hour") +
  geom_point(size = 0.4)+
  theme_light(base_size = 20)+
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black",  size = 16, face = 'bold.italic'),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none')
ggsave(combined_season_habitat, filename = 'plots/means_indices_habitat_seasonal.png', width = 12, height = 20)


#--------------------Boxplots------------------------
box_full = ggplot(long_df, aes(x = Cluster, y = value, fill = Cluster)) +
  geom_boxplot() +
  facet_wrap(~index, scales = "free_y", ncol = 3) +
  theme_light() +
  theme(axis.text.x = element_blank()) +
  labs(x = "",
       y = "Value")
ggsave(filename = paste0('plots/boxplot_indices_and_clusters.png'), plot = box_full, width = 18, height = 12)

#---------------Violin-------------------------------

vio_full = ggplot(long_df, aes(x = Cluster, y = value, fill = Cluster)) +
  geom_violin()+
  facet_wrap(~index, scales = "free_y") +
  theme_light() +
  theme(axis.text.x = element_blank())
vio_full
ggsave(filename = paste0('plots/violinplot_indices_and_clusters.png'), plot = vio_full, width = 18, height = 12)

#----------Means-by-day-----------------

mean_by_day <- merged_df %>%
  group_by(day, Habitat, Cluster) %>%
  summarise(
    mean_AEI = mean(AEI, na.rm = TRUE),
    mean_NDSI = mean(NDSI, na.rm = TRUE),
    mean_ADI = mean(ADI, na.rm = TRUE),
    mean_H = mean(H, na.rm = TRUE),
    mean_BI = mean(BI, na.rm = TRUE),
    mean_ACI = mean(ACI, na.rm = TRUE),
    mean_TQ = mean(TQ, na.rm = TRUE),
    mean_CENT = mean(CENT, na.rm = TRUE),
    mean_nROI = mean(nROI, na.rm = TRUE),
    mean_aROI = mean(aROI, na.rm = TRUE),
    mean_TFSD = mean(TFSD, na.rm = TRUE),
    .groups = "drop"
  )

# Pivot longer to make indices a single column
long_sum <- mean_by_day %>%
  pivot_longer(
    cols = starts_with("mean_"),
    names_to = "Index",
    values_to = "Value"
  )


long_sum$Index <- gsub("mean_", "", long_sum$Index)

# Faceted plot
mean_day = ggplot(long_sum, aes(x = day, y = Value, color = Habitat)) +
  facet_wrap(~ Index, scales = "free_y") +
  theme_light() +
  geom_smooth()+
  labs(x = "Date",
    y = "Mean Value"
  )

ggsave(filename = paste0('plots/means_by_day.png'), plot = mean_day, width = 25, height = 12)

for (idx in indices) {
  
  p <- ggplot(filter(long_sum, Index == idx), 
              aes(x = day, y = Value, color = Habitat, linetype = Cluster)) +
    geom_line(linewidth = 0.6) +
    theme_light() +
    labs(x = "Date", y = idx)
  
  # Save each plot
  ggsave(filename = paste0("plots/individual/mean_", idx, "_by_day.png"), plot = p,width = 10, height = 6)
}


#-------------Heatmap_downtime-------------------------------
samples_per_day <- merged_df %>%
  group_by(day, Device_ID, Cluster) %>%
  summarise(sample_count = n(), .groups = "drop")

samples_per_day <- samples_per_day %>%
  arrange(Cluster, Device_ID) %>%
  mutate(Device_ID = factor(Device_ID, levels = unique(Device_ID)))

Downtime <- ggplot(samples_per_day, aes(x = day, y = Device_ID, fill = sample_count)) +
  geom_tile() +
  facet_grid(rows = vars(Cluster), scales = "free_y", space = "free_y") +  # label clusters on y-axis
  theme_light() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(angle = 0, face = "bold", color = 'black'),      
    strip.background.y = element_blank(),
    panel.spacing.y = unit(0.1, "lines")
  ) +
  scale_fill_viridis_c(option = 'C', name = 'Samples per day', direction = -1)+
  labs(x = "Day", y = "Device (ordered by Cluster)")

ggsave('plots/Downtime_heatmap.png', plot = Downtime, width = 14, height = 8)



#------------Sampling_effort-----------------------
sampling_effort_summary <- merged_df %>%
  group_by(Habitat, diel, season) %>%
  summarise(
    Sampling_Effort = n(),
    .groups = 'drop' # Recommended to drop grouping structure after summarising
  )
write_xlsx(sampling_effort_summary, path = 'analysis/sampling_effort.xlsx')

#---------emmeans_plots-------------
diel_interest = c('BI', 'ACI', 'NDSI')

df_emmeansdielseason <- read_xlsx("analysis/model/emmeans_diel_season.xlsx") %>%
  filter(index %in% indices)

df_emmeanshabitatdiel <- read_xlsx("analysis/model/emmeans_habitat_diel.xlsx") %>%
  filter(index %in% diel_interest)

df_emmeanshabitatseason <- read_xlsx("analysis/model/emmeans_habitat_season.xlsx") %>%
  filter(index %in% indices)

# Diel × Season
plot_emmeans_dielseason <- ggplot(df_emmeansdielseason, aes(x = diel, y = emmean, color = season)) +
  geom_point(position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                position = position_dodge(width = 0.5),
                width = 0.2) +
  facet_wrap(~index, scales = "free_y", nrow = 2) +
  theme_light(base_size = 20) +
  labs(x = "Diel part",
       y = "Estimated marginal means (with CL)",
       color = "Season")+
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black",  size = 16, face = 'bold.italic'),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

ggsave(plot = plot_emmeans_dielseason,
       filename = "plots/emmeans/diel_season.png",
       height = 10, width = 15)

# Diel × Habitat
plot_emmeans_dielhab <- ggplot(df_emmeanshabitatdiel, aes(x = Habitat, y = emmean, color = diel)) +
  geom_point(position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                position = position_dodge(width = 0.5),
                width = 0.2) +
  facet_wrap(~index, scales = "free_y") +
  theme_light() +
  labs(x = "Habitat",
       y = "Estimated marginal means (with CL)",
       color = "Diel part")+
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black",  size = 16, face = 'bold.italic'),
        axis.ticks = element_blank())

ggsave(plot = plot_emmeans_dielhab,
       filename = "plots/emmeans/diel_habitat.png",
       height = 4, width = 10)

# Habitat × Season
plot_emmeans_habseason <- ggplot(df_emmeanshabitatseason, aes(x = Habitat, y = emmean, color = Habitat, shape = season)) +
  geom_point(position = position_dodge(0.5), size = 5) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                position = position_dodge(width = 0.5),
                width = 0.2)  +
  scale_color_manual(values = habitat_colours)+
  facet_wrap(~index, scales = "free_y", ncol = 1) +
  theme_light(base_size = 20) +
  labs(x = "Habitat",
       y = "Estimated marginal means (with CL)",
       shape = "Season")+
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black",  size = 16, face = 'bold.italic'),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

ggsave(plot = plot_emmeans_habseason,
       filename = "plots/emmeans/habitat_season.png",
       height = 18, width = 6)


#----------------patchwork-------------------
patchedplot = combined_season_habitat + 
  plot_emmeans_habseason +
  plot_layout(widths = c(2.5,1)) +
  plot_annotation(tag_levels = 'A')

ggsave(patchedplot, filename = 'plots/composite_habitat_season.png', height = 20, width = 16)


dielpatch = combined_diel / combined_diel_season / plot_emmeans_dielseason+
  plot_annotation(tag_levels = 'A')

ggsave(dielpatch, filename = 'plots/composite_diel_season.png', height = 17, width = 17)




