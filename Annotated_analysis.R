setwd('C:/Users/stijn/Documents/Drive/Earth Sciences/Thesis')
annotated_df = read.csv('Subsampling/annotated_audio_files_extra.csv', sep = ';', header = TRUE)

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(scales)
library(patchwork)
library(tidyverse)

#------------dominance proportion plot--------------
dominance_summary <- annotated_df %>%
  group_by(habitat, season, diel, dominance) %>%
  summarise(n = n(), .groups = "drop")

dominance_colors = c("#39426D", "#44CC44", "#64A8C3", "#93856D")

proportion_dominance = ggplot(dominance_summary, aes(x = diel, y = n, fill = dominance)) +
  geom_bar(stat = "identity", position = "fill") +  
  facet_wrap(habitat ~ season, ncol = 2) +
  scale_fill_manual(values = dominance_colors)+
  labs(y = "Proportion",
    x = "Diel part",
    fill = "Dominant source") +
  theme_light(base_size = 16)+
  theme(strip.background = element_blank(),
        strip.text = element_text(color = "black",  size = 16, face = 'bold.italic'),
        axis.ticks = element_blank())
proportion_dominance

ggsave(proportion_dominance, filename = 'plots/proportion_sound_dominance_habitat.png', height = 12, width = 10)


#-------------Biophony source plot------------------

biophony_colors = c('#99C68E', '#7285A5', '#a87a46')

biophony_n <- annotated_df %>%
  mutate(any_biophony = ifelse(avian + insect + amphibian > 0, 1, 0)) %>%
  group_by(habitat, season, diel) %>%
  summarise(n_files = sum(any_biophony), .groups = "drop")

biophony_summary_prop <- annotated_df %>%
  pivot_longer(
    cols = c(avian, insect, amphibian),
    names_to = "source",
    values_to = "count"
  ) %>%
  group_by(habitat, season, source, diel) %>%
  summarise(total_count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  group_by(habitat, season, diel) %>%
  mutate(proportion = total_count / sum(total_count))

biophony_source <- ggplot(biophony_summary_prop, aes(x = diel, y = proportion, fill = source)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(data = biophony_n,
            aes(x = diel, y = -0.05, label = paste0("n = ", n_files)),
            inherit.aes = FALSE,
            size = 4) +
  scale_fill_manual(values = biophony_colors) +
  facet_wrap(habitat~season, ncol=2) +
  theme_light(base_size = 16) +
  labs(x = "Diel period",
    y = "Proportion",
    fill = "Type of biophony") +
  theme(strip.background = element_blank(), 
        strip.text = element_text(color = "black", size = 16, face = 'bold.italic'),
        legend.position = "right",
        axis.ticks = element_blank())
biophony_source
ggsave(biophony_source, filename = 'plots/biophony_source.png', height = 12, width = 10)

#------------avian source plot-----------------

avian_summary <- annotated_df %>%
  select(habitat, waterbirds, birdsong, other) %>%
  pivot_longer(
    cols = c(waterbirds, birdsong, other),
    names_to = "source",
    values_to = "count"
  ) %>%
  group_by(habitat, source) %>%
  summarise(total_count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  group_by(habitat) %>%
  mutate(proportion = total_count / sum(total_count))

avian_colors <- c(
  "waterbirds" = "#3c5166",
  "birdsong" = "#943E00",
  "other" = "#100C08"
)

avian_plot <- ggplot(avian_summary, aes(x = habitat, y = proportion, fill = source)) +
  geom_bar(stat = "identity", position = "stack") +
#  geom_text(
#    data = subset(avian_summary, proportion > 0),
#    aes(label = percent(proportion, accuracy = 1)),
#    position = position_stack(vjust = 0.5),
#    color = "white",
#    size = 4
#  ) +
  scale_fill_manual(values = avian_colors) +
  guides(fill = guide_legend(reverse = TRUE))+
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Habitat",
    y = "Proportion",
    fill = "Avian Source"
  ) +
  theme_light(base_size = 16) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(color = "black",  size = 16, face = 'bold.italic'),
    legend.position = "right",
    axis.ticks = element_blank()
  )

avian_plot
ggsave(avian_plot, filename = 'plots/avian_biophony_source.png', height = 8, width = 12)


#------------------geophony presence-----------------
habitat_colours = c('#33A02C', '#FF7F00', '#1F78B4')

geophony_summary <- annotated_df %>%
  mutate(geophony_presence = ifelse(geophony != "0" & geophony != 0 & !is.na(geophony), 1, 0)) %>%
  group_by(habitat) %>%
  summarise(
    total = n(),
    presence_count = sum(geophony_presence, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(proportion = presence_count / total)

geophony_plot <- ggplot(geophony_summary, aes(x = habitat, y = proportion, fill = habitat)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = habitat_colours)+
  geom_text(
    aes(label = percent(proportion, accuracy = 1)),
    position = position_dodge(width = 0.9),
    vjust = -0.3,
    size = 4
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    x = "Habitat",
    y = "Recordings with geophony (%)",
    fill = "Habitat"
  ) +
  theme_light(base_size = 16) +
  theme(
    legend.position = "right",
    axis.ticks = element_blank()
  )

geophony_plot

ggsave(geophony_plot, filename = 'plots/geophony_presence.png', height = 6, width = 12)


#---------------------anthropophony presence-------------------------
anthropophony_summary <- annotated_df %>%
  mutate(anthropophony_presence = ifelse(anthropophony != "0" & anthropophony != 0 & !is.na(anthropophony), 1, 0)) %>%
  group_by(habitat) %>%
  summarise(
    total = n(),
    presence_count = sum(anthropophony_presence, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(proportion = presence_count / total)


anthropophony_plot <- ggplot(anthropophony_summary, aes(x = habitat, y = proportion, fill = habitat)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(
    aes(label = percent(proportion, accuracy = 1)),
    position = position_dodge(width = 0.9),
    vjust = -0.3,
    size = 4
  ) +
  scale_fill_manual(values = habitat_colours)+
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    x = "Habitat",
    y = "Recordings with anthropophony (%)",
    fill = "Habitat"
  ) +
  theme_light(base_size = 16) +
  theme(
    axis.ticks = element_blank(),
    legend.position = 'none'
  )

anthropophony_plot

ggsave(anthropophony_plot, filename = 'plots/anthropophony_presence.png', height = 6, width = 12)


#------------composite figure--------------

row1 <- (biophony_source)

row2 <- (avian_plot)

row3 <- (anthropophony_plot | geophony_plot)

# Combine all rows vertically
final_plot <- (row1) / (row2)+ 
  plot_layout(
    heights = c(3, 1)) + 
  plot_annotation(
    tag_levels = 'A')

final_plot

ggsave(final_plot, 
       filename = 'plots/composite_figure_annotation_A4.png', 
       height = 17, 
       width = 12)


appendix_plot = (row3) +
  plot_layout(
    heights = c(1, 1)) + 
  plot_annotation(
    tag_levels = 'A')

ggsave(appendix_plot, 
       filename = 'plots/composite_figure_geoanthro_A4.png', 
       height = 8, 
       width = 12,
       dpi = 500)
