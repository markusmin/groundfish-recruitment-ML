features <- read.csv(here::here("ML_inputs", "features_2025-01-14.csv"))
features %>% 
  dplyr::rename(Year = X) %>% 
  pivot_longer(cols = -Year, names_to = "covariate", values_to = "value") -> features_long

features_ts_plot <- ggplot(subset(features_long, !(covariate %in% c("BEUTI_central", "BEUTI_north"))), aes(x = Year, y = value)) +
  geom_line() +
  facet_wrap(~covariate) +
  theme_minimal()

ggsave(here::here("figures", "feature_extraction", "features_ts_plot.png"), features_ts_plot, height = 8, width = 12)
