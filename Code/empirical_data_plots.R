library(ggplot2)
library(tidyverse)

df_coupling <- read.csv("path to/DatasetS1.csv")

df_coupling$result <- as.factor(df_coupling$Habitat.Coupling.Result..Direction.)
df_coupling$ecosystem <- as.factor(df_coupling$Ecosystem)
df_coupling$stressor <- as.factor(df_coupling$Stressor.Category)
df_coupling$mechanism <- as.factor(df_coupling$Mechanism.of.Shift)


# Anthropogenic Pressures (Figure 2b)

stressor_df <- df_coupling %>%
  count(stressor, result) %>%
  mutate(result = factor(result, levels = c("No change", "Increase", "Decrease")))

stressor_totals <- stressor_df %>%
  group_by(stressor) %>%
  summarise(total_n = sum(n))

stressor_df <- stressor_df %>%
  left_join(stressor_totals, by = "stressor")

stressor_plot <- ggplot(stressor_df, aes(x = reorder(stressor, -total_n), y = n, fill = result)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Decrease" = "steelblue4", "Increase" = "skyblue3", "No change" = "grey71")) +
  theme_classic() +
  labs(fill = "Result", y = "Count", x = "Stressor") +
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Rotate x-axis labels to prevent overlap
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16))

stressor_plot


# Mechanism (Fig. 2c)
mechanism_df <- df_coupling %>%
  group_by(mechanism, result) %>%
  reframe(count = n()) %>%
  na.omit()

mechanism_totals <- mechanism_df %>%
  group_by(mechanism) %>%
  summarise(total_n = sum(count))

mechanism_df <- mechanism_df %>%
  left_join(mechanism_totals, by = "mechanism") %>%
  mutate(mechanism = reorder(mechanism, -total_n))


mechanism_df$result <- factor(mechanism_df$result, levels = c("Increase", "Decrease"))

mechanism_plot <- ggplot(mechanism_df, aes(x = mechanism, y = count, fill = result)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = c("Decrease" = "steelblue4", "Increase" = "skyblue3")) +
  ylim(0, 65) +
  theme_classic()+ 
  labs(fill = "Result", y = "Count", x = "Mechanism")+
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16)) 

mechanism_plot

