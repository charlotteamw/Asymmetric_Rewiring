library(ggplot2)
library(tidyverse)
library(patchwork)
library(knitr)
library(kableExtra)

df_coupling <- read.csv("/Users/charlotteward/Documents/Global Change Rewiring Review/Asymmetric_Rewiring/literature_data/LiteratureTableS1.csv")

df_coupling$result <- as.factor(df_coupling$Result..Direction.)
df_coupling$ecosystem <- as.factor(df_coupling$Ecosystem)
df_coupling$stressor <- as.factor(df_coupling$Stressor.Category)
df_coupling$mechanism <- as.factor(df_coupling$Mechanism.of.Shift)


# mechanism
mechanism_df <- df_coupling %>%
  group_by(mechanism, result) %>%
  reframe(count = n()) %>%
  na.omit()

mechanism_df$result <- factor(mechanism_df$result, levels = c("Increase", "Decrease"))

mechanism_plot <- ggplot(mechanism_df, aes(x = mechanism, y = count, fill = result)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = c("Decrease" = "steelblue4", "Increase" = "skyblue3")) +
  ylim(0, 40) +
  theme_classic()+ 
  labs(fill = "Result", y = "Count", x = "Mechanism")+
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16)) 

mechanism_plot


#ecosystem

ecosystem_df <- df_coupling %>%
  count(ecosystem, result) %>%
  mutate(result = factor(result, levels = c("No change", "Increase", "Decrease")))

ecosystem_totals <- ecosystem_df %>%
  group_by(ecosystem) %>%
  summarise(total_n = sum(n))

ecosystem_df <- ecosystem_df %>%
  left_join(ecosystem_totals, by = "ecosystem")

ecosystem_plot <- ggplot(ecosystem_df, aes(x = reorder(ecosystem, -total_n), y = n, fill = result)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Decrease" = "steelblue4", "Increase" = "skyblue3", "No change" = "grey71")) +
  ylim(0, 20) +
  theme_classic() +
  labs(fill = "Result", y = "Count", x = "Ecosystem") +
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Rotate x-axis labels to prevent overlap
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16))

ecosystem_plot

# stressor

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
  ylim(0, 15) +
  theme_classic() +
  labs(fill = "Result", y = "Count", x = "Stressor") +
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Rotate x-axis labels to prevent overlap
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16))

stressor_plot

# ecosystem & stressor plot

ecosystem_plot <- ecosystem_plot + theme(plot.margin = margin(r = 10, unit = "mm")) # Add right margin to ecosystem plot
stressor_plot <- stressor_plot + theme(plot.margin = margin(l = 10, unit = "mm")) # Add left margin to stressor plot

combined_plot <- ecosystem_plot + stressor_plot + plot_layout(guides = "collect") 

combined_plot



# Calculate the counts for each combination of ecosystem, stressor, and result
count_df <- df_coupling %>%
  group_by(ecosystem, stressor, result) %>%
  summarize(count = n(), .groups = 'drop')

# Calculate the order for ecosystem and stressor by the number of studies
ecosystem_order <- count_df %>%
  group_by(ecosystem) %>%
  summarize(total_count = sum(count)) %>%
  arrange(desc(total_count)) %>%
  mutate(ecosystem = factor(ecosystem, levels = rev(ecosystem)))

stressor_order <- count_df %>%
  group_by(stressor) %>%
  summarize(total_count = sum(count)) %>%
  arrange(desc(total_count)) %>%
  mutate(stressor = factor(stressor, levels = rev(stressor)))

# Join the order information back to the original count_df
count_df <- count_df %>%
  left_join(ecosystem_order, by = "ecosystem") %>%
  left_join(stressor_order, by = "stressor") %>%
  arrange(desc(count))  # Arrange in descending order to plot larger bubbles first

# Create the bubble plot with the ordered factors
bubble_plot <- ggplot(count_df, aes(x = ecosystem, y = stressor, size = count, color = result)) +
  geom_point() +  # removed alpha for non-translucent bubbles
  scale_size(range = c(4, 15)) +  # Adjust the size range as needed
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(angle = 45, vjust = 1)) +
  labs(title = "Bubble Plot of Ecosystem, Stressor, and Results",
       x = "Ecosystem Type",
       y = "Anthropogenic Pressure",
       size = "Count",
       color = "Result") +
  scale_color_manual(values = c("Decrease" = "steelblue4", "Increase" = "skyblue3", "No change" = "grey71")) +
  scale_x_discrete(limits = ecosystem_order$ecosystem) +
  scale_y_discrete(limits = stressor_order$stressor) +
  guides(size = guide_legend(reverse = TRUE))  # Reverse the legend for bubble sizes

# Print the bubble plot
print(bubble_plot)

