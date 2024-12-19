library(ggplot2)


# Exclude the "group" column for PCA (only numeric columns should be used)
pca_data <- raw_data[, c(-1,-2)]  # Remove the 'group' column

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)  # Scale the data for PCA

# Extract PCA results for the first two principal components
original_group<-raw_data$GROUP


modified_group <- ifelse(
  grepl("^C", original_group), "control",original_group)  # Change values starting with 'HA' to 'HA', others remain the same


pca_df <- as.data.frame(pca_result$x)  # PCA scores
pca_df$GROUP <- modified_group  # Add group information for coloring

# Visualize PCA with ellipses using ggplot2
ggplot(pca_df, aes(x = PC1, y = PC2, color = GROUP)) +
  geom_point(size = 3) +  # Scatter plot
  stat_ellipse(aes(fill = GROUP), type = "norm", alpha = 0.2, geom = "polygon") +  # Ellipses
  labs(
    title = "PCA Analysis of identified SNPs with Groups ",
    x = "Principal Component 1",
    y = "Principal Component 2"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")




