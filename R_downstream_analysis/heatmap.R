# Load necessary libraries
library(pheatmap)

# Step 1: Load the data
data <- read.csv("wide_format_data_unique.csv", row.names = 1)  # Use the first column as row names

# Step 2: Extract group information
group_info <- data.frame(Group = data$GROUP)  # Create a data frame for group annotations
data <- data[, c(-1,-2)]  # Remove the 'GROUP' column from the data to keep only numeric values
rownames(group_info) <- rownames(data)  # Set row names to match the data

# Step 3: Replace zeros with NA
  # Replace all 0 values with NA to avoid clutter in the heatmap


#unique_groups <- unique(group_info$Group) 
unique_groups <- c("P1" ,"C_1_70","C_2_70","C_3_70","C_4_70","C_5_70", "C_6_70","HA_1_70","HA_2_70","HA_2_30","HA_3_70","HA_5_70","HA_6_70")
# Get all unique group names
#group_colors <- setNames(
 # colorRampPalette(c("cyan", "white","yellow"))(length(unique_groups)), 
  #unique_groups)
group_colors <- setNames(c("#F0E442",rep("#0072B2",6),rep("#D55E00",6))
  , unique_groups
)

#group_colors<-c("green",rep("#00FFFF",6),rep("yellow",6))


# Step 4: Draw the Heatmap
pheatmap(
  data,
  annotation_row = group_info,       # Add group annotation on the rows
  scale = "row",                     # Scale data by row
  color = colorRampPalette(c("blue", "white", "red"))(50),  # Color gradient
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_colors = list(Group = group_colors),  # Adjust group colors
  show_rownames = FALSE,
  show_colnames = TRUE,
  angle_col = 45  
)
