# Example R dataframe

# Load required library
install.packages("tidyr")
install.packages("dplyr")

# Load required library
library(tidyr)
library(dplyr)

# Read the file
data <- read.table("small_raw03.txt", sep = "\t", header = TRUE)

data_unique <- data %>%
  distinct(SRA_RUN, GENE, .keep_all = TRUE)

# Transform from long to wide format
raw_data <- data_unique %>%
  pivot_wider(names_from = GENE, values_from = QUAL, values_fill = 0)

# Add row names
#rownames(raw_data) <- raw_data$SRA_RUN

# Drop IDENTIFIER column as it is now row names
#raw_data <- raw_data %>% select(-SRA_RUN)


new_dataframe <- raw_data %>%
  select(-SRA_RUN) %>%                             # Remove the SRA_RUN column
  select(-starts_with("A9F03")) %>%                # Remove columns starting with "A9F03"
  group_by(GROUP) %>%                              # Group by the GROUP column
  summarise(across(everything(), sum, na.rm = TRUE)) 

raw_data<-new_dataframe
# Inspect the final data frame
print(raw_data)


write.csv(raw_data, "wide_format_data_unique.csv", row.names = FALSE)

