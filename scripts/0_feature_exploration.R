# FEATURE EXPLORATION FOR COVID STRESS PROJECT
# This script contains code for exploring the potential feature space of 
# a new project through PCoA, regression, and various plots. 

rm(list=ls())
library(tidyverse)
library(vegan)
setwd('~/covidStress/data/metadata')
df <- as.data.frame(read_csv('HYPHYVirtual_StressCandidates.csv'))
rownames(df) <- df$record_id # keep an identifier column

#### LOAD IN DATA ####
# self_pa_v2 has values 1-9 and then "A" "B". Makes no sense. Fix.
df <- df %>% mutate(self_pa_v2= ifelse(self_pa_v2 == "A", 10, 
                                       ifelse(self_pa_v2 == "B", 11, self_pa_v2)))
df$self_pa_v2 <- as.numeric(df$self_pa_v2)

# create a column for bmi categories as factors. 
# Categories by CDC https://www.cdc.gov/bmi/adult-calculator/bmi-categories.html
df <- df %>% mutate(bmi_cat = ifelse(bmi < 18.5, 0, 
                    ifelse(bmi < 25, 1, ifelse(bmi < 30, 2, ifelse(bmi < 35, 3,
                    ifelse(bmi < 40, 4, ifelse(bmi >= 40, 5, 99))))))) # no 99's. Good!

# general data cleanup  
df <- df %>% select(-(1:8)) %>% # first 8 are id columns
  select(where(~ mean(is.na(.)) <= 0.1)) %>% # 23 columns with > 10% na entries
  select(where(is.numeric)) %>% # 4 columns non-numeric
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0.5 * min(., na.rm = TRUE), .))) # coerce NA values

#### Run PCoA. Extract coordinates. ####
pcoa <- capscale(df~1,distance = 'bray', na.rm=TRUE)
pcoa_coords <- as.data.frame(pcoa$CA$u[, 1:2])  # Use first two axes
colnames(pcoa_coords) <- c("PC1", "PC2")

#### MAKE SINGLE FEATURE PCOA PLOTS ####
# Create advanced axis labels with explained variance percentages
xAxisLab <- paste0("PCoA1 (", format(100 * summary(pcoa)$cont$importance[2], digits = 3), "%)")
yAxisLab <- paste0("PCoA2 (", format(100 * summary(pcoa)$cont$importance[5], digits = 3), "%)")

# Prepare for plotting
output_pdf <- "PCoA_individual_feature_plots.pdf"
pdf(file.path('../../plots',output_pdf), width = 8, height = 6)

# Loop through each feature
for (feature in colnames(df)) {
  # Define colors based on the feature values
  feature_colors <- factor(df[[feature]])  # Convert feature to factor
  num_levels <- length(levels(feature_colors))  # Count unique levels
  color_map <- rainbow(num_levels)  # Generate unique colors for each level
  colors <- color_map[as.numeric(feature_colors)]  # Map colors to feat
  
  # Create the ordiplot with custom x and y lims
  pcoa_ordiplot <- ordiplot(pcoa, choices = c(1, 2), type = 'none', 
                     cex.lab = 1.2, xlab = xAxisLab, ylab = yAxisLab, 
                     main = paste("PCoA Colored by", feature))
  
  # Plot the points with custom colors and transparency
  points(pcoa_ordiplot, "sites", col = adjustcolor(colors, alpha.f = 0.2), pch = 16, cex = 2.5)
  
  # Add ellipses for confidence intervals for each level of the feature
  for (i in seq_along(levels(feature_colors))) {
    ordiellipse(pcoa_ordiplot, df[[feature]], kind = 'se', conf = 0.95, 
                lwd = 4, draw = 'lines', col = color_map[i], 
                show.groups = levels(feature_colors)[i], label = TRUE, font = 2, cex = 1)
  }
  
  # Add legend for clarity, using all levels of the feature
  legend('topleft', legend = levels(feature_colors), col = color_map, pch = 19, cex = 1.3)
}
dev.off()

#### REGRESS FEATURES AGAINST PC1 ####
# Extract the first principal component (PC1) from the PCOA result
pc1 <- pcoa$CA$u[, 1]

# Create an empty data frame to store regression results
regression_results <- data.frame(
  feature = character(0),
  adj_r_squared = numeric(0),
  p_value = numeric(0),
  stringsAsFactors = FALSE
)

# Loop through each feature (column) in the numeric data and regress against PC1
for (feature in names(df)) {
  # Skip non-numeric columns
  if (is.numeric(df[[feature]])) {
    # Run the regression: feature ~ PC1
    model <- lm(df[[feature]] ~ pc1)
    
    # Extract coefficient and p-value
    adj_r_squared_value <- summary(model)$adj.r.squared  # Corr coefficient of PC1
    p_val <- summary(model)$coefficients[2, 4]  # p-value of PC1
    
    # Store the results in the data frame
    regression_results <- regression_results %>%
      add_row(feature = feature, adj_r_squared = adj_r_squared_value, p_value = p_val)
  }
}

# Filter for significant features before correction (p-value < 0.05)
significant_results <- regression_results %>%
  filter(p_value < 0.05)

# Print the significant features and their coefficients
print(significant_results)

# Let's make a bar chart!
#write_csv(significant_results,'significant_features.csv') # will manually add categories
significant_results <- read_csv('significant_features.csv')
ggplot(significant_results, aes(x = feature, y = adj_r_squared, fill = category)) +
  geom_bar(stat = "identity") + 
  labs(y = "adjusted R squared") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#### MAKE A LOADINGS PLOT ####
# Calculate loadings for the features
loadings <- scores(pcoa, display = "species")  # Extract feature scores (loadings)

# Select the top contributing features (e.g., top 10 with highest loadings on PC1 and PC2)
top_features <- loadings[order(rowSums(abs(loadings)), decreasing = TRUE), ]
top_n <- min(nrow(top_features), 10)  # Ensure there are at least 10 features
top_features <- top_features[1:top_n, ]

# Scale loadings to fit in the plot with the PCoA scores
scaling_factor <- max(abs(pcoa_coords)) / max(abs(top_features))

# Scale the loadings uniformly by the scaling factor
loadings_scaled <- top_features * scaling_factor

# Create the plot
output_loadings_pdf <- "PCoA_Loadings_Plot.pdf"
pdf(file.path('../../plots', output_loadings_pdf), width = 8, height = 6)

plot(pcoa_coords, pch = 16, col = "darkgrey", 
     xlab = xAxisLab, ylab = yAxisLab, main = "PCoA Loadings Plot")
arrows(0, 0, loadings_scaled[, 1], loadings_scaled[, 2], 
       length = 0.1, col = "grey", lwd = 2)

nudge_x <- c(0.01, -0.02, -0.01, -0.04, 0.02, -0.01, 0.03, -0.02, 0.01, -0.03)
nudge_y <- c(0.02, 0.01, -0.03, 0.04, -0.02, 0.02, -0.01, -0.03, -0.04, -0.01)
text(loadings_scaled[, 1] + nudge_x, 
     loadings_scaled[, 2] + nudge_y, 
     labels = rownames(loadings_scaled), 
     col = "red", cex = 0.9, pos = 3)


# Add a legend for clarity
#legend("topright", legend = "Feature Loadings", col = "red", lty = 1, lwd = 2)

dev.off()
