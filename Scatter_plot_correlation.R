library(Hmisc)
library(tidyverse)
library(truncnorm)
library(car)
library(ggplot2)
library(ggh4x)
library(grid)

#-------- Generate random data ---------
set.seed(123)  # For reproducibility

n_participants <- 400
n_genes <- 10

sample_id <- sample(paste0("Sample_", sprintf("%03d", 1:n_participants)), n_participants, replace = FALSE)
sex <- sample(c("male", "female"), n_participants, replace = TRUE)
Age <- sample(c(40:60), n_participants, replace = TRUE)


#To make weight values  based on sex
Weight <- numeric(n_participants)

for (i in 1:n_participants) {
  if (sex[i] == "male") {
    Weight[i] <- rtruncnorm(1, a = 55, b = 130, mean = 80, sd = 15)
  } else {
    Weight[i] <- rtruncnorm(1, a = 45, b = 100, mean = 65, sd = 12)
  }
}


# Initialize matrix for gene expression
gene_data <- matrix(nrow = n_participants, ncol = n_genes)

# Generate gene expression values based on Sex
for (i in 1:n_participants) {
  if (sex[i] == "male") {
    gene_data[i, ] <- rtruncnorm(n_genes, a = 0, b = 200, mean = 60, sd = 30)
  } else {
    gene_data[i, ] <- rtruncnorm(n_genes, a = 0, b = 200, mean = 40, sd = 20)
  }
}

# Convert to dataframe and add Sex column
colnames(gene_data) <- paste0("Gene_", 1:n_genes)  # Name genes
df <- data.frame(Sampe_ID = sample_id, Sex = sex, Age = Age, 
                 Weight = Weight, gene_data)  # Make the datframe



# Calculate correlation and p-values for each sex and prepare the format
sex_list <- c("female", "male")
# To store the results
corr_p_value_df <- data.frame(sex = character(),
                         correlation = numeric(),
                         P_Value_Base = numeric(), 
                         P_Value_Exponent = numeric(), 
                         stringsAsFactors = FALSE)

for (i in sex_list) {
  # Perform correlation test
  Gene_1_cor <- cor.test(df$Gene_1[df$Sex == i],
                              df$Weight[df$Sex == i], 
                              method = "spearman")
  
  Gene_1_cor_estimate <- Gene_1_cor$estimate
  Gene_1_cor_pvalue <- formatC(Gene_1_cor$p.value, format = "e", digits = 1)
  
  # Split pvalue into base and exponent
  p_value_split <- strsplit(Gene_1_cor_pvalue, "e")[[1]]
  p_value_base <- as.numeric(p_value_split[1])
  p_value_exponent <- as.numeric(p_value_split[2])
  
  # Store results in dataframe
  corr_p_value_df <- rbind(corr_p_value_df, data.frame(Sex = i, Corr_r = Gene_1_cor_estimate,
                                           P_Value_Base = p_value_base, 
                                           P_Value_Exponent = p_value_exponent))
}




p <- ggplot(data = df, aes(x = Weight, y = Gene_1, fill = Sex)) +
  geom_point(aes(colour = Sex), size= 1.1, alpha = 0.8, shape = 19, 
             stroke = 0.4, show.legend = TRUE) +
  geom_smooth(aes(colour = Sex, fill = Sex), linewidth = 1, alpha = 0.4, 
              method="lm", show.legend = FALSE, linetype = "dashed") +
  stat_density_2d(geom = "polygon",
                  aes(fill = Sex, color = Sex), linewidth = 0.1, 
                  bins = 7, alpha =0.1, contour_var = "density") +
  scale_x_continuous(expand = c(0,0), limits = c(40, 110)) +
  scale_y_continuous(expand = c(0,0), limits = c(-5, 130)) +
  scale_color_manual(values = c("dodgerblue","orangered"),
                     labels = c("male" = "Male", "female" = "Female"))+
  scale_fill_manual(values = c("dodgerblue","orangered"))+
    theme(panel.background = element_rect(fill = "white", color = "black", size = 0.2),
          plot.margin = margin(t = 20, r = 5, b = 5, l = 5),
          panel.border = element_rect(linetype = "solid", 
                                      size = 0.1, colour = "black", fill = NA),
          panel.grid.major.x = element_line(colour = "grey90", size = 0.1,),
          panel.grid.major.y = element_line(colour = "grey90", size = 0.1,),
          panel.grid.minor = element_blank(),
          axis.line = element_line(size = 0.1, colour = "black"),
          axis.text.x = element_text(size = 6, colour = "black"),
          axis.text.y = element_text(size = 6, colour = "black"),
          axis.title.y = element_text(size = 6, hjust = 0.5, colour = "black"),
          axis.title.x = element_text(size = 6, hjust = 0.5, colour = "black"),
          plot.title = element_text(size = 6, angle = 0, hjust = 0.5),
          axis.ticks = element_line(size = 0.2, colour = "black"),
          legend.key = element_blank(),
          legend.key.size = unit(0.12, "cm"),)+ #To set the appearance 
    guides(color=guide_legend(title="Sex", override.aes = list(size = 3), direction = "vertical", 
                              position = "right", theme = theme(legend.title.position = "top",
                                                                legend.title = element_text(hjust = 0.5, size = 6),
                                                                legend.text = element_text(size=6),
                                                                legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
                                                                legend.background = element_rect(fill = NA, color = NA))), 
           fill = "none") + #To set the legend
    annotation_custom(xmin = 40, xmax = 90, ymin = 110, ymax = 130, 
                      grob = rectGrob(gp = gpar(fill = "azure", 
                                                lwd = 0.52, col = "black", alpha = 0.5))) + # To set boxes
    annotation_custom(xmin = 40, xmax = 110, ymin = 130, ymax = 150, 
                      grob = rectGrob(gp = gpar(fill = "greenyellow", 
                                                lwd = 0.52, col = "black", alpha = 1))) + # To set boxes
    annotation_custom(xmin = 75, xmax = 75, ymin = 140, ymax = 140, 
                      textGrob(expression("Weight and gene 1 correlation"), 
                               gp = gpar(col = "black",  fontsize = 7))) + # To set title
    annotation_custom(
      grob = textGrob(just = "left",
                      label = bquote(~italic("R:") ~ .(formatC(corr_p_value_df[1,2], digits = 2))~ "; "
                                     ~italic("P") ~ "value:" ~ .(corr_p_value_df[1,3]) ~ "x" ~ 10^.(corr_p_value_df[1,4])),
                      gp = gpar(col = "dodgerblue", fontsize = 6)
      ),
      xmin = 40, xmax = 40, ymin = 125, ymax = 125) + #Statistics For female
    annotation_custom(
      grob = textGrob(just = "left",
                      label = bquote(~italic("R:") ~ .(formatC(corr_p_value_df[2,2], digits = 2))~ "; "
                                     ~italic("P") ~ "value:" ~ .(corr_p_value_df[2,3]) ~ "x" ~ 10^.(corr_p_value_df[2,4])),
                      gp = gpar(col = "orangered", fontsize = 6)
      ),
      xmin = 40, xmax = 40, ymin = 115, ymax = 115) + #Statistics For male
  coord_cartesian(clip = "off") #To make the boxes visible


  ggsave("Annotated_Corr_plot.pdf", plot = p, device = "pdf",   
         width = 7, height = 5.2, units = "cm", ,bg = "transparent")
  
  