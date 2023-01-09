library(promor)
library(VennDiagram)
library(ggplot2)


################################################################################

# Perseus data parsing

################################################################################

# Upload Perseus results
perseus_results <- read.csv("./Perseus_MinProb_noNormalization_DEresults.txt",
  sep = "\t"
)


# Reduce the data frame to significant hits and limit the data frame to
#only those columns we need
de_perseus <- perseus_results[perseus_results$H_vs_L_Significant == "+", c(
  "Majority.protein.IDs",
  "H_vs_L_P.Value",
  "H_vs_L_logFC"
)]


# Add a Protein.IDs column to promor results by extracting the first protein
#from majority_protein_ids
de_perseus$Protein.IDs <- sapply(strsplit(as.character
                                          (de_perseus$Majority.protein.IDs),
                                          ";"), "[", 1)

# remove maj prot id column
de_perseus <- subset(de_perseus, select = -Majority.protein.IDs)

# Add a new column with the name of the method used
de_perseus$method <- "Perseus"

# Let's give both data frames similar column names
colnames(de_perseus) <- c("p_val", "log_fc", "protein", "method")

# Make a list object to build a venn diagram
de_perseus_prot <- de_perseus$protein


################################################################################

# promor analysis

################################################################################

# Create a raw_df object with proteinGroups.txt and exp_design file
raw_df <- create_df(
  prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg1.txt",
  exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed1.txt",
  uniq_pep = 2
)


# Filter out proteins with higher than 40% missing data in either one group
#(in other words - requires 60% valid data in each group to retain the protein)
raw_df_filt <- filterbygroup_na(raw_df, set_na = 0.40, filter_condition = "either")


# Impute missing data
imp_df <- impute_na(raw_df_filt, seed = 327, method = "minProb")

# Find DE proteins
fit_df <- find_dep(imp_df)

# Save results from all DE proteins
fit_df <- find_dep(imp_df,
                   n_top = 1294,
                   save_tophits = TRUE,
                   save_output = TRUE,
                   file_path = "./")

# Upload the TopHits
de_promor <- read.csv("./TopHits.txt", sep = "\t")

# Add a Protein.IDs column to promor results by extracting the first protein
#from majority_protein_ids
de_promor$Protein.IDs <- sapply(strsplit
                                (as.character(de_promor$majority_protein_id),
                                  ";"), "[", 1)

# Extract only those columns we need from de_promor
de_promor <- de_promor[, c("Protein.IDs", "P.Value", "logFC")]

# Add a new column with the method information
de_promor$method <- "promor"

# Let's give both data frames similar column names
colnames(de_promor) <- c("protein", "p_val", "log_fc", "method")

# Make a list object to build a venn diagram
de_promor_prot <- de_promor$protein


################################################################################

# Comparsion - Number of DE proteins

################################################################################


venn.diagram(list("promor" = de_promor_prot, "Perseus" = de_perseus_prot),
  fill = c("#17456B", "#ACF0F2"),
  alpha = c(0.5, 0.5),
  resolution = 400,
  lwd = 5,
  filename = "./venn_diagram.tiff",
  scaled = TRUE,
  ext.pos = 0,
  ext.percent = 0.5,
  fontface = "bold",
  ext.line.lwd = 3
)


################################################################################

# Comparsion - log FC

################################################################################

# combine data from both dataframes into one
df_all <- merge(de_promor, de_perseus, by = "protein")

# Convert non-numeric values to numeric
df_all$log_fc.y <- as.numeric(df_all$log_fc.y)
df_all$p_val.y <- as.numeric(df_all$p_val.y)
attach(df_all)

#calculate and annotate pearson correlation
grob1 <- grobTree(textGrob(paste("Pearson Correlation : ",
                                 round(cor(log_fc.x, log_fc.y), 4)),
  x = 0.5, y = 0.97, hjust = 0,
  gp = gpar(
    col = "black",
    fontsize = 15,
    fontface = "bold"
  )
))

#Make the plot
ggplot(df_all, aes(x = log_fc.x, y = log_fc.y)) +
  geom_point(size = 10, shape = 20, alpha = .4, col = "#17456B") +
  # ggtitle(bquote('promor vs Perseus - log' [2]~ 'fold change')) +
  geom_smooth(method = lm, se = FALSE, lwd = 1) +
  xlab(bquote("promor log"[2] ~ "fold-change")) +
  ylab(bquote("Perseus log"[2] ~ "fold-change")) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = "black", linewidth = 2),
    axis.line.x = element_line(color = "black"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 16)
  ) +
  annotation_custom(grob1)


################################################################################

# Comparsion - p-value

################################################################################

#calculate and annotate pearson correlation
grob2 <- grobTree(textGrob(paste("Pearson Correlation : ",
                                 round(cor(log(p_val.x), log(p_val.y)), 4)),
  x = 0.5, y = 0.97, hjust = 0,
  gp = gpar(
    col = "black",
    fontsize = 15,
    fontface = "bold"
  )
))

#Make the plot
ggplot(df_all, aes(x = log(p_val.x), y = log(p_val.y))) +
  geom_point(size = 10, shape = 20, alpha = .4, col = "#17456B") +
  # ggtitle(bquote('promor vs Perseus - log' [10]~ 'p-value')) +
  geom_smooth(method = lm, se = FALSE, lwd = .5) +
  xlab(bquote("promor log"[10] ~ "p-value")) +
  ylab(bquote("Perseus log"[10] ~ "p-value")) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = "black", linewidth = 2),
    axis.line.x = element_line(color = "black"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 16)
  ) +
  annotation_custom(grob2)

