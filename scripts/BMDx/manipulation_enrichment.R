
#Import the libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ComplexHeatmap)
library(circlize) 
library(writexl)
library(stringr)
library(tibble)

#Set the directory
#dir_path <- "C:/Users/carbi/Desktop/ke_enrichment/"
dir_path <- "C:/Users/Utente/Desktop/denys/results/ke_enrichment/"

ke_df <- read_excel(paste0(dir_path,"ke.xlsx"))
View(ke_df)

types = unique(ke_df$Ke_type)

#################################################################################

## HEATMAP

#Function to create the heatmap
create_heatmap_gg <- function(df) {
  #Order the experiments by time
  experiment_order <- df %>%
    distinct(Experiment, time) %>%
    arrange(time) %>%
    pull(Experiment)
  
  #For each (Ke_description, Experiment) take the minimum BMD_norm
  df_min <- df %>%
    group_by(Ke_description, Experiment) %>%
    summarise(BMD_min = min(BMD_norm, na.rm = TRUE), .groups = "drop")
  
  #Matrix for the clustering
  mat <- df_min %>%
    pivot_wider(names_from = Experiment, values_from = BMD_min) %>%
    column_to_rownames("Ke_description") %>%
    as.matrix()
  
  #Order the columns by time
  mat <- mat[, experiment_order]
  
  #Clustering or the rows
  # Sostituisci gli NA con un valore alto (qui: max*1.5) per non interrompere dist()
  mat_clust <- mat
  mat_clust[is.na(mat_clust)] <- max(mat, na.rm = TRUE) * 1.5
  
  row_order <- hclust(dist(mat_clust, method = "euclidean"), method = "complete")$order
  
  #Data for the plot
  df_plot <- as.data.frame(mat) |>
    rownames_to_column("Ke_description") |>
    pivot_longer(-Ke_description,
                 names_to = "Experiment",
                 values_to = "BMD_norm")
  
  #Factors to have the desired order
  df_plot$Experiment     <- factor(df_plot$Experiment,     levels = experiment_order)
  df_plot$Ke_description <- factor(df_plot$Ke_description, levels = rownames(mat)[row_order])
  
  #Heatmap ggplot
  plt <- ggplot(df_plot,
                aes(x = Experiment, y = Ke_description, fill = BMD_norm)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(
      colours = c("white", "skyblue", "royalblue4"),
      na.value = "grey90",
      name = "BMD_norm") +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      axis.title   = element_blank(),
      panel.grid   = element_blank()) +
    ggtitle("Heatmap BMD_norm")
  
  list(plot = plt,
       nrow = nrow(mat),
       ncol = ncol(mat))
}

for (type in types) {
  df <- ke_df %>%
    filter(Ke_type == type)
  
  ht <- create_heatmap_gg(df)      # crea la heatmap
  p  <- ht$plot                    # oggetto ggplot
  nr <- ht$nrow; nc <- ht$ncol    # n. righe/colonne
  
  # Dimensioni dinamiche (0.30″ per cella, con min/max ragionevoli)
  width  <- pmin(pmax(nc * 0.30,  8), 20)   # tra 8″ e 20″
  height <- pmin(pmax(nr * 0.30, 10), 30)   # tra 10″ e 30″
  
  ggsave(paste0(dir_path, "heatmap_", type, "_bmd_gg.pdf"), plot = p, width = width, height = height, units = "in")
}


###########################################


create_heatmap_gg <- function(df) {
  # Calcola la matrice minima (BMD_norm più basso per combinazione)
  df_min <- df %>%
    group_by(Ke_description, Experiment) %>%
    summarise(BMD_min = min(BMD_norm, na.rm = TRUE), .groups = "drop")
  
  # Matrice wide
  mat <- df_min %>%
    pivot_wider(names_from = Experiment, values_from = BMD_min) %>%
    column_to_rownames("Ke_description") %>%
    as.matrix()
  
  # Gestione NA per clustering (sostituiti temporaneamente con valore alto)
  mat_clust <- mat
  mat_clust[is.na(mat_clust)] <- max(mat, na.rm = TRUE) * 1.5
  
  # Clustering righe
  row_order <- hclust(dist(mat_clust), method = "complete")$order
  # Clustering colonne
  col_order <- hclust(dist(t(mat_clust)), method = "complete")$order
  
  # Ordine finale per i fattori
  row_levels <- rownames(mat)[row_order]
  col_levels <- colnames(mat)[col_order]
  
  # Matrice lunga per ggplot
  df_plot <- as.data.frame(mat) %>%
    rownames_to_column("Ke_description") %>%
    pivot_longer(-Ke_description, names_to = "Experiment", values_to = "BMD_norm")
  
  # Applica i livelli clusterizzati
  df_plot$Ke_description <- factor(df_plot$Ke_description, levels = row_levels)
  df_plot$Experiment     <- factor(df_plot$Experiment, levels = col_levels)
  
  # Heatmap ggplot
  plt <- ggplot(df_plot, aes(x = Experiment, y = Ke_description, fill = BMD_norm)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(
      colours = c("white", "skyblue", "royalblue4"),
      na.value = "grey90",
      name = "BMD_norm") +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title  = element_blank(),
      panel.grid  = element_blank()) +
    ggtitle("Heatmap BMD_norm (automatic clustering)")
  
  return(list(plot = plt, nrow = nrow(mat), ncol = ncol(mat)))
}


for (type in types) {
  df <- ke_df %>%
    filter(Ke_type == type)
  res <- create_heatmap_gg(df)
  ggsave(paste0(dir_path, "heatmap_", type, "_bmd_gg_notimeclust.pdf"), plot = res$plot,
         width  = min(max(res$ncol * 0.3, 8), 20),
         height = min(max(res$nrow * 0.3, 10), 30),
         units = "in")
}


############################################

create_heatmap_gg_horizontal <- function(df) {
  # Order the experiments by time
  experiment_order <- df %>%
    distinct(Experiment, time) %>%
    arrange(time) %>%
    pull(Experiment)
  
  # For each (Ke_description, Experiment) take the minimum BMD_norm
  df_min <- df %>%
    group_by(Ke_description, Experiment) %>%
    summarise(BMD_min = min(BMD_norm, na.rm = TRUE), .groups = "drop")
  
  # Create matrix for clustering
  mat <- df_min %>%
    pivot_wider(names_from = Experiment, values_from = BMD_min) %>%
    column_to_rownames("Ke_description") %>%
    as.matrix()
  
  # Order columns by time
  mat <- mat[, experiment_order]
  
  # Replace NA with high value for clustering so distance calculation doesn't break
  mat_clust <- mat
  mat_clust[is.na(mat_clust)] <- max(mat, na.rm = TRUE) * 1.5
  
  # Cluster rows
  row_order <- hclust(dist(mat_clust, method = "euclidean"), method = "complete")$order
  
  # Prepare data for ggplot
  df_plot <- as.data.frame(mat) %>%
    rownames_to_column("Ke_description") %>%
    pivot_longer(-Ke_description, names_to = "Experiment", values_to = "BMD_norm")
  
  # Set factor levels for horizontal orientation (swap axes)
  df_plot$Ke_description <- factor(df_plot$Ke_description, levels = rownames(mat)[row_order])
  df_plot$Experiment <- factor(df_plot$Experiment, levels = experiment_order)
  
  # Plot with swapped axes: x = Ke_description, y = Experiment
  plt <- ggplot(df_plot, aes(x = Ke_description, y = Experiment, fill = BMD_norm)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(
      colours = c("white", "skyblue", "royalblue4"),
      na.value = "grey90",
      name = "BMD_norm") +
    theme_minimal(base_size = 20) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 20),
      axis.text.y = element_text(angle = 0, hjust = 1, size = 20),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.margin = margin(t = 10, r = 10, b = 40, l = 10)
    ) +
    ggtitle("Heatmap BMD_norm (Horizontal)")
  
  list(plot = plt,
       nrow = nrow(mat),
       ncol = ncol(mat))
}



for (type in types) {
  df <- ke_df %>%
    filter(Ke_type == type)
  
  ht <- create_heatmap_gg_horizontal(df)  # crea heatmap orizzontale
  p  <- ht$plot
  nr <- ht$nrow  # number of KE descriptions (righe)
  nc <- ht$ncol  # number of Experiments (colonne)
  
  # Width basato su nr (KE descriptions), height su nc (Experiments)
  width <- pmax(nr * 0.8, 20)
  height <- pmax(nc * 0.8, 20)
  
  
  ggsave(paste0(dir_path, "heatmap_", type, "_bmd_gg_horizontal.pdf"),
         plot = p, width = width, height = height, units = "in",  limitsize = FALSE)
}


############################################

#Interactive heatmap

library(dplyr)
library(tidyr)
library(plotly)

create_heatmap_plotly <- function(df) {
  # Ordina gli esperimenti per tempo
  experiment_order <- df %>%
    distinct(Experiment, time) %>%
    arrange(time) %>%
    pull(Experiment)
  
  # Calcola il minimo BMD_norm per KE e Esperimento
  df_min <- df %>%
    group_by(Ke_description, Experiment) %>%
    summarise(BMD_min = min(BMD_norm, na.rm = TRUE), .groups = "drop")
  
  # Crea una matrice larga
  mat <- df_min %>%
    pivot_wider(names_from = Experiment, values_from = BMD_min)
  
  # Estrai i nomi di riga e colonna
  row_labels <- mat$Ke_description
  mat <- mat[, -1]  # rimuovi colonna Ke_description
  
  mat_matrix <- as.matrix(mat)
  colnames(mat_matrix) <- experiment_order
  rownames(mat_matrix) <- row_labels
  
  # Crea testo per i tooltip
  hover_text <- expand.grid(
    Ke_description = row_labels,
    Experiment = experiment_order
  ) %>%
    rowwise() %>%
    mutate(
      value = mat_matrix[Ke_description, Experiment],
      label = paste0(
        "KE: ", Ke_description,
        "<br>Experiment: ", Experiment,
        "<br>BMD_norm: ", signif(value, 3)
      )
    ) %>%
    pull(label) %>%
    matrix(nrow = length(row_labels), byrow = TRUE)
  
  # Heatmap interattiva con plotly
  plot_ly(
    x = experiment_order,
    y = row_labels,
    z = mat_matrix,
    type = "heatmap",
    colors = c("white", "skyblue", "royalblue4"),
    text = hover_text,
    hoverinfo = "text",
    showscale = TRUE
  ) %>%
    layout(
      title = "Interactive Heatmap of BMD_norm",
      xaxis = list(title = "Experiment", tickangle = 45),
      yaxis = list(title = "KE Description")
    )
}

for (type in types) {
  df <- ke_df %>%
       filter(Ke_type == type)
 
  ht <- create_heatmap_plotly(df)   
  htmlwidgets::saveWidget(ht, paste0(dir_path, "heatmap_", type, "_prova.html"))
}


############################################


# 1. Count unique KEs per Ke_type, level, and system
ke_counts <- ke_df %>%
  distinct(Ke_description, Ke_type, level, system) %>%  # ensure uniqueness
  group_by(Ke_type, level, system) %>%
  summarise(n_KE = n(), .groups = "drop")

# 2. Boxplot: KE count per Ke_type grouped by LEVEL
boxplot_by_level <- ggplot(ke_counts, aes(x = Ke_type, y = n_KE, fill = Ke_type)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +               # main boxplot
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +             # add individual points
  facet_wrap(~ level, scales = "free_y") +                      # one facet per level
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "Number of KEs per Ke_type by LEVEL",
       y = "Number of KEs", x = "KE Type")

# 2.1. Boxplot: KE count per Ke_type grouped by LEVEL without outliers to have a better visualization
# Function to remove outliers by ke_type and level
remove_outliers <- function(data) {
  data %>%
    group_by(Ke_type, level) %>%
    filter(n_KE >= quantile(n_KE, 0.25, na.rm = TRUE) - 1.5 * IQR(n_KE, na.rm = TRUE),
           n_KE <= quantile(n_KE, 0.75, na.rm = TRUE) + 1.5 * IQR(n_KE, na.rm = TRUE)) %>%
    ungroup()
}

# Apply the function
ke_counts_no_outliers <- remove_outliers(ke_counts)

boxplot_by_level <- ggplot(ke_counts_no_outliers, aes(x = Ke_type, y = n_KE, fill = Ke_type)) +
  geom_boxplot(width = 0.6) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  facet_wrap(~ level, scales = "free_y") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "Number of KEs per Ke_type by LEVEL (Outliers Removed)",
       y = "Number of KEs", x = "KE Type")


# 3. Boxplot: KE count per Ke_type grouped by system
boxplot_by_system <- ggplot(ke_counts, aes(x = Ke_type, y = n_KE, fill = Ke_type)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  facet_wrap(~ system, scales = "free_y") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "Number of KEs per Ke_type by system",
       y = "Number of KEs", x = "KE Type")

ggsave(paste0(dir_path, "boxplot_by_level.pdf"), boxplot_by_level, width = 20, height = 10)
ggsave(paste0(dir_path, "boxplot_by_system.pdf"), boxplot_by_system, width = 20, height = 15)


#####################################
library(RColorBrewer)


## BARPLOT
# Aggregate KE counts by level and Ke_type
ke_sum_level <- ke_counts %>%
  group_by(level, Ke_type) %>%
  summarise(total_n_KE = sum(n_KE), .groups = "drop")

# Plot: X = level, stacked by Ke_type
stacked_bar_by_level <- ggplot(ke_sum_level, aes(x = level, y = total_n_KE, fill = Ke_type)) +
  geom_col(width = 0.7) +
  theme_minimal(base_size = 12) +
  scale_fill_brewer(palette = "Blues", name = "KE Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Total number of KEs per LEVEL (stacked by KE Type)",
       y = "Total number of KEs", x = "LEVEL", fill = "KE Type")

# 1. Separare i sistemi multipli
ke_expanded <- ke_counts %>%
  filter(!is.na(system)) %>%
  mutate(system = str_split(system, "/")) %>%  # divide in lista
  unnest(system) %>%                           # espande in righe
  mutate(system = str_trim(system))            # rimuove spazi extra

# 2. Calcolare il totale KE per tipo e sistema
ke_sum_system <- ke_expanded %>%
  group_by(system, Ke_type) %>%
  summarise(total_n_KE = sum(n_KE), .groups = "drop")

# 3. Calcolare il totale per sistema per ordinarli
system_order <- ke_sum_system %>%
  group_by(system) %>%
  summarise(system_total = sum(total_n_KE)) %>%
  arrange(desc(system_total)) %>%
  pull(system)

# 3. Set factor levels in decreasing order of total KE
ke_sum_system$system <- factor(ke_sum_system$system, levels = system_order)

# 4. Plot: stacked barplot with log scale
stacked_bar_by_system <- ggplot(ke_sum_system, aes(x = system, y = total_n_KE, fill = Ke_type)) +
  geom_col(width = 0.7) +
  theme_minimal(base_size = 12) +
  scale_fill_brewer(palette = "Blues", name = "KE Type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Total number of KEs per SYSTEM (stacked by KE Type)",
       y = "Total number of KEs", x = "System", fill = "KE Type")

ggsave(paste0(dir_path, "stacked_bar_by_level.pdf"), stacked_bar_by_level, width = 6, height = 4)
ggsave(paste0(dir_path, "stacked_bar_by_system.pdf"), stacked_bar_by_system, width = 10, height = 6)


################################################################################

#AO ANALYSIS

ke_desc_freq <- read_excel(paste0(dir_path,"ke_nonunique.xlsx")) %>%
  filter(Ke_type == "AdverseOutcome") %>%
  distinct(TermID, Experiment, Dataset, .keep_all = TRUE) %>%
  count(TermID, Ke_description, sort = TRUE)

print(ke_desc_freq)

ke_freq_plot <- ke_desc_freq %>%
  filter(n > 2)  

# Plot
p <- ggplot(ke_desc_freq, aes(x = reorder(Ke_description, n), y = n, fill = n)) +
  geom_col() +
  coord_flip() +
  labs(title = "Frequency of AdverseOutcome",
       x = "AO_description",
       y = "Frequency") +
  scale_y_continuous(breaks = seq(0, max(ke_freq_plot$n), by = 1)) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave(paste0(dir_path,"frequency_ao_plot.pdf"), plot = p, width = 15, height = 15)

p2 <- ggplot(ke_freq_plot, aes(x = reorder(Ke_description, n), y = n, fill = n)) +
  geom_col() +
  coord_flip() +
  labs(title = "Frequency of AdverseOutcome (from 3 occurrences)",
       x = "AO_description",
       y = "Frequency") +
  scale_y_continuous(breaks = seq(0, max(ke_freq_plot$n), by = 1)) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave(paste0(dir_path,"frequency_from_3_ao_plot.pdf"), plot = p2, width = 15, height = 10)

################ check aop for each ao ########################################

ao <- read_excel(paste0(dir_path,"ke_nonunique.xlsx")) %>% 
  filter(Ke_type == "AdverseOutcome") %>%
  distinct(TermID, Experiment, Dataset, .keep_all = TRUE) %>% 
  group_by(TermID, Ke_description) %>%
  summarise( n = n(),
             Aop = paste(unique(Aop), collapse = ", "),
             .groups = "drop") %>%
  arrange(desc(n))
ao$n <- as.numeric(ao$n)

################################################################################

#PLOT BMD-FREQUENCY AO
library(plotly) 

ke_df_nonunique <- read_excel(paste0(dir_path,"ke_nonunique.xlsx"))

min_mean_bmd <- ke_df_nonunique %>%
  group_by(Ke_description) %>%
  summarise(
    min_BMD  = min(BMD[BMD > 0], na.rm = TRUE),
    mean_BMD = mean(BMD, na.rm = TRUE),
    chemicals = paste(unique(chemical), collapse=", "),
    .groups   = "drop"
  )

df_plot <- ke_desc_freq %>%
  inner_join(min_mean_bmd, by = "Ke_description")

min <- plot_ly(
  data = df_plot,
  x = ~min_BMD,
  y = ~n,
  type = 'scatter',
  mode = 'markers',
  text = ~paste("KE: ", Ke_description,
                "<br>Min BMD: ", signif(min_BMD, 3),
                "<br>Freq: ", n,
                "<br>Chemicals: ", chemicals),
  hoverinfo = 'text',
  marker = list(size = 10, color = 'navy')
) %>%
  layout(
    title = "Interactive dot plot of KE by Frequency and Min BMD",
    xaxis = list(title = "Min BMD", type = "log"),
    yaxis = list(title = "Frequency KE")
  )

mean <- plot_ly(
  data = df_plot,
  x = ~mean_BMD,
  y = ~n,
  type = 'scatter',
  mode = 'markers',
  text = ~paste("KE: ", Ke_description,
                "<br>Mean BMD: ", signif(mean_BMD, 3),
                "<br>Freq: ", n,
                "<br>Chemicals: ", chemicals),
  hoverinfo = 'text',
  marker = list(size = 10, color = 'navy')
) %>%
  layout(
    title = "Interactive dot plot of AO by Frequency and Average BMD",
    xaxis = list(title = "Average BMD"),
    yaxis = list(title = "Frequency AO")
  )

min
mean


library(htmlwidgets)
saveWidget(min, paste0(dir_path,"dotplot_minBMD.html"))
saveWidget(mean, paste0(dir_path,"dotplot_averageBMD.html"))

################################################################################

#PLOT BMD-FREQUENCY AO NORMALIZED PER EXPERIMENT
library(plotly) 
ke_df_nonunique <- read_excel(paste0(dir_path,"ke_nonunique.xlsx"))

min_mean_bmd <- ke_df_nonunique %>%
  group_by(Ke_description) %>%
  summarise(
    min_BMD  = min(BMD_norm, na.rm = TRUE),
    mean_BMD = mean(BMD_norm, na.rm = TRUE),
    chemicals = paste(unique(chemical), collapse=", "),
    .groups   = "drop"
  )

df_plot <- ke_desc_freq %>%
  inner_join(min_mean_bmd, by = "Ke_description")

min <- plot_ly(
  data = df_plot,
  x = ~min_BMD,
  y = ~n,
  type = 'scatter',
  mode = 'markers',
  text = ~paste("KE: ", Ke_description,
                "<br>Min BMD: ", signif(min_BMD, 3),
                "<br>Freq: ", n,
                "<br>Chemicals: ", chemicals),
  hoverinfo = 'text',
  marker = list(size = 10, color = 'navy')
) %>%
  layout(
    title = "Interactive dot plot of KE by Frequency and Min BMD (normalized per experiment)",
    xaxis = list(title = "Min BMD", type = "log"),
    yaxis = list(title = "Frequency KE")
  )

mean <- plot_ly(
  data = df_plot,
  x = ~mean_BMD,
  y = ~n,
  type = 'scatter',
  mode = 'markers',
  text = ~paste("KE: ", Ke_description,
                "<br>Mean BMD: ", signif(mean_BMD, 3),
                "<br>Freq: ", n,
                "<br>Chemicals: ", chemicals),
  hoverinfo = 'text',
  marker = list(size = 10, color = 'navy')
) %>%
  layout(
    title = "Interactive dot plot of AO by Frequency and Average BMD (normalized per experiment)",
    xaxis = list(title = "Average BMD"),
    yaxis = list(title = "Frequency AO")
  )

min
mean


library(htmlwidgets)
saveWidget(min, paste0(dir_path,"dotplot_minBMD_norm_EXP.html"))
saveWidget(mean, paste0(dir_path,"dotplot_averageBMD_norm_EXP.html"))


################################################################################

#PLOT BMD-FREQUENCY AO NORMALIZED PER CHEMICAL
library(plotly) 

#Function to normalize the BMD 
normalize_bmd_chemical <- function(df){
  #Force the column as numeric
  df <- df %>%
    mutate(BMD = as.numeric(BMD))
  
  #The function: groups by chemical, applies a min max normalization and then ungroup the dataframe
  #The min max normalization returns a new values that is a rescaling of the original one so that every value is between 0 and 1
  #In general: normalized_value = (value-min)/(max-min)   -> so we remove the min value to the actual value and we divide for the range
  #In this case all the values will be between 1 and 2 to not have 0 that can be confused with the controls so we add 1 to every value.
  df_normalized <- df %>%
    group_by(chemical) %>%
    mutate(BMD_range = max(BMD, na.rm = TRUE) - min(BMD, na.rm = TRUE),
           BMD_norm = ifelse(BMD_range == 0,
                             1.5,  # return 1.5 if there is only 1 BMD value for that Event
                             1 + (BMD - min(BMD, na.rm = TRUE)) / BMD_range)) %>%
    select(-BMD_range) %>%
    ungroup()
  return(df_normalized)
}


combine_nounique_norm_exp <- function(pattern, path = dir_path) {
  files <- list.files(path = path, pattern = paste0("^", pattern), full.names = TRUE)
  df_list <- lapply(files, function(file) {
    df <- read_excel(file)
    df$file_id <- basename(file)
    #Call the function to merge the dataframe with the Biological_system_annotations
    df <- merge_annotation(df)
    return(df)
  })
  combined_df <- bind_rows(df_list)

  return(combined_df)
}

ke_df_not_norm <- combine_nounique_norm_exp("ke_enrichment_results") %>%
  separate(Experiment, into = c("chemical", "time"), sep = "_", remove = FALSE) %>%
  mutate(time = as.numeric(time))

ke_df_norm_chem <- normalize_bmd_chemical(ke_df_not_norm)

min_mean_bmd <- ke_df_norm_chem %>%
  group_by(Ke_description) %>%
  summarise(
    min_BMD  = min(BMD_norm, na.rm = TRUE),
    mean_BMD = mean(BMD_norm, na.rm = TRUE),
    chemicals = paste(unique(chemical), collapse=", "),
    .groups   = "drop"
  )

df_plot <- ke_desc_freq %>%
  inner_join(min_mean_bmd, by = "Ke_description")

min <- plot_ly(
  data = df_plot,
  x = ~min_BMD,
  y = ~n,
  type = 'scatter',
  mode = 'markers',
  text = ~paste("KE: ", Ke_description,
                "<br>Min BMD: ", signif(min_BMD, 3),
                "<br>Freq: ", n,
                "<br>Chemicals: ", chemicals),
  hoverinfo = 'text',
  marker = list(size = 10, color = 'navy')
) %>%
  layout(
    title = "Interactive dot plot of KE by Frequency and Min BMD (normalized per chemical)",
    xaxis = list(title = "Min BMD", type = "log"),
    yaxis = list(title = "Frequency KE")
  )

mean <- plot_ly(
  data = df_plot,
  x = ~mean_BMD,
  y = ~n,
  type = 'scatter',
  mode = 'markers',
  text = ~paste("KE: ", Ke_description,
                "<br>Mean BMD: ", signif(mean_BMD, 3),
                "<br>Freq: ", n,
                "<br>Chemicals: ", chemicals),
  hoverinfo = 'text',
  marker = list(size = 10, color = 'navy')
) %>%
  layout(
    title = "Interactive dot plot of AO by Frequency and Average BMD (normalized per chemical)",
    xaxis = list(title = "Average BMD"),
    yaxis = list(title = "Frequency AO")
  )

min
mean


library(htmlwidgets)
saveWidget(min, paste0(dir_path,"dotplot_minBMD_norm_CHEM.html"))
saveWidget(mean, paste0(dir_path,"dotplot_averageBMD_norm_CHEM.html"))

