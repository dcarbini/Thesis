#Import the libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(circlize) 
library(writexl)
library(ggbreak)
library(stringr)
library(AOPfingerprintR)

#Set the directory
#dir_path <- "C:/Users/carbi/Desktop/ke_enrichment/"
dir_path <- "C:/Users/Utente/Desktop/denys/results/ke_enrichment/"

#Function to take only the unique Events in a dataset
unique_events <- function(df) {
  df <- df %>% distinct(TermID, Experiment, file_id, .keep_all = TRUE) 
  return(df)
}

#Function to map the dataframe with the Biological_system_annotations
merge_annotation <- function(df){
  #Not needed if used aopfingerprint package
  #load("C:/Users/carbi/Downloads/Biological_system_annotations.rda")
  #load("C:/Users/Utente/Desktop/denys/scripts/bmdx_noshiny-master/Biological_system_annotations.rda")
  annot <- Biological_system_annotations
  # Rename one of the 2 columns 
  names(annot)[names(annot) == "ke"] <- "TermID"
  
  # Merging on the dataframes on the column "TermId"
  final_df <- merge(df, annot, by = "TermID")
  return(final_df)
}

#Function to normalize the BMD 
normalize_bmd <- function(df){
  #Force the column as numeric
  df <- df %>%
    mutate(BMD = as.numeric(BMD))
  
  #The function: groups by experiment, applies a min max normalization and then ungroup the dataframe
  #The min max normalization returns a new values that is a rescaling of the original one so that every value is between 0 and 1
  #In general: normalized_value = (value-min)/(max-min)   -> so we remove the min value to the actual value and we divide for the range
  #In this case all the values will be between 1 and 2 to not have 0 that can be confused with the controls so we add 1 to every value.
  df_normalized <- df %>%
    group_by(Experiment) %>%
    mutate(BMD_range = max(BMD, na.rm = TRUE) - min(BMD, na.rm = TRUE),
           BMD_norm = ifelse(BMD_range == 0,
                             1.5,  # return 1.5 if there is only 1 BMD value for that Event
                             1 + (BMD - min(BMD, na.rm = TRUE)) / BMD_range)) %>%
    select(-BMD_range) %>%
    ungroup()
  return(df_normalized)
}


#Function to read and combine excel files with the same name
combine_excels <- function(pattern, path = dir_path) {
  files <- list.files(path = path, pattern = paste0("^", pattern), full.names = TRUE)
  df_list <- lapply(files, function(file) {
    df <- read_excel(file)
    df$file_id <- basename(file)
    #Call the function to take only the unique events
    df <- unique_events(df)
    #Call the function to merge the dataframe with the Biological_system_annotations
    df <- merge_annotation(df)
    #Normalize the BMD values
    df <- normalize_bmd(df)
    return(df)
  })
  combined_df <- bind_rows(df_list)
  return(combined_df)
}


# Create the dataframes
ke <- combine_excels("ke_enrichment_results")

#Add 2 columns to the KE dataframe
ke <- ke %>% 
  separate(Experiment, into = c("chemical", "time"), sep = "_", remove = FALSE) %>%
  mutate(time = as.numeric(time))

#Checks
View(ke)

write_xlsx(ke, paste0(dir_path,"ke.xlsx"))

#####################################
library(RColorBrewer)

# 1. Count unique KEs per Ke_type, level, system and organ_tissue
ke_counts <- ke %>%
  distinct(Ke_description, Ke_type, level, system, organ_tissue) %>%
  group_by(Ke_type, level, system, organ_tissue) %>%
  summarise(n_KE = n(), .groups = "drop")

type <- aop_ke_table_hure %>%
  group_by(Ke) %>%
  summarise(
    across(
      .cols = everything(),
      .fns  = ~ paste(unique(.x), collapse = "/")
    ),
    .groups = "drop"
  ) %>%
  rename(ke = Ke)

bio <- Biological_system_annotations %>% 
  left_join(type, by = "ke", suffix = c("", "_type")) 



## BARPLOT

##LEVEL
# Aggregate KE counts by level and Ke_type
ke_sum_level <- ke_counts %>%
  group_by(level, Ke_type) %>%
  summarise(total_n_KE = sum(n_KE), .groups = "drop")

# Plot: X = level, stacked by Ke_type
#Order the labels
ke_sum_level$level <- factor(
  ke_sum_level$level,
  levels = c("Molecular", "Cellular", "Tissue", "Organ", "Individual")
)
#Plot
stacked_bar_by_level <- ggplot(ke_sum_level, aes(x = level, y = total_n_KE, fill = Ke_type)) +
  geom_col(width = 0.7) +
  theme_minimal(base_size = 12) +
  scale_fill_brewer(palette = "Blues", name = "KE Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Total number of KEs per LEVEL (stacked by KE Type)",
       y = "Total number of KEs", x = "LEVEL", fill = "KE Type")

ggsave("C:/Users/Utente/Desktop/denys/immagini/multiDose/stacked_bar_by_level.pdf", stacked_bar_by_level, width = 6, height = 4)
ggsave(paste0(dir_path, "stacked_bar_by_level.pdf"), stacked_bar_by_level, width = 6, height = 4)


##SYSTEM

# Separate multiple systems
ke_expanded <- ke_counts %>%
  filter(!is.na(system)) %>%
  mutate(system = str_split(system, "/")) %>%  # divide in list
  unnest(system) %>%                           # expand the rows
  mutate(system = str_trim(system))            # remove extra spaces

# Aggregate KE counts by system and Ke_type
ke_sum_system <- ke_expanded %>%
  group_by(system, Ke_type) %>%
  summarise(total_n_KE = sum(n_KE), .groups = "drop")

# Calculate the total per system to reorder them 
system_order <- ke_sum_system %>%
  group_by(system) %>%
  summarise(system_total = sum(total_n_KE)) %>%
  arrange(desc(system_total)) %>%
  pull(system)

# Set factor levels in decreasing order of total KE
ke_sum_system$system <- factor(ke_sum_system$system, levels = system_order)

# Plot: stacked barplot
stacked_bar_by_system <- ggplot(ke_sum_system, aes(x = system, y = total_n_KE, fill = Ke_type)) +
  geom_col(width = 0.7) +
  theme_minimal(base_size = 12) +
  scale_fill_brewer(palette = "Blues", name = "KE Type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Total number of KEs per SYSTEM (stacked by KE Type)",
       y = "Total number of KEs", x = "System", fill = "KE Type")

ggsave("C:/Users/Utente/Desktop/denys/immagini/multiDose/stacked_bar_by_system.pdf", stacked_bar_by_system, width = 10, height = 6)
ggsave(paste0(dir_path, "stacked_bar_by_system.pdf"), stacked_bar_by_system, width = 10, height = 6)

total_KE_sys_type <- bio %>% 
  filter(!is.na(system)) %>%
  mutate(Ke_type = str_split(Ke_type, "/")) %>%
  unnest(Ke_type) %>%
  mutate(system = str_split(system, "/")) %>%
  unnest(system) %>%
  mutate(system = str_trim(system))%>%
  group_by(system) %>%
  summarise(total_n_KE = n(), .groups = "drop")

total_KE_sys_notype <- bio %>% 
  filter(!is.na(system)) %>%
  mutate(system = str_split(system, "/")) %>%
  unnest(system) %>%
  mutate(system = str_trim(system))%>%
  group_by(system) %>%
  summarise(n_KE = n(), .groups = "drop")

total_KE_sys <- total_KE_sys_type %>%
  left_join(total_KE_sys_notype, by = "system", suffix = c("", "_notype")) %>%
  mutate(difference = total_n_KE - n_KE)

# Separate multiple systems
ke_expanded <- ke_counts %>%
  filter(!is.na(system)) %>%
  mutate(system = str_split(system, "/")) %>%  # divide in list
  unnest(system) %>%                           # expand the rows
  mutate(system = str_trim(system))            # remove extra spaces

# Aggregate KE counts by system and Ke_type
ke_sum_system <- ke_expanded %>%
  group_by(system, Ke_type) %>%
  summarise(total_n_KE = sum(n_KE), .groups = "drop") %>%
  left_join(total_KE_sys, by = "system", suffix = c("", "_total")) %>%
  mutate(norm_n_KE = total_n_KE / total_n_KE_total)

# Calculate the total per system to reorder them 
system_order <- ke_sum_system %>%
  group_by(system) %>%
  summarise(system_total = sum(total_n_KE)) %>%
  arrange(desc(system_total)) %>%
  pull(system)

# Set factor levels in decreasing order of total KE
ke_sum_system$system <- factor(ke_sum_system$system, levels = system_order)

# Plot: stacked barplot
stacked_bar_by_system <- ggplot(ke_sum_system, aes(x = system, y = norm_n_KE, fill = Ke_type)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = total_n_KE), 
            position = position_stack(vjust = 0.5), 
            size = 3, color = "black") +
  theme_minimal(base_size = 12) +
  scale_fill_brewer(palette = "Blues", name = "KE Type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Normalized number of KEs per SYSTEM (stacked by KE Type)",
       y = "Normalized number of KEs", x = "System", fill = "KE Type")

ggsave("C:/Users/Utente/Desktop/denys/immagini/multiDose/stacked_bar_by_system_norm2.pdf", stacked_bar_by_system, width = 10, height = 6)
ggsave(paste0(dir_path, "stacked_bar_by_system_norm.pdf"), stacked_bar_by_system, width = 10, height = 6)

##ORGAN
# Separate multiple organs
ke_expanded_organ <- ke_counts %>%
  filter(!is.na(organ_tissue)) %>%
  mutate(organ_tissue = str_split(organ_tissue, "/")) %>%  # divide in list
  unnest(organ_tissue) %>%                           # expand the rows
  mutate(organ_tissue = str_trim(organ_tissue))            # remove extra spaces

# Aggregate KE counts by organ and Ke_type
ke_sum_organ <- ke_expanded_organ %>%
  group_by(organ_tissue, Ke_type) %>%
  summarise(total_n_KE = sum(n_KE), .groups = "drop")

# Calculate the total per organ to reorder them 
organ_order <- ke_sum_organ %>%
  group_by(organ_tissue) %>%
  summarise(organ_total = sum(total_n_KE)) %>%
  arrange(desc(organ_total)) %>%
  pull(organ_tissue)

# Set factor levels in decreasing order of total KE
ke_sum_organ$organ_tissue <- factor(ke_sum_organ$organ_tissue, levels = organ_order)

# Plot: stacked barplot with log scale
stacked_bar_by_organ <- ggplot(ke_sum_organ, aes(x = organ_tissue, y = total_n_KE, fill = Ke_type)) +
  geom_col(width = 0.7) +
  theme_minimal(base_size = 12) +
  scale_fill_brewer(palette = "Blues", name = "KE Type") +
  scale_y_break(c(30,70)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(
    title = "Total number of KEs per ORGAN_TISSUE (stacked by KE Type)",
    y = "Total number of KEs", x = "Organ_Tissue", fill = "KE Type"
  )

ggsave("C:/Users/Utente/Desktop/denys/immagini/multiDose/stacked_bar_by_organ.pdf", stacked_bar_by_organ, width = 12, height = 6)
ggsave(paste0(dir_path, "stacked_bar_by_organ.pdf"), stacked_bar_by_organ, width = 12, height = 6)


total_KE_organ_type <- bio %>% 
  filter(!is.na(organ_tissue)) %>%
  mutate(Ke_type = str_split(Ke_type, "/")) %>%
  unnest(Ke_type) %>%
  mutate(organ_tissue = str_split(organ_tissue, "/")) %>%
  unnest(organ_tissue) %>%
  mutate(organ_tissue = str_trim(organ_tissue))%>%
  group_by(organ_tissue) %>%
  summarise(total_n_KE = n(), .groups = "drop")

total_KE_organ_notype <- bio %>% 
  filter(!is.na(organ_tissue)) %>%
  mutate(organ_tissue = str_split(organ_tissue, "/")) %>%
  unnest(organ_tissue) %>%
  mutate(organ_tissue = str_trim(organ_tissue))%>%
  group_by(organ_tissue) %>%
  summarise(n_KE = n(), .groups = "drop")

total_KE_organ <- total_KE_organ_type %>%
  left_join(total_KE_organ_notype, by = "organ_tissue", suffix = c("", "_notype")) %>%
  mutate(difference = total_n_KE - n_KE)

# Separate multiple organs
ke_expanded_organ <- ke_counts %>%
  filter(!is.na(organ_tissue)) %>%
  mutate(organ_tissue = str_split(organ_tissue, "/")) %>%  # divide in list
  unnest(organ_tissue) %>%                           # expand the rows
  mutate(organ_tissue = str_trim(organ_tissue))            # remove extra spaces

# Aggregate KE counts by organ and Ke_type
ke_sum_organ <- ke_expanded_organ %>%
  group_by(organ_tissue, Ke_type) %>%
  summarise(total_n_KE = sum(n_KE), .groups = "drop") %>%
  left_join(total_KE_organ, by = "organ_tissue", suffix = c("", "_total")) %>%
  mutate(norm_n_KE = total_n_KE / total_n_KE_total)

# Calculate the total per organ to reorder them 
organ_order <- ke_sum_organ %>%
  group_by(organ_tissue) %>%
  summarise(organ_total = sum(total_n_KE)) %>%
  arrange(desc(organ_total)) %>%
  pull(organ_tissue)

# Set factor levels in decreasing order of total KE
ke_sum_organ$organ_tissue <- factor(ke_sum_organ$organ_tissue, levels = organ_order)

# Plot: stacked barplot with log scale
stacked_bar_by_organ <- ggplot(ke_sum_organ, aes(x = organ_tissue, y = norm_n_KE, fill = Ke_type)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = total_n_KE), 
            position = position_stack(vjust = 0.5), 
            size = 3, color = "black") +
  theme_minimal(base_size = 12) +
  scale_fill_brewer(palette = "Blues", name = "KE Type") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(
    title = "Normalized number of KEs per ORGAN_TISSUE (stacked by KE Type)",
    y = "Normalized number of KEs", x = "Organ_Tissue", fill = "KE Type"
  )

ggsave("C:/Users/Utente/Desktop/denys/immagini/multiDose/stacked_bar_by_organ_norm2.pdf", stacked_bar_by_organ, width = 12, height = 6)
ggsave(paste0(dir_path, "stacked_bar_by_organ_norm.pdf"), stacked_bar_by_organ, width = 12, height = 6)


################################################################################

#AO ANALYSIS

ke_desc_freq <- as.data.frame(ke_df) %>%
  filter(Ke_type == "AdverseOutcome") %>%
  count(Ke_description, sort = TRUE)

print(ke_desc_freq)

ke_freq_plot <- ke_desc_freq %>%
  filter(n > 2)  

# Plot
p <- ggplot(ke_desc_freq, aes(x = reorder(Ke_description, n), y = n, fill = n)) +
  geom_col() +
  coord_flip() +
  labs(title = "Frequency of AdverseOutcome",
       x = "Ke_description",
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
       x = "Ke_description",
       y = "Frequency") +
  scale_y_continuous(breaks = seq(0, max(ke_freq_plot$n), by = 1)) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave(paste0(dir_path,"frequency_from_3_ao_plot.pdf"), plot = p2, width = 15, height = 10)
