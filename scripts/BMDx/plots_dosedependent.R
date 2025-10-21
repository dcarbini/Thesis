library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(AOPfingerprintR)
library(writexl)
library(readxl)
library(dplyr)
library(tidyr)

folder <- "C:/Users/Utente/Desktop/denys/results/DoseDependentGenes/"

all_data <- read_excel(paste0(folder,"all_doseDependent.xlsx"))

dim(all_data)

##################### List Accepted Chemicals ########################################


nemesis <- c("BPA", "BPF", "BPS", "DEP", "DINP", "Hexamoll DINCH", "PFOA", "PFHxS", "ZEA")
eu <- c("BPA", "BPB", "BPS", "BPAF", "paranonylphenol", "Prochloraz", "Propiconazole", "Ziram", "4-(4-propan-2-yloxyphenyl)sulfonylphenol", "dibutylphthalate", "lithiumchloride", "DEHP", "triclosan", "TBBA")
zenodo <- c('RU486','Clobetasol', 'Dihydroxyvitamin3','Dehydrorepiandrosterone',  'BPA','Genistein', 'ValproicAcid', 'Tamoxifen', 'Progesterone', 'Ketoconazole', 'Phenobarbital', 'Flutamide', 'Trenbolone', 'Griseofulvin', 'Methanol','Estradiol', 'Nicotine','Vinblastine', 'FK506/Cyclosporin A',  'diethylstilbestrol', 'Daidzein', 'Ethynylestradiol', 'paranonylphenol', 'BPS', 'BPF',  'Glyphosate', 'Roundup', 'BPAF', 'BPAP','BPB', 'BPZ', 'PPT', 'benzo[a]pyrene', 'Ethanol','TBBA','pendimethalin','stomp_aqua_pendimethalin', 'Reserpine', 'Thiram', 'Propiconazole','Cyproterone.acetate', 'Bifenthrin', 'Trifloxystrobin', 'PFOA', 'Cycloheximide','Simazine', 'Cyproconazole', 'Simvastatin', 'Pyraclostrobin', 'PFOS', 'Triiodothyronine', 'Vinclozolin', '4.Hydroxytamoxifen', 'Maneb','Tetrac', 'Ziram','4.Cumylphenol', 'Lovastatin', 'Cyanazine', 'Fulvestrant', 'Nilutamide','Cypermethrin', 'Prochloraz', 'gesaprim', 'Imazalil', 'Rotenone','Dexamethasone', 'ZEA','TGSA', 'BADGE')
more <- c("2,4.BPF","2,4.BPS", "17B.estradiol", "4,4.BPF", "BPS.MPE", "Dex", "para.nonylphenol")
accepted_chemicals <- unique(c(nemesis, eu, zenodo, more))

##################### Bar plot DDGs per time ####################################

counts <- all_data %>%
  filter(Experiment %in% accepted_chemicals) %>%
  distinct(Experiment, exposure_time, Feature) %>%
  group_by(exposure_time, Experiment) %>%
  summarise(number = n(), .groups = "drop")

experiment_multiple_times <- counts %>%
  group_by(Experiment) %>%
  summarise(n_times = n_distinct(exposure_time)) %>%
  filter(n_times > 1) %>%
  pull(Experiment)

counts <- counts %>%
  filter(Experiment %in% experiment_multiple_times)

counts$exposure_time <- factor(counts$exposure_time,
                               levels = c(6, 8, 12, 24, 48, 72))
experiment_order <- counts %>%
  group_by(Experiment) %>%
  summarise(total_genes = sum(number), .groups = "drop") %>%
  arrange(desc(total_genes)) %>%
  pull(Experiment)

counts$Experiment <- factor(counts$Experiment, levels = experiment_order)
y_breaks <- c(10, 50, 100, 500, 1000, 5000, 10000)

bar_dodge <- ggplot(counts, aes(x = Experiment, y = number, fill = exposure_time)) +
  geom_col(width = 0.7, position = position_dodge(width = 0.8)) +
  theme_minimal(base_size = 18) +
  scale_y_log10(breaks = y_breaks, labels = y_breaks) + 
  scale_fill_brewer(palette = "GnBu", name = "exposure_time", direction = -1)   +   
  labs(title = "Total number of dose dependent genes per Chemical",
       y = "Total number of dose dependent genes",
       x = "Chemicals",
       fill = "exposure_time") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("C:/Users/Utente/Desktop/denys/results/DoseDependentGenes/bar.pdf", bar_dodge)


#################### NO #########################################

#NO

plot_data <- all_data %>%
  mutate(Dataset = ifelse(str_detect(Dataset, "GSE249377_"),
                          "GSE249377",
                          Dataset))

plot_data$exposure_time <- factor(plot_data$exposure_time,
                                  levels = c(6, 8, 12, 24, 48, 72))

dot_plot <- ggplot(plot_data, aes(x = Experiment, y = exposure_time, color = Dataset)) +
  geom_point(alpha = 0.6, position = position_jitter(width = 0.2, height = 0.2)) +
  theme_minimal(base_size = 12) +
  scale_color_brewer(palette = "Set3", name = "Dataset") +
  labs(title = "Dose dependent genes per Chemical and Exposure time",
       x = "Chemicals",
       y = "Exposure time") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

###################### Number DDGs per dataset ############################################

#all_data$dataset_time <- paste(all_data$Dataset, all_data$exposure_time, sep = "_")

# aggiungiamo colonna combinata
counts <- all_data %>%
  mutate(ChemTime = paste(Experiment, exposure_time, sep = "_")) %>%
  group_by(Dataset, ChemTime) %>%
  summarise(n_features = n_distinct(Feature), .groups = "drop")

# Tutti i dataset tranne GSE249377
counts_other <- counts %>% 
  filter(!grepl("GSE249377", Dataset))%>% 
  filter(!grepl("GSE17624", Dataset))%>% 
  filter(!grepl("GSE153320", Dataset))
  
bar_plot_other <- ggplot(counts_other, aes(x = ChemTime, y = n_features, fill = Dataset)) +
  geom_col(width = 0.7) +
  facet_wrap(~ Dataset, scales = "free") +
  theme_minimal(base_size = 12) +
  scale_y_log10(labels = scales::number_format())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(
    title = "Number of dose dependent genes per Dataset",
    x = "Chemical + Exposure time",
    y = "Number of Dose Dependent Genes",
    fill = "Dataset"
  )

bar_plot_other  
  

counts_other <- counts %>% 
  filter(grepl("GSE17624|GSE153320", Dataset))


counts_other <- counts_other %>%
  mutate(ChemTime = factor(ChemTime, levels = c("BPA_72", "BPF_72", "BPA_8", "BPA_24", "BPA_48", "BPS_72", "TBBA_72")))


bar_plot_other <- ggplot(counts_other, aes(x = ChemTime, y = n_features, fill = Dataset)) +
    geom_col(width = 0.7) +
    facet_wrap(~ Dataset, scales = "free") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(
      title = "Number of dose dependent genes per Dataset",
      x = "Chemical + Exposure time",
      y = "Number of Dose Dependent Genes",
      fill = "Dataset"
    )
  
bar_plot_other  



# Solo GSE249377
gse <- all_data %>% filter(grepl("GSE249377", Dataset))
counts_gse <- gse %>%
  group_by(Dataset, Experiment) %>%
  summarise(n_features = n_distinct(Feature), .groups = "drop")

# barplot per ogni dataset
bar_plot <- ggplot(counts_gse, aes(x = Experiment, y = n_features, fill = Dataset)) +
  geom_col(width = 0.7) +
  facet_wrap(~ Dataset, scales = "free",) + 
  scale_y_log10(labels = scales::number_format())+
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = "Number of dose dependent genes per chemical in GSE249377",
       x = "Chemical",
       y = "Number of Dose Dependent Genes",
       fill = "Dataset")

bar_plot

###################### Comparison number DEGs and DDGs ############################################

deg <- read_excel(paste0(folder,"all_DEGs.xlsx"))

deg <- deg %>% mutate(Category = "DEGs")

all_data <- all_data %>% mutate(Category = "DDGs")

combined_df <- bind_rows(deg, all_data)


counts <- combined_df %>%
  mutate(Experiment = ifelse(Experiment == "paranonylphenol", "para.nonylphenol", Experiment),
         Experiment = ifelse(Experiment == "17B.estradiol", "Estradiol", Experiment),
         Experiment = ifelse(Experiment == "Ethynyl", "Ethynylestradiol", Experiment),
         Experiment = ifelse(Experiment == "Valproic", "ValproicAcid", Experiment))%>%
  filter(Experiment %in% accepted_chemicals) %>%
  distinct(Experiment, Category, Feature) %>%
  group_by(Category, Experiment) %>%
  summarise(number = n(), .groups = "drop")

experiment_order <- counts %>%
  group_by(Experiment) %>%
  summarise(total_genes = sum(number), .groups = "drop") %>%
  arrange(desc(total_genes)) %>%
  pull(Experiment)

counts$Experiment <- factor(counts$Experiment, levels = experiment_order)
y_breaks <- c(10, 50, 100, 500, 1000, 5000, 10000)
col <- c("#ffadad","#a0c4ff")
bar_dodge <- ggplot(counts, aes(x = Experiment, y = number, fill = Category)) +
  geom_col(width = 0.7, position = position_dodge(width = 0.8)) +
  theme_minimal(base_size = 12) +
  scale_y_log10(breaks = y_breaks, labels = y_breaks) + 
  scale_fill_manual(values = col, name = "Category") +    
  labs(title = "Total number of dose dependent genes per Chemical",
       y = "Number of genes",
       x = "Chemicals",
       fill = "Category") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



ggsave("C:/Users/Utente/Desktop/denys/results/DoseDependentGenes/confronto_DEG_DDG.pdf", bar_dodge)

####################### Upset plot comparison DEGs and DDGs per chemical #####################################
combined_df <- combined_df %>%
  mutate(Experiment = ifelse(Experiment == "paranonylphenol", "para.nonylphenol", Experiment),
         Experiment = ifelse(Experiment == "17B.estradiol", "Estradiol", Experiment),
         Experiment = ifelse(Experiment == "Ethynyl", "Ethynylestradiol", Experiment),
         Experiment = ifelse(Experiment == "Valproic", "ValproicAcid", Experiment))%>%
  filter(Experiment %in% accepted_chemicals)

library(UpSetR)

list <- list()

for(chem in unique(combined_df$Experiment)) {
  sub <- combined_df[combined_df$Experiment == chem, ]
  contents <- split(sub$Feature, sub$Category)
  
  png(paste0(folder, "/upset/upset_", chem, ".png"), width = 1600, height = 1200, res = 150)
  upset(fromList(contents),
              sets = names(contents),
              order.by = "freq",
              mainbar.y.label = paste("Intersections - ", chem),
              sets.x.label = "Number of genes",
              keep.order = TRUE)
  dev.off()
  
  # list[[chem]] <- grDevices::recordPlot()
  # 
  # # Pulisci il device per il prossimo plot
  # dev.off()  # <-- ATTENZIONE: chiude la finestra grafica corrente
  # grDevices::dev.new() 
}

print(list)

for (chem in names(list)) {
  png(paste0(folder, "/upset/upset_", chem, ".png"), width = 1600, height = 1200, res = 150)
  replayPlot(list[[chem]])
  dev.off()
}


 
chem <-   "Estradiol"            


sub <- combined_df[combined_df$Experiment == chem, ]
contents <- split(sub$Feature, sub$Category)

png(paste0(folder, "/upset/upset_", chem, ".png"), width = 1600, height = 1200, res = 150)
upset(fromList(contents),
      sets = names(contents),
      order.by = "freq",
      mainbar.y.label = paste("Intersections - ", chem),
      sets.x.label = "Number of genes",
      keep.order = TRUE,
      text.scale = c(2, 1.5, 1.5, 1.5))
dev.off()


library(magick)

# Cartella con le immagini
img_files <- list.files(paste0(folder, "/upset/"), pattern = "\\.png$", full.names = TRUE)

# Carica tutte le immagini come oggetti magick-image
imgs <- lapply(img_files, image_read)

# Trasforma la lista in un vettore di immagini (magick richiede un vettore, non lista)
imgs_vec <- do.call(c, imgs)

# Suddividi in blocchi di 9 immagini per pagina
pages <- split(imgs_vec, ceiling(seq_along(imgs_vec)/9))

# Funzione per combinare immagini in multipanel
combine_images <- function(img_list, ncol = 3) {
  # Crea righe
  rows <- lapply(seq(1, length(img_list), ncol), function(i) {
    row_imgs <- img_list[i:min(i+ncol-1, length(img_list))]
    image_append(row_imgs, stack = FALSE)  # append orizzontale
  })
  # Combina tutte le righe verticalmente
  combined <- image_append(do.call(c, rows), stack = TRUE)
  return(combined)
}

# Salva ogni pagina multipanel
for (i in seq_along(pages)) {
  page_combined <- combine_images(pages[[i]], ncol = 3)
  image_write(page_combined, path = file.path(paste0(folder, "/upset/"), paste0("upset_page_", i, ".png")))
}

###################### Upset plot DDGs shared across chemicals ######################################


all <- all_data %>%
  mutate(Experiment = ifelse(Experiment == "paranonylphenol", "para.nonylphenol", Experiment),
         Experiment = ifelse(Experiment == "17B.estradiol", "Estradiol", Experiment),
         Experiment = ifelse(Experiment == "Ethynyl", "Ethynylestradiol", Experiment),
         Experiment = ifelse(Experiment == "Valproic", "ValproicAcid", Experiment))
  distinct(Experiment, Feature) 


counts <- all %>%
  group_by(Experiment) %>%
  summarise(n_geni = n_distinct(Feature))  # conto i geni unici

# Prendo solo i chemical con almeno 50 geni 
valid_chemicals <- counts %>%
  filter(n_geni >= 50) %>%
  pull(Experiment)

# Filtro il DataFrame originale
all <- all %>%
  filter(Experiment %in% valid_chemicals)

contents <- split(all$Feature, all$Experiment)

upset(
  fromList(contents),
  sets = names(contents),
  order.by = "freq",
  mainbar.y.label = paste("Intersections"),
  sets.x.label = "Number of genes"
)
