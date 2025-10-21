library(AOPfingerprintR)
library(readxl)
library(dplyr)
library(ggplot2)
library(writexl)



dir_path <- "C:/Users/Utente/Desktop/denys/results/"


nemesis <- c("BPA", "BPF", "BPS", "DEP", "DINP", "Hexamoll DINCH", "PFOA", "PFHxS", "ZEA")
eu <- c("BPA", "BPB", "BPS", "BPAF", "paranonylphenol", "Prochloraz", "Propiconazole", "Ziram", "4-(4-propan-2-yloxyphenyl)sulfonylphenol", "dibutylphthalate", "lithiumchloride", "DEHP", "triclosan", "TBBA")
zenodo <- c('RU486','Clobetasol', 'Dihydroxyvitamin3','Dehydrorepiandrosterone',  'BPA','Genistein', 'ValproicAcid', 'Tamoxifen', 'Progesterone', 'Ketoconazole', 'Phenobarbital', 'Flutamide', 'Trenbolone', 'Griseofulvin', 'Methanol','Estradiol', 'Nicotine','Vinblastine', 'FK506/Cyclosporin A',  'diethylstilbestrol', 'Daidzein', 'Ethynylestradiol', 'paranonylphenol', 'BPS', 'BPF',  'Glyphosate', 'Roundup', 'BPAF', 'BPAP','BPB', 'BPZ', 'PPT', 'benzo[a]pyrene', 'Ethanol','TBBA','pendimethalin','stomp_aqua_pendimethalin', 'Reserpine', 'Thiram', 'Propiconazole','Cyproterone.acetate', 'Bifenthrin', 'Trifloxystrobin', 'PFOA', 'Cycloheximide','Simazine', 'Cyproconazole', 'Simvastatin', 'Pyraclostrobin', 'PFOS', 'Triiodothyronine', 'Vinclozolin', '4.Hydroxytamoxifen', 'Maneb','Tetrac', 'Ziram','4.Cumylphenol', 'Lovastatin', 'Cyanazine', 'Fulvestrant', 'Nilutamide','Cypermethrin', 'Prochloraz', 'gesaprim', 'Imazalil', 'Rotenone','Dexamethasone', 'ZEA','TGSA', 'BADGE')
more <- c("2,4.BPF","2,4.BPS", "17B.estradiol", "4,4.BPF", "BPS.MPE", "Dex", "para.nonylphenol")
accepted_chemicals <- unique(c(nemesis, eu, zenodo, more))

ke_nonunique <- read_excel(paste0(dir_path, "ke_enrichment/ke_nonunique.xlsx"))%>%
  mutate(chemical = ifelse(chemical == "paranonylphenol", "para.nonylphenol", chemical),
         chemical = ifelse(chemical == "17B.estradiol", "Estradiol", chemical),
         Experiment = ifelse(Experiment == "paranonylphenol_48", "para.nonylphenol_48", Experiment),
         Experiment = ifelse(Experiment == "17B.estradiol_48", "Estradiol_48", Experiment)) %>%
  filter(chemical %in% accepted_chemicals)


################# Distributions BMD_norm ##################################
##### per chemical

all_DDG <- read_excel(paste0(dir_path, "DoseDependentGenes/all_doseDependent.xlsx")) %>%
  distinct(Experiment, Feature, .keep_all = T)

bin0 <- all_DDG %>%
  filter(BMD != 0, Experiment != "BPA")%>%
  group_by(Experiment) %>%
  filter(n() <= 40) %>%
  ungroup()

bin0_p <- ggplot(bin0, aes(x = BMD_norm, fill = Experiment)) +
  geom_histogram(
    position = "identity",
    alpha = 0.4,
    bins = 100
  ) +
  labs(
    title = "Chemicals with less than 40 dose dependent genes",
    x = "BMD_norm",
    y = "Number of dose dependent genes"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

bin0_p

bin1 <- all_DDG %>%
  filter(BMD != 0, Experiment != "BPA")%>%
  group_by(Experiment) %>%
  filter(n() <= 150, n() > 40) %>%
  ungroup()

bin1_p <- ggplot(bin1, aes(x = BMD_norm, fill = Experiment)) +
  geom_histogram(
    position = "identity",
    alpha = 0.4,
    bins = 100
  ) +
  labs(
    title = "Chemicals with number of dose dependent genes between 41 and 150",
    x = "BMD_norm",
    y = "Number of dose dependent genes"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

bin1_p

bin5 <- all_DDG %>%
  filter(BMD != 0, Experiment != "BPA")%>%
  group_by(Experiment) %>%
  filter(n() <= 400, n() > 150) %>%
  ungroup()

bin5_p <- ggplot(bin5, aes(x = BMD_norm, fill = Experiment)) +
  geom_histogram(
    position = "identity",
    alpha = 0.4,
    bins = 100
  ) +
  labs(
    title = "Chemicals with number of dose dependent genes between 151 and 500",
    x = "BMD_norm",
    y = "Number of dose dependent genes"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

bin5_p

bin2 <- all_DDG %>%
  distinct(Experiment, Feature, .keep_all = T)%>%
  filter(BMD != 0, Experiment != "BPA")%>%
  group_by(Experiment) %>%
  filter(n() > 400, n() <= 1000) %>%
  ungroup()

bin2_p <- ggplot(bin2, aes(x = BMD_norm, fill = Experiment)) +
  geom_histogram(
    position = "identity",
    alpha = 0.4,
    bins = 100
  ) +
  labs(
    title = "Chemicals with number of dose dependent genes between 501 and 1000",
    x = "BMD_norm",
    y = "Number of dose dependent genes"
  ) +
  theme_minimal() +
  theme(legend.position = "right")
bin2_p

bin3 <- all_DDG %>%
  distinct(Experiment, Feature, .keep_all = T)%>%
  filter(BMD != 0, Experiment != "BPA")%>%
  group_by(Experiment) %>%
  filter(n() > 1000, n() <= 2000) %>%
  ungroup()

bin3_p <-ggplot(bin3, aes(x = BMD_norm, fill = Experiment)) +
  geom_histogram(
    position = "identity",
    alpha = 0.4,
    bins = 100
  ) +
  labs(
    title = "Chemicals with number of dose dependent genes between 1001 and 2000",
    x = "BMD_norm",
    y = "Number of dose dependent genes"
  ) +
  theme_minimal() +
  theme(legend.position = "right")
bin3_p


bin4 <- all_DDG %>%
  distinct(Experiment, Feature, .keep_all = T)%>%
  filter(BMD != 0, Experiment != "BPA")%>%
  group_by(Experiment) %>%
  filter(n() > 2000) %>%
  ungroup()

bin4_p <-ggplot(bin4, aes(x = BMD_norm, fill = Experiment)) +
  geom_histogram(
    position = "identity",
    alpha = 0.4,
    bins = 100
  ) +
  labs(
    title = "Chemicals with more than 2000 dose dependent genes",
    x = "BMD_norm",
    y = "Number of dose dependent genes"
  ) +
  theme_minimal() +
  theme(legend.position = "right")
bin4_p

bin <- all_DDG %>%
  filter(Experiment == "BPA")

BPA_p <-ggplot(bin, aes(x = BMD_norm, fill = Experiment)) +
  geom_histogram(
    position = "identity",
    alpha = 0.4,
    bins = 100
  ) +
  labs(
    title = "BPA",
    x = "BMD_norm",
    y = "Number of dose dependent genes"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
BPA_p

library(patchwork)

combined_plot <- (bin1_p | bin5_p) / (bin2_p |bin3_p) / (bin4_p | BPA_p)

combined_plot + 
  plot_annotation(
    title = "Distribution of number of genes per BMD_norm",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

