##### BMDX ANALYSIS

#renv::install("/project/omics/public_data/endocrine_disruptors/bmdx_package")
library(bmdx)
library(ggplot2)

# READ DATA
working_dir = "/project/omics/public_data/endocrine_disruptors/GSE153320/results/bmdx/"
setwd(working_dir)
gse = "GSE153320"

# load the data 
metadata = bmdx::read_excel_allsheets(filename =  paste0(working_dir, "log_nocyt/phenodata_bmdx_noCyt.xlsx"), check_numeric = F)
experimental_data = bmdx::read_excel_allsheets(filename =  paste0(working_dir,"log_nocyt/vst_counts_bmdx_noCyt.xlsx"), first_col_as_rownames = TRUE,is_rnaseq_raw_count = FALSE,check_numeric = FALSE)

### log2 doses

log_transform = function(x){
  log2(x+1)
}
inverse_transform = function(log_x){
  (2^log_x) - 1
}

doses = metadata$TBBA$dose_amount
log_doses = log_transform(doses)
inverse_log_doses = inverse_transform(log_doses)
all.equal(doses,inverse_log_doses)

# log transform doses 

for(i in 1:length(metadata)){
  treatment = names(metadata)[i]
  metadata[[treatment]][,"dose_amount"] = log_transform(metadata[[treatment]][,"dose_amount"])
}

for(i in names(metadata)){
  print(all.equal(colnames(experimental_data[[i]]),metadata[[i]]$geo_accession))
}


# CONVERT DATA FOR MODELLING
x = "dose"
y = "expr"
sample_id_col = "geo_accession"
dose_id_col = "dose_amount"
time_col_id = "exposure_time"
other_variables_id_col = NULL

data_dictionary = create_data_structure(experimental_data,
                                        metadata,
                                        sample_id_col = sample_id_col,
                                        dose_id_col = dose_id_col,
                                        other_variables_id_col = c(time_col_id,other_variables_id_col) ,
                                        x = x,
                                        y = y)

# BUILD LIST OF MODELS FOR FITTING
model_list = build_models(model_names = c("linear","poly2","hill","exp2","power"), max_iter = 1024, data_type = "continuous", x = x, y = y)

deviation_type = "standard"
rl = 1.349
variance_type = "constant"
confidence_interval = 0.95
significance_level = 0.05 # change the description in the manual
nCores = 1
is_parallel = FALSE

all_fitted_models = fitting_list(data_dictionary,
                                 model_list,
                                 deviation_type = deviation_type,
                                 rl = rl,
                                 confidence_interval = confidence_interval,
                                 variance_type = variance_type,
                                 significance_level = significance_level,
                                 nCores = nCores,
                                 is_parallel = is_parallel)

all_stats = compute_model_statistics(fitted_models = all_fitted_models,
                                     other_variables_id_col = c(time_col_id,other_variables_id_col),
                                     nCores = 1,is_parallel = F)

loofth = lower_bound_th = upper_bound_th =  0.1
bmd_bmdl_th = bmdu_bmd_th = 20
bmdu_bmdl_th = 40
r2_th = 0.6
filter_by_monotonicity = F
bmd_na_filter = bmdl_na_filter = bmdu_na_filter = T
ic50_na_filter = F
filter_by_lack_of_fit = F
r2_filter = T
ratio_filter = filter_lower_bound = filter_upper_bound  = F

filtered_models = model_filtering(fitted_models = all_fitted_models,
                                  loofth = loofth,
                                  lower_bound_th = lower_bound_th,
                                  upper_bound_th = upper_bound_th,
                                  bmd_bmdl_th = bmd_bmdl_th,
                                  bmdu_bmd_th = bmdu_bmd_th,
                                  bmdu_bmdl_th = bmdu_bmdl_th,
                                  filter_lower_bound = filter_lower_bound,
                                  filter_upper_bound = filter_upper_bound,
                                  filter_by_lack_of_fit = filter_by_lack_of_fit,
                                  ratio_filter = ratio_filter,
                                  bmd_na_filter = bmd_na_filter,
                                  bmdl_na_filter = bmdl_na_filter,
                                  bmdu_na_filter = bmdu_na_filter,
                                  ic50_na_filter = ic50_na_filter,
                                  r2_filter = r2_filter,
                                  r2_th = r2_th,
                                  filter_by_monotonicity = filter_by_monotonicity)


all_stats_filtered = compute_model_statistics(filtered_models,
                                              other_variables_id_col = c(time_col_id,other_variables_id_col),
                                              nCores = 1)


res = select_optimal_models(filtered_models,
                            method = "AIC",
                            time_col_id = "exposure_time",
                            optional_col_ids = NULL,
                            nCores = 1)

optimal_models = res$optimal_models
optimal_models_stats = res$BMD_tab_optimal
dim(optimal_models_stats)

optimal_models_stats$BMD = inverse_transform(optimal_models_stats$BMD)
optimal_models_stats$BMDL = inverse_transform(optimal_models_stats$BMDL)
optimal_models_stats$BMDU = inverse_transform(optimal_models_stats$BMDU)
optimal_models_stats$AC50 = inverse_transform(optimal_models_stats$AC50)

save(optimal_models_stats, file = paste0(working_dir,"log_nocyt/optimal_models_stats.RData"))

# Load the package
library(writexl)
# Save as an Excel file
write_xlsx(optimal_models_stats, paste0(working_dir,"log_nocyt/optimal_models_stats.xlsx"))
write.csv(optimal_models_stats, file = paste0(working_dir,"log_nocyt/optimal_models_stats.csv"), row.names = FALSE)




#### PLOTTING

# Load the package
library(writexl)
# Save as an Excel file
write_xlsx(optimal_models_stats, paste0(working_dir,"log_nocyt/optimal_models_stats.xlsx"))
write.csv(optimal_models_stats, file = paste0(working_dir,"log_nocyt/optimal_models_stats.csv"), row.names = FALSE)




#### PLOTTING BY CHEMICALS
library(dplyr)

current_stat = optimal_models_stats
current_stat$Experiment = as.factor(current_stat$Experiment)

p = plot_histogram(current_stat, y_val = "BMD",
                   color_by = "Experiment",
                   group_by = "Experiment",
                   group_by2 = NULL,
                   filter_column = NULL,
                   filter_by = list(c("24")))+ ggplot2::theme(
                     axis.title.x = ggplot2::element_text(size = 16),
                     axis.title.y = ggplot2::element_text(size = 16),
                     axis.text.x = ggplot2::element_text(size = 10),
                     axis.text.y = ggplot2::element_text(size = 10),
                     strip.text = ggplot2::element_text(size = 16),  # group_by label size
                     legend.title = ggplot2::element_text(size = 16), # legend title
                     legend.text = ggplot2::element_text(size = 14) 
                   )

p

ggsave(p, file = paste0(working_dir,"log_nocyt/BMD_distribution.pdf"), width = 20)

p = plot_histogram(current_stat, y_val = "BMD",
                   color_by = "Experiment",
                   group_by = "exposure_time",
                   group_by2 = NULL,
                   filter_column = NULL,
                   filter_by = list(c("24")))+ ggplot2::theme(
                     axis.title.x = ggplot2::element_text(size = 16),
                     axis.title.y = ggplot2::element_text(size = 16),
                     axis.text.x = ggplot2::element_text(size = 10),
                     axis.text.y = ggplot2::element_text(size = 10),
                     strip.text = ggplot2::element_text(size = 16),  # group_by label size
                     legend.title = ggplot2::element_text(size = 16), # legend title
                     legend.text = ggplot2::element_text(size = 14) 
                   )

p

ggsave(p, file = paste0(working_dir,"log_nocyt/BMD_distribution_onePlot.pdf"), width = 20)

# p = p +   ggplot2::scale_fill_manual(values = c("#235888", "#F7BD03", "#357D8A"))
# p = p + ggplot2::xlim(c(0,30)) + ggplot2::ylim(c(0,0.3))
# p
# 
# ggsave(p, file = paste0(working_dir,"log_nocyt/BMD_distribution_blue_colors.pdf"))

current_stat2 = current_stat
current_stat2$exposure_time  = as.factor(current_stat2$exposure_time)

p=plot_scatter(current_stat2, x_val = "BMDL",y_val = "BMD",
               color_by = "Model",
               group_by = "Experiment",
               group_by2 = NULL,
               filter_column = NULL,
               filter_by = list(c("exposure_time")))
p

ggsave(p, file = paste0(working_dir,"log_nocyt/BMDL_BMD_ratio.pdf"), width = 20)


p=plot_pie_chart(current_stat, category = "Model",
                 group_by = time_col_id,
                 group_by2 = "Experiment",
                 filter_column = NULL,
                 filter_by = list(c("24")))+ 
  ggplot2::facet_wrap(~ get(time_col_id) + Experiment, ncol = 2) +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(size = 16),
    axis.title.y = ggplot2::element_text(size = 16),
    axis.text.x = ggplot2::element_text(size = 10),
    axis.text.y = ggplot2::element_text(size = 10),
    strip.text = ggplot2::element_text(size = 16),  # group_by label size
    legend.title = ggplot2::element_text(size = 16), # legend title
    legend.text = ggplot2::element_text(size = 14) 
  )
p
ggsave(p, file = paste0(working_dir,"log_nocyt/model_distribution.pdf"))


p = aggregate_rows_time(current_stat,
                        gen_feat = "Feature",
                        first_feat = "exposure_time",
                        group_by = "Experiment",
                        filter_column = NULL,
                        filter_by = NULL) + ggplot2::theme(
                          axis.title.x = ggplot2::element_text(size = 16),
                          axis.title.y = ggplot2::element_text(size = 16),
                          axis.text.x = ggplot2::element_text(size = 16),
                          axis.text.y = ggplot2::element_text(size = 16)
                        )
p

ggsave(p, file = paste0(working_dir,"log_nocyt/n_DD_genes.pdf"))


p = aggregate_rows_time(current_stat,
                        gen_feat = "Feature",
                        first_feat = time_col_id,
                        group_by = c("Experiment","Adverse_direction"),
                        filter_column = NULL,
                        filter_by = NULL) + 
  ggplot2::facet_grid(Experiment ~ Adverse_direction, scales = "free_y") +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(size = 16),
    axis.title.y = ggplot2::element_text(size = 16),
    axis.text.x = ggplot2::element_text(size = 15),
    axis.text.y = ggplot2::element_text(size = 16),
    strip.text = ggplot2::element_text(size = 16),  # group_by label size
    legend.title = ggplot2::element_text(size = 16), # legend title
    legend.text = ggplot2::element_text(size = 14) 
  )
p

ggsave(p, file = paste0(working_dir,"log_nocyt/n_DD_genes_by_direction.pdf"), width = 10)

# To see the resulting dataframe:

p = ecdf_plots (mod_stats = current_stat,
                rel_variable = "Experiment",
                group_by = "exposure_time",
                is_group_by_numeric = FALSE,
                other_variables = NULL,
                number_of_column = 2,
                scaling = F, # to use when compounds with different doses are analysed
                filter_column = NULL,
                plot_type = "ecdf",linewidth = 1.2)+ ggplot2::theme(
                  axis.title.x = ggplot2::element_text(size = 16),
                  axis.title.y = ggplot2::element_text(size = 16),
                  axis.text.x = ggplot2::element_text(size = 16),
                  axis.text.y = ggplot2::element_text(size = 16)
                )
p

ggsave(p, file = paste0(working_dir,"log_nocyt/ECDF.pdf"), width = 10)

p =  upset_plot(mod_stats = current_stat,
                rel_variable = "Experiment",
                group_by = NULL,
                other_variables = NULL,
                filter_column = NULL,
                filter_by = list(),
                nintersects = 10,
                group.by = "degree",
                order.by ="degree",text.scale = 2)

pdf(paste0(working_dir,"log_nocyt/upset.pdf"), width = 10, height = 6)
p
dev.off()


tb <- table(optimal_models_stats$Experiment)
write.table(tb, file = paste0(working_dir,"log_nocyt/num_models.tsv"), sep = "\t", col.names =TRUE)

#BPA

print(paste("BPA number of optimal models:", length(unique(optimal_models_stats$Feature[optimal_models_stats$Experiment == "BPA"]))))
#to_print = c("", "", "", "", "", "")
to_print = c("ENSG00000162618$poly2", "ENSG00000166206$power", "ENSG00000101546$linear", "ENSG00000121297$exp2", "ENSG00000215908$hill", "ENSG00000161638$poly2")

for (i in seq_along(to_print)) {
  entry = to_print[i]
  print(entry)
  parts <- strsplit(entry, "\\$")[[1]]
  gene_id <- parts[1]
  model_type <- parts[2]
  
  if(!(gene_id %in% optimal_models_stats$Feature[optimal_models_stats$Experiment == "BPA"])){
    message(paste("Gene not found:", gene_id))
    next
  }
  
  model <- optimal_models[[paste0("BPA_72_", gene_id)]][model_type]
  model <- model[[1]]
  
  if (is.null(model)){
    cat("CHANGED MODEL")
    model <- optimal_models[[paste0("BPA_72_", gene_id)]][[1]]
  }
  
  title_text <- paste("BPA", gene_id, model_type, "R2", model$R2)
  
  p <- plot_bmdx(
    model = model,
    cex = 16,
    plot_ic50 = FALSE,
    title_label = title_text
  )
  
  # p<- plot_bmdx(
  #   model = optimal_models$BPA_72_ENSG00000089847$linear,
  #   cex = 16,
  #   plot_ic50 = F,
  #   title_label = title_text
  # )
  # 
  # p<- plot_bmdx(
  #   model = optimal_models[["BPA_72_ENSG00000089847"]]$linear,
  #   cex = 16,
  #   plot_ic50 = F,
  #   title_label = title_text
  # )
  
  print(p)
  assign(paste0("p", i), p)
}

pdf(paste0(working_dir,"log_nocyt/BPA_genes.pdf"),width=20)
easyGgplot2::ggplot2.multiplot(p1, p2, p3, p4, p5, p6, cols = 3)
dev.off()


#TBBA

print(paste("TBBA number of optimal models:", length(unique(optimal_models_stats$Feature[optimal_models_stats$Experiment == "TBBA"]))))
to_print = c("ENSG00000100604$poly2", "ENSG00000121057$linear", "ENSG00000175197$exp2", "ENSG00000073578$poly2", "ENSG00000273521$power", "ENSG00000116661$linear")

for (i in seq_along(to_print)) {
  entry = to_print[i]
  print(entry)
  parts <- strsplit(entry, "\\$")[[1]]
  gene_id <- parts[1]
  model_type <- parts[2]
  
  if(!(gene_id %in% optimal_models_stats$Feature[optimal_models_stats$Experiment == "TBBA"])){
    message(paste("Gene not found:", gene_id))
    next
  }
  
  model <- optimal_models[[paste0("TBBA_72_", gene_id)]][model_type]
  model <- model[[1]]
  
  if (is.null(model)){
    cat("CHANGED MODEL")
    model <- optimal_models[[paste0("TBBA_72_", gene_id)]][[1]]
  }
  
  title_text <- paste("TBBA", gene_id, model_type, "R2", model$R2)
  
  p <- plot_bmdx(
    model = model,
    cex = 16,
    plot_ic50 = FALSE,
    title_label = title_text
  )
  
  print(p)
  assign(paste0("p", i), p)
}

pdf(paste0(working_dir,"log_nocyt/TBBA_genes.pdf"),width=20)
easyGgplot2::ggplot2.multiplot(p1, p2, p3, p4, p5, p6, cols = 3)
dev.off()

#BPS

print(paste("BPS number of optimal models:", length(unique(optimal_models_stats$Feature[optimal_models_stats$Experiment == "BPS"]))))
to_print = c("ENSG00000115457$poly2", "ENSG00000280707$power", "ENSG00000169981$hill", "ENSG00000115232$linear", "ENSG00000000457$poly2", "ENSG00000005889$exp2")

for (i in seq_along(to_print)) {
  entry = to_print[i]
  print(entry)
  parts <- strsplit(entry, "\\$")[[1]]
  gene_id <- parts[1]
  model_type <- parts[2]
  
  if(!(gene_id %in% optimal_models_stats$Feature[optimal_models_stats$Experiment == "BPS"])){
    message(paste("Gene not found:", gene_id))
    next
  }
  
  model <- optimal_models[[paste0("BPS_72_", gene_id)]][model_type]
  model <- model[[1]]
  
  if (is.null(model)){
    cat("CHANGED MODEL")
    model <- optimal_models[[paste0("BPS_72_", gene_id)]][[1]]
  }
  
  title_text <- paste("BPS", gene_id, model_type, "R2", model$R2)
  
  p <- plot_bmdx(
    model = model,
    cex = 16,
    plot_ic50 = FALSE,
    title_label = title_text
  )
  
  print(p)
  assign(paste0("p", i), p)
}

pdf(paste0(working_dir,"log_nocyt/BPS_genes.pdf"),width=20)
easyGgplot2::ggplot2.multiplot(p1, p2, p3, p4, p5, p6, cols = 3)
dev.off()

#BPF

print(paste("BPF number of optimal models:", length(unique(optimal_models_stats$Feature[optimal_models_stats$Experiment == "BPF"]))))
to_print = c("ENSG00000117868$poly2", "ENSG00000142303$power", "ENSG00000137413$hill", "ENSG00000105329$exp2", "ENSG00000134690$linear", "ENSG00000253320$poly2")

for (i in seq_along(to_print)) {
  entry = to_print[i]
  print(entry)
  parts <- strsplit(entry, "\\$")[[1]]
  gene_id <- parts[1]
  model_type <- parts[2]
  
  if(!(gene_id %in% optimal_models_stats$Feature[optimal_models_stats$Experiment == "BPF"])){
    message(paste("Gene not found:", gene_id))
    next
  }
  
  model <- optimal_models[[paste0("BPF_72_", gene_id)]][model_type]
  model <- model[[1]]
  
  if (is.null(model)){
    cat("CHANGED MODEL")
    model <- optimal_models[[paste0("BPF_72_", gene_id)]][[1]]
  }
  
  title_text <- paste("BPF", gene_id, model_type, "R2", model$R2)
  
  p <- plot_bmdx(
    model = model,
    cex = 16,
    plot_ic50 = FALSE,
    title_label = title_text
  )
  
  print(p)
  assign(paste0("p", i), p)
}

pdf(paste0(working_dir,"log_nocyt/BPF_genes.pdf"),width=20)
easyGgplot2::ggplot2.multiplot(p1, p2, p3, p4, p5, p6, cols = 3)
dev.off()



####TPOD -> transcription point of departure (when all the genes are activated somehow)
#There are different way to choose it: minimum dose, mean dose, dose where there is the max increment of response 

library(scam)

# Parameters for tPOD computation
pod_value <- "BMD" #or "BMDL", "BMDU"
percentile <- 0.20 # a number between 0 and 1
lowest_method = "lowest" #or "LCRD"

model_stats_list <- setNames(
  split(optimal_models_stats, interaction(optimal_models_stats[["Experiment"]], 
                                          optimal_models_stats[["exposure_time"]], 
                                          drop = TRUE, sep = "_")),
  levels(interaction(optimal_models_stats[["Experiment"]], 
                     optimal_models_stats[["exposure_time"]], 
                     drop = TRUE, sep = "_"))
)

tpod_methods_list <- c("percentile", "first_mode", "lowest", "accumulation")

# Compute tPOD for each experimental condition and method
result_tPOD <- lapply(model_stats_list, function(model_stats) {
  res = apply_tpod_methods(model_stats = model_stats,
                           tpod_methods_list = tpod_methods_list,
                           pod_value = "BMD",
                           percentile = percentile,
                           lowest_method = lowest_method)
  rownames(res) = res$Method
  return(res)
})


p1 = tpod_plot(pod_vector = model_stats_list$BPA_72$BMD,
               tpod_method = "accumulation",
               pod_value = "BMD",
               tPOD = as.numeric(result_tPOD$BPA_72["accumulation","tPOD"]),
               xlog = T,
               subtitle = "tPOD BPA exposure 72h") + theme(legend.position = "bottom")

p2 = tpod_plot(model_stats_list$TBBA_72$BMD,
               tpod_method = "accumulation",
               pod_value = "BMD",
               tPOD = as.numeric(result_tPOD$TBBA_72["accumulation","tPOD"]),
               xlog = T,
               subtitle = "tPOD TBBA exposure 72h")+ theme(legend.position = "bottom")

p3 = tpod_plot(model_stats_list$BPS_72$BMD,
               tpod_method = "accumulation",
               pod_value = "BMD",
               tPOD = as.numeric(result_tPOD$BPS_72["accumulation","tPOD"]),
               xlog = T,
               subtitle = "tPOD BPS exposure 72h")+ theme(legend.position = "bottom")

p4 = tpod_plot(model_stats_list$BPF_72$BMD,
               tpod_method = "accumulation",
               pod_value = "BMD",
               tPOD = as.numeric(result_tPOD$BPF_72["accumulation","tPOD"]),
               xlog = T,
               subtitle = "tPOD BPF exposure 72h")+ theme(legend.position = "bottom")

pdf(paste0(working_dir,"log_nocyt/accumulation.pdf"),width = 20)
easyGgplot2::ggplot2.multiplot(p1,p2,p3,p4, cols = 4)
dev.off()


# Plot density for BMD
library(ggplot2)

# Define your data frame with the BMD values and statistics
#BPA
BMD_values <- model_stats_list$BPA_72$BMD

lowest = as.numeric(result_tPOD$BPA_72["lowest","tPOD"])
percentile = as.numeric(result_tPOD$BPA_72["percentile","tPOD"])
mean_value = as.numeric(result_tPOD$BPA_72["mean","tPOD"])
accumulation = as.numeric(result_tPOD$BPA_72["accumulation","tPOD"])
first_mode = as.numeric(result_tPOD$BPA_72["first_mode","tPOD"])

p = plot_BMD_tPOD_density(BMD_values, lowest = lowest,
                          percentile = percentile,
                          mean_value = mean_value,
                          accumulation = accumulation,
                          first_mode = first_mode)
print(p)
ggsave(p, file = paste0(working_dir,"log_nocyt/tpod_on_distribution_BPA.pdf"))

plot_BMD_tPOD_density(BMD_values, lowest = lowest,
                      percentile = percentile,
                      mean_value = mean_value,
                      accumulation = accumulation,
                      first_mode = first_mode,xrange = c(0,5))


plot_BMD_tPOD_density(BMD_values, lowest = lowest,
                      percentile = percentile,
                      mean_value = NA,
                      accumulation = accumulation,
                      first_mode = first_mode)

plot_BMD_tPOD_density(BMD_values = BMD_values, 
                      lowest = NA,
                      percentile = NA,
                      mean_value = NA,
                      accumulation = NA,
                      first_mode = NA,
                      xrange = NULL)


#BPS
BMD_values <- model_stats_list$BPS_72$BMD

lowest = as.numeric(result_tPOD$BPS_72["lowest","tPOD"])
percentile = as.numeric(result_tPOD$BPS_72["percentile","tPOD"])
mean_value = as.numeric(result_tPOD$BPS_72["mean","tPOD"])
accumulation = as.numeric(result_tPOD$BPS_72["accumulation","tPOD"])
first_mode = as.numeric(result_tPOD$BPS_72["first_mode","tPOD"])

p = plot_BMD_tPOD_density(BMD_values, lowest = lowest,
                          percentile = percentile,
                          mean_value = mean_value,
                          accumulation = accumulation,
                          first_mode = first_mode)
print(p)
ggsave(p, file = paste0(working_dir,"log_nocyt/tpod_on_distribution_BPS.pdf"))

plot_BMD_tPOD_density(BMD_values, lowest = lowest,
                      percentile = percentile,
                      mean_value = mean_value,
                      accumulation = accumulation,
                      first_mode = first_mode,xrange = c(0,5))


plot_BMD_tPOD_density(BMD_values, lowest = lowest,
                      percentile = percentile,
                      mean_value = NA,
                      accumulation = accumulation,
                      first_mode = first_mode)

plot_BMD_tPOD_density(BMD_values = BMD_values, 
                      lowest = NA,
                      percentile = NA,
                      mean_value = NA,
                      accumulation = NA,
                      first_mode = NA,
                      xrange = NULL)

#BPF
BMD_values <- model_stats_list$BPF_72$BMD

lowest = as.numeric(result_tPOD$BPF_72["lowest","tPOD"])
percentile = as.numeric(result_tPOD$BPF_72["percentile","tPOD"])
mean_value = as.numeric(result_tPOD$BPF_72["mean","tPOD"])
accumulation = as.numeric(result_tPOD$BPF_72["accumulation","tPOD"])
first_mode = as.numeric(result_tPOD$BPF_72["first_mode","tPOD"])

p = plot_BMD_tPOD_density(BMD_values, lowest = lowest,
                          percentile = percentile,
                          mean_value = mean_value,
                          accumulation = accumulation,
                          first_mode = first_mode)
print(p)
ggsave(p, file = paste0(working_dir,"log_nocyt/tpod_on_distribution_BPF.pdf"))

plot_BMD_tPOD_density(BMD_values, lowest = lowest,
                      percentile = percentile,
                      mean_value = mean_value,
                      accumulation = accumulation,
                      first_mode = first_mode,xrange = c(0,5))


plot_BMD_tPOD_density(BMD_values, lowest = lowest,
                      percentile = percentile,
                      mean_value = NA,
                      accumulation = accumulation,
                      first_mode = first_mode)

plot_BMD_tPOD_density(BMD_values = BMD_values, 
                      lowest = NA,
                      percentile = NA,
                      mean_value = NA,
                      accumulation = NA,
                      first_mode = NA,
                      xrange = NULL)

#TBBA
BMD_values <- model_stats_list$TBBA_72$BMD

lowest = as.numeric(result_tPOD$TBBA_72["lowest","tPOD"])
percentile = as.numeric(result_tPOD$TBBA_72["percentile","tPOD"])
mean_value = as.numeric(result_tPOD$TBBA_72["mean","tPOD"])
accumulation = as.numeric(result_tPOD$TBBA_72["accumulation","tPOD"])
first_mode = as.numeric(result_tPOD$TBBA_72["first_mode","tPOD"])

p = plot_BMD_tPOD_density(BMD_values, lowest = lowest,
                          percentile = percentile,
                          mean_value = mean_value,
                          accumulation = accumulation,
                          first_mode = first_mode)
print(p)
ggsave(p, file = paste0(working_dir,"log_nocyt/tpod_on_distribution_TBBA.pdf"))

plot_BMD_tPOD_density(BMD_values, lowest = lowest,
                      percentile = percentile,
                      mean_value = mean_value,
                      accumulation = accumulation,
                      first_mode = first_mode,xrange = c(0,5))


plot_BMD_tPOD_density(BMD_values, lowest = lowest,
                      percentile = percentile,
                      mean_value = NA,
                      accumulation = accumulation,
                      first_mode = first_mode)

plot_BMD_tPOD_density(BMD_values = BMD_values, 
                      lowest = NA,
                      percentile = NA,
                      mean_value = NA,
                      accumulation = NA,
                      first_mode = NA,
                      xrange = NULL)



#before and after tpod accumulation
source("/project/omics/public_data/endocrine_disruptors/GSE153320/scripts/bmdx_noshiny-master/supporting_functions/enrichment_optimalstats_filtering.R")
library(scam)
opt_tpod_acc <- filter_tpod_acc(optimal_models_stats, timepoints_col_variable = "exposure_time")
ke_enrich_tpod_acc <- ke_enrichment_optimalstats_filtering(optimal_models_filtered = opt_tpod_acc,parameter = "tpod acc",time_var = "exposure_time")
# gost_enrich_tpod_acc<- gost_enrichment_optimalstats_filtering(optimal_models_filtered = opt_tpod_acc)

p1 = ke_enrich_tpod_acc$p_ngene + theme(
  axis.title.x = ggplot2::element_text(size = 16),
  axis.title.y = ggplot2::element_text(size = 16),
  axis.text.x = ggplot2::element_text(size = 15),
  axis.text.y = ggplot2::element_text(size = 16),
  strip.text = ggplot2::element_text(size = 16),  # group_by label size
  legend.title = ggplot2::element_text(size = 16), # legend title
  legend.text = ggplot2::element_text(size = 14) 
) + ggplot2::scale_fill_manual(values = c("#235888", "#F7BD03"))+
  # theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  ggtitle(label  = "Number of genes below accumulation tPOD")


p2 = ke_enrich_tpod_acc$p_nke + theme(
  axis.title.x = ggplot2::element_text(size = 16),
  axis.title.y = ggplot2::element_text(size = 16),
  axis.text.x = ggplot2::element_text(size = 15),
  axis.text.y = ggplot2::element_text(size = 16),
  strip.text = ggplot2::element_text(size = 16),  # group_by label size
  legend.title = ggplot2::element_text(size = 16), # legend title
  legend.text = ggplot2::element_text(size = 14) 
) + ggplot2::scale_fill_manual(values = c("#235888", "#F7BD03"))+
  ggtitle(label  = "Number of Key Events below accumulation tPOD")


##Before and after antimode
opt_antimode <- filter_antimode(optimal_models_stats,timepoints_col_variable = "exposure_time")
ke_enrich_antimode <- ke_enrichment_optimalstats_filtering(optimal_models_filtered = opt_antimode,parameter = "antimode")

p3 = ke_enrich_antimode$p_ngene + theme(
  axis.title.x = ggplot2::element_text(size = 16),
  axis.title.y = ggplot2::element_text(size = 16),
  axis.text.x = ggplot2::element_text(size = 15),
  axis.text.y = ggplot2::element_text(size = 16),
  strip.text = ggplot2::element_text(size = 16),  # group_by label size
  legend.title = ggplot2::element_text(size = 16), # legend title
  legend.text = ggplot2::element_text(size = 14) 
) + ggplot2::scale_fill_manual(values = c("#235888", "#F7BD03"))+
  # theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  ggtitle(label  = "Number of genes below antimode")

p4 = ke_enrich_antimode$p_nke + theme(
  axis.title.x = ggplot2::element_text(size = 16),
  axis.title.y = ggplot2::element_text(size = 16),
  axis.text.x = ggplot2::element_text(size = 15),
  axis.text.y = ggplot2::element_text(size = 16),
  strip.text = ggplot2::element_text(size = 16),  # group_by label size
  legend.title = ggplot2::element_text(size = 16), # legend title
  legend.text = ggplot2::element_text(size = 14) 
) + ggplot2::scale_fill_manual(values = c("#235888", "#F7BD03"))+
  ggtitle(label  = "Number of Key Events below antimode")

easyGgplot2::ggplot2.multiplot(p1,p2,p3,p4,cols = 2)

library(patchwork)

# Arrange plots with varying column widths
layout <- (p1 + p2 + plot_layout(widths = c(1, 2))) / 
  (p3 + p4 + plot_layout(widths = c(1, 2)))

# Display the combined plot
layout

pdf(paste0(working_dir, "log_nocyt/tpod_summary.pdf"),width = 15,height = 10)
easyGgplot2::ggplot2.multiplot(p1,p2,p3,p4,cols = 2)
dev.off()


##########################################################

### AOP enrichment

# Sys.setenv(GITHUB_PAT = "il tuo token github")
# renv::install("https://github.com/fhaive/AOPfingerprintR")

library(AOPfingerprintR)
only_significant = T
pval_th = 0.05
adj.method = "fdr"

experiment_var = "Experiment"
time_var = "exposure_time"
other_variables_id_col = c(time_var)

BMD_TAB = optimal_models_stats

experiments = unique(BMD_TAB$Experiment)
tp = unique(BMD_TAB$exposure_time)
GList = list()

# Made a list for each chemical with gene and doses (BMD, BMDL e BMDU) saved
for (expi in experiments) {
  for (tpi in tp) {
    idx = intersect(which(BMD_TAB$Experiment == expi), which(BMD_TAB$exposure_time == tpi))
    GList[[paste(expi,tpi,sep = "_")]] = BMD_TAB[idx,]
  }
}


# For each list make KE and AOP enrichment 
original_gene_set = rownames(experimental_data$BPA)
# For each GSE you have to define which are the genes in the platform, that is your background distribution, because they are the only ones that have a probability to be enriched 
background = original_gene_set 
ke_enrichment_results = enrich_KEs_AOPs(GList = GList,
                                        list_gene_sets = human_ens_clusters,
                                        only_significant = T,
                                        pval_th = 0.05,
                                        adj.method = "fdr",
                                        merge_by = "Ke",
                                        numerical_properties = c("BMD","BMDL","BMDU"),background = background)

source("/project/omics/public_data/endocrine_disruptors/GSE153320/scripts/bmdx_noshiny-master/supporting_functions/aop_bmd_functions.R")

res = plot_bmd_boxplots_by_ke_type(ke_enrichment_results,
                                   fill_colors = c(MIE = "#235888", KE = "#F7BD03", AO = "#357D8A"),
                                   text_size = 16,
                                   facet_scales = "free",
                                   angle_x = 0)

pbmdl = res$BMDL
pbmd= res$BMD
pbmdu = res$BMDU

pdf(paste0(working_dir, "log_nocyt/boxplot_BMD_by_MIE_KE_AO.pdf"),width = 15,height = 10)
easyGgplot2::ggplot2.multiplot(pbmdl,pbmd, pbmdu, cols = 1)
dev.off()

source("/project/omics/public_data/endocrine_disruptors/GSE153320/scripts/bmdx_noshiny-master/supporting_functions/aop_bmd_functions.R")

p <- plot_bmd_by_ke_and_level(
  ke_enrichment_results,
  Biological_system_annotations,
  text_size = 14,
  angle_x = 30,
  facet_scales = "free_y"
)

print(p)
ggsave(p, file = paste0(working_dir, "log_nocyt/boxplot_BMD_by_ke_type.pdf"), width = 10, height = 5)

source("/project/omics/public_data/endocrine_disruptors/GSE153320/scripts/bmdx_noshiny-master/supporting_functions/aop_bmd_functions.R")
p <- ke_bmd_distribution_by_system(ke_enrichment_results, 
                                   Biological_system_annotations,
                                   c("#235888", "#F7BD03", "#357D8A", "#FFDD00", 
                                     "#68A5D4", "#FFB600", "#D9D7CD", "#641DFF",
                                     "#0022d2", "#505a74", "#bf212f"))
print(p)
ggsave(p, file = paste0(working_dir, "log_nocyt/boxplot_BMD_by_ke_system.pdf"), width = 20, height = 10)

background = original_gene_set
aop_enrichment_results = AOPfingerprintR::enrich_KEs_AOPs(GList = GList,
                                                          list_gene_sets = human_ens_aop,
                                                          only_significant = T,
                                                          pval_th = 0.05,
                                                          adj.method = "fdr",
                                                          merge_by = "Aop",
                                                          numerical_properties = c("BMD","BMDL","BMDU"),
                                                          background = background)

p <- plot_aop_distribution_over_ssbd_categories(
  aop_enrichment_results,
  Annotate_AOPs,
  fill_colors = c("#235888", "#F7BD03", "#357D8A", "#bf212f"),
  text_size = 14
)

print(p)
ggsave(p, file = paste0(working_dir,"log_nocyt/aop_by_jrc_category.pdf"), width = 15, height = 8)

library(dplyr)
library(tidyr)
p <- plot_average_bmd_in_aop_of_a_ssbd_category(
  aop_enrichment_results,
  Annotate_AOPs,
  ssbd_category = "Carcinogenicity",
  low_color = "#235888",
  high_color = "#F7BD03",
  text_size = 14
)
print(p)
ggsave(p, file = paste0(working_dir,"log_nocyt/aop_carcinogenicity.pdf"), width = 10, height = 8)

p <- plot_average_bmd_in_aop_of_a_ssbd_category(
  aop_enrichment_results,
  Annotate_AOPs,
  ssbd_category = "STO_tox_Lung",
  low_color = "#235888",
  high_color = "#F7BD03",
  text_size = 14
)
print(p)
ggsave(p, file = paste0(working_dir,"log_nocyt/aop_STO_tox_Lung.pdf"), width = 10, height = 6)

p <- plot_average_bmd_in_aop_of_a_ssbd_category(
  aop_enrichment_results,
  Annotate_AOPs,
  ssbd_category = "STO_tox_Kidney",
  low_color = "#235888",
  high_color = "#F7BD03",
  text_size = 10
)
print(p)
ggsave(p, file = paste0(working_dir,"log_nocyt/aop_STO_tox_Kidney.pdf"))

p <- plot_average_bmd_in_aop_of_a_ssbd_category(
  aop_enrichment_results,
  Annotate_AOPs,
  ssbd_category = "Endocrine disruption (human health)",
  low_color = "#235888",
  high_color = "#F7BD03",
  text_size = 10
)
print(p)
ggsave(p, file = paste0(working_dir,"log_nocyt/aop_ED.pdf"))

# library(stringr)
# Annotate_AOPs2 <- Annotate_AOPs
# Annotate_AOPs2$SSbD_category <- str_wrap(Annotate_AOPs$SSbD_category, width = 20)

p <- boxplot_bmd_distribution_across_ssbd_categories(
  aop_enrichment_results,
  Annotate_AOPs,
  text_size = 14,
  angle_x = 45,
  fill_colors = c("#235888", "#F7BD03", "#D9D7CD",
                  "navyblue", "brown", "orange",
                  "#357D8A", "royalblue", "#68A5D4","#006f3c",
                  "#999999", "#E69F00", "#56B4E9",
                  "#0022d2", "#505a74", "#bf212f")
)
print(p)
ggsave(p, file = paste0(working_dir,"log_nocyt/BMD_distribution_by_ssbd_category.pdf"), width=30, height=18)

# Save as an Excel file
library(writexl)
write_xlsx(ke_enrichment_results, paste0(working_dir,"log_nocyt/ke_enrichment_results.xlsx"))
write_xlsx(aop_enrichment_results, paste0(working_dir,"log_nocyt/aop_enrichment_results.xlsx"))



#range plot
group_by = "Experiment"
group_by2 = NULL
filter_column = "Experiment"
filter_by = c("BPA_72", "TBBA_72", "BPS_72", "BPF_72")
is_group_by_numeric = F


source("/project/omics/public_data/endocrine_disruptors/GSE153320/scripts/bmdx_noshiny-master/supporting_functions/render_range_plot.R")
p = render_range_plot(enrichement_data = ke_enrichment_results,
                      group_by, group_by2,
                      filter_column=NULL,
                      filter_by,
                      is_group_by_numeric,
                      display = "Ke")
p
ggsave(p, file = paste0(working_dir,"log_nocyt/render_range_plot.pdf"), width=30, height=18)


#filtering: aop almeno 5 KE e almeno il 60% dei KE è enriched
res = AOPfingerprintR::build_aop_for_aop_fingeprints(aop_enrichment_results,
                                                     ke_enrichment_results,
                                                     min_aop_length = 5,
                                                     percentage_enriched_ke = 0.33)
#percentage_enriched_ke = 0.66 se i risultati con 33 sono troppi
#percentage_enriched_ke = 0.33 se i risultati con 66 sono tropo pochi --> only 2 results

enrichement_data = res$detailed_results_only_enriched
write_xlsx(enrichement_data, paste0(working_dir,"log_nocyt/aopfingerprint_enrichment_results.xlsx"))

source("/project/omics/public_data/endocrine_disruptors/GSE153320/scripts/bmdx_noshiny-master/supporting_functions/aop_bmd_functions.R")
# plot barplot of frequecny of AOP over ssbd categories
p <- plot_aop_distribution_over_ssbd_categories(
  enrichement_data,
  Annotate_AOPs,
  fill_colors = c("#235888", "#F7BD03", "#357D8A", "#bf212f"),
  text_size = 14
)
print(p)
ggsave(p, file = "log_nocyt/aop_fingerprint_by_jrc_category.pdf", width = 15, height = 6)

idx = which(enrichement_data$SSbD_category == "Specific target organ toxicity")
enrichement_data$Organ[is.na(enrichement_data$Organ)] = "Unspecified"
enrichement_data$SSbD_category[idx] = paste("STO_tox",enrichement_data$Organ[idx], sep = "_")
categories_count = table(enrichement_data$SSbD_category)

text_cex = 15
x_axis_var = "Experiment"
library(ggplot2)
p = AOPfingerprintR::render_aop_fingerprint_bubble_plot(enrichement_data = enrichement_data,
                                                        group_by = NULL,
                                                        group_by2 = NULL,
                                                        x_axis_var = x_axis_var,
                                                        y_axis_var = "AOP_name",
                                                        filter_column = NULL,
                                                        filter_by = "STO_tox_Lung",
                                                        is_group_by_numeric = F,
                                                        threshold_proportion = 0.33,
                                                        text_cex = text_cex,
                                                        group_AOPs = NULL)
p

ggsave(p, file = paste0(working_dir,"log_nocyt/aop_fingerprint.pdf"), width = 20, height = 12)


detailed_results = ke_enrichment_results #KE enrichment only

ke_id = "TermID"
numerical_variables = c("BMDL","BMD","BMDU")
pval_variable = "padj"
gene_variable = "Genes"
experiment = "BPA_72"
enlarge_ke_selection = T
group_by = "ssbd"

nodes_edges = AOPfingerprintR::make_visNetwork(detailed_results,
                                               experiment = experiment,
                                               enlarge_ke_selection = enlarge_ke_selection,
                                               ke_id, numerical_variables = numerical_variables,
                                               pval_variable,
                                               gene_variable,
                                               max_path_length = 7,
                                               n_AOs = 7,
                                               n_MIEs = 4)


vn = AOPfingerprintR::plot_visNetwork(nodes = nodes_edges$nodes,
                                      edges = nodes_edges$edges,
                                      group_by = group_by,
                                      numerical_variables = numerical_variables)

vn
saveRDS(vn,file=paste0(working_dir,"log_nocyt/network_",experiment,".rds"))

experiment = "BPS_72"
nodes_edges = AOPfingerprintR::make_visNetwork(detailed_results,
                                               experiment = experiment,
                                               enlarge_ke_selection = enlarge_ke_selection,
                                               ke_id, numerical_variables = numerical_variables,
                                               pval_variable,
                                               gene_variable,
                                               max_path_length = 7,
                                               n_AOs = 7,
                                               n_MIEs = 4)


vn = AOPfingerprintR::plot_visNetwork(nodes = nodes_edges$nodes,
                                      edges = nodes_edges$edges,
                                      group_by = group_by,
                                      numerical_variables = numerical_variables)

vn
saveRDS(vn,file=paste0(working_dir,"log_nocyt/network_",experiment,".rds"))

experiment = "BPF_72"
nodes_edges = AOPfingerprintR::make_visNetwork(detailed_results,
                                               experiment = experiment,
                                               enlarge_ke_selection = enlarge_ke_selection,
                                               ke_id, numerical_variables = numerical_variables,
                                               pval_variable,
                                               gene_variable,
                                               max_path_length = 7,
                                               n_AOs = 7,
                                               n_MIEs = 4)


vn = AOPfingerprintR::plot_visNetwork(nodes = nodes_edges$nodes,
                                      edges = nodes_edges$edges,
                                      group_by = group_by,
                                      numerical_variables = numerical_variables)

vn
saveRDS(vn,file=paste0(working_dir,"log_nocyt/network_",experiment,".rds"))

experiment = "TBBA_72"
nodes_edges = AOPfingerprintR::make_visNetwork(detailed_results,
                                               experiment = experiment,
                                               enlarge_ke_selection = enlarge_ke_selection,
                                               ke_id, numerical_variables = numerical_variables,
                                               pval_variable,
                                               gene_variable,
                                               max_path_length = 7,
                                               n_AOs = 7,
                                               n_MIEs = 4)


vn = AOPfingerprintR::plot_visNetwork(nodes = nodes_edges$nodes,
                                      edges = nodes_edges$edges,
                                      group_by = group_by,
                                      numerical_variables = numerical_variables)

vn
saveRDS(vn,file=paste0(working_dir,"log_nocyt/network_",experiment,".rds"))
