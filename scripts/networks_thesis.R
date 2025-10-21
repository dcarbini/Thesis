library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ComplexHeatmap)
library(circlize) 
library(writexl)
library(stringr)
library(AOPfingerprintR)
library(igraph)
library(visNetwork)

dir_path <- "C:/Users/Utente/Desktop/denys/results/networks/tesi/"

#upload the file with all the KE enriched from all the experiments
data <- as.data.frame(read_xlsx("C:/Users/Utente/Desktop/denys/results/ke_enrichment/ke_nonunique.xlsx"))
data$BMDL <- as.numeric(data$BMDL)
data$BMDU <- as.numeric(data$BMDU)
data$BMD_norm <- as.numeric(data$BMD_norm)

##### Create the general network #####

agg_fun <- function(column) {
  if (is.numeric(column)) {
    return(median(column, na.rm = TRUE))
  } else {
    return(paste(unique(column[!is.na(column)]), collapse = ","))
  }
}

data_net <- data %>% distinct(TermID, .keep_all = TRUE) 
data_net$Experiment <- "all"

#Create a unique network
ke_id = "TermID"
numerical_variables = c("BMD_norm")
pval_variable = "padj"
gene_variable = "Genes"
enlarge_ke_selection = F #only the KE that are present in my file

nodes_edges = AOPfingerprintR::make_visNetwork(as.data.frame(data_net),
                                               experiment="all",
                                               enlarge_ke_selection = enlarge_ke_selection,
                                               ke_id, numerical_variables = numerical_variables,
                                               pval_variable,
                                               gene_variable,
                                               max_path_length=3,
                                               n_AOs = 2,
                                               n_MIEs = 2)

df_grouped <- data %>%
  mutate(chemical = gsub(",", "-", chemical)) %>%
  group_by(TermID) %>%
  summarise(across(everything(), agg_fun), .groups = "drop")

nodes_edges$nodes <- nodes_edges$nodes%>% 
  select(id) %>% 
  left_join(df_grouped, by = c(id = "TermID"))


##### FUNCTIONS #####
#Function to groupby aop and colored by number of aops
plot_grouped_aop <- function(nodes,
                             edges,
                             numerical_variables,
                             chemicals_column = "chemicals",
                             aop_column = "Aop") {
  
  if (!requireNamespace("visNetwork", quietly = TRUE))
    stop("Package 'visNetwork' is required but not installed.")
  
  # -- COUNT CHEMICALS PER NODE
  chem_counts <- sapply(
    strsplit(as.character(nodes[[chemicals_column]]), ","),
    function(x) length(unique(trimws(x[x != ""])))
  )
  
  # -- COUNT AOP PER NODE
  nodes_aop_lists <- strsplit(as.character(nodes[[aop_column]]), ",")
  nodes_aop_lists <- lapply(nodes_aop_lists, function(x) trimws(x[x != ""]))
  
  aop_counts <- sapply(nodes_aop_lists, function(x) length(unique(trimws(x[x != ""]))))
  
  
  # -- COLOR BY AOP cCOUNT
  aop_counts <- sapply(nodes_aop_lists, length)
  nodes$color <- ifelse(
    aop_counts == 1, "#caffbf",
    ifelse(aop_counts == 2, "#a5d8ff",
           ifelse(aop_counts == 3, "#ffd6a5",
                  "#A9A9A9"))
  )
  
  
  # -- SHAPE BY KE TYPE
  nd <- nodes$id
  ke_type <- character(length(nd))
  for (n in seq_along(nd)) {
    ke <- nd[n]
    ke_annotation <- aop_ke_table_hure$Ke_type[aop_ke_table_hure$Ke == ke]
    if (!"MolecularInitiatingEvent" %in% ke_annotation &
        !"AdverseOutcome"          %in% ke_annotation) {
      ke_type[n] <- "KeyEvent"
    } else {
      tab <- table(ke_annotation)
      if ("KeyEvent" %in% names(tab)) tab <- tab[names(tab) != "KeyEvent"]
      ke_type[n] <- names(sort(tab, decreasing = TRUE))[1]
    }
  }
  my_shapes <- c("MolecularInitiatingEvent" = "triangle",
                 "KeyEvent"                = "square",
                 "AdverseOutcome"          = "star")
  nodes$shape <- my_shapes[ke_type]
  
  # -- GROUP BY AOP
  nodes$group <- sapply(
    strsplit(as.character(nodes[[aop_column]]), ","),
    function(x) paste(trimws(x), collapse = ", ")
  )
  
  # -- NODE SIZE
  #nodes$size <- chem_counts * 2
  
  # -- ADD LABEL
  nodes$label <- nodes$Ke_description
  
  # -- TOOLTIP
  nodes$title <- paste0(
    "<p>Id: ", nodes$id, "</p>",
    "<p>Description: ", nodes$Ke_description, "</p>"
  )
  for (nv in numerical_variables) {
    nodes$title <- paste0(nodes$title, "<p>", nv, ": ", nodes[[nv]], "</p>")
  }
  nodes$title <- paste0(
    nodes$title,
    "<p># Chemicals: ", chem_counts, "</p>",
    "<p>Chemicals: ", nodes[[chemicals_column]], "</p>",
    "<p># AOPs: ", aop_counts, "</p>",
    "<p>AOPs: ", nodes$group, "</p>"
  )
  
  # -- LEGEND
  legend_nodes <- list(
    list(label = "1 AOP",  color = "#caffbf"),
    list(label = "2 AOP",  color = "#a5d8ff"),
    list(label = "3 AOP",  color = "#ffd6a5"),
    list(label = ">3 AOP", color = "#A9A9A9"),
    list(label = "MolecularInitiatingEvent", shape = "triangle"),
    list(label = "KeyEvent",                shape = "square"),
    list(label = "AdverseOutcome",          shape = "star")
  )
  
  # -- PLOT
  visNetwork::visNetwork(nodes, edges) %>% 
    visNetwork::visLegend(position = "right",
                          addNodes  = legend_nodes,
                          useGroups = FALSE) %>% 
    visNetwork::visOptions(selectedBy = list(variable = "group",
                                             multiple = TRUE))
}

#Function to show the network and color it by number of AOP, grouping by chemical and the dimension of the node is the number of chemicals that enriched that KE. 
try_bmd <- function(nodes,
                    edges,
                    numerical_variables,
                    time_column      = "time",
                    chemicals_column = "chemicals",
                    aop_column = "Aop") {
  
  if (!requireNamespace("visNetwork", quietly = TRUE))
    stop("Package 'visNetwork' is required but not installed.")
  
  # -- COUNT CHEMICALS PER NODE 
  chem_counts <- sapply(
    strsplit(as.character(nodes[[chemicals_column]]), ","),
    function(x) length(unique(trimws(x[x != ""])))
  )
  
  # -- COUNT AOP PER NODE
  aop_counts <- sapply(
    strsplit(as.character(nodes[[aop_column]]), ","),
    function(x) length(unique(trimws(x[x != ""])))
  )
  
  # -- COLOR BY NUMBER OF CHEMICALS
  nodes$color <- ifelse(
    chem_counts == 1,                      "#ffd6a5",
    ifelse(chem_counts <= 5,               "#f28b82",
           ifelse(chem_counts <= 10,            "#caffbf",
                  ifelse(chem_counts <= 15,          "#a5d8ff",
                         "#ffc6ff")
           )
    )
  )
  
  # -- SHAPE BY KE TYPE
  nd <- nodes$id
  ke_type <- character(length(nd))
  for (n in seq_along(nd)) {
    ke <- nd[n]
    ke_annotation <- aop_ke_table_hure$Ke_type[aop_ke_table_hure$Ke == ke]
    if (!"MolecularInitiatingEvent" %in% ke_annotation &
        !"AdverseOutcome"          %in% ke_annotation) {
      ke_type[n] <- "KeyEvent"
    } else {
      tab <- table(ke_annotation)
      if ("KeyEvent" %in% names(tab)) tab <- tab[names(tab) != "KeyEvent"]
      ke_type[n] <- names(sort(tab, decreasing = TRUE))[1]
    }
  }
  my_shapes <- c("MolecularInitiatingEvent" = "triangle",
                 "KeyEvent"                = "square",
                 "AdverseOutcome"          = "star")
  nodes$shape <- my_shapes[ke_type]
  
  # -- GROUP BY CHEMICALS
  nodes$group <- sapply(
    strsplit(as.character(nodes[[chemicals_column]]), ","),
    function(x) paste(trimws(x), collapse = ", ")
  )
  
  # -- NODE SIZE
  nodes$size <- (exp((nodes$BMD_norm)*3))
  
  # -- ADD LABEL
  nodes$label <- nodes$Ke_description
  
  # -- TOOLTIP
  nodes$title <- paste0(
    "<p>Id: ", nodes$id, "</p>",
    "<p>Description: ", nodes$Ke_description, "</p>"
  )
  for (nv in numerical_variables) {
    nodes$title <- paste0(nodes$title, "<p>", nv, ": ", nodes[[nv]], "</p>")
  }
  nodes$title <- paste0(
    nodes$title,
    "<p># Chemicals: ", chem_counts, "</p>",
    "<p>Chemicals: ", nodes$group, "</p>",
    "<p># AOPs: ", aop_counts, "</p>",
    "<p>AOPs: ", nodes[[aop_column]], "</p>"
  )
  
  # -- LEGEND
  legend_nodes <- list(
    list(label = "1 chemical",       color = "#ffd6a5"),
    list(label = "2–5 chemicals",    color = "#f28b82"),
    list(label = "6–10 chemicals",   color = "#caffbf"),
    list(label = "11–15 chemicals",  color = "#a5d8ff"),
    list(label = ">15 chemicals",    color = "#ffc6ff"),
    list(label = "MolecularInitiatingEvent", shape = "triangle"),
    list(label = "KeyEvent",               shape = "square"),
    list(label = "AdverseOutcome",         shape = "star")
  )
  
  # -- PLOT
  visNetwork::visNetwork(nodes, edges) %>% 
    visNetwork::visLegend(position = "right",
                          addNodes  = legend_nodes,
                          useGroups = FALSE) %>% 
    visNetwork::visOptions(selectedBy = list(variable = "group",
                                             multiple = TRUE))
}

#Function to show the network and color it by number of AOP, grouping by chemical and the dimension of the node is the number of chemicals that enriched that KE. 
plot_aop <- function(nodes,
                     edges,
                     numerical_variables,
                     time_column      = "time",
                     chemicals_column = "chemicals",
                     aop_column = "Aop") {
  
  if (!requireNamespace("visNetwork", quietly = TRUE))
    stop("Package 'visNetwork' is required but not installed.")
  
  #COUNT CHEMICALS PER NODE
  chem_counts <- sapply(
    strsplit(as.character(nodes[[chemicals_column]]), ","),
    function(x) length(unique(trimws(x[x != ""])))
  )
  
  # COUNT AOP PER NODE 
  aop_counts <- sapply(
    strsplit(as.character(nodes[[aop_column]]), ","),
    function(x) length(unique(trimws(x[x != ""])))
  )
  
  # COLOR BY NUMBER OF AOP 
  nodes$color <- ifelse(
    aop_counts == 1, "#caffbf",
    ifelse(aop_counts == 2, "#a5d8ff",
           ifelse(aop_counts == 3, "#ffd6a5",
                  "#A9A9A9"))
  )
  
  # SHAPE BY KE TYPE
  nd <- nodes$id
  ke_type <- character(length(nd))
  for (n in seq_along(nd)) {
    ke <- nd[n]
    ke_annotation <- aop_ke_table_hure$Ke_type[aop_ke_table_hure$Ke == ke]
    if (!"MolecularInitiatingEvent" %in% ke_annotation &
        !"AdverseOutcome"          %in% ke_annotation) {
      ke_type[n] <- "KeyEvent"
    } else {
      tab <- table(ke_annotation)
      if ("KeyEvent" %in% names(tab)) tab <- tab[names(tab) != "KeyEvent"]
      ke_type[n] <- names(sort(tab, decreasing = TRUE))[1]
    }
  }
  my_shapes <- c("MolecularInitiatingEvent" = "triangle",
                 "KeyEvent"                = "square",
                 "AdverseOutcome"          = "star")
  nodes$shape <- my_shapes[ke_type]
  
  # -- GROUP BY CHEMICALS
  nodes$group <- sapply(
    strsplit(as.character(nodes[[chemicals_column]]), ","),
    function(x) paste(trimws(x), collapse = ", ")
  )
  
  # -- NODE SIZE
  #nodes$size <- chem_counts * 2
  
  # -- ADD LABEL
  nodes$label <- nodes$Ke_description
  
  # -- TOOLTIP
  nodes$title <- paste0(
    "<p>Id: ", nodes$id, "</p>",
    "<p>Description: ", nodes$Ke_description, "</p>"
  )
  for (nv in numerical_variables) {
    nodes$title <- paste0(nodes$title, "<p>", nv, ": ", nodes[[nv]], "</p>")
  }
  nodes$title <- paste0(
    nodes$title,
    "<p># Chemicals: ", chem_counts, "</p>",
    "<p>Chemicals: ", nodes$group, "</p>",
    "<p># AOPs: ", aop_counts, "</p>",
    "<p>AOPs: ", nodes[[aop_column]], "</p>"
  )
  
  # -- LEGEND
  legend_nodes <- list(
    list(label = "1 AOP",  color = "#caffbf"),
    list(label = "2 AOP",  color = "#a5d8ff"),
    list(label = "3 AOP",  color = "#ffd6a5"),
    list(label = ">3 AOP", color = "#A9A9A9"),
    list(label = "MolecularInitiatingEvent", shape = "triangle"),
    list(label = "KeyEvent",                shape = "square"),
    list(label = "AdverseOutcome",          shape = "star")
  )
  
  # -- PLOT
  visNetwork::visNetwork(nodes, edges) %>% 
    visNetwork::visLegend(position = "right",
                          addNodes  = legend_nodes,
                          useGroups = FALSE) %>% 
    visNetwork::visOptions(selectedBy = list(variable = "group",
                                             multiple = TRUE))
}


#Function to plot only complete aop
plot_complete_aop <- function(nodes,
                              edges,
                              numerical_variables,
                              time_column      = "time",
                              chemicals_column = "chemicals",
                              aop_column = "Aop") {
  
  #FIND COMPLETE AOPs 
  aops_all <- unique(aop_ke_table_hure$Aop)
  
  complete_aops <- sapply(aops_all, function(aop) {
    ke_expected <- unique(aop_ke_table_hure$Ke[aop_ke_table_hure$Aop == aop])
    #ke_present  <- unique(nodes$id[grepl(aop, nodes$Aop)])
    ke_present  <- unique(nodes$id[grepl(aop, nodes$Aop)])
    all(ke_expected %in% ke_present)
  })
  
  complete_aops <- aops_all[complete_aops]
  
  if (length(complete_aops) == 0) {
    stop("No complete AOPs found.")
  }
  
  message("AOP complete found: ", paste(complete_aops, collapse = ", "))
  
  #FILTER NODES that belongs to the complete AOPs 
  nodes$Aop_list <- strsplit(as.character(nodes[[aop_column]]), ",")
  
  nodes$has_complete_aop <- sapply(nodes$Aop_list, function(aops) {
    any(trimws(aops) %in% complete_aops)
  })
  
  nodes_sub <- nodes[nodes$has_complete_aop, ]
  
  #FILTER EDGES: both the nodes have to belog to the AOPs
  edges_sub <- edges[edges$from %in% nodes_sub$id & edges$to %in% nodes_sub$id, ]
  
  if (nrow(nodes_sub) == 0) {
    stop("No nodes found belonging to the complete AOPs.")
  }
  
  #USE THE BASE FUNCTION TO PLOT THE RESULTS
  plot_aop(nodes = nodes_sub,
           edges = edges_sub,
           numerical_variables = numerical_variables,
           time_column = time_column,
           chemicals_column = chemicals_column,
           aop_column = aop_column)
  
  # AOPfingerprintR::plot_visNetwork(nodes = nodes_sub,
  #                 edges = edges_sub,
  #                 group_by = "ssbd",
  #                 numerical_variables = numerical_variables)
}

#Function to show only the component connected to an AOP, colored by num of chemicals and grouping by chemical.
num_chem <- function(nodes,
                     edges,
                     numerical_variables,
                     time_column      = "time",
                     chemicals_column = "chemicals",
                     aop_column = "Aop") {
  
  # -- COLOR NODES BY CHEMICAL COUNT
  chem_counts <- sapply(
    strsplit(as.character(nodes[[chemicals_column]]), ","),
    function(x) length(unique(trimws(x[x != ""])))
  )
  
  # -- COUNT AOP PER NODE
  aop_counts <- sapply(
    strsplit(as.character(nodes[[aop_column]]), ","),
    function(x) length(unique(trimws(x[x != ""])))
  )
  
  nodes$color <- ifelse(
    chem_counts == 1,                      "#ffd6a5",
    ifelse(chem_counts <= 5,               "#f28b82",
           ifelse(chem_counts <= 10,            "#caffbf",
                  ifelse(chem_counts <= 15,          "#a5d8ff",
                         "#ffc6ff")
           )
    )
  )
  
  # -- SHAPE BY KE TYPE
  nd <- nodes$id
  ke_type <- character(length(nd))
  for (n in seq_along(nd)) {
    ke <- nd[n]
    ke_annotation <- aop_ke_table_hure$Ke_type[aop_ke_table_hure$Ke == ke]
    if (!"MolecularInitiatingEvent" %in% ke_annotation &
        !"AdverseOutcome"          %in% ke_annotation) {
      ke_type[n] <- "KeyEvent"
    } else {
      tab <- table(ke_annotation)
      if ("KeyEvent" %in% names(tab)) tab <- tab[names(tab) != "KeyEvent"]
      ke_type[n] <- names(sort(tab, decreasing = TRUE))[1]
    }
  }
  my_shapes <- c("MolecularInitiatingEvent" = "triangle",
                 "KeyEvent"                = "square",
                 "AdverseOutcome"          = "star")
  nodes$shape <- my_shapes[ke_type]
  
  # -- GROUP BY CHEMICALS
  nodes$group <- sapply(
    strsplit(as.character(nodes[[chemicals_column]]), ","),
    function(x) paste(trimws(x), collapse = ", ")
  )
  # -- ADD LABEL
  nodes$label <- nodes$Ke_description
  # -- TOOLTIP
  nodes$title <- paste0(
    "<p>Id: ", nodes$id, "</p>",
    "<p>Description: ", nodes$Ke_description, "</p>"
  )
  for (nv in numerical_variables) {
    nodes$title <- paste0(nodes$title, "<p>", nv, ": ", nodes[[nv]], "</p>")
  }
  nodes$title <- paste0(
    nodes$title,
    "<p># Chemicals: ", chem_counts, "</p>",
    "<p>Chemicals: ", nodes$group, "</p>",
    "<p># AOPs: ", aop_counts, "</p>",
    "<p>AOPs: ", nodes[[aop_column]], "</p>"
  )
  
  # -- LEGEND
  legend_nodes <- list(
    list(label = "1 chemical",       color = "#ffd6a5"),
    list(label = "2–5 chemicals",    color = "#f28b82"),
    list(label = "6–10 chemicals",   color = "#caffbf"),
    list(label = "11–15 chemicals",  color = "#a5d8ff"),
    list(label = ">15 chemicals",    color = "#ffc6ff"),
    list(label = "MolecularInitiatingEvent", shape = "triangle"),
    list(label = "KeyEvent",               shape = "square"),
    list(label = "AdverseOutcome",         shape = "star")
  )
  
  # -- PLOT
  visNetwork::visNetwork(nodes, edges) %>% 
    visNetwork::visLegend(position = "right",
                          addNodes  = legend_nodes,
                          useGroups = FALSE) %>% 
    visNetwork::visOptions(selectedBy = list(variable = "group",
                                             multiple = TRUE))
}


#Function to show the full network and color it by the AOP of interest, grouping by chemical. the color is grey if it is part of more than one aop
plot_colored_aop <- function(nodes,
                             edges,
                             numerical_variables,
                             chemicals_column = "chemicals",
                             aop_column = "Aop",
                             aop_list = NULL) {
  
  if (!requireNamespace("visNetwork", quietly = TRUE))
    stop("Package 'visNetwork' is required but not installed.")
  
  # -- COUNT CHEMICALS PER NODE
  chem_counts <- sapply(
    strsplit(as.character(nodes[[chemicals_column]]), ","),
    function(x) length(unique(trimws(x[x != ""])))
  )
  
  # -- COUNT AOP PER NODE
  nodes_aop_lists <- strsplit(as.character(nodes[[aop_column]]), ",")
  nodes_aop_lists <- lapply(nodes_aop_lists, function(x) trimws(x[x != ""]))
  
  aop_counts <- sapply(nodes_aop_lists, function(x) length(unique(trimws(x[x != ""]))))
  
  
  # -- COLOR BY AOP
  if (is.null(aop_list) || length(aop_list) == 0) {
    #if aop_list is empty, color by number of aop
    aop_counts <- sapply(nodes_aop_lists, length)
    nodes$color <- ifelse(
      aop_counts == 1, "#caffbf",
      ifelse(aop_counts == 2, "#a5d8ff",
             ifelse(aop_counts == 3, "#ffd6a5",
                    "#A9A9A9"))
    )
    mapping <- NULL
  } else {
    max_colors <- 6
    base_palette <- c(
      "#bdb2ff", # verde
      "#ffadad", # arancione
      "#a5d8ff", # blu
      "#ffc6ff", # rosa
      "#caffbf", # verde chiaro
      "#ffd6a5"  # giallo ocra
    )
    
    if (!is.null(aop_list) && length(aop_list) > max_colors) {
      stop("Too many AOP in the list, the maximum number supported is 6.")
    }
    
    mapping <- data.frame(
      AOP = aop_list,
      Color = base_palette[seq_along(aop_list)],
      stringsAsFactors = FALSE
    )
    
    nodes$color <- sapply(nodes_aop_lists, function(aops) {
      in_mapping <- intersect(aops, mapping$AOP)
      if (length(in_mapping) == 0) {
        "#A9A9A980" # doesn't belong to any aop in the list
      } else if (length(in_mapping) == 1) {
        mapping$Color[mapping$AOP == in_mapping]
      } else {
        "#1b9e77" # belongs to more than one aop in the list
      }
    })
  }
  
  # -- SHAPE BY KE TYPE
  nd <- nodes$id
  ke_type <- character(length(nd))
  for (n in seq_along(nd)) {
    ke <- nd[n]
    ke_annotation <- aop_ke_table_hure$Ke_type[aop_ke_table_hure$Ke == ke]
    if (!"MolecularInitiatingEvent" %in% ke_annotation &
        !"AdverseOutcome"          %in% ke_annotation) {
      ke_type[n] <- "KeyEvent"
    } else {
      tab <- table(ke_annotation)
      if ("KeyEvent" %in% names(tab)) tab <- tab[names(tab) != "KeyEvent"]
      ke_type[n] <- names(sort(tab, decreasing = TRUE))[1]
    }
  }
  my_shapes <- c("MolecularInitiatingEvent" = "triangle",
                 "KeyEvent"                = "square",
                 "AdverseOutcome"          = "star")
  nodes$shape <- my_shapes[ke_type]
  
  # -- GROUP BY CHEMICALS
  nodes$group <- sapply(
    strsplit(as.character(nodes[[chemicals_column]]), ","),
    function(x) paste(trimws(x), collapse = ", ")
  )
  
  # -- NODE SIZE
  #nodes$size <- chem_counts * 2
  
  # -- ADD LABEL
  has_aop <- vapply(
    nodes_aop_lists,
    function(aops) {
      aops <- aops[!is.na(aops) & nzchar(aops)]
      length(intersect(aops, mapping$AOP)) > 0
    },
    logical(1)
  )
  
  nodes$label <- ifelse(has_aop, nodes$Ke_description, "")
  
  #nodes$label <- nodes$Ke_description
  
  # -- TOOLTIP
  nodes$title <- paste0(
    "<p>Id: ", nodes$id, "</p>",
    "<p>Description: ", nodes$Ke_description, "</p>"
  )
  for (nv in numerical_variables) {
    nodes$title <- paste0(nodes$title, "<p>", nv, ": ", nodes[[nv]], "</p>")
  }
  nodes$title <- paste0(
    nodes$title,
    "<p># Chemicals: ", chem_counts, "</p>",
    "<p>Chemicals: ", nodes$group, "</p>",
    "<p># AOPs: ", aop_counts, "</p>",
    "<p>AOPs: ", nodes[[aop_column]], "</p>"
  )
  
  # -- LEGEND
  legend_nodes <- list(
    list(label = "Multiple AOPs (custom)", color = "#1b9e77"),
    list(label = "No AOP in list", color = "#A9A9A980"),
    list(label = "MolecularInitiatingEvent", shape = "triangle"),
    list(label = "KeyEvent",                shape = "square"),
    list(label = "AdverseOutcome",          shape = "star")
  )
  
  if (is.null(mapping)) {
    legend_nodes <- list(
      list(label = "1 AOP",  color = "#caffbf"),
      list(label = "2 AOP",  color = "#a5d8ff"),
      list(label = "3 AOP",  color = "#ffd6a5"),
      list(label = ">3 AOP", color = "#A9A9A9"),
      list(label = "MolecularInitiatingEvent", shape = "triangle"),
      list(label = "KeyEvent",                shape = "square"),
      list(label = "AdverseOutcome",          shape = "star")
    )
  } else {
    for (i in seq_len(nrow(mapping))) {
      legend_nodes[[length(legend_nodes) + 1]] <- list(
        label = mapping$AOP[i],
        color = mapping$Color[i]
      )
    }
  }
  
  colored_nodes <- nodes$id[nodes$color != "#A9A9A980"]
  edges_filtered <- edges[edges$from %in% colored_nodes | edges$to %in% colored_nodes, ]
  connected_nodes <- unique(c(edges_filtered$from, edges_filtered$to, colored_nodes))
  nodes_filtered <- nodes[nodes$id %in% connected_nodes, ]
  
  # -- PLOT
  visNetwork::visNetwork(nodes_filtered, edges_filtered) %>% 
    visNetwork::visLegend(position = "right",
                          addNodes  = legend_nodes,
                          useGroups = FALSE) %>% 
    visNetwork::visOptions(selectedBy = list(variable = "group",
                                             multiple = TRUE))
}

#Function to color by categories
categories <- function(nodes,
                       edges,
                       numerical_variables,
                       chemicals_column = "chemicals",
                       aop_column = "Aop",
                       aop_list = NULL,
                       show_grey) {
  
  if (!requireNamespace("visNetwork", quietly = TRUE))
    stop("Package 'visNetwork' is requi#f28b82 but not installed.")
  
  # -- COUNT AOP PER NODE 
  nodes_aop_lists <- strsplit(as.character(nodes[[aop_column]]), ",")
  nodes_aop_lists <- lapply(nodes_aop_lists, function(x) trimws(x[x != ""]))
  
  # -- COLOR BY SSbD CATEGORY
  nodes$color <- rep("#d9d9d9", nrow(nodes))  # default grey
  
  category_aop_map <- list(
    Endometrium = c("Aop:167", "Aop:112"),
    Reproductive_glands      = c("Aop:165", "Aop:440", "Aop:170", "Aop:168"),
    Liver      = c("Aop:1", "Aop:118", "Aop:37", "Aop:108", "Aop:107", "Aop:117", "Aop:46"),
    Kidney_UrinaryTract      = c("Aop:116", "Aop:105", "Aop:121"),
    Endocrine_Glands      = c("Aop:162", "Aop:110", "Aop:119", "Aop:169"),
    Lung_SquamousEpithelium     = c("Aop:114", "Aop:115", "Aop:109")
  )
  
  palette_colors <- c(
    Endometrium = "#ffadad",  
    Reproductive_glands      = "#ffd6a5",  
    Liver     = "#bdb2ff",  
    Kidney_UrinaryTract       = "#caffbf",  
    Endocrine_Glands    = "#9bf6ff",  
    Lung_SquamousEpithelium   = "#1b9e77"   
  )

  # Ciclo sui nodi e assegna colori in base alla categoria
  for (i in seq_len(nrow(nodes))) {
    aops <- nodes_aop_lists[[i]]
    matched_category <- NULL
    
    # Controlla a quale categoria appartiene il nodo
    for (cat in names(category_aop_map)) {
      if (any(aops %in% category_aop_map[[cat]])) {
        matched_category <- cat
        break
      }
    }
    if (!is.null(matched_category)) {
      nodes$color[i] <- palette_colors[matched_category]
    }
  }
  
  # -- COUNT CHEMICALS / AOP 
  chem_counts <- sapply(
    strsplit(as.character(nodes[[chemicals_column]]), ","),
    function(x) length(unique(trimws(x[x != ""])))
  )
  
  aop_counts <- sapply(strsplit(as.character(nodes[[aop_column]]), ","),
                       function(x) length(unique(trimws(x[x != ""]))))
  
  # -- SHAPE BY KE TYPE
  nd <- nodes$id
  ke_type <- character(length(nd))
  for (n in seq_along(nd)) {
    ke <- nd[n]
    ke_annotation <- aop_ke_table_hure$Ke_type[aop_ke_table_hure$Ke == ke]
    if (!"MolecularInitiatingEvent" %in% ke_annotation &
        !"AdverseOutcome"          %in% ke_annotation) {
      ke_type[n] <- "KeyEvent"
    } else {
      tab <- table(ke_annotation)
      if ("KeyEvent" %in% names(tab)) tab <- tab[names(tab) != "KeyEvent"]
      ke_type[n] <- names(sort(tab, decreasing = TRUE))[1]
    }
  }
  my_shapes <- c("MolecularInitiatingEvent" = "triangle",
                 "KeyEvent"                = "square",
                 "AdverseOutcome"          = "star")
  nodes$shape <- my_shapes[ke_type]
  
  # -- GROUP BY CHEMICALS
  nodes$group <- sapply(
    strsplit(as.character(nodes[[chemicals_column]]), ","),
    function(x) paste(trimws(x), collapse = ", ")
  )
  
  # -- LABEL
  nodes$label <- ifelse((nodes$color != "#d9d9d9"), nodes$Ke_description, "")
  #nodes$label <- nodes$Ke_description
  
  # -- TOOLTIP
  nodes$title <- paste0(
    "<p>Id: ", nodes$id, "</p>",
    "<p>Description: ", nodes$Ke_description, "</p>"
  )
  for (nv in numerical_variables) {
    nodes$title <- paste0(nodes$title, "<p>", nv, ": ", nodes[[nv]], "</p>")
  }
  nodes$title <- paste0(
    nodes$title,
    "<p># Chemicals: ", chem_counts, "</p>",
    "<p>Chemicals: ", nodes$group, "</p>",
    "<p># AOPs: ", aop_counts, "</p>",
    "<p>AOPs: ", nodes[[aop_column]], "</p>"
  )
  
  # -- LEGEND
  legend_nodes <- list(
    list(label = "Endometrium",          color = "#ffadad"), 
    list(label = "Reproductive glands",     color = "#ffd6a5"),
    list(label = "Liver", color = "#bdb2ff"),
    list(label = "Kidney and Urinary Tract",      color = "#caffbf"),
    list(label = "Endocrine Glands",      color = "#9bf6ff"),
    list(label = "Lung and Squamous Epithelium",      color = "#1b9e77")
  )
  legend_nodes <- c(legend_nodes,
                    list(list(label = "MolecularInitiatingEvent", shape = "triangle"),
                         list(label = "KeyEvent",                 shape = "square"),
                         list(label = "AdverseOutcome",           shape = "star")))
  
  if (!show_grey){
    colored_nodes <- nodes$id[nodes$color != "#d9d9d9"]
    edges_filtered <- edges[edges$from %in% colored_nodes | edges$to %in% colored_nodes, ]
    connected_nodes <- unique(c(edges_filtered$from, edges_filtered$to, colored_nodes))
    nodes_filtered <- nodes[nodes$id %in% connected_nodes, ]
  } else {
    edges_filtered <- edges
    nodes_filtered <- nodes
  }
  
  # -- PLOT
  visNetwork::visNetwork(nodes_filtered, edges_filtered) %>% 
    visNetwork::visLegend(position = "right",
                          addNodes  = legend_nodes,
                          useGroups = FALSE) %>% 
    visNetwork::visOptions(selectedBy = list(variable = "group",
                                             multiple = TRUE))
}

find_nodes <- function(nodes, edges) {
  #Given in input a list of nodes and a dataframe of edges it returns the list of the nodes that are part of the connected components
  
  #find the edges that have a node in nodes
  relevant_edges <- edges[edges$from %in% nodes | edges$to %in% nodes, ]
  
  #take all the nodes that are involved in the found edges
  new_nodes <- unique(c(relevant_edges$from, relevant_edges$to))
  
  #if the list is not changed end the recursion
  if (setequal(new_nodes, nodes)) {
    return(new_nodes)
  } else {
    #if the list are not equal means that there are others nodes that are connected so recall the function with the new list of nodes
    find_nodes(new_nodes, edges)
  }
}

##### General Network grouped by AOP #####
# AIM: make the network with all the enriched KE, grouped by AOP and colored by the num of aop

grouped_aop = plot_grouped_aop(nodes = nodes_edges$nodes,
                               edges = nodes_edges$edges,
                               numerical_variables = numerical_variables,
                               chemicals_column = "chemical")

grouped_aop
#visSave(grouped_aop, file = paste0(dir_path, "net_grouped_aop_all.html"))

##### General Network with BMD dimension #####
# AIM: make the network with all the enriched KE, grouped by chemical, colored by the num of chemical and dimension of nodes proportional to BMD_norm

bmd_net = try_bmd(nodes = nodes_edges$nodes,
             edges = nodes_edges$edges,
             numerical_variables = numerical_variables,
             chemicals_column = "chemical")

bmd_net
#visSave(bmd_net, file = paste0(dir_path, "bmd_norm.html"))

##### Complete AOP #####
## AIM: make a network in which are present only the enriched KE that are part of AOPs for which all the KE are enriched.

complete_aop = plot_complete_aop(nodes = nodes_edges$nodes,
                                 edges = nodes_edges$edges,
                                 numerical_variables = numerical_variables,
                                 time_column = "time",
                                 chemicals_column = "chemical")

complete_aop
#visSave(complete_aop, file = paste0(dir_path,"net_complete_aop.html"))


##### Connected to an AO #####
#only components connected to an AO

ao_only <- as.data.frame(data) %>% filter(Ke_type == "AdverseOutcome")
ao_list <- unique(ao_only$TermID)

new_nodes_list <- find_nodes(ao_list, nodes_edges$edges)

new_nodes <- nodes_edges$nodes %>% filter(id %in% new_nodes_list)
new_nodes2 <- new_nodes %>% 
  select(id) %>% 
  left_join(df_grouped, by = c(id = "TermID"))

new_edges <- nodes_edges$edges %>% filter(from %in% new_nodes_list | to %in% new_nodes_list)

new_nodes_edges <- list(nodes = new_nodes2, edges = new_edges)


vn = num_chem(nodes = new_nodes_edges$nodes,
              edges = new_nodes_edges$edges,
              numerical_variables = numerical_variables,
              time_column = "time",
              chemicals_column = "chemical")

vn
#visSave(vn, file = paste0(dir_path,"aop_of_ao.html"))


##### Personalized colored AOP #####
## AIM: make the network with all the enriched KE and color only the events that are part of the aop of interest 
aop_lists <- list(
  gse17624_1 = c("Aop:57", "Aop:281", "Aop:437", "Aop:276", "Aop:420", "Aop:416"),
  gse17624_2 = c("Aop:166", "Aop:19", "Aop:306", "Aop:288"),
  fingerprint_1 = c("Aop:37", "Aop:41", "Aop:46", "Aop:439", "Aop:393"),
  fingerprint_2 = c("Aop416", "Aop:451", "Aop:372", "Aop:106", "Aop:322"),
  general_moa = c("Aop:296", "Aop:423", "Aop:293", "Aop:294", "Aop:205"),
  lung = c("Aop:109", "Aop:272", "Aop:303", "Aop:417", "Aop:4168", "Aop:451"),
  others = c("Aop:115", "Aop:114", "Aop:316", "Aop:202", "Aop:171", "Aop:409")
)

for (name in names(aop_lists)){
  aopList <- aop_lists[[name]]
  colored_aop = plot_colored_aop(nodes = nodes_edges$nodes,
                               edges = nodes_edges$edges,
                               numerical_variables = numerical_variables,
                               chemicals_column = "chemical",
                               aop_list = aopList)


  visSave(colored_aop, file = paste0(dir_path, "specific/", name,".html"))
}

p = categories(nodes = nodes_edges$nodes,
               edges = nodes_edges$edges,
               numerical_variables = numerical_variables,
               chemicals_column = "chemical",
               show_grey = TRUE)

visSave(p, file = paste0(dir_path,"specific/Adenomas_Carcinomas.html"))

p_no = categories(nodes = nodes_edges$nodes,
               edges = nodes_edges$edges,
               numerical_variables = numerical_variables,
               chemicals_column = "chemical",
               show_grey = FALSE)

visSave(p_no, file = paste0(dir_path,"specific/Adenomas_Carcinomas_nogrey.html"))
