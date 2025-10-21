# Take the arguments from the command line
args <- commandArgs(trailingOnly = TRUE)  

if (length(args) == 0) {
  stop("Missing argument!")
}

GSE <- args[1]  


# Set the path
directory_path <- paste0("/data/omics/public_data/endocrine_disruptors/", GSE, "/processed_data/unique_sorted/")

if (!dir.exists(directory_path)) {
  stop("The path does not exist!")
}

# List all .tsv files in the directory
bam_files_list <- list.files(path = directory_path, pattern = "\\.bam$", full.names = TRUE)

gtf_file <- "/data/omics/references/hsapiens/gtf/Homo_sapiens.GRCh38.110.gtf"


# Raw_files directory
raw_data_dir <- paste0("/data/omics/public_data/endocrine_disruptors/", GSE, "/raw_data_GSM_merged/")

# If there are paired end files (*_1.fastq.gz), set paired.end = TRUE, otherwise set it at FALSE
if (length(list.files(raw_data_dir, pattern = "_1.fastq.gz", full.names = TRUE)) > 0) {
  paired <- TRUE
} else {
  paired <- FALSE
}


# function

get_raw_counts <- function(bam.files, annotation, paired.end, ncores=20) {

  if(is.null(bam.files)){stop("Error: please provide one or more .bam files!")}
  if(is.null(annotation)){stop("Error: please provide an annotation file in .gtf or .gff format!")}
  if(!class(paired.end)=="logical"){stop("Error: indicate whether the reads are paired end")}
  if(!class(ncores)=="numeric"){stop("Error: please input a numeric value!")}

  raw.counts1 <- Rsubread::featureCounts(bam.files, annot.ext=annotation, isPairedEnd=paired.end, isGTFAnnotationFile = TRUE, nthreads = ncores)

  total.counts.matrix <- raw.counts1[[1]]
  genelengths <- raw.counts1$annotation$Length
  names(genelengths) <- raw.counts1$annotation$GeneID
  return(list(total.counts.matrix, genelengths))
}

# do the job

output <- capture.output({
  counts_list <- get_raw_counts(bam.files = bam_files_list, annotation = gtf_file, paired.end = paired)
})

# Transform the list into a matrix
count_matrix <- unlist(counts_list[[1]]) # extracts the counts which are in list element 1


# remove .unique.sorted.bam from colnames
colnames(count_matrix) <- sub("_unique_sorted.bam", "" ,colnames(count_matrix))


# save results

write.table(count_matrix, paste0("/data/omics/public_data/endocrine_disruptors/",GSE,"/processed_data/counts/count_matrix.tsv"),
            sep = "\t", col.names = T, row.names = T, quote = F)
save(count_matrix, file = paste0("/data/omics/public_data/endocrine_disruptors/",GSE,"/processed_data/counts/count_matrix.RData"))

# Write the captured output to a file
writeLines(output, paste0("/data/omics/public_data/endocrine_disruptors/",GSE,"/processed_data/counts/output_text_featureCounts.txt"))

