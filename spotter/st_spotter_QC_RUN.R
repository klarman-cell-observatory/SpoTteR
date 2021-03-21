# Load libraries 
suppressMessages(require(stringr))
suppressMessages(require(imager))
suppressMessages(require(vioplot))
suppressMessages(require(rmarkdown))
suppressMessages(require(wesanderson))

# Expect command line args at the end 
args = commandArgs(trailingOnly = TRUE)

# Source the functions file (make sure path is valid)
script_path=args[4]
source(paste0(script_path,"/", "st_spotter_QC_FUNCTIONS.R"))

# Read in expression matrix from ST pipeline v1.3.1
all_barcodes <- read.delim(args[1], sep="\t",quote="", row.names = 1)

# Exchange gene names in expression tsv file
#all_barcodes = convert.gene.ids(all_barcodes, "mouse")
#print(colnames(all_barcodes))

# Read in batch number 
batch = args[2]

# Read in folder with aligned tsv files
folder = args[3]

# Set names
sample = sapply(strsplit(basename(args[1]),split="_"),"[[",2)
well = sapply(strsplit(basename(args[1]),split="_"),"[[",3)
array = sapply(strsplit(basename(args[1]),split="_"),"[[",1)

# Send to stout
print(paste0("Processing sample: ", paste(folder, "/", array, "_", sample, "_", batch, ".", well, "-Spot000001_inferred_points.tsv", sep="")))

# Read in tsv file as output from SpoTter
spots = read.delim(paste0("../output/spots/",  array, sample, "_", well, "_stdata_under_tissue_spots.txt"), header = T)
img = load.image(paste0("../processing/", array, sample, "_", well, "_HE.jpg"))

# Test and output permeabilization
dir.create("../output/QC_plots", recursive = T)
print("checnk bellow")
setwd("../output/QC_plots")
test.plot(all_barcodes, spots, img)

# Get expression matrix for selected spots (as defined by SpoTter)
exp.values = exp.under.tissue(all_barcodes, spots) # matrix under tissue
exp.not.values = exp.not.under.tissue(all_barcodes, spots) # matrix outside tissue

# Make QC plots of top genes
top.genes.plot(all_barcodes, spots)

# Make QC plots of interesting genes
interesting.genes.plot(all_barcodes, spots, genes=c("Snap25", "Penk", "Ppp3ca", "Kctd12", "Fabp7", "Mbp", "Epha6", "Shisa9")) # MOB
# interesting.genes.plot(all_barcodes, spots, genes=c("Slc17a6", "Ier5", "Nes", "Kctd4", "Pcdh20", "Dlk1", "Nnat", "Hap1")) # Cortex
# interesting.genes.plot(all_barcodes, spots, genes=c("Cd3e", "Adgre1", "Cd68", "Cd163", "Cd19", "Cd4", "Cd8a", "Siglec1")) #Proteins
# interesting.genes.plot(all_barcodes, spots, genes=c("Cd8a", "Cd4", "Itgam", "Ptprc", "Pdcd1", "Tcf7", "Adgre1", "Cd3e")) # MC38

# Write file with correct coordinates as a matrix (tab-delimited)
dir.create("../aligned_counts", recursive = T)
setwd("../aligned_counts")
write.coor.matrix(exp.values)
write.not.coor.matrix(exp.not.values)

# Write spots files for celery 
## this refines to tae only spots that have received a value on sequencing
dir.create("../spots", recursive = T)
setwd("../spots")
write.coor.spots(exp.values)
### Extra QC stuff needed for the Metadata collection
inside.stdata.file = exp.values
dir.create("../QC_logs", recursive = T)
setwd("../QC_logs")
file.create(paste(array,sample,well,"QC_data.txt",sep="_"))
# write("Expression levels for the top five genes under the tissue:\t", paste(array,sample,well,"QC_data.txt",sep="_"))
# write(names(head(sort(rowSums(inside.stdata.file), decreasing = T), n=5)), paste(array,sample,well,"QC_data.txt",sep="_"), ncolumns=5, append=TRUE, sep="\t")
# write(head(sort(rowSums(inside.stdata.file), decreasing = T), n=5),paste(array,sample,well,"QC_data.txt",sep="_"), ncolumns=5, append=TRUE, sep="\t")
write(paste("Mean number of transcripts per feature under tissue:", round(mean(colSums(inside.stdata.file)), digits = 2), sep = "\t"),paste(array,sample,well,"QC_data.txt",sep="_"),append=TRUE)
write(paste("Mean number of genes per feature under tissue:", round(mean(colSums(inside.stdata.file>0)), digits = 2), sep = "\t"), paste(array,sample,well,"QC_data.txt",sep="_"), append=TRUE)
write(paste("Number of features covered:", dim(inside.stdata.file)[2], sep = "\t"), paste(array,sample,well,"QC_data.txt",sep="_"), append = TRUE)

# Be sure to remove files
rm(list = ls())
