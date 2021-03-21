### Runner script for automatic ST spot detection from HE images ### 

# Expect command line args at the end
args = commandArgs(trailingOnly = TRUE)

### Important: size of input images needs to be around 500x500 pixels ###
image_file=args[1] # small image
script_path=args[2]
output=args[3] # st batch number

# Load libraries
suppressMessages(library(imager))
suppressMessages(library(ggplot2));theme_set(theme_bw())
suppressMessages(library(data.table))
ggpl_settings=list(geom_raster(aes(fill=value)),
scale_y_continuous(trans=scales::reverse_trans()),
scale_fill_gradient(low="black",high="white"),
guides(fill=FALSE),
coord_fixed())

# Source f(x) files
source(paste0(script_path,"/", "spotter_functions_hex.R"))

# change wd to export tsv with aligned spots
setwd(paste0(output))

# Set you presumed array dimensions (incuding the frame if any)
exp_cols = 66
exp_rows = 63

#-------PROCESS one input image---------#
process_HE(image_file)





