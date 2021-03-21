# load libraries
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(jsonlite))

# Expect command line args at the end 
args = commandArgs(trailingOnly = TRUE)

# Read in batch number 
batch = args[2]

# Read in folder with aligned tsv files
folder = args[3]

# Read in imaging objective
obj = basename(args[4])
print(paste0("arg 4 detecting objective: ", obj))

# Set names
sample = sapply(strsplit(args[1],split="_"),"[[",2)
well = sapply(strsplit(args[1],split="_"),"[[",3)
array = sapply(strsplit(args[1],split="_"),"[[",1)

# Send to stout
print(paste0("Processing sample: ", paste0(array, "_", sample, "_", well)))

# Read in tsv file as output from SpoTter
coors_small = read.delim(paste(folder, "/", array, "_", sample, "_", batch, ".", well, "-Spot000001_inferred_points.tsv", sep=""))

# get new coordiantes that match: Note: 2%=downsamling small images for spotter and 85%=downsapling for GUI
tissue_small = coors_small

# check if envoking imaging at 10x or 20x
dm = 0
if (obj == "spotter_10X.sh"){
    tissue_small$x = tissue_small$x*(100/2)*(85/100) # 85 to match app final file size and 2 match initial file size used for detection of spots
    tissue_small$y = tissue_small$y*(100/2)*(85/100)
    dm = 300
    fact = 50
} else {
    tissue_small$x = tissue_small$x*(100/1)*(42.5/100)
    tissue_small$y = tissue_small$y*(100/1)*(42.5/100)
    dm = 600
    fact = 100
}

# get x,y, w and h and set up some controls
x1min=min(tissue_small$x_id)
y1min=min(tissue_small$y_id)
x2max=max(tissue_small$x_id)
y2max=max(tissue_small$y_id)

if ((x1min == "1") & (y1min == "1") & (x2max == "33") & (y2max == "35")){

  x1 = tissue_small[tissue_small$bc == "1x35",]$x
  y1 = tissue_small[tissue_small$bc == "1x35",]$y
  h_x2 = tissue_small[tissue_small$bc == "1x1",]$x
  h_y2 = tissue_small[tissue_small$bc == "1x1",]$y
  w_x2 = tissue_small[tissue_small$bc == "33x35",]$x # input as x
  w_y2 = tissue_small[tissue_small$bc == "33x35",]$y # input as y
  hg = sqrt((h_x2-x1)^2+(h_y2-y1)^2) # input as width
  wd = sqrt((w_x2-x1)^2+(w_y2-y1)^2) # input as height
  
  # write cropping file for imagemagick
  write.table(t(c(round(wd), round(hg), round(w_x2), round(w_y2))), sep =",", col.names = F, row.names = F, file = paste0(folder,"/", array, "_", sample, "_", batch, ".", well, "-Spot000001_cropping_points.tsv", sep=""), quote = F)

  # print out coordinates for celery GUI
  ut = coors_small[coors_small$feature_spot == TRUE,]$bc
  ut = str_replace(string = ut, pattern = "x", replacement = "_")
  dir.create("../output/spots/", recursive = T)
  write.table(t(ut), file=paste("../output/spots/", array, sample, "_", well,   "_stdata_under_tissue_spots.txt",sep=""), sep = "\t", quote=F, col.names=F, row.names=F)
  
  # save transformation matrix from original image to spots 
  trans1s = (coors_small[coors_small$x_id == 33 & coors_small$y_id == 1,]$x_precise*fact-coors_small[coors_small$x_id == 1 & coors_small$y_id == 1,]$x_precise*fact)/32
  trans5s = (coors_small[coors_small$x_id == 1 & coors_small$y_id == 35,]$y_precise*fact-coors_small[coors_small$x_id == 1 & coors_small$y_id == 1,]$y_precise*fact)/34
  
  trans7s = coors_small[coors_small$x_id == 1 & coors_small$y_id == 1,]$x_precise*fact - trans1s*1
  trans8s = coors_small[coors_small$x_id == 1 & coors_small$y_id == 1,]$y_precise*fact - trans5s*1
  
  dir.create("../output/transformation_matrix/", recursive = T)
  write.table(t(c(trans1s, 0, 0, 0, trans5s, 0, trans7s, trans8s, 1)), file = paste0("../output/transformation_matrix/", array, sample, "_", well,"_transformation_matrix.tsv"),sep="\t",row.names=FALSE,quote=FALSE, col.names = FALSE)
  
  # prepare data for json file format to read into spaceranger
  #coors_small = read.delim("/Users/svickovi/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SpoTter/STB46/processing/10015_CN108_STB46.C1-Spot000001_inferred_points.tsv")
  
  flag =c("2_2", "2_3", "2_4", "2_5",
                  "3_2", "3_3", "3_4", "3_5",
                  "4_2", "4_3", "4_4", "4_5",
                  "5_2", "5_3", "5_4", "5_5")
  ftoligo = rbind(coors_small[paste(coors_small$x_id, coors_small$y_id, sep ="_") %in% flag,], coors_small[coors_small$x_id == 1, ], coors_small[coors_small$x_id == 33, ],coors_small[coors_small$y_id == 1, ], coors_small[coors_small$y_id == 35, ] )
  fiducial = data.table(ftoligo$x, ftoligo$y, ftoligo$y_id, ftoligo$x_id, as.numeric(as.character(dm)), "frame", ftoligo$x_precise*fact ,ftoligo$y_precise*fact)
  colnames(fiducial) = c("x", "y", "row", "col", "dia" ,"fidName", "imageX", "imageY")
  
  coors_small = coors_small[!paste(coors_small$y_id, coors_small$x_id, sep ="_") %in% paste(fiducial$row, fiducial$col, sep ="_"),]
  oligos = data.table(coors_small$x, coors_small$y, coors_small$y_id, coors_small$x_id, as.numeric(as.character(dm)), coors_small$x_precise*fact, coors_small$y_precise*fact, coors_small$feature_spot)
  colnames(oligos) = c("x", "y", "row", "col", "dia" ,"imageX", "imageY", "tissue")

  transform = data.table(c(trans1s, 0, trans7s) ,c(0, trans5s, trans8s), c(0, 0, 1))
  colnames(transform) = NULL
  
  serialNumber = paste0(array, sample)
  area = paste0(well)
  
  # transform to start at zero
  oligos$row = oligos$row-1
  oligos$col = oligos$col-1
  fiducial$row = fiducial$row-1
  fiducial$col = fiducial$col-1
  
  #make json
  ls = c(list(oligos, fiducial, transform, serialNumber,area))
  names(ls) = c("oligo", "fiducial", "transform", "serialNumber", "area")

  json = toJSON((ls))
  dir.create("../output/json/", recursive = T)
  cat(json, file = "tmp.json")
  system(paste0(paste0("cat tmp.json | sed 's/,\"tissue\":false//g' | sed 's/\"serialNumber\":\\[/\"serialNumber\":/g' | sed 's/\\],\"area\":\\[/,\"area\":/'g  | sed 's/\\]}/}/g' > ", paste0("../output/json/", serialNumber,"_", area, ".json"))))
  file.remove("tmp.json")


} else {print("Array not aligned correctly.")}

