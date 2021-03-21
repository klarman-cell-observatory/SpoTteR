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
#coors_small = read.delim("/Users/svickovi/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SpoTter/RA/processing/2018Diamond_CNO-254_RA.E2-Spot000001_inferred_points.tsv")
# get new coordiantes that match: Note: 2%=downsamling small images for spotter and 85%=downsapling for GUI
tissue_small = coors_small

#tissue_small$x = tissue_small$x*(100/2)*(85/100)
#tissue_small$y = tissue_small$y*(100/2)*(85/100)

# check if envoking imaging at 10x or 20x
if (obj == "spotter_10X_hex.sh"){
    tissue_small$x_prec = tissue_small$x_prec*(100/2)*(85/100) # 85 to match app final file size and 2 match initial file size used for detection of spots
    tissue_small$y_prec = tissue_small$y_prec*(100/2)*(85/100)
    dm = 300
    fact = 50
} else {
    tissue_small$x_prec = tissue_small$x_prec*(100/1)*(42.5/100)
    tissue_small$y_prec = tissue_small$y_prec*(100/1)*(42.5/100)
    dm = 600
    fact = 100
}

# get x,y, w and h and set up some controls
x1min=min(tissue_small$x_id)
y1min=min(tissue_small$y_id)
x2max=max(tissue_small$x_id)
y2max=max(tissue_small$y_id)

print(x1min)
print(y1min)
print(x2max)
print(y2max)

if ((x1min == "1") & (y1min == "1") & (x2max == "66") & (y2max == "63")){
  

  x1 = tissue_small[tissue_small$bc == "1x62",]$x_prec
  y1 = tissue_small[tissue_small$bc == "2x63",]$y_prec
  print(x1)
  print(y1)
  
  h_x2 = tissue_small[tissue_small$bc == "1x2",]$x_prec
  h_y2 = tissue_small[tissue_small$bc == "2x1",]$y_prec
  print(h_x2)
  print(h_y2)
  
  w_x2 = tissue_small[tissue_small$bc == "66x63",]$x_prec # input as x
  w_y2 = tissue_small[tissue_small$bc == "66x63",]$y_prec # input as y
  print(w_x2)
  print(w_y2)
  
  hg = sqrt((h_x2-x1)^2+(h_y2-y1)^2) # input as width
  wd = sqrt((w_x2-x1)^2+(w_y2-y1)^2) # input as height
  print(hg)
  print(wd)
  
  # write cropping file for imagemagick
  write.table(t(c(round(wd), round(hg), round(w_x2), round(w_y2))), sep =",", col.names = F, row.names = F, file = paste0(folder,"/", array, "_", sample, "_", batch, ".", well, "-Spot000001_cropping_points.tsv", sep=""), quote = F)

  # print out coordinates for celery GUI
  library(data.table)
  coors_small = data.table(coors_small)
  coors_small[, final_y:= if (x_id %% 2==0) y_id=seq(1,63,by=2), by = x_id]
  coors_small[, final_y:= if (x_id %% 2!=0) y_id=seq(2,62,by=2), by = x_id]
  coors_small[,final_x:=(max(x_id)+1)-x_id,]
  coors_small[,bc:=paste0(final_x,"x",final_y),]
  
  # output spots 
  ut = coors_small[coors_small$feature_spot == TRUE,]$bc
  ut = str_replace(string = ut, pattern = "x", replacement = "_")
  dir.create("../output/spots/", recursive = T)
  write.table(t(ut), file=paste("../output/spots/", array, sample, "_", well,   "_stdata_under_tissue_spots.txt",sep=""), sep = "\t", quote=F, col.names=F, row.names=F)
  
  # save transformation matrix from original image to spots 
  trans1s = (coors_small[coors_small$x_id == 1 & coors_small$new_y == 62,]$x_precise*fact-coors_small[coors_small$x_id == 1 & coors_small$new_y == 2,]$x_precise*fact)/62
  trans5s = (coors_small[coors_small$x_id == 2 & coors_small$new_y == 63,]$y_precise*fact-coors_small[coors_small$x_id == 2 & coors_small$new_y == 1,]$y_precise*fact)/65
  
  trans7s = coors_small[coors_small$x_id == 1 & coors_small$new_y == 2,]$x_precise*fact - trans1s*1
  trans8s = coors_small[coors_small$x_id == 2 & coors_small$new_y == 1,]$y_precise*fact - trans5s*1
  
  dir.create("../output/transformation_matrix/", recursive = T)
  write.table(t(c(trans1s, 0, 0, 0, trans5s, 0, trans7s, trans8s, 1)), file = paste0("../output/transformation_matrix/", array, sample, "_", well,"_transformation_matrix.tsv"),sep="\t",row.names=FALSE,quote=FALSE, col.names = FALSE)
  
  # prepare data for json file format to read into spaceranger
  #coors_small = read.delim("/Users/svickovi/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SpoTter/STB46/processing/10015_CN108_STB46.C1-Spot000001_inferred_points.tsv")
  
  flag = c("2_2", "2_4", "3_3", "3_5", "4_2", "4_4", "6_2", "6_4", "5_3", "5_5", "8_2", "8_4",
           "7_3", "7_5", "10_2", "10_4", "9_3", "9_5", "62_62", "64_62")
  frame = c(paste(1, seq(1,63, by = 2), sep ="_"), paste(66, seq(2,63, by = 2), sep ="_"),
            paste(seq(1,66, by = 2), 1, sep ="_"), paste(seq(1,63, by = 2),63, sep ="_"))
  
  ftoligo = rbind(coors_small[paste(coors_small$x_id, coors_small$y_id, sep ="_") %in% flag,], coors_small[paste(coors_small$x_id, coors_small$y_id, sep ="_") %in% frame,])
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

