# Functions for the st_spotter
convert.gene.ids = function(my.data){

  # Load genenames from input stdata.tsv matrix
  # my.data = all_barcodes
  id = colnames(my.data)
  coors = row.names(my.data) # save coordinate names
  
  # check which species ID and check if ensembl id
  organism = ""
  if (grep("ENSMUSG", id[1], value = F) == 1) 
    {organism = "mouse"} else {organism = "human"}

  # Load organism
  if (organism == "mouse") {df=read.delim(paste0(script_path, '/data/Gene_names_mm.txt'), row.names = 1, header = T)}
  if (organism == "human")  {df=read.delim(paste0(script_path, '/data/Gene_names_hg.txt'), row.names = 1, header = T)}
  
  # Load organism
  #if (organism == "mouse") {df=read.delim('/Users/svickovi/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SpoTter/STaut/spotter/data/Gene_names_mm.txt', row.names = 1, header = T)}
  #if (organism == "human")  {df=read.delim('/Users/svickovi/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SpoTter/STaut/spotter/data/Gene_names_hg19.txt', row.names = 1, header = T)}
  
  #if (organism == "mouse") {df=read.delim('/Users/sanjavickovic/Desktop/SpoTter/STaut/spotter/data/Gene_names_mm.txt', row.names = 1, header = T)}
  #if (organism == "human")  {df=read.delim('/Users/sanjavickovic/Desktop/SpoTter/STaut/spotter/data/Gene_names_hg19.txt', row.names = 1, header = T)}
  
  
  # clean up dots in gene names 
  df[,1] = sapply(strsplit(row.names(df), split="_"), "[[",2)
  row.names(df) = sapply(strsplit(row.names(df), split="_"), "[[",1)
  row.names(df) = gsub("\\..*","", row.names(df))
  colnames(df) = "nm_refseq"
  
  # subset to name same genes as mm10
  colnames(my.data) = gsub("\\..*","", colnames(my.data))
  my.data = my.data[,colnames(my.data) %in% row.names(df)]
  id = colnames(my.data)
  coors = row.names(my.data) # save coordinate names

  # Check for name pairs
  ssg= ""
  for (genen in id){
      if (!is.na(df[genen,])){
          ssg = c(ssg, as.character(df[genen,]))
      }
      else {ssg = c(ssg, genen)}
  }
  ssg = ssg[-1]
  
  # Rename things in original matrix
  colnames(my.data) = ssg
  row.names(my.data) = coors
  
  # Save modified matrices with correct gene names
  return(my.data)
}

mean.region.data <- function(region.table){
  barcode.sums <- apply(region.table, 1, sum) 
  total.sum <- median(barcode.sums)
  norm.region.data <- (region.table/barcode.sums)*total.sum
  return(norm.region.data)
}

tpm.region.data <- function(region.table){
  barcode.sums <- apply(region.table, 1, sum) 
  total.sum <- 1000000
  norm.region.data <- (region.table/barcode.sums)*total.sum
  return(norm.region.data)
}

test.plot = function(all_barcodes, spots, img){
  
    # all_barcodes = read.delim("/Users/svickovi/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SpoTter/STB49/input/10015_CN107_D2_stdata.tsv", sep="\t",quote="", row.names = 1)
    
    #all_barcodes = read.delim("/Users/sanjavickovic/Desktop/SpoTter/STB49/input/10015_CN98_C1_stdata.tsv", sep="\t",quote="", row.names = 1)
  
    A = all_barcodes
    A = t(A)
    A = mean.region.data(A)
    B = as.data.frame(A[rowSums(A) != 0,])

    xa1 = colnames(B)
    xa = as.numeric(sapply(strsplit(xa1,split="x"),"[[",1))
    ya = as.numeric(sapply(strsplit(xa1,split="x"),"[[",2))
    za = colSums(B)
    za = log2(za+1)
    b.max = max(za)
  
    # spots = read.delim("/Users/svickovi/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SpoTter/STB49/output/spots/10015CN107_D2_stdata_under_tissue_spots.txt", header = T)
    #spots  = read.delim("/Users/sanjavickovic/Desktop/SpoTter/STB49/output/spots/10015CN98_C1_stdata_under_tissue_spots.txt", header = T)
    #library(stringr)
    spots = str_remove(colnames(spots), "X")
    xs = sapply(strsplit(spots, split = "_"),"[[",1)
    ys = sapply(strsplit(spots, split = "_"),"[[",2)
    
    #library(wesanderson)
    breaks = seq(0, max(za), by = (max(za)/254))
    cl.scale.all = wes_palette("Zissou1", 256, type = "continuous")
    dat.col = cl.scale.all[as.numeric(cut(za, breaks = breaks))]
    dat.col[is.na(dat.col)] <- cl.scale.all[1]  
    
    # img=load.image("/Users/svickovi/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SpoTter/STB30/output/rotated_images/10005CN11_C1_HE_small.jpg")
    # library(imager)
    # img = load.image("/Users/sanjavickovic/Desktop/SpoTter/STB37/output/rotated_images/10005CN47_C1_HE_small.jpg")
    out.file = paste("ST_Report_", array,"_", sample, "_", well, ".pdf", sep = "")
    pdf(out.file, width = 21, height = 7)
    # pdf("test.pdf", width = 21, height = 7)
    
    # plot HE image of tissue from experiment
    par(mar = c(1,1,3,1) , pty = "s", mfrow=c(1,6))
    par(mar = c(1,1,3.5,1) , pty = "s", cex.main = 2)
    HE=img
    plot(HE, main="ST HE:\n Tissue Image", axes=FALSE)
    par(new=T)
    rect(1, 20, nrow(HE), ncol(HE)-20, lwd = 10)
    par(new=T)
    rect(1,20,nrow(HE)/8,nrow(HE)/7, lwd = 10, density = 100)

    # plot frame and all ST points
    par(mar = c(1,1,3,0.5) , pty = "s")
    plot(xa, ya, main = "ST Heatmap:\n All points ", cex.main=2, cex = 2, ylim=c(35, 1), xlim=c(1, 33), ylab = '', xlab = '', xaxt='n', col = dat.col, yaxt='n', bty = 'n', pch = 15)
    points(xs, ys, cex = 0.45, ylim=c(35, 1), xlim=c(1, 33), ylab = '', xlab = '', xaxt='n', col = "black", yaxt='n', bty = 'n', pch = 4)
    par(new=T)
    rect(1, 1, 33, 35, lwd = 10)
    par(new=T)
    rect(1,1,5,5, lwd = 10, density = 100)

    # subset matrix
    A = all_barcodes
    A = t(A)
    B = A[, colnames(A) %in% str_replace(spots, "_", "x")]
    B = as.data.frame(A[rowSums(A) != 0,])
    B = mean.region.data(A)
    za = log2(colSums(B)+1)
    breaks = seq(0, max(za), by = (max(za)/254))
    cl.scale = wes_palette("Zissou1", 256, type = "continuous")
    dat.col = cl.scale[as.numeric(cut(za,breaks = breaks))]
    dat.col[is.na(dat.col)] <- cl.scale[1] 

    # plot UMIs per spot only under tissue
    par(mar = c(1,1,3,0.5) , pty = "s")
    plot(xa, ya, main = "ST Heatmap:\n Tissue points ",cex.main=2, cex = 2, ylim=c(35, 1), xlim=c(1, 33), ylab = '', xlab = '', xaxt='n', col = cl.scale[1], yaxt='n', bty = 'n', pch = 15)
    points(xs, ys, cex = 2, ylim=c(35, 1), xlim=c(1, 33), ylab = '', xlab = '', xaxt='n', col = dat.col, yaxt='n', bty = 'n', pch = 15)
    
    par(new=T)
    rect(1, 1, 33, 35, lwd = 10)
    par(new=T)
    rect(1,1,5,5, lwd = 10, density = 100)
    
    # plot scale
    xl <- 0
    yb <- 15
    xr <- 0.1
    yt <- 30
    par(mar = c(0,0,0,0))
    plot(NA,type="n",ann=FALSE,xlim=c(0,1),ylim=c(0,37),xaxt="n",yaxt="n",bty="n")
    rect(
      xl, head(seq(yb,yt,(yt-yb)/255),-1),
      xr, tail(seq(yb,yt,(yt-yb)/255),-1),
      col=cl.scale.all,
      border=cl.scale.all
    )
    text(adj = 0, labels = seq(round(0, digits = 2), round(b.max, digits = 2), by = (round(b.max, digits = 2)-round(0, digits = 2))/5), x = 0.12, y = seq(yb, yt, by=(yt-yb)/5), las=1, cex=2, col=c("black"))
    text(labels = "median norm", cex = 2, x = (xr-xl)/2, y = (yt-yb)+(yt-yb)/2, adj=0.5, srt = 90)
    text(labels = "tissue spots", cex = 2, x = 1.9*xr, y = (yt-yb)/2)
    for (i in c(1:6)){
      lines(x = c(0.1,0.12), y = rep(seq(yb,yt, by=(yt-yb)/5)[i],2), type = "l", pch=22, cex = 1)
    }
    points(x = 4*xr, y = (yt-yb)/2, pch = 7, cex = 3)
    par(new=T)
    rect(xl, yb, xr, yt, lty = 1, lwd = 1, border = "black") 
    dev.off()
    
    # output into new pdf file
    out.file = paste("ST_Report_Violins_", array,"_", sample, "_", well, ".pdf", sep = "")
    pdf(out.file, width = 21, height = 7)
    # pdf("test.pdf", width = 21, height = 7)
    
    # plot violins under/outside tissue from experiment
    par(mar = c(2,3,3,1) , pty = "s", mfrow=c(1,6))
    # plot avg reads inside/outside tissue boundaries
    A = all_barcodes
    A = t(A)
    B = A[, colnames(A) %in% str_replace(spots, "_", "x")]
    C = A[, !colnames(A) %in% str_replace(spots, "_", "x")]
    
    # violing plot for UMIs
    inside <- data.frame(colSums(B), "Density")
    names(inside) = c("UMIS", "Sample")
    outside <- data.frame(colSums(C), "Density")
    names(outside) = c("UMIS", "Sample")
    
    #library(vioplot)
    par(mar = c(1,3,2,0.5) , pty = "s")
    vioplot(UMIS~Sample, data=outside, ylim = c(0, 0.8*max(inside$UMIS)), axes=F, cex.axis = 2, col = "#239D8B", plotCentre = "point", side = "left", pchMed = 21, colMed = "black", colMed2 = "black")
    vioplot(UMIS~Sample, data=inside, axes=F,  cex.axis = 2, col = "#B89942", plotCentre = "point", side = "right", pchMed = 21, colMed = "black", colMed2 = "black", add = T) 
    points(1:1, median(outside$UMIS), pch = 21, col = "#239D8B", bg = "#239D8B")
    points(1:1, median(inside$UMIS), pch = 21, col = "#B89942", bg = "#B89942")
    title(main = "UMIs", cex.main=2)
    legend("topleft", fill = c("#239D8B", "#B89942"), legend = c("Outside tissue", "Inside tissue"), box.lwd = "n")

    # violing plot for genes
    b = apply(B,2,function(x)length(x[x>0]))
    c = apply(C,2,function(x)length(x[x>0]))
    
    inside <- data.frame(b,  "Density")
    names(inside) = c("Genes", "Sample")
    outside <- data.frame(c, "Density")
    names(outside) = c("Genes", "Sample")
    
    par(mar = c(1,3,2,0.5) , pty = "s")
    vioplot(Genes~Sample, ylim = c(0, 0.8*max(inside$Genes)), data=outside,  cex.axis = 2, axes=F, col = "#239D8B", plotCentre = "point", side = "left", pchMed = 21, colMed = "black", colMed2 = "black")
    vioplot(Genes~Sample, data=inside,  cex.axis = 2, axes=F, col = "#B89942", plotCentre = "point", side = "right", pchMed = 21, colMed = "black", colMed2 = "black", add = T) 
    points(1:1, median(outside$Genes), pch = 21, col = "#239D8B", bg = "#239D8B")
    points(1:1, median(inside$Genes), pch = 21, col = "#B89942", bg = "#B89942")
    title(main = "Genes", cex.main = 2)
    legend("topleft", fill = c("#239D8B", "#B89942"), legend = c("Outside tissue", "Inside tissue"), box.lwd = "n")
    
    dev.off()
    
}



exp.under.tissue = function(exp.values, spots){
  
  #clean up spots names
  # spots = read.delim("/Users/svickovi/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SpoTter/STB30/output/spots/10005CN11_C1_stdata_under_tissue_spots.txt", header = T)
  spots = str_remove(colnames(spots), "X")
  
  # clean up column names
  A = all_barcodes
  A = t(A)
  colnames(A) = str_replace(colnames(A), "x", "_")
  
  # select barcodes outside tissue
  selected.barcodes = A[,colnames(A) %in% spots]
  
  #return matrix
  return(selected.barcodes)
}

exp.not.under.tissue = function(exp.values, spots){
  
  #clean up spots names
  # spots = read.delim("/Users/svickovi/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SpoTter/STB30/output/spots/10005CN11_C1_stdata_under_tissue_spots.txt", header = T)
  spots = str_remove(colnames(spots), "X")
  
  # clean up column names
  A = all_barcodes
  A = t(A)
  colnames(A) = str_replace(colnames(A), "x", "_")
  
  # select barcodes outside tissue
  selected.barcodes = A[,!colnames(A) %in% spots]

  #return matrix
  return(selected.barcodes)
}


write.coor.matrix = function(exp.values) {

  write.table(exp.values, file=paste(array, sample, "_", well,   "_stdata_under_tissue_IDs.txt",sep=""), sep = "\t", quote=F, col.names=T, row.names=T)

}

write.coor.spots = function(exp.values) {
  
  write.table(t(colnames(exp.values)), file=paste(array, sample, "_", well,   "_stdata_under_tissue_spots.txt",sep=""), sep = "\t", quote=F, col.names=F, row.names=F)
  
}

write.not.coor.matrix = function(exp.values) {
  
  write.table(exp.values, file=paste(array, sample, "_", well,   "_stdata_outside_tissue_IDs.txt",sep=""), sep = "\t", quote=F, col.names=T, row.names=T)
  
}

top.genes.plot = function(all_barcodes, spots){
  
  #all_barcodes= read.delim("/Users/svickovi/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SpoTter/STB48/input/10015_CN96_C1_stdata.tsv", header=T, row.names=1)
  #spots = read.delim("/Users/svickovi/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SpoTter/STB48/output/spots/10015CN96_C1_stdata_under_tissue_spots.txt", header = T)
  spots = str_remove(colnames(spots), "X")
  
  # clean up column names
  A = all_barcodes
  A = convert.gene.ids(A)
  A = t(A)

  colnames(A) = str_replace(colnames(A), "x", "_")
  xa = sapply(strsplit(colnames(A), split = "_"),"[[",1)
  ya = sapply(strsplit(colnames(A), split = "_"),"[[",2)
  
  # select barcodes inside tissue
  A = A[,colnames(A) %in% spots]
  nmtop = names(head(sort(rowSums(A), decreasing = T), n=5))
  A = mean.region.data(A)
  A = log2(A+1)
  A = A[nmtop,]
  
  xs = sapply(strsplit(colnames(A), split = "_"),"[[",1)
  ys = sapply(strsplit(colnames(A), split = "_"),"[[",2)
  sm = A
  sm[is.na(sm)] <- 0
  mincol = min(sm) 
  maxcol = max(sm)
  #library(wesanderson)
  
  out.file = paste("ST_Report_Top_Genes_", array,"_", sample, "_", well, ".pdf", sep = "")
  pdf(out.file, width = 21, height = 7)
  # pdf("test.pdf", width = 21, height = 7)
  
  # plot all five top genes in one plot
  par(mar = c(0,1,0,1) , pty = "s", mfrow=c(1,6))

  for (i in row.names(A)){
    
    # plot frame and first gene ST points under tissue
    # subset matrix
    B = A[i,]

    # get color scale so it is the same for all sampels
    za = log2(B+1)
    breaks = seq(0, maxcol, by = (maxcol-0)/256)
    cl.scale = wes_palette("Zissou1", 255, type = "continuous")
    dat.col = cl.scale[as.numeric(cut(as.numeric(za), breaks = breaks))]
    dat.col[is.na(dat.col)] <- cl.scale[1] 

    par(mar = c(1,1,4,0.5) , pty = "s")
    plot(xa, ya, main = paste0("ST Heatmap:\n ", i),  cex.main=2, cex = 2, ylim=c(35, 1), xlim=c(1, 33), ylab = '', xlab = '', xaxt='n', col = cl.scale[1], yaxt='n', bty = 'n', pch = 15)
    points(xs, ys, cex = 2, ylim=c(35, 1), xlim=c(1, 33), ylab = '', xlab = '', xaxt='n', col = dat.col, yaxt='n', bty = 'n', pch = 15)
    par(new=T)
    rect(1, 1, 33, 35, lwd = 10)
    par(new=T)
    rect(1,1,5,5, lwd = 10, density = 100)
  
  }

  # plot scale
  xl <- 0
  yb <- 15
  xr <- 0.1
  yt <- 30
  par(mar = c(0,0,0,0))
  plot(NA,type="n",ann=FALSE,xlim=c(0,1),ylim=c(0,37),xaxt="n",yaxt="n",bty="n")
  rect(
    xl, head(seq(yb,yt,(yt-yb)/255),-1),
    xr, tail(seq(yb,yt,(yt-yb)/255),-1),
    col=cl.scale,
    border=cl.scale
  )
  text(adj = 0, labels = seq(round(0, digits = 2), round(max(za), digits = 2), by = (round(max(za), digits = 2)-round(0, digits = 2))/5), x = 0.12, y = seq(yb,yt, by=(yt-yb)/5),las=1,cex=2, col=c("black"))
  text(labels = "mean exp",cex = 2,  x = (xr-xl)/2, y = (yt-yb)+(yt-yb)/2, adj=0.5, srt = 90)
  for (i in c(1:6)){
    lines(x = c(0.1,0.12), y = rep(seq(yb,yt, by=(yt-yb)/5)[i],2), type = "l", pch=22, cex = 1)
  }
  par(new=T)
  rect(xl, yb, xr, yt, lty = 1, lwd = 1, border = "black") 
  
  dev.off()

}

interesting.genes.plot = function(all_barcodes, spots, genes){
  
  # all_barcodes = read.delim("/Users/sanjavickovic/Desktop/SpoTter/STB37/input/10005_CN47_C1_stdata.tsv", sep="\t",quote="", row.names = 1)
  # spots = read.delim("/Users/sanjavickovic/Desktop/SpoTter/STB37/output/spots/10005CN47_C1_stdata_under_tissue_spots.txt", header = T)
  spots = str_remove(colnames(spots), "X")
  
  # clean up column names
  A = all_barcodes
  A = convert.gene.ids(A)
  #head(A)
  A = t(A)

  colnames(A) = str_replace(colnames(A), "x", "_")
  xa = sapply(strsplit(colnames(A), split = "_"),"[[",1)
  ya = sapply(strsplit(colnames(A), split = "_"),"[[",2)
  
  # select barcodes inside tissue
  A = A[,colnames(A) %in% spots]
  A = mean.region.data(A)
  A = log2(A+1)
  A = A[rowSums(A)>0,]

  # check if interesting genes exhist in the dataset; otherwise pick 6 random genes
  if (nrow(A[row.names(A) %in% genes,]) < 8){
    nmtop = names(head(sort(rowSums(A), decreasing = T), n=50))
    ext = sample(nmtop[!nmtop %in% genes], 8-nrow(A[row.names(A) %in% genes,]))
    gene = c(genes, ext)
  } else gene=genes
  
  # subset based on gene
  A = A[row.names(A) %in% gene,]
  xs = sapply(strsplit(colnames(A), split = "_"),"[[",1)
  ys = sapply(strsplit(colnames(A), split = "_"),"[[",2)
  
  sm = A
  sm[is.na(sm)] <- 0
  mincol = min(sm)
  maxcol = max(sm)
  
  out.file = paste("ST_Report_Interesting_Genes_", array,"_", sample, "_", well, ".pdf", sep = "")
  pdf(out.file, width = 21, height = 7)
  # pdf("test.pdf", width = 21, height = 7)
  
  # plot all five top genes in one plot
  par(mar = c(0,1,0,1) , pty = "s", mfrow=c(2,8))
  
  for (i in rownames(A)){
    
    # plot frame and first gene ST points under tissue
    # subset matrix
    B = sm[i,]
    
    # get color scale so it is the same for all sampels
    za = B
    mincol = min(za)
    maxcol = max(za)
    breaks = seq(mincol, maxcol, by = (maxcol-mincol)/254)
    cl.scale = wes_palette("Zissou1", 256, type = "continuous")
    dat.col = cl.scale[as.numeric(cut(as.numeric(za), breaks = breaks))]
    dat.col[is.na(dat.col)] <- cl.scale[1] 
    
    par(mar = c(1,1,1,0.5) , pty = "s")
    plot(xa, ya, main = paste0(i),  cex.main=2, cex = 2, ylim=c(35, 1), xlim=c(1, 33), ylab = '', xlab = '', xaxt='n', col = cl.scale[1], yaxt='n', bty = 'n', pch = 15)
    points(xs, ys, cex = 1.5, ylim=c(35, 1), xlim=c(1, 33), ylab = '', xlab = '', xaxt='n', col = dat.col, yaxt='n', bty = 'n', pch = 15)
    par(new=T)
    rect(1, 1, 33, 35, lwd = 10)
    par(new=T)
    rect(1,1,5,5, lwd = 10, density = 100)
    
    # plot scale
    xl <- 0
    yb <- 15
    xr <- 0.1
    yt <- 30
    par(mar = c(0,0,0,0))
    plot(NA,type="n",ann=FALSE,xlim=c(0,1),ylim=c(0,37),xaxt="n",yaxt="n",bty="n")
    rect(
      xl, head(seq(yb,yt,(yt-yb)/255),-1),
      xr, tail(seq(yb,yt,(yt-yb)/255),-1),
      col=cl.scale,
      border=cl.scale
    )
    text(adj = 0, labels = seq(round(min(za), digits = 2), round(max(za), digits = 2), by = (round(max(za), digits = 2)-round(min(za), digits = 2))/5), x = 0.12, y = seq(yb,yt, by=(yt-yb)/5),las=1,cex=2, col=c("black"))
    text(labels = "mean exp", cex = 2, x = (xr-xl)/2, y = (yt-yb)+(yt-yb)/2, adj=0.5, srt = 90)
    for (i in c(1:6)){
      lines(x = c(0.1,0.12), y = rep(seq(yb,yt, by=(yt-yb)/5)[i],2), type = "l", pch=22, cex = 1)
    }
    par(new=T)
    rect(xl, yb, xr, yt, lty = 1, lwd = 1, border = "black") 
  }
  
  dev.off()
  
}



