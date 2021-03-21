#-------FUNCTIONS---------#
#find good point candidates based on initial blob detection
find_candidates = function(sub_data_table,n_nn,alpha) {
    dist_threshold = 0.03 # 3% relative error in the distance
    phi_threshold = 0.03 # 3%*pi absolute error in the angle
    inner_combis=as.data.table(t(combn(n_nn,2)))
    inner_combis[,dist_quality:=abs(sub_data_table[V1]$dist-sub_data_table[V2]$dist)/(sub_data_table[V1]$dist+sub_data_table[V2]$dist),]
    inner_combis[,linear_phi_quality:=abs((((sub_data_table[V1]$phi-sub_data_table[V2]$phi)/pi+4.0)%%2.0)-1.0),]
    inner_combis[,ortho_phi_quality:=abs((((sub_data_table[V1]$phi-sub_data_table[V2]$phi)/pi+4.0)%%1.0)-0.5),]
    inner_combis[,dist_candidate:=alpha*tanh(1/alpha*dist_threshold/dist_quality),]
    inner_combis[,linear_phi_candidate:=alpha*tanh(1/alpha*phi_threshold/linear_phi_quality),]
    inner_combis[,ortho_phi_candidate:=alpha*tanh(1/alpha*phi_threshold/ortho_phi_quality),]
    inner_combis[,candidate:=dist_candidate+linear_phi_candidate+ortho_phi_candidate,]
    
    return(inner_combis[candidate>1.8*alpha,.(V2=sub_data_table[V1]$V2,V3=sub_data_table[V2]$V2,dist_V1V2=sub_data_table[V1]$dist,dist_V1V3=sub_data_table[V2]$dist,phi_V1V2=sub_data_table[V1]$phi,phi_V1V3=sub_data_table[V2]$phi,linear_phi_quality=linear_phi_quality,ortho_phi_quality=ortho_phi_quality,dist_quality=dist_quality,dist_candidate=dist_candidate,linear_phi_candidate=linear_phi_candidate,ortho_phi_candidate=ortho_phi_candidate,candidate=candidate),])
}


#objective function for fitter
modified_chi2_lines=function(parameters, angle, points){
    length = parameters[1]
    origin = parameters[2]
    cos_angle = cos(angle)
    sin_angle = sin(angle)
    # Create the points on the grid.
    points[,rot_x:=cos_angle*x-sin_angle*y,]
    points[,rot_y:=sin_angle*x+cos_angle*y,]
    points[,diff:=((rot_y-origin+length/2)%%length)-length/2,]
    points[,penalty:=-cos(diff/length*2*pi),]
    total_penalty = sum(points$penalty)
    return(total_penalty)
}

# Self brewed global minimizer
sbgm = function(points, range_dy, range_oy, angle){
    # global minimum for dy and oy
    eval_grid=as.data.table(expand.grid(dy = range_dy, oy = range_oy))
    sub_fun=function(x){ modified_chi2_lines(x, points = points, angle = angle) }
    eval_grid$eval=apply(eval_grid, 1, sub_fun)
    best_guess = eval_grid[which.min(eval)]
    fitter <- nlminb(
    start=c(best_guess$dy, best_guess$oy),
    objective=modified_chi2_lines, points=points, angle = angle
    )
    return(fitter)
}

#create grid
create_grid=function(angles,length,origin,nrows,ncols){
    # Create the points on the grid.
    basis <- rbind(cos(angles * pi/180), sin(angles * pi/180)) %*% diag(length)
    x <- as.vector(outer(1:ncols, 1:nrows, FUN=function(x,y) x-1))
    y <- as.vector(outer(1:ncols, 1:nrows, FUN=function(x,y) y-1))
    grid <- as.data.table(t(basis %*% rbind(x, y) + origin))
    setnames(grid, names(grid),c("x","y"))
    return(grid)
}

#remove point that don't fit the grid well in order to perform another rould of fitting based on an improved set of points
remove_worst_fitting_points = function (pts, f_chi, relative_share_to_remove=0.002) {
    better_pts = pts
    n_rows = nrow(better_pts)
    better_pts[,f_chi:=f_chi(as.data.table(.SD)),by=1:n_rows]
    n_rows_to_keep = ceiling(n_rows * (1.0 - relative_share_to_remove))
    #  if(n_rows_to_keep==n_rows){n_rows_to_keep=n_rows_to_keep-1} #always remove at least 1 point
    #better_pts=head(better_pts[,.SD,keyby=f_chi], n_rows_to_keep)
    better_pts=head(better_pts[order(f_chi)], n_rows_to_keep)
    return(better_pts[,c("x","y","f_chi"),])
}

pink=TRUE
bubble=TRUE

#-------ANALYSIS---------#

process_HE=function(image_file,pink=FALSE,bubble=TRUE){
    #-------LOAD DATA---------#
    HE_image_name=gsub(".jpg","",image_file)
    HE=load.image(paste0(HE_image_name,".jpg")) #%>%imresize(scale=0.01)
    HE_bw=grayscale(HE)
    p0=plot(HE, xlim = c(0, nrow(HE)), ylim = c(0, ncol(HE)))
    p0
    
    #-------ANALYSIS---------#
    #remove non-grid features (i.e. the tissue sections)
    #HE2=copy(HE_bw)
    HE_spl=imsplit(HE,"c")
    HE2=HE_spl$`c = 1`
    cutoff=quantile(HE2,0.1)
    HE2[,,,1][HE2[,,,1]<cutoff]=min(HE2[,,,1][HE2[,,,1]>cutoff])
    p1=ggplot(as.data.frame(HE2),aes(x,y))+ggpl_settings+ggtitle("Removed non-grid features")
    p1
    
    #detect blobs (grid points)
    ib=isoblur(crop.borders(HE2,nPix = 4), 3)
    Hdet <- with(imhessian(ib),(xx+yy))
    hc.clean=Hdet ; hc.clean[dilate_square(hc.clean,3) != hc.clean] <- 0
    hc.clean=pad(hc.clean,4,pos=-1,"xy")
    
    #store points in data.table
    hc_dt=as.data.table(as.data.frame(hc.clean))
    hc_dt_red=hc_dt[value>0.0001&x<dim(HE)[1]-5&y<dim(HE)[2]-5&x>5&y>5]
    sc_pts=hc_dt_red[,c("x","y"),]
    sc_pts[,xy:=paste0(x,"_",y)]
    
    #calculate all x and y distances between points (to select those points that likely are part of the grid)
    combis=as.data.table(t(cbind(combn(nrow(sc_pts),2),combn(nrow(sc_pts),2)[c(2,1),]))) # duplicate the combn to get (i,j) and (j,i)
    combis[,diff_x:=sc_pts[V1]$x-sc_pts[V2]$x,]
    combis[,diff_y:=sc_pts[V1]$y-sc_pts[V2]$y,]
    combis[,dist:=sqrt(diff_x*diff_x+diff_y*diff_y),]
    combis[,phi:=atan2(diff_y,diff_x),]
    #keep only the nearest neighbours
    n_nn = 8 # number of nearest neighbours to find
    nn_combis=combis[,head(.SD[,.SD,keyby=dist],n_nn),by=V1]
    
    #___________________#
    
    alpha = 0.7 # alpha*tanh(1/alpha*.) to cutoff extremely large values for the candidate scores
    candidates=nn_combis[,find_candidates(.SD,n_nn,alpha),by=V1]
    
    candidates_long=melt(candidates[,c("V1","V2","V3","candidate"),],id.vars = "candidate")
    selected_candidates=candidates_long[,list(combined_candidate=sum(candidate)),by="value"]
    
    min_number_of_candidates = 3
    selected_indices = selected_candidates[combined_candidate > min_number_of_candidates * 2*alpha]$value
    length(selected_indices)
    
    p2=ggplot()+geom_point(data=sc_pts,aes(x,y),pch=1)+
    geom_point(data=sc_pts[selected_indices],aes(x,y),pch=3)+
    coord_fixed()+
    scale_y_reverse()+
    ggtitle("Selected blobs to start")
    p2
    
    
    #create feature mask
    if (pink == TRUE){
        HE3=HE_spl$`c = 3`
    }else{
        HE3=HE_spl$`c = 1`
    }
    #remove grey artefacts (e.g. bubbles)
    if (bubble == TRUE){
        c1=HE_spl$`c = 1`[,,,1]
        c2=HE_spl$`c = 2`[,,,1]
        c3=HE_spl$`c = 3`[,,,1]
        
        sl=0.065#0.045
        bubbles=c1<=c2+sl&c1>=c2-sl&c1<=c3+sl&c1>=c3-sl&c3<=c2+sl&c3>=c2-sl
        bubbles=as.numeric(bubbles)
        bubble_mask=HE3
        bubble_mask[,,,1]=bubbles
    }else{
        bubble_mask=HE3
        bubble_mask[,,,1]=0
    }
    
    #bulid feature mask
    start_coords=which(HE3[,,,1] > 0.9,arr.ind = TRUE)[1,]
    mean_bg=mean(HE3[,,,1][1:dim(HE3)[1],c(0:20,dim(HE3)[2]-20:dim(HE3)[2])])
    sd_bg=sd(HE3[,,,1][1:dim(HE3)[1],c(0:20,dim(HE3)[2]-20:dim(HE3)[2])])
    msk= px.flood(HE3,y=start_coords["col"],x=start_coords["row"],sigma=.1) %>% as.cimg
    HE3[,,,1][HE3[,,,1]>mean_bg-min(sd_bg/2,0.05)]=1
    HE3[,,,1][HE3[,,,1]<=mean_bg-min(sd_bg/2,0.05)]=0
    HE_msk=(msk|HE3)|bubble_mask
    HE4=isoblur(HE_msk,10)
    HE4[,,,1][HE4[,,,1]>0.9|HE_msk[,,,1]==1]=1
    HE5=isoblur(HE4,1)
    feature_dt=as.data.table(as.data.frame(HE5))
    feature_dt[,xy:=paste0(x,"_",y)]
    
    #now fit grid parameters to selected candidate points
    good_pts=sc_pts[selected_indices][!xy%in%feature_dt[as.double(value)<0.999]$xy]
    dy_fitter=sbgm(good_pts, seq(10,13,0.1), seq(0,13,1.0), 0)
    dx_fitter=sbgm(good_pts, seq(10,13,0.1), seq(0,13,1.0), pi/2)
    
    #rough_grid=create_grid(angles = c(0,90),length = c(dx_fitter$par[1],dy_fitter$par[1]),origin = c(dx_fitter$par[2]-1*dx_fitter$par[1],dy_fitter$par[2]-1*dy_fitter$par[1]),ncols = 42,nrows = 42)
    rough_grid=create_grid(angles = c(0,90),length = c(dx_fitter$par[1],dy_fitter$par[1]),origin = c(dx_fitter$par[2]-1*dx_fitter$par[1],dy_fitter$par[2]-1*dy_fitter$par[1]),ncols = 42,nrows = 42)
    
    p3=ggplot(as.data.frame(HE2),aes(x,y))+ggpl_settings+
    annotate(geom = "point",x = good_pts$x,y = good_pts$y,pch=1,color="red",size=1)+
    annotate(geom = "point",x=rough_grid$x,y=rough_grid$y,pch=3, col="Blue",size=1)+
    ggtitle("Initial fit")
    p3
    
    #plot(HE2,main="Initial fit")
    #points(good_pts[,c("x","y")])
    #points(create_grid(angles = c(0,90),length = c(dx_fitter$par[1],dy_fitter$par[1]),origin = c(dx_fitter$par[2]-100*dx_fitter$par[1],dy_fitter$par[2]-100*dy_fitter$par[1]),ncols = 300,nrows = 300), asp=1, pch=3, col="Blue")
    
    
    # find and remove_worst_fitting_points
    better_pts = remove_worst_fitting_points(good_pts,f_chi=function(pts) { return (
        modified_chi2_lines(dx_fitter$par, pi/2, pts)
        + modified_chi2_lines(dy_fitter$par,    0, pts)
        ) },relative_share_to_remove = 0.25)
    dy_fitter=sbgm(better_pts, seq(10,13,0.1), seq(0,13,1.0), 0)
    dx_fitter=sbgm(better_pts, seq(10,13,0.1), seq(0,13,1.0), pi/2)
    # get the global offset by checking the minimal bins where the points fall into
    better_pts[,bin_x:=floor((x-dx_fitter$par[2])/dx_fitter$par[1]+0.5)]
    better_pts[,bin_y:=floor((y-dy_fitter$par[2])/dy_fitter$par[1]+0.5)]
    better_pts[,bin_x_count:=.N,by="bin_x"]
    better_pts[,bin_y_count:=.N,by="bin_y"]
    offset_x = min(better_pts$bin_x) * dx_fitter$par[1] + dx_fitter$par[2]
    offset_y = min(better_pts$bin_y) * dy_fitter$par[1] + dy_fitter$par[2]
    grid_size_x = max(better_pts$bin_x)-min(better_pts$bin_x)+1
    grid_size_y = max(better_pts$bin_y)-min(better_pts$bin_y)+1
    
    #fitted_grid=create_grid(angles = c(0, 90),length = c(dx_fitter$par[1],dy_fitter$par[1]),origin = c(offset_x,offset_y),ncols = grid_size_x,nrows = grid_size_y)
    fitted_grid=create_grid(angles = c(0, 90),length = c(dx_fitter$par[1],dy_fitter$par[1]),origin = c(offset_x,offset_y),ncols = grid_size_x,nrows = grid_size_y)
    
    improvement_counter=1
    while (((grid_size_x>exp_cols)|(grid_size_y>exp_rows))&improvement_counter<=50){
        print(paste0("improvement loop # ", improvement_counter))
        
        ##remove bad points regardless of position
        better_pts = remove_worst_fitting_points(better_pts, function(pts) { return (
            modified_chi2_lines(dx_fitter$par, pi/2, pts)
            + modified_chi2_lines(dy_fitter$par,    0, pts)
            ) } )
        dy_fitter=sbgm(better_pts, seq(10,13,0.1), seq(0,13,1.0), 0)
        dx_fitter=sbgm(better_pts, seq(10,13,0.1), seq(0,13,1.0), pi/2)
        print(c(0, dx_fitter$par[1], dy_fitter$par[1], dx_fitter$par[2], dy_fitter$par[2]))
        
        # get the global offset by checking the minimal bins where the points fall into
        better_pts[,bin_x:=floor((x-dx_fitter$par[2])/dx_fitter$par[1]+0.5)]
        better_pts[,bin_y:=floor((y-dy_fitter$par[2])/dy_fitter$par[1]+0.5)]
        # number of points per bin
        better_pts[,bin_x_count:=.N,by="bin_x"]
        better_pts[,bin_y_count:=.N,by="bin_y"]
        # combined f_chi per bin
        better_pts[,bin_x_f_chi:=sum(-f_chi),by="bin_x"]
        better_pts[,bin_y_f_chi:=sum(-f_chi),by="bin_y"]
        
        # get the grid size by checking the minimal and maximal bins where the points fall into
        grid_size_x = max(better_pts$bin_x)-min(better_pts$bin_x)+1
        grid_size_y = max(better_pts$bin_y)-min(better_pts$bin_y)+1
        
        ##if there are still too many rows or columns remove the weakest ones and re-fit
        if(grid_size_x>exp_cols){
            #remove=better_pts[bin_x%in%c(min(bin_x),max(bin_x)),bin_x[which.min(bin_x_count)],]
            remove=better_pts[bin_x%in%c(min(bin_x),max(bin_x)),bin_x[which.min(bin_x_f_chi)],]
            better_pts=better_pts[bin_x!=remove]
        }
        
        if(grid_size_y>exp_rows){
            #remove=better_pts[bin_y%in%c(min(bin_y),max(bin_y)),bin_y[which.min(bin_y_count)],]
            remove=better_pts[bin_y%in%c(min(bin_y),max(bin_y)),bin_y[which.min(bin_y_f_chi)],]
            better_pts=better_pts[bin_y!=remove]
        }
        #re-fit
        dy_fitter=sbgm(better_pts, seq(10,13,0.1), seq(0,13,1.0), 0)
        dx_fitter=sbgm(better_pts, seq(10,13,0.1), seq(0,13,1.0), pi/2)
        print(c(0, dx_fitter$par[1], dy_fitter$par[1], dx_fitter$par[2], dy_fitter$par[2]))
        # get the global offset by checking the minimal bins where the points fall into
        better_pts[,bin_x:=floor((x-dx_fitter$par[2])/dx_fitter$par[1]+0.5)]
        better_pts[,bin_y:=floor((y-dy_fitter$par[2])/dy_fitter$par[1]+0.5)]
        # number of points per bin
        better_pts[,bin_x_count:=.N,by="bin_x"]
        better_pts[,bin_y_count:=.N,by="bin_y"]
        # get the grid size by checking the minimal and maximal bins where the points fall into
        grid_size_x = max(better_pts$bin_x)-min(better_pts$bin_x)+1
        grid_size_y = max(better_pts$bin_y)-min(better_pts$bin_y)+1
        
        
        offset_x = min(better_pts$bin_x) * dx_fitter$par[1] + dx_fitter$par[2]
        offset_y = min(better_pts$bin_y) * dy_fitter$par[1] + dy_fitter$par[2]
        
        #fitted_grid=create_grid(angles = c(0, 90),length = c(dx_fitter$par[1],dy_fitter$par[1]),origin = c(offset_x,offset_y),ncols = grid_size_x,nrows = grid_size_y)
        fitted_grid=create_grid(angles = c(0, 90),length = c(dx_fitter$par[1],dy_fitter$par[1]),origin = c(offset_x,offset_y),ncols = grid_size_x,nrows = grid_size_y)
        improvement_counter=improvement_counter+1
        
        #plot(HE_bw^3)
        #points(better_pts)
        message(grid_size_x)
        message(grid_size_y)
    }
    
    
    p4=ggplot(as.data.frame(HE2),aes(x,y))+ggpl_settings+
    annotate(geom = "point",x = better_pts$x,y = better_pts$y,pch=1,color="red",size=1)+
    annotate(geom = "point",x=fitted_grid$x,y=fitted_grid$y,pch=3, col="Blue",size=1)+
    ggtitle("Final fit")
    p4
    #p4=plot(HE2,main="Final fit")
    #p4=p4+points(better_pts[,c("x","y")], col="red")
    #p4=p4+points(fitted_grid[,c("x","y"),], pch=3, col="Blue")
    
    #keep precise values for scaling later
    setnames(fitted_grid,c("x","y"),c("x_precise","y_precise"))
    fitted_grid[,x:=round(x_precise,0),]
    fitted_grid[,y:=round(y_precise,0),]
    
    
    #now select gridpoints that overlap with the features
    #merge
    feature_spots=merge(fitted_grid,feature_dt[as.double(value)<0.999],all.x = TRUE,by=c("x","y"))
    feature_spots[,x_id:=.GRP,b="x"]
    feature_spots[,y_id:=.GRP,b="y"]
    
    #re-orient grid point IDs
    feature_spots[,x_id:=(exp_cols+1)-x_id,]
    feature_spots[,y_id:=(exp_rows+1)-y_id,]
    feature_spots[,bc:=paste0(x_id,"x",y_id),]
    
    feature_spots[,feature_spot:=!is.na(value),]
    feature_spots[,value:=NULL,]
    write.table(feature_spots,paste0(HE_image_name,"_inferred_points.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
    
    #sanity checks --> create warning, if any indication of problems while infering the grid
    #Nomber of columns or rows divert from expected
    col_row_check=""
    if (grid_size_x!=exp_cols|grid_size_y!=exp_rows) {col_row_check="\nError: Didn't find the expected number of colums or rows."}
    
    #Number of dots based on which the grid was inferred is low
    dot_number_check=""
    if (nrow(better_pts)<0.2*exp_rows*exp_cols){dot_number_check="\nWarning: Low number of grid points detected. Grid might be inaccurate"}
    
    #Low support of border rows or columns --> might indicate wrong origin
    border_support_check=""
    if (better_pts[c(which.min(bin_x),which.max(bin_x)),any(bin_x_count<3)]|better_pts[c(which.min(bin_y),which.max(bin_y)),any(bin_y_count<3)]){border_support_check="\nWarning: Low number of border points detected. Origin might be inaccurate."}
    
    p5=ggplot(as.data.frame(HE_bw),aes(x,y))+ggpl_settings+
    annotate(geom = "point",x = feature_spots$x,y = feature_spots$y,pch=3,color="blue",size=1)+
    annotate(geom = "point",x=feature_spots[feature_spot==TRUE]$x,y=feature_spots[feature_spot==TRUE]$y,pch=16, col="red",size=1)+
    ggtitle(paste0("Final grid and feature points: Grid pts=",max(feature_spots$x_id),"x",max(feature_spots$y_id),", feature pts=",sum(feature_spots$feature_spot),col_row_check,dot_number_check,border_support_check))
    
    #p5=plot(HE_bw^10,main="Final grid and feature points")+points(x=feature_spots$x,y=feature_spots$y,col="blue",pch=3)
    #p5=p5+points(x=feature_spots$x,y=feature_spots$y,col="blue",pch=3)
    #p5=p5+points(x=feature_spots[feature_spot==TRUE]$x,y=feature_spots[feature_spot==TRUE]$y,col="red",pch=16)
    
    if(col_row_check!=""|dot_number_check!=""|border_support_check!=""){
        p5=p5+theme(plot.title = element_text(color="red", size=10, face="bold.italic"))  }
    p5
    
    #now make report
    pdf(paste0(HE_image_name,"_fit_report.pdf"),height=8,width=6, useDingbats = F)
    plot(p5)
    plot(p4)
    plot(p1)
    plot(p3)
    plot(p2)
    plot(p0)
    
    
    dev.off()
}
