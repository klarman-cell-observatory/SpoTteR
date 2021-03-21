# SpoTteR
This is [code](./spotter) for automatically detecting (x,y) ST coordinates, H&E stained tissue and creating QC metrics from the corresponding ST or SM-Omics experiment (as given by ST Pipeline ([v.1.7.6](https://github.com/SpatialTranscriptomicsResearch/st_pipeline/releases/tag/1.7.6)). 

Note: for SpoTteR to detect ST coordinates and H&E tissue, recommended image size for the spot detection step is 500x500 px. Currenlty, SpoTteR enables you to input raw 10X (25kx25k) or 20X (40kx40k) jpg images and runs the image downsampling, spot detection, alignments and QC for you. Installation, should only take a few seconds on your destop computer.

### File Structure Overview
The code presumed following folder structure: 

    . 
    ├── STaut                                 # Your github repo 
    │   ├── spotter 
    │   │   ├── cropping_spotter.R            # runs cropping 
    │   │   ├── cropping_spotter_hex.R        # runs cropping  (hex matrix design)
    │   │   ├── cropping.sh                   # runs cropping 
    │   │   ├── spotter_functions.R           # f(x) for spotter 
    │   │   ├── spotter_10X.sh                # runs whole spotter tool with 10X imaging data
    │   │   ├── spotter_10X_hex.sh            # runs whole spotter tool with 10X imaging data (hex matrix design)
    │   │   ├── spotter_20X.sh                # runs whole spotter tool with 20X imaging data
    │   │   ├── spotter_20X_hex.sh            # runs whole spotter tool with 20X imaging data (hex matrix design)
    │   │   ├── spotter_testing.sh            # runs whole spotter tool with any imaging data and does only alignment (no QC)
    │   │   ├── st_spotter_QC_FUNCTIONS.R     # f(x) for QC 
    │   │   ├── st_spotter_QC_RUN.R           # runs QC 
    │   │   ├── st_spotter_RUN.R              # runs spotter f(x) 
    │   │   ├── st_spotter_QC_FUNCTIONS_hex.R # f(x) for QC  (hex matrix design)
    │   │   ├── st_spotter_QC_RUN_hex.R       # runs QC (hex matrix design)
    │   │   ├── st_spotter_RUN_hex.R          # runs spotter f(x) (hex matrix design)
    │   │   ├── st_spotter_QC_report.Rmd      # creates html QC report     
    ├── STBX                                  # Directory with all ST files from corresponding batch/experiment
    │   ├── input                             # tsv expression files from ST Pipeline v1.7.6
    │   ├── original_images                   # Original jpg images from the imaging platform used
    │   ├── pipeline_logs                     # pipeline logs generated from ST Pipeline v1.7.6 
    │   ├── output                            # This directory is created with SpoTter
    │   │   ├── aligned_counts                # This directory has output matrix files
    │   │   │   ├── *stdata_outside_tissue_IDs.txt*    # expression matrix with ST data outside the detected tissue boundaries
    │   │   │   ├── *stdata_under_tissue_IDs.txt*      # expression matrix with ST data inside the detected tissue boundaries
    │   │   ├── QC_logs           # This directory is created with SpoTter
    │   │   │   ├── *QC_data.txt* # This file holds mean UMIs, genes, # spots covered with tissue and top expressed genes
    │   │   ├── QC_plots          # This directory is created with SpoTter
    │   │   │   ├── *Alignment.pdf* # Spatial heatmap image of total detected UMIs per ST spot on the whole ST array area
    │   │   ├── rotated_images    # This directory has cropped and rotated images # to GUI
    │   │   ├── spots             # This directory has a single row file with (x_y) under tissue coordiantes # to GUI
    │   │   ├── json              # creates a json file with spots coordinates that can be preloaded to spaceranger pipeline   
    │   │   ├── reports           # html QC report (final output report from SpoTteR)
    │   └── processing            # Unit tests
    │   │   ├── *cropping_points.tsv*    # This file has 4 integers reffering to (x1,y1), width and height of the ST frame 
    │   │   ├── *fit_report.pdf*    # This file has initial and final spot and tissue alignment plots
    │   │   ├── *inferred_points.tsv*    # This file has px and (x,y) ST coordinates of each spot (TRUE=under tissue)
    └── ...

### Example [data](./spotter/data)
images: 10005_CN48_STB37.C1-Spot000001.jpg  
input stdata files: 10005_CN48_C1_stdata.tsv  
pipeline_log: 10005CN48_C1_log.txt  

### Running SpoTteR
`bash spotter_10X.sh STB39 10005_CN52_STB39.C1-Spot000001.jpg`

### Expected QC report
Please find example QC report generated with SpoTteR [here](./spotter/data/ST_QC_Report_10005_CN48_C1.html). 

### Renv
Env availabe [here](./spotter/data/sessionInfo_run.txt).
