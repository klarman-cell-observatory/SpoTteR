#!/usr/bin/env bash
# bash script for running st_spotter

# invoking help by user
if [ "$1" == "-h" ]; then
echo "Usage: `basename $0` -batch (ex. STB30) -image (ex. 10005_CN11_STB30.C1-Spot000001.jpg)"
exit 0
fi

# input ST batch number
batch=${1?"Error: no ST batch name given"}

# go to folder with images
ABSPATH=`pwd`

# input jpg file
img=${2?"Error: no ST image file given"}
img_file_nm=$(echo ${img##*/})

# make variables
array=$(echo $img_file_nm | cut -f 1 -d "_")
cn=$(echo $img_file_nm | cut -f 2 -d "_")
batch=$(echo $img_file_nm | cut -f 3 -d "_" | cut -f 1 -d ".")
cde=$(echo $img_file_nm | cut -f 3 -d "_" | cut -f 2 -d "." | cut -d "-" -f1)

# make smaller images at 500x500 px
# echo ${ABSPATH%/*}
sm=$(echo "${ABSPATH%/*}/processing") # output dir for small images and inferred points tsv
mkdir "$(echo "$sm")"
now=$(date +"%T")
echo "Current time : $now"

mogrify -path "$(echo "$sm"/)" -define jpeg:size=$(echo 500x500) -resize 16% -quality 100%% $img_file_nm # to run locally
now=$(date +"%T")
echo "Current time : $now"
# get spotter path
sp="$(echo ${ABSPATH%/*/*}/STaut/spotter)"

# run spotter R script on small image
RScript "$(echo $sp/st_spotter_RUN_hex.R)" "$(echo "$sm/$img")" "$(echo "$sp")" "$(echo "$sm")"
now=$(date +"%T")
echo "Current time : $now"
#rm small image that is not needed anymore
rm "$(echo "$sm/$img")"

# read in file from ST pipeline (tsv file)
st="$(echo ${ABSPATH%/*})/input/$(echo "$array"_"$cn"_"$cde"_stdata.tsv)"

# run cropping script
bash "$(echo $sp/cropping_hex.sh)" "$(echo ${ABSPATH%/*}/processing/${img_file_nm/.jpg/}_inferred_points.tsv)" "$(echo basename $0)"

# run spotter QC R script on spotter output from previous step
RScript "$(echo $sp/st_spotter_QC_RUN_hex.R)" "$(echo "$st")" $batch "$(echo "$sm")" "$(echo "$sp")"
rm "$(echo $sm/"$array""$cn"_"$cde"_HE.jpg)"

# trim pdf QC reports # in QC_plots dir already
cd "$(dirname "$sm")/output/QC_plots"
convert -density 300 ST_Report_Interesting_Genes_"$array"_"$cn"_"$cde".pdf -size 150%x150% ST_Report_Interesting_Genes_"$array"_"$cn"_"$cde".pdf
convert -density 300 ST_Report_Top_Genes_"$array"_"$cn"_"$cde".pdf -shave 10x400 -size 150%x150% -gravity East -chop 750x0 ST_Report_Top_Genes_"$array"_"$cn"_"$cde".pdf
convert -density 300 ST_Report_Violins_"$array"_"$cn"_"$cde".pdf -shave 10x400 -size 150%x150% -gravity East -chop 2500x0 ST_Report_Violins_"$array"_"$cn"_"$cde".pdf
convert -density 300 ST_Report_"$array"_"$cn"_"$cde".pdf -shave 10x400 -size 150%x150% -gravity East -chop 750x0 ST_Report_"$array"_"$cn"_"$cde".pdf

# run R Markdown QC spotter script
mkdir "$(dirname "$sm")/output/reports"
cd "$(dirname "$sm")/output/reports"
rpt_name=ST_QC_Report_"$array"_"$cn"_"$cde".Rmd
cp "$sp"/st_spotter_QC_report_hex.Rmd ./
mv st_spotter_QC_report_hex.Rmd "$rpt_name"
R -e rmarkdown::render"('$rpt_name')" --args "$(echo "${ABSPATH%/*}/output/QC_plots")" "$array"_"$cn"_"$cde"
rm "$rpt_name"
