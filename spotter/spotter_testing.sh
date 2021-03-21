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
echo "${img_file_nm}"
sm=$(echo "${ABSPATH%/*}/processing") # output dir for small images and inferred points tsv
mkdir "$(echo "$sm")"
now=$(date +"%T")
echo "Current time : $now"
# mogrify -path "$(echo "$sm")" -limit thread 2 -resize 1% -quality 100%% $img_file_nm # to run locally
mogrify -path "$(echo "$sm"/)" -define jpeg:size=$(echo 500x500) -resize 8% -quality 100%% $img_file_nm # to run locally
now=$(date +"%T")
echo "Current time : $now"
# get spotter path
sp="$(echo ${ABSPATH%/*/*}/STaut/spotter)"

# run spotter R script on small image
RScript "$(echo $sp/st_spotter_RUN.R)" "$(echo "$sm/$img")" "$(echo "$sp")" "$(echo "$sm")"
now=$(date +"%T")
echo "Current time : $now"
