#!/usr/bin/env bash
# bash script for croping images
# works with input inffered points tsv file

# read input
d=${1?"Error: no input points.tsv name given"}
p=$(echo ${d##*/})

# make variables
array=$(echo $p | cut -f 1 -d "_")
cn=$(echo $p | cut -f 2 -d "_")
batch=$(echo $p | cut -f 3 -d "_" | cut -f 1 -d ".")
cde=$(echo $p | cut -f 3 -d "_" | cut -f 2 -d "." | cut -d "-" -f1)

echo "Running ${array}_${cn}_${cde}"

#move to processing directory
ABSPATH=`pwd`
sm=$(echo "${ABSPATH%/*}/processing") # output dir for small images and inferred points tsv
cd "$(echo "$sm")"

# get spotter path
sp="$(echo ${ABSPATH%/*/*}/STaut/spotter)"

# Run cropping script
RScript "$(echo $sp/cropping_spotter_hex.R)" "$(echo "${array}_${cn}_${cde}")" $batch "$(echo "$sm")" "$(echo "$2")"

# read in crop file name
crop_file="$(echo "${p/_inferred_points.tsv/}_cropping_points.tsv")"
# echo "$crop_file"

# read in variables for cropping
width=$(cut -d "," -f1 $crop_file)
height=$(cut -d "," -f2 $crop_file)
xa=$(cut -d "," -f3 $crop_file)
ya=$(cut -d "," -f4 $crop_file)

# if array not aligned correctly
if [ -z "$ya" ]
then
    echo "ERROR: Exiting."
    exit 1 # terminate and indicate error
fi

# read in original image file name
img="$(echo "${ABSPATH%/*}/original_images/${p/_inferred_points.tsv/}.jpg")"

# Make 85% sized images from original images
if [ "${2##*/}" = "spotter_10X_hex.sh" ]; then
# This is for 10X
    mogrify -path ./ -resize 85% -quality 42.5%% "$(echo $img)"
# echo 10xres
else
# This is for 20X
    mogrify -path ./ -resize 42.5% -quality 85%% "$(echo $img)"
# echo 20xress
fi

# rename img variable to 85%
img="$(echo "$sm/${p/_inferred_points.tsv/}.jpg")"
echo "$img"

# check you are dealing with same cropping file and same image file
img_file_nm="$(echo ${img##*/})"

# crop
convert "$(echo $img)" -crop $(echo $width"x"$height+$xa+$ya) "$(echo "$array""$cn"_"$cde"_HE.jpg)" # otherise _HE_crop.jpg
#echo "starting to flip flop"

# mogrify -flop -flip "$(echo "$array""$cn"_"$cde"_HE.jpg)"
#convert -flip -flop "$(echo "$array""$cn"_"$cde"_HE_crop.jpg)" "$(echo "$array""$cn"_"$cde"_HE.jpg)"
#echo "finished flip flopping"

# copy croppped image to correct folder
mkdir ../output/rotated_images
cp "$(echo "$array""$cn"_"$cde"_HE.jpg)" ../output/rotated_images

# Make smaller rotated img for QC plotting
mogrify -path ./ -resize 5% -quality 100%% "$(echo "$array""$cn"_"$cde"_HE.jpg)"

# rm downscalled un-cropped and cropped images
#rm "$(echo $img)"
rm "$(echo Rplots.pdf)"
rm "$(echo "$array""$cn"_"$cde"_HE_crop.jpg)"
