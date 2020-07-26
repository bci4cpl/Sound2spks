#!/bin/bash
# removes wav and associated mat file, should be excuted from /datasets/spiking/Sound2spks/wav_data

set -e
# Any subsequent(*) commands which fail will cause the shell script to exit immediately

input=$1

#filename=$(basename -- "$fullfile")
#extension="${filename##*.}"
#filename="${filename%.*}"

filename="${input%%.*}"
#extension="${input#*.}"

#echo $filename
#echo $extension

WAV_FILE=$filename.wav
MAT_FILE=$filename.mat

if [ ! -f $WAV_FILE ]; then
    echo "ERROR: '$WAV_FILE' not found!"
    exit 1
fi

rm $WAV_FILE
rm ../result_mats/$MAT_FILE