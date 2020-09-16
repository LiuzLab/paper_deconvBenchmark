#!/bin/bash

cibersort_path="/path/to/CIBERSORT"
input_path=$1 
output_path=$2
mix_prefix=$3
ref_prefix=$4
nComp=6


# initialize ref file, this would help to keep the output file in the same format

for((i=1;i<=nComp;i++))
do
java -jar -Xmx3g -Xms3g -jar "$cibersort_path"/CIBERSORT.jar -M "$input_path"/"$mix_prefix""_C$i.txt" -P "$input_path"/"$ref_prefix""_C$i.txt" -c "$input_path"/"sim3_ref_anno_cibersort_C$i.txt" > "$output_path"/"Results_""$mix_prefix""_C$i.txt"
done


# formally run the code 
for((i=1;i<=nComp;i++))
do
java -jar -Xmx3g -Xms3g -jar "$cibersort_path"/CIBERSORT.jar -M "$input_path"/"$mix_prefix""_C$i.txt" -P "$input_path"/"$ref_prefix""_C$i.txt" -c "$input_path"/"sim3_ref_anno_cibersort_C$i.txt" > "$output_path"/"Results_""$mix_prefix""_C$i.txt"
done
