#!/bin/bash

cibersort_path="/path/to/CIBERSORT"
input_path=$1 
output_path=$2
mix_prefix=$3
ref_prefix=$4
nMarker=3
nDataset=3
nGrid=10


# initialize ref file, this would help to keep the output file in the same format
for((m=1;m<=nMarker;m++))
do
for((i=1;i<=1;i++))
do
for((j=1;j<=1;j++))
do
java -jar -Xmx3g -Xms3g -jar "$cibersort_path"/CIBERSORT.jar -M "$input_path"/"$mix_prefix""_D$i""P$j.txt" -P "$input_path"/"$ref_prefix""_D$m.txt" -c "$input_path"/"sim1_ref_anno_cibersort_D$m.txt" > "$output_path"/"Results_""$mix_prefix""_D$i""P$j""M$m.txt"
done
done
done

# formally run the code 
for((m=1;m<=nMarker;m++))
do
for((i=1;i<=nDataset;i++))
do
for((j=1;j<=nGrid;j++))
do
java -jar -Xmx3g -Xms3g -jar "$cibersort_path"/CIBERSORT.jar -M "$input_path"/"$mix_prefix""_D$i""P$j.txt" -P "$input_path"/"$ref_prefix""_D$m.txt" -c "$input_path"/"sim1_ref_anno_cibersort_D$m.txt" > "$output_path"/"Results_""$mix_prefix""_D$i""P$j""M$m.txt"
done
done
done