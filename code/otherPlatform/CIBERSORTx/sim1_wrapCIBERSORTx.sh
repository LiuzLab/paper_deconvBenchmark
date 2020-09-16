#!/bin/bash

input_path=$1 
output_path=$2
mix_prefix=$3
ref_prefix=$4
nMarker=3
nDataset=3
nGrid=10


for((m=1;m<=nMarker;m++))
do
for((i=1;i<=nDataset;i++))
do
for((j=1;j<=nGrid;j++))
do
docker run -v "$input_path":/src/data -v "$output_path""/$mix_prefix""_D$i""P$j""M$m":/src/outdir cibersortx/fractions --username your_username --token your_token --refsample "$ref_prefix""_D$m.txt" --phenoclasses "sim1_ref_anno_cibersort_D$m.txt" --mixture "$mix_prefix""_D$i""P$j.txt"
cp "$output_path""/$mix_prefix""_D$i""P$j""M$m"/CIBERSORTx_Results.txt "$output_path"/"Results_""$mix_prefix""_D$i""P$j""M$m.txt"
done
done
done

