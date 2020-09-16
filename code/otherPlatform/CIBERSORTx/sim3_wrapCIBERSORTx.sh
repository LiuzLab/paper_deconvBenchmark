#!/bin/bash


input_path=$1 
output_path=$2
mix_prefix=$3
ref_prefix=$4
nComp=6

for((i=1;i<=nComp;i++))
do
docker run -v "$input_path":/src/data -v "$output_path""/$mix_prefix""_C$i":/src/outdir cibersortx/fractions --username your_username --token your_token --refsample "$ref_prefix""_C$i.txt" --phenoclasses "sim3_ref_anno_cibersort_C$i.txt" --mixture "$mix_prefix""_C$i.txt"
cp "$output_path""/$mix_prefix""_C$i"/CIBERSORTx_Results.txt "$output_path"/"Results_""$mix_prefix""_C$i.txt"
done
