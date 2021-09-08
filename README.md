# Mutect2 --force-active has an issue

This directory contains input files and code needed to reproduce GATK4-MUTECT2
issue described in here. The description is also available here.

## Folder content description
1. `S6.bam` bam file containing **tumor** reads mapped to **GRCm38**. To decrease
demonstration run time of Mutect2 only reads mapped to 
**chr2:124652644-124669380** were kept.
2. `S6_GL.bam` bam file containing **germline** reads mapped to **GRCm38**. To decrease
demonstration run time of Mutect2 only reads mapped to 
**chr2:124652644-124669380** were kept.
3. `1_downsample_consequent.sh` bash script to perform downsampling of `S6.bam`
to 25%, 50%, 75% and 90%. **For the sake of reproducibility results of this 
script are deposited in downsampled_bams folder** because downsampling can be 
stochastic.
3. `2_mutect2_downsample.sh` bash script which runs Mutect2 on all downsampled bams.

## What is not included
Please ensure that you have `samtools` (preferably v1.9) and `singularity` available.

## How to reproduce the issue:
1. Clone this repository
```
git clone https://github.com/McGranahanLab/mutect2_forceActive_behaviour
```
2. Download reference genome files and unzip them:
```
wget https://www.dropbox.com/s/f8s9sxdsd1msam6/reference_genome.zip
unzip reference_genome.zip
```
3. Pull singularity image containing GATK:
```
singularity pull library://marialitovchenko/default/gatk3_and_4:latest
```
4. Run `2_mutect2_downsample.sh`
```
chmod +x 2_mutect2_downsample.sh
./2_mutect2_downsample.sh
```