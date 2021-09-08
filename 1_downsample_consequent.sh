#!/bin/bash

# -----------------------------------------------------------------------------
# DESCRIPTION
# Script to perform CONSEQUITIVE downsampling of bam files to 25%, 50%, 75%, 
# 90% of total reads. Consequitive downsampling means that to create, for 
# example, 75% downsampled bam, we will downsample 90% bam, rather than 
# original 100% bam. This techniques removes a factor of random downsampling 
# and therefore presence/absence of some variants due to change of sampling or
# not sampling alternative allele. Downsampling is performed with samtools.
#
# HOW TO RUN: chmod +x 1_downsample_consequent.sh
#             ./1_downsample_consequent.sh
#
# INPUT FILE: S6.bam
# INPUT VALUES: FRACT - fractions (percentages) to downsample bam to and SEED
#               seed to be used in samtools
#
# OUTPUT FILE(S): folder downsampled_bams will be created. It will contain 
#                 downsampled bam files. Names fill follow pattern:
#                 S6.downsampling_percentage.bam
#
# ATTENTION: MAKE SURE SAMTOOLS v1.9 IS IN YOUR PATH
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# LOAD ALL NESSECARY SOFTWARE
# -----------------------------------------------------------------------------
# MAKE SURE SAMTOOLS v1.9 IS IN YOUR PATH

# -----------------------------------------------------------------------------
# Inputs
# -----------------------------------------------------------------------------
BASE_DIR=`pwd`

FULL_BAM='S6.bam'

# fractions we want to downsample to. Always in descending order! Always 
# include 100 as first element
FRACT=( 100 90 75 50 25 )
# seed to be used in samtools
SEED=42

# -----------------------------------------------------------------------------
# Outputs
# -----------------------------------------------------------------------------
DOWNSAMP_DIR=$BASE_DIR"/downsampled_bams/"
mkdir -p $DOWNSAMP_DIR

# -----------------------------------------------------------------------------
# Perform CONSEQUITIVE downsampling
# -----------------------------------------------------------------------------
sample=`basename $FULL_BAM`
sample=`echo $sample | sed 's/.bam//g'`

for i in "${!FRACT[@]}"; do 
	if [ "$i" -gt 0 ]; then
		fract="${FRACT[$i]}"

		# calculate updated fraction: a % which is needed to take from bam
		# from the previous step in order to achieve original desired fraction
		prevIdx=$(expr $i - 1)
		prevFract="${FRACT[$prevIdx]}"
		updFract=$(expr "${FRACT[$i]}"  \* 100 / $prevFract)
		updFract=$(bc -l <<<$updFract)

		# get bam to downsample from: either original, or 
		if [ "$i" -eq 1 ]; then
			bamToDowns=$FULL_BAM
		else
			bamToDowns=$DOWNSAMP_DIR$sample"."$prevFract".bam"
		fi

		downsampledBAM=$DOWNSAMP_DIR$sample"."$fract".bam"

		CURR_TIME=`date`
		echo '[' $CURR_TIME '] Started downsampling to '$fract' from '$bamToDowns
		echo '[' $CURR_TIME ']     Updated fraction: '$updFract
		samtools view -bs $SEED"."$updFract $bamToDowns > $downsampledBAM
		samtools index $downsampledBAM
		CURR_TIME=`date`
		echo '[' $CURR_TIME '] Finished downsampling to '$fract' from '$bamToDowns
	fi
done

# copy full bam to downsampling folder to have everything in one place
cp $FULL_BAM $DOWNSAMP_DIR$sample'.100.bam'
samtools index $DOWNSAMP_DIR$sample'.100.bam'