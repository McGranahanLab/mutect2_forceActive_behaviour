#!/bin/bash

# -----------------------------------------------------------------------------
# DESCRIPTION
# Script to perform calling of SOMATIC tumor variants on downsampled bams.
#
# HOW TO RUN: chmod +x 2_mutect2_downsample.sh
#             ./2_mutect2_downsample.sh
#
# INPUT FILES: Folder downsampled_bams containing downsampled to various 
#              percentages Bams. Names of the files should follow a pattern:
#              S6.<downsampling_percentage>.bam
# INPUT VALUES: FRACT - fractions (percentages) to downsample bam
#               TUMOR_ID - ID of tumor sample
#               GERMLINE_ID - ID of germline sample
#               refGenDir - path to directory with indexed reference genome
#               refGenName - name of the fasta file containing reference
#                            genome
#               PROBES - Path to Agilent probed Bed file
# 
# IN ORDER TO EMULATE ABCENSE OF -L OPTION: change probes.bed to roi.bed
#
# OUTPUT FILE(S): folder mutect2_on_downsampled or 
#                 mutect2_on_downsampled_noProbe will be created containing
#                 results of mutect2: vcf and bam.
# 
#
# ATTENTION: MAKE SURE SAMTOOLS v1.9 IS IN YOUR PATH
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Inputs
# -----------------------------------------------------------------------------
BASE_DIR=`pwd`

DOWNSAMPLE_BAM_DIR=$BASE_DIR"/downsampled_bams"

# downsample fractions
FRACT=( 25 50 75 90 100 )
# tumor and germline ID's = file names without extention
TUMOR_ID=S6
GERMLINE_ID=S6_GL

# reference genome & probes paths
refGenDir=$BASE_DIR"/reference_genome/"
refGenName=GRCm38_68_chr2_ROI.fa
PROBES=$BASE_DIR"/reference_genome/probes.bed"

# -----------------------------------------------------------------------------
# LOAD ALL NESSECARY SOFTWARE
# -----------------------------------------------------------------------------
gatkSif='gatk3_and_4_latest.sif'

# -----------------------------------------------------------------------------
# Outputs
# -----------------------------------------------------------------------------
if [[ $PROBES == *"probes.bed"* ]]; then
    MUTECT2_DIR=$BASE_DIR"/mutect2_on_downsampled"
    mkdir -p $MUTECT2_DIR
fi

if [[ $PROBES == *"roi.bed"* ]]; then
    MUTECT2_DIR=$BASE_DIR"/mutect2_on_downsampled_noProbe"
    mkdir -p $MUTECT2_DIR
fi

# -----------------------------------------------------------------------------
# Perform calling with MUTECT2
# -----------------------------------------------------------------------------
for i in "${!FRACT[@]}"; do 
    fract="${FRACT[$i]}"
    tumorBAM=$DOWNSAMPLE_BAM_DIR"/"$TUMOR_ID"."$fract".bam"
    glBAM=$BASE_DIR"/"$GERMLINE_ID".bam"

    CURR_TIME=`date`
    echo '[' $CURR_TIME '] Started MUTECT2 on '$tumorBAM
    singularity exec --bind $BASE_DIR:$BASE_DIR $gatkSif gatk Mutect2 \
                     --java-options "-Xmx4g -Xms4g -XX:ParallelGCThreads=1" \
                     -R $refGenDir$refGenName -L $PROBES \
                     -I $tumorBAM -I $glBAM \
                     --tumor-sample $TUMOR_ID --normal-sample $GERMLINE_ID \
                     --tumor-lod-to-emit 0.0 --force-active true \
                     -bamout $MUTECT2_DIR'/'$TUMOR_ID"."$fract".m2.bam" \
                     -O $MUTECT2_DIR'/'$TUMOR_ID"."$fract".m2.vcf" 
    CURR_TIME=`date`
    echo '[' $CURR_TIME '] Finished MUTECT2 on '$tumorBAM

    samtools index $MUTECT2_DIR'/'$TUMOR_ID"."$fract".m2.bam"
    
    # Use FilterMutectCalls, provided in the GATK package, to remove probable 
    # technical or germline artifacts.
    CURR_TIME=`date`
    echo '[' $CURR_TIME '] Started FilterMutectCalls on '$tumorBAM
    singularity exec --bind $BASE_DIR:$BASE_DIR $gatkSif gatk FilterMutectCalls \
                     --java-options "-Xmx4g -Xms4g -XX:ParallelGCThreads=1" \
                     --reference $refGenDir$refGenName \
                     -V $MUTECT2_DIR'/'$TUMOR_ID"."$fract".m2.vcf" \
                     -O $MUTECT2_DIR'/'$TUMOR_ID"."$fract".m2.filt.vcf"
    CURR_TIME=`date`
    echo '[' $CURR_TIME '] Finished FilterMutectCalls on '$tumorBAM

    grep '^#' $MUTECT2_DIR'/'$TUMOR_ID"."$fract".m2.filt.vcf" > $MUTECT2_DIR'/'$TUMOR_ID"."$fract".m2.passFilt.vcf"
    grep -w PASS $MUTECT2_DIR'/'$TUMOR_ID"."$fract".m2.filt.vcf" | grep -v '^#' >> $MUTECT2_DIR'/'$TUMOR_ID"."$fract".m2.passFilt.vcf"

done
