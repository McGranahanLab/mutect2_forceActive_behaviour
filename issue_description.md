# Mutect2 issue with --force-active

I've been working on detection of tumor somatic variants on mouse WES data (tumor – germline pairs)
and noticed an unexpected behaviour of Mutect2. Originally, I was asked by the collaborators to 
perform a read depth check of the samples and see if saturation was reached and samples' read depth
could be reduced to cut the costs. To answer that question I sequentially downsampled tumor bam to 
25%, 50%, 75% and 90% using samtools. By sequential downsampling I mean that I first downsampled 
full bam to 90% and then I downsampled obtained 90% bam to  83.3% which would result into 75% coverage
from 100% bam. This was done to avoid a situation then some reads, just by chance are present in, for 
example, 50% downsampled bam but are not present in 75% bam. Next, I run Mutect2 on all of my 
downsampled bams. Knowing that Mutect2 can be rather conservative, I used set up arguments 
`–tumor-lod-to-emit` to `0.0` and `--force-active` to `true`. I expected that sets of mutations 
detected in downsampled bams will be nested, i.e. all variants detected in 25% bam will be found in 50%
bam, all variants detected in 50% bam will be found in 75% bam, etc. However, this is not what I saw 
during manual exception of raw Mutect2 calls:

![Image1](https://github.com/McGranahanLab/mutect2_forceActive_behaviour/blob/master/IGVsnapshots/0_rawvcf.png)

if one zoomes in to `chr2:124,659,282-124,662,971`:

![Image2](https://github.com/McGranahanLab/mutect2_forceActive_behaviour/blob/master/IGVsnapshots/0_vcf_zoomI.png)

**There are some variants which are detected in 50% downsampled bam, but not in 100% bam!** This situation can 
also be found in relation to other percentages of downsampling. Consequently I examined mutect2 bams:

![Image3](https://github.com/McGranahanLab/mutect2_forceActive_behaviour/blob/master/IGVsnapshots/2_mutect2bams_overview.png)

So, **variants were not detected because it seems that some of regions which have significant coverage 
in 75% mutect2 bam (or other % of downsampling) are not covered at all at 100% mutect2 bam**. This is 
especially unexpected  taking into account that I set `--force-active` to `true`. Basically, `--force-active` true 
somehow does not mean that all regions will be set to active in contrary to the option description. One may
argue that the issue stems from downsampling and that reads covering those regions just were not included, 
but it's not true as we can see from IGV view of downsampled bams before Mutect2:

![Image4](https://github.com/McGranahanLab/mutect2_forceActive_behaviour/blob/master/IGVsnapshots/1_downsampled_bams.png)

Just to clarify, my **IGV was set not to downsample bams while loading.**

Since I have WES I naturally used `-L` option in Mutect2 and gave it probes locations. I decided to try
to run Mutect2 without this option:

![Image5](https://github.com/McGranahanLab/mutect2_forceActive_behaviour/blob/master/IGVsnapshots/5_mutect2bamsnoProbes_overview.png)

As you can see, the situation somewhat improved, but also not, as regions with coverage at lower percentage 
of downsampling yet with absence of it at higher percentage still exist:

![Image6](https://github.com/McGranahanLab/mutect2_forceActive_behaviour/blob/master/IGVsnapshots/5_mutect2bamsnoProbes_zoomIn.png)

I decided to track the destiny of some reads to double check that they indeed disappear from 
higher % of downsampling  mutect2 bam. I took a read named `A01065:54:HMHGLDRXX:1:1175:8169:3583` which 
maps originally with CIGAR string 98M3S (in input downsampled bams) to `124664309`. This read is present 
in all downsampled bams: 

![Image7](https://github.com/McGranahanLab/mutect2_forceActive_behaviour/blob/master/IGVsnapshots/Screenshot%202021-09-08%20at%2015.24.31.png)

Then I checked if it's present in mutect2 bams:

![Image8](https://github.com/McGranahanLab/mutect2_forceActive_behaviour/blob/master/IGVsnapshots/Screenshot%202021-09-08%20at%2015.24.31.png)

**It was only found in 50% downsampled bam and 75% downsampled bam, but not in the rest!**

This issue is important because if my estimations are correct and if this behaviour is indeed 
erroneous it could be that we only detect about 50% of somatic variants from what we could have detected.

## Short summary of the used software:
* GATK v4.1.8.1
* openjdk version 1.8.0_265
* Preprocessing: reads were trimmed with cutadapt to minimum quality of 20 and minimum length of 20, then mapped with BWA mem, duplicates were removed with sambamba, BaseRecalibrator followed by ApplyBQSR was applied.
