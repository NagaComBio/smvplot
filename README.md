# smvplot
smvplot is a cmd line python tool to generate IGV-like screenshots.

### Install
Install the package via pip

```
 pip install smvplot
```

### Usage
```
$ smvplot --help
usage: visualize.py [-h] --bam_paths STR --bam_names STR --ref FILE [--exclude_flag INT] [--map_quality INT] [--base_quality INT] [--max_depth_plot INT] [--vaf] [--for_gSmVs] [--vcf FILE] [--bed FILE] [--annotations FILE] [--annotation_names STR] [--prefix PREFIX]
                    [--window N] [--samtoolsbin N] [--tabixbin N] [--plot_dir DIR] [--out_format STR] [--ref_base STR] [--alt_base STR]
                    [region]

This script generates a png file for each entry in a vcf file, a bed file or a manually specified region.

positional arguments:
  region                syntax either 'chr:start-end' or 'chr:center', use --vcf or --bed for more convenience

optional arguments:
  -h, --help            show this help message and exit
  --bam_paths STR       input list of bam files separated by comma. Maximum 3 BAM files
  --bam_names STR       input list of names separated by comma. Same length as BAM files
  --ref FILE            input reference genome file (fastq format)
  --exclude_flag INT    Exclude the reads with corresponding SAM flags, [default = 3840]
  --map_quality INT     Minimum mapping quality for the reads, [default = 20]
  --base_quality INT    Minimum base quality for the variant, [default = 13]
  --max_depth_plot INT  Maximum read depth used to plot the high coverage region, [default = 500]
  --vaf                 Include the VAF of the central position in the plot title. Requires reference genome
  --for_gSmVs           For the gSmVs workflow used internally in DKFZ, the VAFs are directly sourced from the input VCF.
  --vcf FILE            input vcf file ( as an alternative use --bed )
  --bed FILE            input bed file ( as an alternative use --vcf )
  --annotations FILE    Annotation track in bed format is indexed with a tabix. The fourth column could contain the annotation text for the segments. A comma can separate multiple files.
  --annotation_names STR
                        annotation names separated by comma. Same length as annotation files
  --prefix PREFIX       target directory and file name prefix for generated output files, [default = smvplot]
  --window N            the output file for position X will show the region [X-window,X+window], [default = 100]
  --samtoolsbin N       the path to the samtools binary, [default = samtools]
  --tabixbin N          the path to the tabix binary, [default = tabix]
  --plot_dir DIR        subfolder for the plots
  --out_format STR      Output format of the plot, [default = pdf]
  --ref_base STR        Reference base for the variant entry, [default = ]
  --alt_base STR        Alternate base for the variant entry, [default = ]
```

### Example plots
On the GIAB samples

1. For a single variant from a single BAM
```
smvplot 
  --bam_paths HG001_merged.mdup.bam \
  --bam_names HG001 \
  --ref GRCh38_decoy_ebv_phiX_alt_hla_chr.fa \
  --plot_dir ~/smvplot_test \
  --prefix giab_HG001 \
  --out_format png
  chr1:3339544
```
![](examples/giab_HG001_chr1_3339544.png)

2. For a single variant from a TRIO (3 BAMs)

```
smvplot \
  --bam_paths HG002_merged.mdup.bam,HG003_merged.mdup.bam,HG004_merged.mdup.bam \
  --bam_names HG002_Son,HG003_Father,HG004_Mother \
  --ref GRCh38_decoy_ebv_phiX_alt_hla_chr.fa \
  --plot_dir ~/smvplot_test \
  --prefix giab_HG00234 \
  --out_format png \
  chr1:783175
```
![](examples/giab_HG00234_chr1_783175.png)

3. For multiple variants from a VCF/BED file

```
smvplot \
  --bam_paths HG002_merged.mdup.bam,HG003_merged.mdup.bam,HG004_merged.mdup.bam \
  --bam_names HG002_Son,HG003_Father,HG004_Mother \
  --ref GRCh38_decoy_ebv_phiX_alt_hla_chr.fa \
  --plot_dir ~/smvplot_test \
  --prefix giab_HG00234 \
  --out_format png \
  --vcf giab_benchmark_variants.vcf # --bed giab_benchmark_variants.vcf
```
### Changelog

**0.1.0**
- Minor: Read sorting by input ref/alt based via cmd parameter

**0.0.5.2**
- Bug fix in plot_region argument

**0.0.5**
- Patch: Update VAF calculation

**0.0.4.14**
- Container generated via github actions and pushed to dockerhub

**0.0.4**
- Added container

**0.0.3.2**
- Minor bug fixes

**0.0.3**
- Removed the underhand issue in the RNAseq histograms
- Limit the VAF float decimals
- In a multi-BAM settings, ignore the BAMs if the path does not exist

**0.0.2**
- Add VAF to the title via pysamstats

**0.0.1**
- Inital version upload to the PIP 

### Acknowledgements
The `visualize.py` was originally written for the [DKFZ somatic indel workflow](https://github.com/DKFZ-ODCF/IndelCallingWorkflow) by Philip Ginsbach and Ivo Buchhalter. Here I have updated script to a python package and added a possibility of a third BAM and a RNAseq BAM file. And also generalized the BAM inputs.
