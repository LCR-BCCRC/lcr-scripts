# Augment SSM/MAF

These scripts allow for searching an index .bam for read support covering variants found in additional samples.

## Versions: 
### 1.0
-is the original version, designed for integration directly into another snakemake file.

### 2.0
-contains a streamlined version of 1.0, and code has been wrapped in a custom class, and so can be imported into other scripts
or called via the cli that has also been added. This should make it more flexible and easier to integrate into various workflows. 

Additionally, an argument was added to set a min alt read limit for variants being added. 

It also adds a new column called "variant_source" which will gain a value of "additional_maf" if this script adds a variant.

### Example usage

example CLI usage
```
python augmentMAF.py --sample_id sample1 --threads 24 --index_maf index.maf --index_bam index.bam --add_maf_files add1.maf add2.maf --genome_build hg38 --alt_count_min 3 --output augmented.maf

```
example importing it into your own python script

```
from augmentMAF import AugmentMAF as AM

AM_object = AM(sample_id="sample1", index_maf="index.maf", index_bam="index.bam",
                        add_maf_files=["maf1.maf", "maf2.maf"], genome_build="hg38",
                        output="output.maf", threads=6, min_alt_count=3)

```