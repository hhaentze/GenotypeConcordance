# GenotypeConcordance
Compare a VCF file with a given ground truth and calculate genotype concordane, recall, precision and F1-score.
True positives are calculated on haplotype level.

# Installation
To install `GenotypeConcordane`,
```sh
git clone https://github.com/hhaentze/GenotypeConcordance.git
cd GenotypeConcordance; make
```

# Usage
Call set and ground truth must have the same format and contain the same variants. 
See bcftools (e.g. sort, norm, isec).

Calculate genotype concordance:

    GenotypeConcordance -c <callset.vcf> -t <ground_truth.vcf>

Calculate genotype concordance and write into file:

    GenotypeConcordance -c <callset.vcf> -t <ground_truth.vcf> -O > results.txt