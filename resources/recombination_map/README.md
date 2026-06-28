# Recombination Map

Put recombination map file to this directory. Recombination map should be a text file, and it should include the following columns:

- First column: Chromosome (This is not used in the simulation)
- Second column: Position(bp)
- Third column: Rate(cM/Mb)

The elements should be delimited by tab, '\t'.

We are planning to use the recombination map that can be obtained from the following link, https://stdpopsim.s3.us-west-2.amazonaws.com/genetic_maps/HomSap/HapMapII_GRCh38.tar.gz, where this is taken from stdpopsim GitHub repository, https://github.com/popsim-consortium/stdpopsim/blob/0d746bd58753af6b410ffe78e88c56520a3393da/stdpopsim/catalog/HomSap/genetic_maps.py#L56.

The shell script to download this recombination map is in `obtain_recombination_map.sh`.

Explanation (taken from stdpopsim's GitHub repository):

This genetic map is from the Phase II Hapmap project and based on 3.1 million genotyped SNPs from 270 individuals across four populations (YRI, CEU, CHB and JPT). Genome wide recombination rates were estimated using LDHat. This version is lifted over to GRCh38 using liftover from the HapMap Phase II map previously lifted over to GRCh37. Liftover was performed using the liftOver_catalog.py script from stdpopsim/maintainance. Exact command used is as follows:

`python <path_to_stdpopsim>/stdpopsim/maintenance/liftOver_catalog.py --species HomSap --map HapMapII_GRCh37 --chainFile <path_to_chainfiles>/chainfiles/hg19ToHg38.over.chain.gz --validationChain <path_to_chainfiles>/chainfiles/hg38ToHg19.over.chain.gz --winLen 1000 --useAdjacentAvg --retainIntermediates --gapThresh 1000000`