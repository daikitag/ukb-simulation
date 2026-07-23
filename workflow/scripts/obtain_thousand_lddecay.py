import os
import string
import subprocess
import sys
import tempfile

import pandas as pd

# This script is used to output a subset of VCF
bcftools_script = r"""
bcftools view \
    -S $individual_id \
    -m2 -M2 -v snps \
    -Oz \
    -r $region \
    -o $output \
    $thousand_genomes_vcf
"""


def main():
    sys.stderr = open(snakemake.log[0], "w", buffering=1)

    chromosome = int(snakemake.params.chromosome)
    arm = snakemake.params.arm

    position_df = pd.read_csv(snakemake.input.position_file)

    chromosome_arm = str(chromosome) + arm
    recom_position = position_df[position_df.chromosome == chromosome_arm]

    region = (
        f"chr{chromosome}:{recom_position.left.item()}-{recom_position.right.item()}"
    )

    for pop in ["ceu", "yri", "chb", "jpt"]:
        with tempfile.TemporaryDirectory() as temp_dir:
            vcf_file_name = os.path.join(temp_dir, "analysis.vcf.gz")

            bcftools_command = string.Template(bcftools_script).substitute(
                individual_id=snakemake.params[f"{pop}_founder"],
                region=region,
                thousand_genomes_vcf=snakemake.input.thousand_genomes,
                output=vcf_file_name,
            )
            subprocess.run(bcftools_command, shell=True, check=True)

            lddecay_command = [
                snakemake.params.poplddecay,
                "-InVCF",
                vcf_file_name,
                "-MAF",
                "0.05",
                "-OutStat",
                snakemake.output[f"{pop}_ld_decay"],
            ]

            subprocess.run(lddecay_command, check=True)


if __name__ == "__main__":
    main()
