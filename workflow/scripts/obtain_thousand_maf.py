import string
import subprocess
import sys

import pandas as pd

bcftools_script = r"""
bcftools view \
    -S $individual_id \
    -m2 -M2 -v snps \
    -Ou \
    -r $region \
    $thousand_genomes_vcf \
| bcftools +fill-tags -Ou -- -t AF \
| bcftools query -f '%CHROM\t%POS\t%INFO/AF\n' \
> $output
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
        bcftools_command = string.Template(bcftools_script).substitute(
            individual_id=snakemake.params[f"{pop}_founder"],
            region=region,
            thousand_genomes_vcf=snakemake.input.thousand_genomes,
            output=snakemake.output[f"{pop}_maf"],
        )
        subprocess.run(bcftools_command, shell=True, check=True)


if __name__ == "__main__":
    main()
