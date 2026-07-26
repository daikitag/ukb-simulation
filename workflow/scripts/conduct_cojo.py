import os
import subprocess
import sys
import tempfile

import pandas as pd


def main():
    sys.stderr = open(snakemake.log[0], "w", buffering=1)

    gwas_result = pd.read_csv(snakemake.input.gwas_result, sep="\t")
    gwas_result = gwas_result[gwas_result["TEST"] == "ADD"]

    cojo_df = pd.DataFrame(
        {
            "SNP": gwas_result["ID"],
            "A1": gwas_result["ALT"],
            "A2": gwas_result["REF"],
            "freq": gwas_result["A1_FREQ"],
            "b": gwas_result["BETA"],
            "se": gwas_result["SE"],
            "p": gwas_result["P"],
            "N": gwas_result["OBS_CT"],
        }
    )
    cojo_df.to_csv(snakemake.output.cojo_file, sep="\t", index=False)

    with tempfile.TemporaryDirectory() as temp_dir:
        subset_plink_file_name = os.path.join(temp_dir, "subset_gwas")

        plink_command = [
            "plink2",
            "--bed",
            snakemake.input.bed,
            "--bim",
            snakemake.input.bim,
            "--fam",
            snakemake.input.fam,
            "--keep",
            snakemake.input.individual_id,
            "--make-bed",
            "--out",
            subset_plink_file_name,
        ]

        subprocess.run(plink_command, check=True, stdout=sys.stderr, stderr=sys.stderr)

        subprocess.run(
            [
                snakemake.params.manc_cojo,
                "--bfile",
                subset_plink_file_name,
                "--cojo-file",
                snakemake.output.cojo_file,
                "--out",
                snakemake.output.cojo_result.removesuffix(".sumstat.jma.cojo"),
                "--cojo-p",
                "1e-6",
                "--cojo-slct",
            ],
            check=True,
            stdout=sys.stderr,
            stderr=sys.stderr,
        )


if __name__ == "__main__":
    main()
