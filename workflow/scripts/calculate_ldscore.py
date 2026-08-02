import os
import subprocess
import sys
import tempfile


def main():
    sys.stderr = open(snakemake.log[0], "w", buffering=1)

    with tempfile.TemporaryDirectory() as temp_dir:
        subset_plink_file_name = os.path.join(temp_dir, "subset_ceu")

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
            "--maf",
            str(snakemake.params.maf),
            "--make-bed",
            "--out",
            subset_plink_file_name,
        ]

        subprocess.run(plink_command, check=True, stdout=sys.stderr, stderr=sys.stderr)

        subprocess.run(
            [
                snakemake.params.cojo,
                "--bfile",
                subset_plink_file_name,
                "--ld-score",
                "--ld-wind",
                "1000",
                "--out",
                snakemake.output.ld_score.removesuffix(".score.ld"),
                "--thread-num",
                str(snakemake.threads),
            ],
            check=True,
            stdout=sys.stderr,
            stderr=sys.stderr,
        )


if __name__ == "__main__":
    main()
