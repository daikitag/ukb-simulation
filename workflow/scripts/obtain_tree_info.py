import os
import sys

import tszip


def main():
    sys.stderr = open(snakemake.log[0], "w", buffering=1)

    ts = tszip.load(snakemake.input.tsz)

    file_size_bytes = os.path.getsize(snakemake.input.tsz)
    file_size_gb = file_size_bytes / (1024**3)

    with open(snakemake.output.tree_info, "w") as f:
        print(ts, file=f)
        print(f"\nTszip File Size (GB): {file_size_gb}", file=f)


if __name__ == "__main__":
    main()
