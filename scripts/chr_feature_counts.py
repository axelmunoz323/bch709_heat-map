#!/usr/bin/env python3

from __future__ import annotations

import csv
import gzip
from pathlib import Path

def load_chrom_sizes(path: Path) -> list[tuple[str, int]]:
    rows: list[tuple[str, int]] = []
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for fields in reader:
            if not fields:
                continue
            chrom = fields[0]
            chrom_length_bp = int(fields[1])
            rows.append((chrom, chrom_length_bp))
    return rows


def parse_gff_counts(gff_path: Path, valid_chroms: set[str]) -> tuple[dict[str, dict[str, int]], set[str], int]:
    counts: dict[str, dict[str, int]] = {}
    exon_seen: dict[str, set[tuple[int, int, str]]] = {}

    dropped_seqids: set[str] = set()
    dropped_feature_lines = 0

    with gzip.open(gff_path, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for fields in reader:
            if not fields:
                continue
            if fields[0].startswith("#"):
                continue
            if len(fields) < 9:
                continue

            seqid = fields[0]
            feature_type = fields[2]
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                continue
            strand = fields[6]

            if seqid not in valid_chroms:
                dropped_seqids.add(seqid)
                dropped_feature_lines += 1
                continue

            if seqid not in counts:
                counts[seqid] = {
                    "n_gene": 0,
                    "n_exon_unique": 0,
                    "n_tRNA": 0,
                    "n_snoRNA": 0,
                }
                exon_seen[seqid] = set()

            if feature_type == "gene":
                counts[seqid]["n_gene"] += 1
            elif feature_type == "exon":
                exon_key = (start, end, strand)
                if exon_key not in exon_seen[seqid]:
                    exon_seen[seqid].add(exon_key)
                    counts[seqid]["n_exon_unique"] += 1
            elif feature_type == "tRNA":
                counts[seqid]["n_tRNA"] += 1
            elif feature_type == "snoRNA":
                counts[seqid]["n_snoRNA"] += 1

    return counts, dropped_seqids, dropped_feature_lines


def main() -> None:
    chrom_sizes_path = Path("data/chrom.sizes")
    gff_path = Path("data/saccharomyces_cerevisiae.gff.gz")
    results_dir = Path("results")
    results_dir.mkdir(parents=True, exist_ok=True)

    output_counts_path = results_dir / "chr_feature_counts.tsv"
    output_dropped_path = results_dir / "dropped_seqids.txt"

    chrom_rows = load_chrom_sizes(chrom_sizes_path)
    valid_chroms = {chrom for chrom, _ in chrom_rows}

    counts, dropped_seqids, dropped_feature_lines = parse_gff_counts(gff_path, valid_chroms)

    result_rows: list[dict[str, str | int | float]] = []
    for chrom, chrom_length_bp in chrom_rows:
        row_counts = counts.get(
            chrom,
            {
                "n_gene": 0,
                "n_exon_unique": 0,
                "n_tRNA": 0,
                "n_snoRNA": 0,
            },
        )
        result_rows.append(
            {
                "chrom": chrom,
                "chrom_length_bp": chrom_length_bp,
                "n_gene": row_counts["n_gene"],
                "n_exon_unique": row_counts["n_exon_unique"],
                "n_tRNA": row_counts["n_tRNA"],
                "n_snoRNA": row_counts["n_snoRNA"],
                "gene_per_Mb": round(row_counts["n_gene"] / (chrom_length_bp / 1e6), 4),
                "exon_unique_per_Mb": round(row_counts["n_exon_unique"] / (chrom_length_bp / 1e6), 4),
                "tRNA_per_Mb": round(row_counts["n_tRNA"] / (chrom_length_bp / 1e6), 4),
                "snoRNA_per_Mb": round(row_counts["n_snoRNA"] / (chrom_length_bp / 1e6), 4),
            }
        )

    result_rows.sort(key=lambda row: float(row["gene_per_Mb"]), reverse=True)

    with output_counts_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "chrom",
                "chrom_length_bp",
                "n_gene",
                "n_exon_unique",
                "n_tRNA",
                "n_snoRNA",
                "gene_per_Mb",
                "exon_unique_per_Mb",
                "tRNA_per_Mb",
                "snoRNA_per_Mb",
            ]
        )
        for row in result_rows:
            writer.writerow(
                [
                    row["chrom"],
                    row["chrom_length_bp"],
                    row["n_gene"],
                    row["n_exon_unique"],
                    row["n_tRNA"],
                    row["n_snoRNA"],
                    f"{float(row['gene_per_Mb']):.4f}",
                    f"{float(row['exon_unique_per_Mb']):.4f}",
                    f"{float(row['tRNA_per_Mb']):.4f}",
                    f"{float(row['snoRNA_per_Mb']):.4f}",
                ]
            )

    with output_dropped_path.open("w", encoding="utf-8") as handle:
        for seqid in sorted(dropped_seqids):
            handle.write(f"{seqid}\n")

    print(f"Dropped seqids: {len(dropped_seqids)}")
    print(f"Excluded feature lines: {dropped_feature_lines}")
    print("\nTop 5 rows:")
    headers = [
        "chrom",
        "chrom_length_bp",
        "n_gene",
        "n_exon_unique",
        "n_tRNA",
        "n_snoRNA",
        "gene_per_Mb",
        "exon_unique_per_Mb",
        "tRNA_per_Mb",
        "snoRNA_per_Mb",
    ]
    print("\t".join(headers))
    for row in result_rows[:5]:
        print(
            "\t".join(
                [
                    str(row["chrom"]),
                    str(row["chrom_length_bp"]),
                    str(row["n_gene"]),
                    str(row["n_exon_unique"]),
                    str(row["n_tRNA"]),
                    str(row["n_snoRNA"]),
                    f"{float(row['gene_per_Mb']):.4f}",
                    f"{float(row['exon_unique_per_Mb']):.4f}",
                    f"{float(row['tRNA_per_Mb']):.4f}",
                    f"{float(row['snoRNA_per_Mb']):.4f}",
                ]
            )
        )


if __name__ == "__main__":
    main()
