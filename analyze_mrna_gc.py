#!/usr/bin/env python3
import gzip
import math
import os
from dataclasses import dataclass
from typing import Iterator, Tuple, List

import numpy as np
import matplotlib.pyplot as plt


@dataclass
class Record:
    accession: str
    length: int
    gc_content: float


def parse_fasta_gz(path: str) -> Iterator[Tuple[str, str]]:
    header = None
    seq_chunks: List[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if header is not None:
        yield header, "".join(seq_chunks)


def accession_from_header(header: str) -> str:
    return header.split()[0] if header else ""


def gc_content(seq: str) -> float:
    s = seq.upper()
    n = len(s)
    if n == 0:
        return float("nan")
    gc = s.count("G") + s.count("C")
    return gc / n


def gaussian_kde_1d(x: np.ndarray, grid: np.ndarray) -> np.ndarray:
    x = x[np.isfinite(x)]
    n = x.size
    if n == 0:
        return np.zeros_like(grid)

    sd = float(np.std(x, ddof=1)) if n > 1 else 0.0
    bw = 1e-3 if sd == 0.0 else 1.06 * sd * (n ** (-1 / 5))
    bw = max(bw, 1e-4)

    diff = (grid[:, None] - x[None, :]) / bw
    dens = np.exp(-0.5 * diff * diff).sum(axis=1) / (n * bw * math.sqrt(2 * math.pi))
    return dens


def main() -> None:
    in_path = "data/mrna.fa.gz"
    out_tsv = "results/mrna_metrics.tsv"
    out_png = "results/gc_content_distribution.png"
    os.makedirs("results", exist_ok=True)

    records: List[Record] = []
    for header, seq in parse_fasta_gz(in_path):
        acc = accession_from_header(header)
        if not acc:
            continue
        gc = gc_content(seq)
        if not np.isfinite(gc):
            continue
        records.append(Record(acc, len(seq), gc))

    records.sort(key=lambda r: r.gc_content, reverse=True)

    with open(out_tsv, "w", encoding="utf-8") as out:
        out.write("accession\tlength\tgc_content\n")
        for r in records:
            out.write(f"{r.accession}\t{r.length}\t{r.gc_content:.4f}\n")

    gc_vals = np.array([r.gc_content for r in records], dtype=float)
    n = gc_vals.size
    mean = float(np.mean(gc_vals)) if n else float("nan")
    median = float(np.median(gc_vals)) if n else float("nan")
    sd = float(np.std(gc_vals, ddof=1)) if n > 1 else 0.0

    plt.figure(figsize=(8, 4.5), dpi=200)  # 1600x900 px
    plt.hist(gc_vals, bins=50, range=(0, 1), density=True, alpha=0.6, edgecolor="black", linewidth=0.3)

    grid = np.linspace(0, 1, 500)
    dens = gaussian_kde_1d(gc_vals, grid)
    plt.plot(grid, dens, linewidth=1.5)

    plt.axvline(mean, linestyle="--", linewidth=1.2)
    plt.axvline(median, linestyle="--", linewidth=1.2)

    plt.xlim(0, 1)
    plt.xlabel("GC content (0-1)")
    plt.ylabel("Density")
    plt.title("Yeast mRNA GC content distribution")

    caption = f"n={n}, mean={mean:.4f}, median={median:.4f}, sd={sd:.4f}"
    plt.figtext(0.5, 0.01, caption, ha="center", va="bottom", fontsize=9)

    plt.tight_layout(rect=[0, 0.04, 1, 1])
    plt.savefig(out_png, dpi=200)
    plt.close()

    print("Wrote:", out_tsv)
    print("Wrote:", out_png)


if __name__ == "__main__":
    main()
