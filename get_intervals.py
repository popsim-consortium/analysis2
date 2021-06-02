#!/usr/bin/env python3
from collections.abc import Iterable
import operator
import pathlib
from typing import List, Tuple

import allel
import numpy as np
import stdpopsim

GFF_URL = (
    "ftp://ftp.ensembl.org/pub/release-104/gff3/homo_sapiens/"
    "Homo_sapiens.GRCh38.104.gff3.gz"
)
GFF_SHA256 = "313ad46bd4af78b45b9f5d8407bbcbd3f87f4be0747060e84b3b5eb931530ec1"
OUTPUT_DIRECTORY = "./intervals/HomSap"
CHROM_IDS = [chrom.id for chrom in stdpopsim.get_species("HomSap").genome.chromosomes]


def merged(
    intervals: Iterable[Tuple[int, int]], *, closed: bool
) -> List[Tuple[int, int]]:
    """
    Merge overlapping and adjacent intervals.

    :param intervals: An iterable of (start, end) coordinates.
    :param bool closed: If True, [start, end] coordinates are closed,
        so [1, 2] and [3, 4] are adjacent intervals and will be merged.
        If False, [start, end) coordinates are half-open,
        so [1, 2) and [3, 4) are not adjacent and will not be merged.
    """

    def iter_merged(
        intervals: Iterable[Tuple[int, int]], *, closed: bool
    ) -> Iterable[Tuple[int, int]]:
        """
        Generate tuples of (start, end) coordinates for merged intervals.
        """
        intervals = sorted(intervals, key=operator.itemgetter(0))
        if len(intervals) == 0:
            return
        start, end = intervals[0]
        for a, b in intervals[1:]:
            assert a <= b
            if a > end + closed:
                # No intersection with the current interval.
                yield start, end
                start, end = a, b
            else:
                # Intersects, or is contiguous with, the current interval.
                end = max(end, b)
        yield start, end

    return list(iter_merged(intervals, closed=closed))


def test_merged():
    rng = np.random.default_rng(1234)

    for closed in (True, False):
        assert merged([], closed=closed) == []
        assert merged([(10, 20), (15, 30)], closed=closed) == [(10, 30)]
        assert merged([(10, 20), (20, 30)], closed=closed) == [(10, 30)]
        assert merged([(10, 20), (22, 30)], closed=closed) == [(10, 20), (22, 30)]
        assert merged([(10, 20), (12, 16)], closed=closed) == [(10, 20)]
        assert merged([(12, 16), (10, 20), (13, 15)], closed=closed) == [(10, 20)]

        # Check merging is idempotent.
        for _ in range(100):
            starts = rng.integers(1, 1000, size=100)
            ends = starts + rng.integers(1, 100, size=len(starts))
            merged_intervals = merged(zip(starts, ends), closed=closed)
            assert merged_intervals == merged(merged_intervals, closed=closed)

    assert merged([(10, 20), (21, 30)], closed=True) == [(10, 30)]
    assert merged([(10, 20), (21, 30)], closed=False) == [(10, 20), (21, 30)]


def gff_recarray_to_stdpopsim_intervals(gff):
    """
    Merge overlapping intervals and convert coordinates. GFF intervals are
    1-based [i,j], but stdpopsim intervals are 0-based [i-1,j).
    """
    intervals = np.array(merged(zip(gff.start, gff.end), closed=True))
    intervals[:, 0] = intervals[:, 0] - 1
    return intervals


def get_gff_recarray(url, sha256):
    local_path = pathlib.Path(url).name

    if not pathlib.Path(local_path).exists():
        print(f"downloading {url}")
        stdpopsim.utils.download(url, local_path)

    print("checking sha256")
    local_sha256 = stdpopsim.utils.sha256(local_path)
    if local_sha256 != sha256:
        print(
            f"{local_path}: sha256: expected {sha256}, but found {local_sha256}. "
            "Delete the file to download it again."
        )
        exit(1)

    print(f"loading {local_path} into numpy recarray")
    gff = allel.gff3_to_recarray(local_path)
    return gff


if __name__ == "__main__":
    gff = get_gff_recarray(GFF_URL, GFF_SHA256)

    print("extracting exons")
    exons = gff[
        np.where(np.logical_and(gff.source == "ensembl_havana", gff.type == "exon"))
    ]

    out_dir = pathlib.Path(OUTPUT_DIRECTORY)
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"merging overlapping regions and dumping to {out_dir}/")

    for chrom_id in CHROM_IDS:
        chrom_exons = exons[np.where(exons.seqid == chrom_id)]
        if len(chrom_exons) == 0:
            continue
        intervals = gff_recarray_to_stdpopsim_intervals(chrom_exons)
        # double check that the intervals can be used in stdpopsim
        stdpopsim.utils.check_intervals_validity(intervals)

        out_file = out_dir / f"ensembl_havana_exons_{chrom_id}.txt"
        np.savetxt(out_file, intervals, fmt="%d")
