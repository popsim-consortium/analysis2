"""
Utilities functions for masking
"""
import numpy as np
import pandas as pd
import tskit
import stdpopsim
import operator


def get_mask_from_file(mask_file, chromID):
    """
    Get the mask for a specific chromosome from the specified file.

    :mask_file: path to the mask file
    :chromID: chromosome ID
    """
    mask_table = pd.read_csv(mask_file, sep="\t", header=None)
    # get the mask for the specified chromosome
    mask = mask_table[mask_table[0] == chromID]
    # turn into a numpy array
    mask = np.array(mask.values[:, 1:3])
    return mask


def get_mask_from_chrom_annotation(speciesID, chrom_annotation, chromID):
    """
    Get annotation intervals from the specified species/chrom/annotation.

    :speciesID: species ID (e.g., "HomSap")
    :chrom_annotation: chromosome annotation name (e.g., "ensembl_havana_104_CDS")
    :chromID: chromosome ID (e.g., "chr1")
    """
    species = stdpopsim.get_species(speciesID)
    a = species.get_annotations(chrom_annotation)
    mask = a.get_chromosome_annotations(chromID)
    return mask


def merged(intervals, *, closed: bool):
    """
    Merge overlapping and adjacent intervals.

    :param intervals: An iterable of (start, end) coordinates.
    :param bool closed: If True, [start, end] coordinates are closed,
        so [1, 2] and [3, 4] are adjacent intervals and will be merged.
        If False, [start, end) coordinates are half-open,
        so [1, 2) and [3, 4) are not adjacent and will not be merged.
    """

    def iter_merged(intervals, *, closed: bool):
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


def get_combined_masks(species, mask_file, chromID, chrom_annotation=None):
    """
    Get the combined mask for the specified species/chromosome/mask_file/chrom_annotation

    :species: species ID (e.g., "HomSap")
    :mask_file: path to the mask file
    :chromID: chromosome ID (e.g., "chr1")
    :chrom_annotation: chromosome annotation name (e.g., "ensembl_havana_104_CDS"),
    if set to None, the mask will be obtained from the mask_file only
    """
    if mask_file is None and chrom_annotation is None:
        return None
    # get the mask for the specified chromosome
    if mask_file is not None:
        mask = get_mask_from_file(mask_file, chromID)
    else:
        mask = np.array([], dtype=int).reshape(0, 2)
    if chrom_annotation is not None:
        an_mask = get_mask_from_chrom_annotation(species, chrom_annotation, chromID)
        mask = np.concatenate((mask, an_mask))
    # merge overlapping and adjacent intervals
    mask = merged(mask, closed=True)
    # turn into a numpy array
    mask = np.array(mask)
    return mask 
