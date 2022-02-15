"""
Functions to prepare CWAS annotation step
"""
from __future__ import annotations

import multiprocessing as mp
from collections import defaultdict
from functools import partial
from pathlib import Path
from typing import Generator

import cwas.utils.log as log
import numpy as np
from cwas.core.preparation.bedreader import BedReader


def merge_bed_files(
    out_merge_bed: Path,
    bed_file_and_keys: list[tuple[Path, str]],
    num_proc: int = 1,
    force_overwrite: int = 0,
):
    """ Merge all information of coordinates from all the BED files
    into one file.

    For example, assume that following three coordinates exist
    in bed file A, B and C, respectively.
    - chr1  100 300 (from A BED file)
    - chr1  200 400 (from B BED file)
    - chr1  250 350 (from C BED file)

    If these coordinates are merged,
    than following coordinates will be newly produced.
    - chr1  100 200 annotation A
    - chr1  200 250 annotation A, B
    - chr1  250 300 annotation A, B, C
    - chr1  300 350 annotation B, C
    - chr1  350 400 annotation B

    The annotation information for each coordinate will be represented
    as an annotation integer, which is a one-hot encoding of annotation
    information. In this case, A and C are represented as a least significant
    bit and a most significant bit, respectively.

    For example, if all the BED files are A, B, and C,
    then an annotation integer of the second coordinate above is 3(b'011).

    So here is a final result.
    - chr1  100 200 1 (b'001)
    - chr1  200 250 3 (b'011)
    - chr1  250 300 7 (b'111)
    - chr1  300 350 6 (b'110)
    - chr1  350 400 2 (b'010)

    """
    bed_gz_path = Path(str(out_merge_bed) + ".gz")
    if not force_overwrite and (out_merge_bed.exists() or bed_gz_path.exists()):
        log.print_warn("Merged BED file already exists. Skip merging.")
        return

    tmp_dir = out_merge_bed.parent / "tmp"
    tmp_dir.mkdir(exist_ok=True)

    # Make temporary merged BED paths for each chromosome
    chroms = [f"chr{n}" for n in range(1, 23)]
    merge_bed_paths = {}

    for chrom in chroms:
        merge_bed_name = out_merge_bed.name.replace(".bed", f".{chrom}.bed")
        merge_bed_paths[chrom] = tmp_dir / merge_bed_name

    # Make a merged BED file for each chromosome
    bed_paths, bed_ids = zip(*bed_file_and_keys)
    try:
        if num_proc == 1:
            for chrom in chroms:
                merge_bed_files_by_chrom(
                    merge_bed_paths[chrom], chrom, bed_paths, force_overwrite
                )
        else:
            pool = mp.Pool(num_proc)
            pool.starmap(
                partial(
                    merge_bed_files_by_chrom,
                    bed_file_paths=bed_paths,
                    force_overwrite=force_overwrite,
                ),
                [(merge_bed_paths[chrom], chrom) for chrom in chroms],
            )
            pool.close()
            pool.join()
    except Exception:
        log.print_err(
            "Merging BED files has failed because some error occurred."
        )
        # Remove the temporary directory if empty
        if not list(tmp_dir.glob("*")):
            log.print_progress(f'Remove "{tmp_dir}"')
            tmp_dir.rmdir()
        raise

    with out_merge_bed.open("w") as outfile:
        annot_key_str = "|".join(bed_ids)
        print(f"#ANNOT={annot_key_str}", file=outfile)
        print("#chrom", "start", "end", "annot_int", sep="\t", file=outfile)

        for chrom in chroms:
            merge_bed_path = merge_bed_paths[chrom]
            with merge_bed_path.open() as merge_bed_file:
                for line in merge_bed_file:
                    print(line, end="", file=outfile)
            merge_bed_path.unlink()

    tmp_dir.rmdir()


def merge_bed_files_by_chrom(
    out_merge_bed: Path,
    chrom: str,
    bed_file_paths: list[Path],
    force_overwrite: int = 0,
):
    if not force_overwrite and out_merge_bed.exists():
        log.print_warn(
            f'"{out_merge_bed}" already exists. '
            f"Skip merging BED files for {chrom}."
        )
        return
    try:
        log.print_progress(f"Merge BED files for {chrom}")
        _merge_bed_files(out_merge_bed, chrom, bed_file_paths)
    except Exception:
        log.print_err(
            f"Some error occurred during merging BED files for {chrom}."
        )
        if out_merge_bed.exists():
            # Clean the incomplete output
            log.print_progress(f'Remove "{out_merge_bed}"')
            out_merge_bed.unlink()
        raise


def bed_reader_iter(
    bed_file_paths: list[Path], chrom: str
) -> Generator(BedReader):
    """Generate the BedReader instance for each BED file paths"""
    for bed_file_path in bed_file_paths:
        bedreader_inst = BedReader(bed_file_path)
        bedreader_inst.set_contig(chrom)
        yield bedreader_inst


def _merge_bed_files(
    out_merge_bed: Path, chrom: str, bed_file_paths: list[Path],
):
    """ Merge annotation information of all BED coordinates of one chromosome
    from all the BED files """
    start_to_key_idx = defaultdict(list)
    end_to_key_idx = defaultdict(list)
    pos_set = set()
    bed_readers = bed_reader_iter(bed_file_paths, chrom)

    # Read coordinates from each annotation bed file
    for i, bed_reader in enumerate(bed_readers):
        for _, start, end, *_ in bed_reader:
            start_to_key_idx[start].append(i)
            end_to_key_idx[end].append(i)
            pos_set.add(start)
            pos_set.add(end)

    # Create a BED file listing new coordinates
    # with merged annotation information
    pos_list = sorted(pos_set)
    annot_cnt = np.zeros(len(bed_file_paths))
    prev_pos = -1
    n_bed = 0

    with out_merge_bed.open("w") as outfile:
        for pos in pos_list:
            if n_bed > 0:  # Make a new coordinate
                one_hot = np.vectorize(lambda x: 1 if x else 0)(annot_cnt)
                annot_int = _one_hot_to_int(one_hot)
                bed_entry = (chrom, prev_pos, pos, annot_int)
                print(*bed_entry, sep="\t", file=outfile)

            end_key_ind = end_to_key_idx.get(pos)
            start_key_ind = start_to_key_idx.get(pos)

            # This position is an end of at least one bed coordinate.
            if end_key_ind is not None:
                for end_key_idx in end_key_ind:
                    annot_cnt[end_key_idx] -= 1
                n_bed -= len(end_key_ind)

            # This position is a start of at least one bed coordinate.
            if start_key_ind is not None:
                for start_key_idx in start_key_ind:
                    annot_cnt[start_key_idx] += 1
                n_bed += len(start_key_ind)

            prev_pos = pos


def _one_hot_to_int(one_hot: np.ndarray) -> int:
    """
    Note: Index 0 is the least significant bit.
    e.g.
    [1, 1, 0, 1] -> 0'b1011 = 11
    """
    n = 0

    for i in range(len(one_hot)):
        if one_hot[i]:
            n += 2 ** i

    return n
