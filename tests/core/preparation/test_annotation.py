"""
Tests cwas.core.prepartion.annotation
"""
from __future__ import annotations

import random

import cwas.core.preparation.annotation as annotation
import pytest
from cwas.utils.cmd import compress_bed_file, index_bed_file

NUM_CHROM = random.randint(1, 22)


def create_bed_file(bed_path_txt: pathlib.Path, coordinates: list[tuple]):
    with bed_path_txt.open("w") as bed_file:
        for coordinate in coordinates:
            print(*coordinate, sep="\t", file=bed_file)

    bed_path = compress_bed_file(bed_path_txt)
    _ = index_bed_file(bed_path)

    return bed_path


@pytest.fixture(scope="module")
def input_coordinates():
    """Make input coordinates without a contig (chromosome)
    """
    input_coordinates = []

    # Number of appends = Number of BED files
    input_coordinates.append(
        [(100, 300, 1), (1000, 1100, 1),]
    )
    input_coordinates.append(
        [(200, 400, 1), (1150, 1200, 1),]
    )
    input_coordinates.append(
        [(250, 350, 1), (1250, 1350, 1),]
    )

    return input_coordinates


@pytest.fixture(scope="module")
def output_coordinates():
    output_coordinates = [
        (100, 200, 1),
        (200, 250, 3),
        (250, 300, 7),
        (300, 350, 6),
        (350, 400, 2),
        (1000, 1100, 1),
        (1150, 1200, 2),
        (1250, 1350, 4),
    ]
    return output_coordinates


@pytest.fixture(scope="module")
def bed_file_paths(cwas_workspace, input_coordinates):
    bed_file_path_txts = [
        cwas_workspace / f"test{i + 1}.bed"
        for i in range(len(input_coordinates))
    ]
    chroms = [f"chr{i + 1}" for i in range(NUM_CHROM)]
    bed_file_paths = []

    for bed_file_path_txt, coordinates_per_file in zip(
        bed_file_path_txts, input_coordinates
    ):
        bed_coordinates = []

        for chrom in chroms:
            for coordinate in coordinates_per_file:
                bed_coordinates.append((chrom, *coordinate))

        bed_file_path = create_bed_file(bed_file_path_txt, bed_coordinates)
        bed_file_paths.append(bed_file_path)

    yield bed_file_paths

    for bed_file_path in bed_file_paths:
        bed_file_path.unlink()


def test_merge_bed_file_by_chrom(
    cwas_workspace, bed_file_paths, output_coordinates
):
    chrom = f"chr{random.randint(1, NUM_CHROM)}"
    result_bed_path = cwas_workspace / "result.bed"
    expected = [(chrom, *coordinate) for coordinate in output_coordinates]
    annotation.merge_bed_files_by_chrom(result_bed_path, chrom, bed_file_paths)

    with result_bed_path.open("r") as result_bed_file:
        for i, line in enumerate(result_bed_file):
            expected_line = "\t".join([str(field) for field in expected[i]])
            assert expected_line == line.strip()

    result_bed_path.unlink()


def test_merge_bed_file_by_chrom_already_exists(cwas_workspace, bed_file_paths):
    chrom = f"chr{random.randint(1, NUM_CHROM)}"
    result_bed_path = cwas_workspace / "result_empty.bed"
    result_bed_path.touch()

    # Expect the merging process will not run.
    annotation.merge_bed_files_by_chrom(result_bed_path, chrom, bed_file_paths)

    with result_bed_path.open("r") as result_bed_file:
        assert not result_bed_file.read()  # Expect this file is empty.

    result_bed_path.unlink()


def test_merge_bed_file_by_chrom_force_overwrite(
    cwas_workspace, bed_file_paths, output_coordinates
):
    chrom = f"chr{random.randint(1, NUM_CHROM)}"
    result_bed_path = cwas_workspace / "result_overwritten.bed"
    result_bed_path.touch()

    with result_bed_path.open("r") as result_bed_file:
        assert not result_bed_file.read()  # Expect this file is empty.

    expected = [(chrom, *coordinate) for coordinate in output_coordinates]
    annotation.merge_bed_files_by_chrom(
        result_bed_path, chrom, bed_file_paths, 1
    )

    with result_bed_path.open("r") as result_bed_file:
        for i, line in enumerate(result_bed_file):
            expected_line = "\t".join([str(field) for field in expected[i]])
            assert expected_line == line.strip()

    result_bed_path.unlink()


def test_merge_bed_files(cwas_workspace, bed_file_paths, output_coordinates):
    result_bed_path = cwas_workspace / "result.bed"
    bed_file_and_keys = [
        (bed_file_path, f"BED{i}")
        for i, bed_file_path in enumerate(bed_file_paths)
    ]
    expected_annot_key = "|".join(
        [f"BED{i}" for i in range(len(bed_file_paths))]
    )
    expected = [
        (f"#ANNOT={expected_annot_key}",),
        ("#chrom", "start", "end", "annot_int"),
    ]

    for i in range(3):
        chrom = f"chrom{i + 1}"
        for coordinate in output_coordinates:
            expected.append((chrom, *coordinate))

    annotation.merge_bed_files(result_bed_path, bed_file_and_keys)
    with result_bed_path.open("r") as result_bed_file:
        for i, line in enumerate(result_bed_file):
            expected_line = "\t".join([str(field) for field in expected[i]])
            assert expected_line == line.strip()

    result_bed_path.unlink()


def test_merge_bed_files_multiprocessing(
    cwas_workspace, bed_file_paths, output_coordinates
):
    result_bed_path = cwas_workspace / "result.bed"
    bed_file_and_keys = [
        (bed_file_path, f"BED{i}")
        for i, bed_file_path in enumerate(bed_file_paths)
    ]
    expected_annot_key = "|".join(
        [f"BED{i}" for i in range(len(bed_file_paths))]
    )
    expected = [
        (f"#ANNOT={expected_annot_key}",),
        ("#chrom", "start", "end", "annot_int"),
    ]

    for i in range(NUM_CHROM):
        chrom = f"chrom{i + 1}"
        for coordinate in output_coordinates:
            expected.append((chrom, *coordinate))

    num_proc = random.randint(1, NUM_CHROM)
    annotation.merge_bed_files(result_bed_path, bed_file_and_keys, num_proc)
    with result_bed_path.open("r") as result_bed_file:
        for i, line in enumerate(result_bed_file):
            expected_line = "\t".join([str(field) for field in expected[i]])
            assert expected_line == line.strip()

    result_bed_path.unlink()
