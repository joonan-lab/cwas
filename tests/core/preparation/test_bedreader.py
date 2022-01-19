"""
Test the 'BedReader' class in the cwas.core.preparation.bedreader module
"""
from pathlib import Path

import pytest
from cwas.core.preparation.bedreader import BedReader
from cwas.utils.cmd import compress_using_bgzip, index_using_tabix


@pytest.fixture(scope="module")
def bed_coordinates():
    coordinates = [
        ("chr1", 1, 2, "A"),
        ("chr2", 20, 30, "B"),
        ("chr3", 100, 200, "C"),
    ]
    return coordinates


@pytest.fixture(scope="module")
def bed_txt_path(cwas_workspace):
    return cwas_workspace / ".test.bed"


@pytest.fixture(scope="module")
def bed_gz_path(bed_txt_path):
    return Path(str(bed_txt_path) + ".gz")


def create_bed_file(bed_txt_path, bed_coordinates):
    with bed_txt_path.open("w") as bed_file:
        for bed_coordinate in bed_coordinates:
            print(*bed_coordinate, sep="\t", file=bed_file)

    bed_gz_path = compress_using_bgzip(bed_txt_path)
    _ = index_using_tabix(bed_gz_path)


@pytest.fixture(scope="module", autouse=True)
def setup(cwas_workspace, bed_txt_path, bed_coordinates):
    cwas_workspace.mkdir()
    create_bed_file(bed_txt_path, bed_coordinates)


def create_bed_file(bed_txt_path, bed_coordinates):
    with bed_txt_path.open("w") as bed_file:
        for bed_coordinate in bed_coordinates:
            print(*bed_coordinate, sep="\t", file=bed_file)

    bed_gz_path = compress_using_bgzip(bed_txt_path)
    _ = index_using_tabix(bed_gz_path)


@pytest.fixture(scope="module", autouse=True)
def teardown(cwas_workspace):
    yield
    remove_workspace(cwas_workspace)


def remove_workspace(cwas_workspace):
    for f in cwas_workspace.glob("*"):
        f.unlink()
    cwas_workspace.rmdir()


def test_init_bedreader(bed_gz_path):
    bedreader = BedReader(bed_gz_path)
    assert isinstance(bedreader, BedReader)


def test_fail_init_bedreader(cwas_workspace):
    invalid_bed_path = cwas_workspace / ".not_exist.bed"
    with pytest.raises(OSError):
        BedReader(invalid_bed_path)


def test_read_coordinates(bed_gz_path, bed_coordinates):
    bedreader = BedReader(bed_gz_path)
    expected = bed_coordinates
    result = list(bedreader)
    assert expected == result


def test_read_coordinates_with_set_contig(bed_gz_path, bed_coordinates):
    contig = "chr1"
    bedreader = BedReader(bed_gz_path)
    bedreader.set_contig(contig)
    expected = [
        bed_coordinate
        for bed_coordinate in bed_coordinates
        if bed_coordinate[0] == contig
    ]
    result = list(bedreader)
    assert expected == result


def test_read_coordinates_with_invalid_contig(bed_gz_path):
    invalid_contig = "test"
    bedreader = BedReader(bed_gz_path)
    bedreader.set_contig(invalid_contig)
    assert list(bedreader) == []
