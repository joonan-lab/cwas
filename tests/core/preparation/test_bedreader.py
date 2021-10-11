"""
Test the 'BedReader' class in the cwas.core.preparation.bedreader module
"""
import pytest
from cwas.core.preparation.bedreader import BedReader
from cwas.core.preparation.utils import compress_bed_file, index_bed_file


@pytest.fixture(scope="module")
def bed_coordinates():
    coordinates = [
        ("chr1", 1, 2),
        ("chr2", 20, 30),
        ("chr3", 100, 200),
    ]
    return coordinates


@pytest.fixture(scope="module")
def tmp_bed_path(cwas_workspace, bed_coordinates):
    txt_bed_path = cwas_workspace / ".test.bed"
    with txt_bed_path.open("w") as bed_file:
        for bed_coordinate in bed_coordinates:
            print(*bed_coordinate, sep="\t", file=bed_file)

    bed_path = compress_bed_file(txt_bed_path)
    _ = index_bed_file(bed_path)

    return bed_path


def test_init_bedreader(tmp_bed_path):
    bedreader = BedReader(tmp_bed_path)
    assert isinstance(bedreader, BedReader)


def test_fail_init_bedreader(cwas_workspace):
    invalid_bed_path = cwas_workspace / ".not_exist.bed"
    with pytest.raises(OSError):
        BedReader(invalid_bed_path)


def test_read_coordinates(tmp_bed_path, bed_coordinates):
    bedreader = BedReader(tmp_bed_path)
    expected = bed_coordinates
    result = list(bedreader)
    assert expected == result


def test_read_coordinates_with_set_contig(tmp_bed_path, bed_coordinates):
    contig = "chr1"
    bedreader = BedReader(tmp_bed_path)
    bedreader.set_contig(contig)
    expected = [
        bed_coordinate
        for bed_coordinate in bed_coordinates
        if bed_coordinate[0] == contig
    ]
    result = list(bedreader)
    assert expected == result
