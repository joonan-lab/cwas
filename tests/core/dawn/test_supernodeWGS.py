"""
Tests for cwas.core.dawn.supernodeWGS module
"""
import os

import numpy as np
import pandas as pd
import pytest

from cwas.core.dawn.supernodeWGS import data_collection, supernodeWGS_func


# --- supernodeWGS_func tests ---

@pytest.fixture
def sample_corr_mat():
    n = 10
    rng = np.random.default_rng(42)
    mat = np.eye(n) + rng.normal(0, 0.1, (n, n))
    mat = (mat + mat.T) / 2
    return mat


@pytest.fixture
def sample_fit_res():
    return pd.DataFrame({
        "annotation": [f"cat_{i}" for i in range(10)],
        "cluster": [1, 1, 1, 2, 2, 2, 3, 3, 3, 3],
    })


@pytest.fixture
def supernode_inst(sample_corr_mat, sample_fit_res, tmp_path):
    return supernodeWGS_func(
        corr_mat=sample_corr_mat,
        fit_res=sample_fit_res,
        max_cluster=3,
        cores=1,
        output_dir_path=str(tmp_path),
        tag="test",
        seed=42,
    )


def test_supernodeWGS_init(supernode_inst):
    assert supernode_inst.max_cluster == 3
    assert supernode_inst.cores == 1
    assert supernode_inst.tag == "test"
    assert supernode_inst.seed == 42
    assert supernode_inst.verbose is True


def test_index_to_pair(supernode_inst):
    assert supernode_inst._index_to_pair(1, max_val=3) == (1, 1)
    assert supernode_inst._index_to_pair(2, max_val=3) == (1, 2)
    assert supernode_inst._index_to_pair(3, max_val=3) == (1, 3)
    assert supernode_inst._index_to_pair(4, max_val=3) == (2, 2)
    assert supernode_inst._index_to_pair(5, max_val=3) == (2, 3)
    assert supernode_inst._index_to_pair(6, max_val=3) == (3, 3)


def test_index_to_pair_invalid(supernode_inst):
    with pytest.raises(AssertionError):
        supernode_inst._index_to_pair(0, max_val=3)
    with pytest.raises(AssertionError):
        supernode_inst._index_to_pair(7, max_val=3)


def test_matching(supernode_inst):
    name1 = np.array(["a", "b", "c"])
    name2 = np.array(["c", "a", "b"])
    result = supernode_inst._matching_(name1, name2)
    assert result[0] == 2  # "a" is at index 2 in name2 (1-based)
    assert result[1] == 3  # "b" is at index 3
    assert result[2] == 1  # "c" is at index 1


def test_matching_safety_assert(supernode_inst):
    name1 = np.array(["a", "d"])
    name2 = np.array(["a", "b", "c"])
    with pytest.raises(AssertionError):
        supernode_inst._matching_(name1, name2, safety=True)


def test_cluster_size(supernode_inst):
    vec = np.array([0.5, 0.3, np.nan, 0.1])
    clustering = [1, 1, 2, 2]
    flag_vec = np.array([False, False, False, False])
    result = supernode_inst._cluster_size_(vec, clustering, flag_vec)
    assert result[0] == 2  # cluster 1: 2 valid
    assert result[1] == 1  # cluster 2: 1 valid (nan excluded)


def test_cluster_size_with_flags(supernode_inst):
    vec = np.array([0.5, 0.3, 0.1, 0.2])
    clustering = [1, 1, 2, 2]
    flag_vec = np.array([True, False, False, True])
    result = supernode_inst._cluster_size_(vec, clustering, flag_vec)
    assert result[0] == 1  # cluster 1: 1 valid (first flagged)
    assert result[1] == 1  # cluster 2: 1 valid (last flagged)


def test_clusters_property(supernode_inst):
    clusters = supernode_inst.clusters
    assert clusters == [1, 1, 1, 2, 2, 2, 3, 3, 3, 3]


def test_supernode_dir_creation(supernode_inst):
    sdir = supernode_inst.supernodeDir
    assert os.path.exists(sdir)
    assert "supernodeWGS_results" in sdir
    assert "blocks_test" in sdir


def test_report_results(supernode_inst):
    vec = ["cat_a", "cat_b", "cat_c"]
    posterior = [0.9, 0.5, 0.1]
    pvalue = [0.01, 0.05, 0.5]
    iupdate = [1, 1, 0]
    result = supernode_inst.report_results(vec, posterior, pvalue, iupdate)
    assert isinstance(result, pd.DataFrame)
    assert list(result.columns) == ["Name", "p.value", "FDR", "indicator"]
    assert len(result) == 3


def test_report_results_length_mismatch(supernode_inst):
    with pytest.raises(AssertionError):
        supernode_inst.report_results(["a"], [0.9, 0.5], [0.01], [1])


def test_value_to_color(supernode_inst):
    cmap = np.array([[1, 0, 0, 1], [0, 1, 0, 1], [0, 0, 1, 1]])
    color = supernode_inst._value_to_color(cmap, maxz=10, minz=0, x=5)
    assert len(color) == 4  # RGBA


def test_node_color(supernode_inst):
    color = supernode_inst._node_color_(maxz=10, minz=0, x=5)
    assert isinstance(color, str)
    assert color.startswith("#")


def test_term_freq(supernode_inst):
    x = ["A_B_C", "A_B_D", "A_E_F"]
    result = supernode_inst._term_freq(x)
    assert "A" in result  # "A" appears 3 times


# --- data_collection tests ---

@pytest.fixture
def data_coll(tmp_path):
    return data_collection(
        path=str(tmp_path),
        cores=1,
        max_cluster=3,
        seed=42,
    )


def test_data_collection_init(data_coll):
    assert data_coll.cores == 1
    assert data_coll.max_cluster == 3
    assert data_coll.seed == 42
    assert data_coll.verbose is True


def test_pair_to_index(data_coll):
    assert data_coll._pair_to_index_(1, 1, max_val=3) == 1
    assert data_coll._pair_to_index_(1, 2, max_val=3) == 2
    assert data_coll._pair_to_index_(1, 3, max_val=3) == 3
    assert data_coll._pair_to_index_(2, 2, max_val=3) == 4
    assert data_coll._pair_to_index_(2, 3, max_val=3) == 5
    assert data_coll._pair_to_index_(3, 3, max_val=3) == 6


def test_pair_to_index_symmetry(data_coll):
    assert data_coll._pair_to_index_(1, 3, max_val=5) == data_coll._pair_to_index_(3, 1, max_val=5)
    assert data_coll._pair_to_index_(2, 4, max_val=5) == data_coll._pair_to_index_(4, 2, max_val=5)


def test_pair_to_index_invalid(data_coll):
    with pytest.raises(AssertionError):
        data_coll._pair_to_index_(0, 1, max_val=3)
    with pytest.raises(AssertionError):
        data_coll._pair_to_index_(4, 1, max_val=3)


def test_soft(data_coll):
    x = np.array([3.0, -2.0, 0.5, -0.5])
    d = 1.0
    result = data_coll.soft(x, d)
    np.testing.assert_array_equal(result, [2.0, -1.0, 0.0, 0.0])


def test_l2n(data_coll):
    vec = np.array([3.0, 4.0])
    assert data_coll.l2n(vec) == 5.0


def test_l2n_zero_vector(data_coll):
    vec = np.array([0.0, 0.0])
    assert data_coll.l2n(vec) == 0.05  # fallback value


def test_determine_sign_positive(data_coll):
    vec = np.array([1.0, 2.0, 3.0, -0.5])
    assert data_coll._determine_sign_(vec) == 1


def test_determine_sign_negative(data_coll):
    vec = np.array([-1.0, -2.0, -3.0, 0.5])
    assert data_coll._determine_sign_(vec) == -1


def test_safesvd(data_coll):
    x = np.array([[1, 2], [3, 4], [5, 6]], dtype=float)
    U, s, Vt = data_coll.safesvd(x)
    assert U.shape[0] == 3
    assert len(s) == 2


def test_binary_search_zero(data_coll):
    argu = np.array([0.0, 0.0, 0.0])
    result = data_coll.BinarySearch(argu, sumabs=1.0)
    assert result == 0


def test_binary_search_already_satisfies(data_coll):
    argu = np.array([0.5, 0.3])
    result = data_coll.BinarySearch(argu, sumabs=10.0)
    assert result == 0


def test_extract_eigenvector(data_coll):
    mat = np.array([[2.0, 1.0], [1.0, 2.0]])
    eig = data_coll._extract_eigenvector_(mat, sparse=False)
    assert len(eig) == 2
    assert not np.any(np.isnan(eig))
