"""
Tests for cwas.core.dawn.clustering module

Note: kmeans_cluster requires R via rpy2, so tests are skipped
if rpy2/R is not available.
"""
import numpy as np
import pandas as pd
import pytest

try:
    from cwas.core.dawn.clustering import kmeans_cluster
    HAS_RPY2 = True
except (ImportError, OSError):
    HAS_RPY2 = False


@pytest.fixture
def sample_tsne_data():
    rng = np.random.default_rng(42)
    return pd.DataFrame(
        rng.normal(size=(50, 2)),
        columns=["t-SNE1", "t-SNE2"],
    )


@pytest.mark.skipif(not HAS_RPY2, reason="rpy2/R not available")
class TestKmeansCluster:
    def test_init(self, sample_tsne_data):
        km = kmeans_cluster(sample_tsne_data, seed=42)
        assert km.seed == 42
        assert km.tsne_out is sample_tsne_data

    def test_tsne_out_property(self, sample_tsne_data):
        km = kmeans_cluster(sample_tsne_data, seed=42)
        pd.testing.assert_frame_equal(km.tsne_out, sample_tsne_data)

    def test_seed_property(self, sample_tsne_data):
        km = kmeans_cluster(sample_tsne_data, seed=123)
        assert km.seed == 123

    def test_optimal_k(self, sample_tsne_data, tmp_path):
        km = kmeans_cluster(sample_tsne_data, seed=42)
        output_path = str(tmp_path / "silhouette.png")
        opt_k = km.optimal_k("2,5", output_path)
        assert isinstance(opt_k, (int, np.integer))
        assert 2 <= opt_k <= 5

    def test_center_init(self, sample_tsne_data):
        km = kmeans_cluster(sample_tsne_data, seed=42)
        init_pts = km.center_init(k=3)
        assert len(init_pts) == 3
        for pt in init_pts:
            assert 0 <= pt < len(sample_tsne_data)
