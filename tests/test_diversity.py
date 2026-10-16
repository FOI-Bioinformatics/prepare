"""
Tests for diversity analysis functions.
"""

import numpy as np
import pandas as pd
import pytest
from scipy.spatial.distance import pdist, squareform

from prepare.core.exceptions import ValidationError
from prepare.processing.diversity import (
    _permanova_test,
    calculate_alpha_diversity,
    calculate_beta_diversity,
    perform_permanova,
)


class TestDiversityMetrics:
    """Test diversity metric calculations."""

    def setup_method(self):
        """Set up test data."""
        # Create test abundance matrix (5 samples, 4 species)
        self.abundance_data = np.array(
            [
                [10, 5, 2, 0],  # Sample 1
                [0, 8, 3, 1],  # Sample 2
                [15, 0, 0, 5],  # Sample 3
                [5, 5, 5, 5],  # Sample 4
                [2, 2, 2, 14],  # Sample 5
            ]
        )

        self.abundance_df = pd.DataFrame(
            self.abundance_data,
            index=[f"Sample_{i}" for i in range(1, 6)],
            columns=[f"Species_{i}" for i in range(1, 5)],
        )

        self.metadata = pd.DataFrame(
            {"treatment": ["A", "A", "B", "B", "B"], "site": ["X", "Y", "X", "Y", "X"]}, index=self.abundance_df.index
        )

    def test_calculate_alpha_diversity_default(self):
        """Test alpha diversity calculation with default metrics."""
        result = calculate_alpha_diversity(self.abundance_df)

        # Check structure
        assert isinstance(result, pd.DataFrame)
        assert result.shape[0] == 5  # 5 samples
        assert "shannon" in result.columns
        assert "simpson" in result.columns
        assert "chao1" in result.columns
        assert "observed_species" in result.columns

        # Check values are reasonable
        assert all(result["shannon"] >= 0)
        assert all(result["simpson"] >= 0)
        assert all(result["simpson"] <= 1)
        assert all(result["chao1"] >= result["observed_species"])

    def test_calculate_alpha_diversity_single_metric(self):
        """Test alpha diversity with single metric."""
        result = calculate_alpha_diversity(self.abundance_df, metrics=["shannon"])

        assert result.shape[1] == 1
        assert "shannon" in result.columns

    def test_calculate_alpha_diversity_numpy_input(self):
        """Test alpha diversity with numpy array input."""
        result = calculate_alpha_diversity(self.abundance_data)

        assert isinstance(result, pd.DataFrame)
        assert result.shape == (5, 4)  # Default metrics

    def test_calculate_beta_diversity_bray_curtis(self):
        """Test beta diversity with Bray-Curtis distance."""
        result = calculate_beta_diversity(self.abundance_df, metric="bray_curtis")

        # Check structure
        assert isinstance(result, pd.DataFrame)
        assert result.shape == (5, 5)  # Square matrix
        assert np.allclose(np.diag(result), 0)  # Diagonal should be zero
        assert np.allclose(result.values, result.values.T)  # Should be symmetric

        # Check values are in valid range [0, 1]
        assert np.all(result.values >= 0)
        assert np.all(result.values <= 1)

    def test_calculate_beta_diversity_numpy_output(self):
        """Test beta diversity with numpy array output."""
        result = calculate_beta_diversity(self.abundance_df, metric="bray_curtis", return_dataframe=False)

        assert isinstance(result, np.ndarray)
        assert result.shape == (5, 5)

    def test_perform_permanova(self):
        """Test PERMANOVA analysis."""
        # Calculate beta diversity first
        distances = calculate_beta_diversity(self.abundance_df, metric="bray_curtis", return_dataframe=True)

        # Run PERMANOVA
        result = perform_permanova(
            distances, self.metadata, formula="treatment", n_permutations=99  # Use fewer permutations for speed
        )

        # Check structure
        assert isinstance(result, dict)
        assert "treatment" in result
        assert "F_statistic" in result["treatment"]
        assert "p_value" in result["treatment"]
        assert "significant" in result["treatment"]

        # Check values
        assert result["treatment"]["F_statistic"] >= 0
        assert 0 <= result["treatment"]["p_value"] <= 1
        assert isinstance(result["treatment"]["significant"], bool)

    def test_invalid_metric(self):
        """Test handling of invalid metrics."""
        # Should not crash, just warn and skip unknown metrics
        result = calculate_alpha_diversity(self.abundance_df, metrics=["shannon", "invalid_metric"])

        assert "shannon" in result.columns
        assert "invalid_metric" not in result.columns


class TestDiversityEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_matrix(self):
        """Test with empty abundance matrix."""
        empty_df = pd.DataFrame()

        result = calculate_alpha_diversity(empty_df)

        # Should return empty DataFrame
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0

    def test_single_sample(self):
        """Test with single sample."""
        single_sample = pd.DataFrame([[1, 2, 3, 4]], columns=["A", "B", "C", "D"])

        result = calculate_alpha_diversity(single_sample)
        assert result.shape == (1, 4)  # Default metrics

    def test_zero_abundances(self):
        """Test with all-zero abundances."""
        zero_df = pd.DataFrame([[0, 0, 0, 0], [0, 0, 0, 0]])

        result = calculate_alpha_diversity(zero_df)

        # Shannon should be 0 for no species
        assert result["shannon"].iloc[0] == 0
        assert result["observed_species"].iloc[0] == 0

    def test_permanova_missing_variable(self):
        """Test PERMANOVA with missing variable."""
        distances = pd.DataFrame([[0, 1], [1, 0]])
        metadata = pd.DataFrame({"other_var": ["A", "B"]})

        with pytest.raises(Exception):  # Should raise ProcessingError
            perform_permanova(distances, metadata, "missing_var")


def _reference_f(d, labels):
    """Plain-loop PERMANOVA pseudo-F, Anderson 2001 eq. 4, written independently of the code."""
    n = len(labels)
    groups = sorted(set(labels))
    a = len(groups)
    ss_t = sum(d[i][j] ** 2 for i in range(n) for j in range(i + 1, n)) / n
    ss_w = 0.0
    for g in groups:
        idx = [i for i in range(n) if labels[i] == g]
        ss_w += sum(d[i][j] ** 2 for i in idx for j in idx if i < j) / len(idx)
    return ((ss_t - ss_w) / (a - 1)) / (ss_w / (n - a))


def test_permanova_matches_the_reference_formula():
    rng = np.random.default_rng(0)
    points = np.vstack([rng.normal(0, 1, (5, 3)), rng.normal(2, 1, (5, 3))])
    d = squareform(pdist(points))
    labels = ["a"] * 5 + ["b"] * 5
    f, p = _permanova_test(d, pd.Series(labels), 199, rng=np.random.default_rng(1))
    assert f == pytest.approx(_reference_f(d, labels))
    assert 0 < p <= 1 and p >= 1 / 200


def test_permanova_two_identical_clusters_has_f_near_one_and_high_p():
    """Labels assigned at random within one cluster carry no group structure."""
    rng = np.random.default_rng(2)
    points = rng.normal(0, 1, (12, 3))
    d = squareform(pdist(points))
    labels = ["a", "b"] * 6
    f, p = _permanova_test(d, pd.Series(labels), 199, rng=np.random.default_rng(3))
    assert f == pytest.approx(_reference_f(d, labels))
    assert 0.2 < f < 3.0
    assert p > 0.2


def test_permanova_nan_labels_are_dropped_and_single_group_is_nan():
    d = np.ones((4, 4)) - np.eye(4)
    f, p = _permanova_test(d, pd.Series(["a", "a", np.nan, np.nan]), 10)
    assert np.isnan(f) and np.isnan(p)


def test_permanova_nan_labels_are_dropped_before_the_statistic():
    rng = np.random.default_rng(4)
    points = np.vstack([rng.normal(0, 1, (4, 2)), rng.normal(3, 1, (4, 2)), rng.normal(1, 1, (2, 2))])
    d = squareform(pdist(points))
    labels = ["a"] * 4 + ["b"] * 4 + [np.nan, np.nan]
    f, _ = _permanova_test(d, pd.Series(labels), 9, rng=np.random.default_rng(5))
    assert f == pytest.approx(_reference_f(d[:8, :8], labels[:8]))


def test_perform_permanova_seed_is_reproducible():
    rng = np.random.default_rng(6)
    points = np.vstack([rng.normal(0, 1, (4, 2)), rng.normal(1, 1, (4, 2))])
    names = [f"S{i}" for i in range(8)]
    d = pd.DataFrame(squareform(pdist(points)), index=names, columns=names)
    metadata = pd.DataFrame({"treatment": ["x"] * 4 + ["y"] * 4}, index=names)
    first = perform_permanova(d, metadata, "treatment", n_permutations=99, seed=11)
    second = perform_permanova(d, metadata, "treatment", n_permutations=99, seed=11)
    assert first == second


def test_perform_permanova_names_samples_missing_from_metadata():
    names = ["S1", "S2", "S3", "S4"]
    d = pd.DataFrame(np.ones((4, 4)) - np.eye(4), index=names, columns=names)
    metadata = pd.DataFrame({"treatment": ["x", "y"]}, index=["S1", "S2"])
    with pytest.raises(ValidationError, match="S3.*S4"):
        perform_permanova(d, metadata, "treatment", n_permutations=9)


if __name__ == "__main__":
    pytest.main([__file__])


def test_permanova_keeps_the_first_row_of_a_duplicated_sample():
    """A sample listed twice in the metadata is read from its first row, as the plots do."""
    names = [f"S{i}" for i in range(1, 6)]
    abundance = pd.DataFrame(
        [[10, 5, 2, 0], [0, 8, 3, 1], [15, 0, 0, 5], [5, 5, 5, 5], [2, 2, 2, 14]],
        index=names,
        columns=["a", "b", "c", "d"],
    )
    distances = calculate_beta_diversity(abundance, metric="bray_curtis", return_dataframe=True)
    metadata = pd.DataFrame({"treatment": ["A", "A", "B", "B", "B"]}, index=names)
    duplicated = pd.concat([metadata, pd.DataFrame({"treatment": ["A"]}, index=["S5"])])

    expected = perform_permanova(distances, metadata, "treatment", n_permutations=9, seed=1)
    result = perform_permanova(distances, duplicated, "treatment", n_permutations=9, seed=1)

    assert result["treatment"]["F_statistic"] == pytest.approx(expected["treatment"]["F_statistic"])
