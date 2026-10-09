"""
Tests for diversity analysis functions.
"""

import numpy as np
import pandas as pd
import pytest

from prepare.processing.diversity import (
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


if __name__ == "__main__":
    pytest.main([__file__])
