"""
Test processing.statistics module functionality.

Tests for statistical analysis functions including enrichment and distance calculations.
"""

from unittest.mock import patch
import numpy as np
import pandas as pd
import importlib.util
import pytest

from prepare.processing.statistics import (
    calculate_enrichment,
    _validate_distance_method,
    _create_binary_presence_matrix,
    _calculate_jaccard_distance,
    _calculate_dice_distance,
    _calculate_cosine_distance,
    calculate_distance_matrix,
    perform_hypergeometric_test,
)
from prepare.core.exceptions import ProcessingError

requires_analysis = pytest.mark.skipif(
    not all(importlib.util.find_spec(m) for m in ("sklearn", "scipy")), reason="needs prepare[analysis]"
)


@requires_analysis
class TestCalculateEnrichment:
    """Test calculate_enrichment function."""

    def test_calculate_enrichment_basic(self):
        """Test basic enrichment calculation."""
        observed = {"gene_A": 10, "gene_B": 5, "gene_C": 15}
        expected = {"gene_A": 5, "gene_B": 10, "gene_C": 15}

        result = calculate_enrichment(observed, expected, min_count=1)

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 3
        assert "fold_enrichment" in result.columns
        assert "p_value" in result.columns
        assert "adjusted_p_value" in result.columns

        # Check enrichment values
        gene_a_row = result[result["key"] == "gene_A"].iloc[0]
        assert gene_a_row["fold_enrichment"] == 2.0  # 10/5

        gene_b_row = result[result["key"] == "gene_B"].iloc[0]
        assert gene_b_row["fold_enrichment"] == 0.5  # 5/10

    def test_calculate_enrichment_no_common_keys(self):
        """Test with no common keys between observed and expected."""
        observed = {"gene_A": 10}
        expected = {"gene_B": 5}

        result = calculate_enrichment(observed, expected)

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0

    def test_calculate_enrichment_min_count_filter(self):
        """Test min_count filtering."""
        observed = {"gene_A": 10, "gene_B": 2}
        expected = {"gene_A": 5, "gene_B": 1}

        result = calculate_enrichment(observed, expected, min_count=5)

        assert len(result) == 1
        assert result.iloc[0]["key"] == "gene_A"

    def test_calculate_enrichment_zero_expected(self):
        """Test with zero expected count."""
        observed = {"gene_A": 10}
        expected = {"gene_A": 0}

        result = calculate_enrichment(observed, expected, min_count=1)

        assert len(result) == 1
        assert result.iloc[0]["fold_enrichment"] == float("inf")

    def test_calculate_enrichment_zero_observed(self):
        """Test with zero observed count."""
        observed = {"gene_A": 0}
        expected = {"gene_A": 10}

        result = calculate_enrichment(observed, expected, min_count=0)

        assert len(result) == 1
        assert result.iloc[0]["fold_enrichment"] == 0.0

    def test_calculate_enrichment_single_result(self):
        """Test with single result (no multiple testing correction)."""
        observed = {"gene_A": 10}
        expected = {"gene_A": 5}

        result = calculate_enrichment(observed, expected, min_count=1)

        assert len(result) == 1
        # With single result, adjusted p-value should equal p-value
        assert result.iloc[0]["adjusted_p_value"] == result.iloc[0]["p_value"]

    def test_calculate_enrichment_empty_results(self):
        """Test when all results are filtered out."""
        observed = {"gene_A": 1}
        expected = {"gene_A": 1}

        result = calculate_enrichment(observed, expected, min_count=10)

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0

    def test_calculate_enrichment_exception_handling(self):
        """Test exception handling."""
        # Invalid input type should raise ProcessingError
        with pytest.raises(ProcessingError, match="Error calculating enrichment"):
            calculate_enrichment("invalid", {"gene_A": 5})


class TestValidateDistanceMethod:
    """Test _validate_distance_method function."""

    def test_validate_distance_method_valid(self):
        """Test with valid methods."""
        valid_methods = ["jaccard", "dice", "cosine"]

        for method in valid_methods:
            # Should not raise exception
            _validate_distance_method(method)

    def test_validate_distance_method_invalid(self):
        """Test with invalid method."""
        with pytest.raises(ProcessingError, match="Invalid distance method"):
            _validate_distance_method("invalid_method")


class TestCreateBinaryPresenceMatrix:
    """Test _create_binary_presence_matrix function."""

    def test_create_binary_presence_matrix(self):
        """Test binary presence matrix creation."""
        summary_df = pd.DataFrame(
            {"GCF_001": [0.5, 0.1, 0.8], "GCF_002": [0.2, 0.6, 0.0], "other_col": [1, 2, 3]},
            index=["sample1", "sample2", "sample3"],
        )

        genome_columns = ["GCF_001", "GCF_002"]
        threshold = 0.3

        result = _create_binary_presence_matrix(summary_df, genome_columns, threshold)

        assert isinstance(result, pd.DataFrame)
        assert list(result.columns) == genome_columns
        assert list(result.index) == list(summary_df.index)

        # Check binary values
        assert result.loc["sample1", "GCF_001"] == 1  # 0.5 > 0.3
        assert result.loc["sample2", "GCF_001"] == 0  # 0.1 < 0.3
        assert result.loc["sample2", "GCF_002"] == 1  # 0.6 > 0.3


class TestDistanceFunctions:
    """Test distance calculation functions."""

    def test_calculate_jaccard_distance(self):
        """Test Jaccard distance calculation."""
        vec1 = np.array([1, 1, 0, 1])
        vec2 = np.array([1, 0, 1, 1])

        distance = _calculate_jaccard_distance(vec1, vec2)

        # Intersection: 2 (positions 0 and 3)
        # Union: 4 (all positions have at least one 1)
        # Jaccard similarity: 2/4 = 0.5
        # Jaccard distance: 1 - 0.5 = 0.5
        assert distance == 0.5

    def test_calculate_jaccard_distance_identical(self):
        """Test Jaccard distance with identical vectors."""
        vec1 = np.array([1, 0, 1, 0])
        vec2 = np.array([1, 0, 1, 0])

        distance = _calculate_jaccard_distance(vec1, vec2)

        assert distance == 0.0

    def test_calculate_jaccard_distance_no_union(self):
        """Test Jaccard distance with no union (all zeros)."""
        vec1 = np.array([0, 0, 0])
        vec2 = np.array([0, 0, 0])

        distance = _calculate_jaccard_distance(vec1, vec2)

        assert distance == 1.0

    def test_calculate_dice_distance(self):
        """Test Dice distance calculation."""
        vec1 = np.array([1, 1, 0, 1])
        vec2 = np.array([1, 0, 1, 1])

        distance = _calculate_dice_distance(vec1, vec2)

        # Intersection: 2, sum1: 3, sum2: 3
        # Dice similarity: 2*2/(3+3) = 4/6 = 0.667
        # Dice distance: 1 - 0.667 = 0.333
        assert abs(distance - (1 - 4 / 6)) < 1e-6

    def test_calculate_dice_distance_no_sum(self):
        """Test Dice distance with zero sums."""
        vec1 = np.array([0, 0, 0])
        vec2 = np.array([0, 0, 0])

        distance = _calculate_dice_distance(vec1, vec2)

        assert distance == 1.0

    def test_calculate_cosine_distance(self):
        """Test Cosine distance calculation."""
        vec1 = np.array([1, 1, 0], dtype=float)
        vec2 = np.array([1, 0, 1], dtype=float)

        distance = _calculate_cosine_distance(vec1, vec2)

        # Dot product: 1, norm1: sqrt(2), norm2: sqrt(2)
        # Cosine similarity: 1/(sqrt(2)*sqrt(2)) = 1/2 = 0.5
        # Cosine distance: 1 - 0.5 = 0.5
        assert abs(distance - 0.5) < 1e-10

    def test_calculate_cosine_distance_zero_norm(self):
        """Test Cosine distance with zero norm."""
        vec1 = np.array([0, 0, 0])
        vec2 = np.array([1, 0, 1])

        distance = _calculate_cosine_distance(vec1, vec2)

        assert distance == 1.0


class TestCalculateDistanceMatrix:
    """Test calculate_distance_matrix function."""

    def test_calculate_distance_matrix_jaccard(self, tmp_path):
        """Test distance matrix calculation with Jaccard method."""
        # Create test file
        summary_file = tmp_path / "summary.txt"
        summary_data = pd.DataFrame(
            {
                "GCF_001": [0.5, 0.1, 0.8],
                "GCF_002": [0.2, 0.6, 0.0],
            },
            index=["sample1", "sample2", "sample3"],
        )
        summary_data.to_csv(summary_file, sep="\t", index_label="sample")

        result = calculate_distance_matrix(summary_file, threshold=0.3, method="jaccard")

        assert isinstance(result, pd.DataFrame)
        assert result.shape == (3, 3)
        assert list(result.index) == list(result.columns)

        # Diagonal should be zeros
        np.testing.assert_array_equal(np.diag(result), [0.0, 0.0, 0.0])

        # Matrix should be symmetric
        assert result.loc["sample1", "sample2"] == result.loc["sample2", "sample1"]

    def test_calculate_distance_matrix_different_methods(self, tmp_path):
        """Test distance matrix with different methods."""
        # Create test file
        summary_file = tmp_path / "summary.txt"
        summary_data = pd.DataFrame(
            {
                "GCF_001": [0.5, 0.1],
                "GCF_002": [0.2, 0.6],
            },
            index=["sample1", "sample2"],
        )
        summary_data.to_csv(summary_file, sep="\t", index_label="sample")

        for method in ["jaccard", "dice", "cosine"]:
            result = calculate_distance_matrix(summary_file, method=method)
            assert isinstance(result, pd.DataFrame)
            assert result.shape == (2, 2)

    def test_calculate_distance_matrix_invalid_method(self, tmp_path):
        """Test with invalid distance method."""
        summary_file = tmp_path / "summary.txt"
        summary_data = pd.DataFrame(
            {
                "GCF_001": [0.5, 0.1],
            },
            index=["sample1", "sample2"],
        )
        summary_data.to_csv(summary_file, sep="\t", index_label="sample")

        with pytest.raises(ProcessingError, match="Invalid distance method"):
            calculate_distance_matrix(summary_file, method="invalid")

    def test_calculate_distance_matrix_no_genome_columns(self, tmp_path):
        """Test with no genome columns."""
        summary_file = tmp_path / "summary.txt"
        summary_data = pd.DataFrame(
            {
                "max_containment": [0.5, 0.1],
                "max_containment_annotation": ["A", "B"],
            },
            index=["sample1", "sample2"],
        )
        summary_data.to_csv(summary_file, sep="\t", index_label="sample")

        with pytest.raises(ProcessingError, match="No genome columns found"):
            calculate_distance_matrix(summary_file)

    def test_calculate_distance_matrix_file_not_found(self):
        """Test with non-existent file."""
        with pytest.raises(ProcessingError, match="Error calculating distance matrix"):
            calculate_distance_matrix("nonexistent_file.txt")


@requires_analysis
class TestPerformHypergeometricTest:
    """Test perform_hypergeometric_test function."""

    def test_perform_hypergeometric_test_basic(self):
        """Test basic hypergeometric test."""
        fold_enrichment, p_value = perform_hypergeometric_test(
            success_in_sample=5, sample_size=100, success_in_population=10, population_size=1000
        )

        # Expected = (10/1000) * 100 = 1
        # Fold enrichment = 5/1 = 5
        assert fold_enrichment == 5.0
        assert isinstance(p_value, float)
        assert 0 <= p_value <= 1

    def test_perform_hypergeometric_test_no_enrichment(self):
        """Test with no enrichment."""
        fold_enrichment, p_value = perform_hypergeometric_test(
            success_in_sample=1, sample_size=100, success_in_population=10, population_size=1000
        )

        # Expected = (10/1000) * 100 = 1
        # Fold enrichment = 1/1 = 1
        assert fold_enrichment == 1.0

    def test_perform_hypergeometric_test_zero_expected(self):
        """Test with zero expected count."""
        fold_enrichment, p_value = perform_hypergeometric_test(
            success_in_sample=5, sample_size=100, success_in_population=0, population_size=1000
        )

        # Expected = 0, fold_enrichment should be inf
        assert fold_enrichment == float("inf")

    def test_perform_hypergeometric_test_zero_observed(self):
        """Test with zero observed count but nonzero expected."""
        fold_enrichment, p_value = perform_hypergeometric_test(
            success_in_sample=0, sample_size=100, success_in_population=10, population_size=1000
        )

        # Expected = 1, observed = 0
        assert fold_enrichment == 0.0

    def test_perform_hypergeometric_test_zero_expected_zero_observed(self):
        """Test with both zero expected and observed."""
        fold_enrichment, p_value = perform_hypergeometric_test(
            success_in_sample=0, sample_size=100, success_in_population=0, population_size=1000
        )

        # Both zero, should return 1.0
        assert fold_enrichment == 1.0

    def test_perform_hypergeometric_test_invalid_sample_size(self):
        """Test with invalid sample size."""
        with pytest.raises(ProcessingError, match="Sample and population sizes must be positive"):
            perform_hypergeometric_test(5, 0, 10, 1000)

    def test_perform_hypergeometric_test_invalid_population_size(self):
        """Test with invalid population size."""
        with pytest.raises(ProcessingError, match="Sample and population sizes must be positive"):
            perform_hypergeometric_test(5, 100, 10, 0)

    def test_perform_hypergeometric_test_success_exceeds_sample(self):
        """Test with success count exceeding sample size."""
        with pytest.raises(ProcessingError, match="Success count cannot exceed sample size"):
            perform_hypergeometric_test(150, 100, 10, 1000)

    def test_perform_hypergeometric_test_success_exceeds_population(self):
        """Test with success count exceeding population size."""
        with pytest.raises(ProcessingError, match="Success count cannot exceed population size"):
            perform_hypergeometric_test(5, 100, 1500, 1000)

    def test_perform_hypergeometric_test_edge_case_equal_sizes(self):
        """Test edge case where sample equals population."""
        fold_enrichment, p_value = perform_hypergeometric_test(
            success_in_sample=10, sample_size=100, success_in_population=10, population_size=100
        )

        # Expected = (10/100) * 100 = 10
        # Fold enrichment = 10/10 = 1
        assert fold_enrichment == 1.0


class TestStatisticsIntegration:
    """Integration tests for statistics module."""

    @requires_analysis
    def test_full_enrichment_workflow(self):
        """Test complete enrichment analysis workflow."""
        # Simulate realistic data
        observed = {"metabolism": 50, "transport": 30, "regulation": 20, "unknown": 10}
        expected = {"metabolism": 40, "transport": 40, "regulation": 15, "unknown": 5}

        result = calculate_enrichment(observed, expected, min_count=5)

        assert len(result) == 4
        assert all(
            col in result.columns
            for col in ["key", "observed", "expected", "fold_enrichment", "p_value", "adjusted_p_value"]
        )

        # Check that results are sorted by fold enrichment
        assert result["fold_enrichment"].is_monotonic_decreasing

    def test_distance_matrix_workflow(self, tmp_path):
        """Test complete distance matrix workflow."""
        # Create realistic test data
        summary_file = tmp_path / "summary.txt"
        np.random.seed(42)  # For reproducible results

        # Create data with some patterns
        data = {
            "GCF_001": [0.8, 0.1, 0.9, 0.0],
            "GCF_002": [0.2, 0.9, 0.1, 0.8],
            "GCF_003": [0.7, 0.0, 0.8, 0.1],
        }
        summary_data = pd.DataFrame(data, index=["sample1", "sample2", "sample3", "sample4"])
        summary_data.to_csv(summary_file, sep="\t", index_label="sample")

        # Test all methods
        for method in ["jaccard", "dice", "cosine"]:
            result = calculate_distance_matrix(summary_file, threshold=0.5, method=method)

            # Basic validation
            assert result.shape == (4, 4)

            # Convert to float to avoid type issues
            result_numeric = result.astype(float)
            assert np.allclose(np.diag(result_numeric.values), 0)  # Diagonal should be zero
            assert np.allclose(result_numeric.values, result_numeric.values.T)  # Should be symmetric


class TestEdgeCasesAndErrorHandling:
    """Test edge cases and error handling not covered elsewhere."""

    @requires_analysis
    def test_calculate_enrichment_both_zero_counts(self):
        """Test enrichment when both observed and expected are zero."""
        observed = {"gene_A": 0}
        expected = {"gene_A": 0}

        result = calculate_enrichment(observed, expected, min_count=0)

        assert len(result) == 1
        assert result.iloc[0]["fold_enrichment"] == 1.0

    @requires_analysis
    def test_calculate_enrichment_large_numbers(self):
        """Test enrichment with very large numbers."""
        observed = {"gene_A": 1000000}
        expected = {"gene_A": 500000}

        result = calculate_enrichment(observed, expected)

        assert len(result) == 1
        assert result.iloc[0]["fold_enrichment"] == 2.0

    @requires_analysis
    @patch("scipy.stats.fisher_exact")
    @patch("scipy.stats.false_discovery_control")
    def test_calculate_enrichment_statistical_functions(self, mock_fdr, mock_fisher):
        """Test that statistical functions are called correctly."""
        mock_fisher.return_value = (2.0, 0.05)
        mock_fdr.return_value = [0.1, 0.15]

        observed = {"gene_A": 10, "gene_B": 5}
        expected = {"gene_A": 5, "gene_B": 8}

        calculate_enrichment(observed, expected)

        # Should call fisher_exact twice (once per gene)
        assert mock_fisher.call_count == 2
        # Should call false_discovery_control once (multiple results)
        mock_fdr.assert_called_once()

    def test_create_binary_presence_matrix_empty_input(self):
        """Test binary matrix creation with empty DataFrame."""
        empty_df = pd.DataFrame(index=["sample1"])

        result = _create_binary_presence_matrix(empty_df, [], 0.5)

        assert len(result.columns) == 0
        assert list(result.index) == ["sample1"]

    def test_create_binary_presence_matrix_extreme_threshold(self):
        """Test binary matrix with extreme threshold values."""
        summary_df = pd.DataFrame(
            {
                "GCF_001": [0.5, 0.1, 0.8],
            },
            index=["sample1", "sample2", "sample3"],
        )

        # Very low threshold - all should be 1
        result_low = _create_binary_presence_matrix(summary_df, ["GCF_001"], 0.0)
        assert result_low["GCF_001"].sum() == 3

        # Very high threshold - all should be 0
        result_high = _create_binary_presence_matrix(summary_df, ["GCF_001"], 1.0)
        assert result_high["GCF_001"].sum() == 0

    def test_distance_functions_with_different_array_types(self):
        """Test distance functions with different numpy array types."""
        # Test with different dtypes
        vec1_int = np.array([1, 0, 1, 0], dtype=int)
        vec2_int = np.array([1, 1, 0, 0], dtype=int)

        vec1_bool = np.array([True, False, True, False], dtype=bool)
        vec2_bool = np.array([True, True, False, False], dtype=bool)

        # Results should be the same regardless of dtype
        jaccard_int = _calculate_jaccard_distance(vec1_int, vec2_int)
        jaccard_bool = _calculate_jaccard_distance(vec1_bool, vec2_bool)
        assert jaccard_int == jaccard_bool

    def test_distance_functions_with_very_long_vectors(self):
        """Test distance functions with very long vectors."""
        # Test performance and correctness with large vectors
        size = 10000
        vec1 = np.random.choice([0, 1], size=size, p=[0.7, 0.3])
        vec2 = np.random.choice([0, 1], size=size, p=[0.7, 0.3])

        # Should not raise exceptions and should return valid distances
        jaccard_dist = _calculate_jaccard_distance(vec1, vec2)
        dice_dist = _calculate_dice_distance(vec1, vec2)
        cosine_dist = _calculate_cosine_distance(vec1.astype(float), vec2.astype(float))

        assert 0 <= jaccard_dist <= 1
        assert 0 <= dice_dist <= 1
        assert 0 <= cosine_dist <= 1

    def test_calculate_distance_matrix_single_sample(self, tmp_path):
        """Test distance matrix with only one sample."""
        summary_file = tmp_path / "summary.txt"
        summary_data = pd.DataFrame(
            {
                "GCF_001": [0.8],
            },
            index=["sample1"],
        )
        summary_data.to_csv(summary_file, sep="\t", index_label="sample")

        result = calculate_distance_matrix(summary_file)

        assert result.shape == (1, 1)
        assert result.iloc[0, 0] == 0.0

    def test_calculate_distance_matrix_all_identical_samples(self, tmp_path):
        """Test distance matrix where all samples are identical."""
        summary_file = tmp_path / "summary.txt"
        summary_data = pd.DataFrame(
            {
                "GCF_001": [0.8, 0.8, 0.8],
                "GCF_002": [0.2, 0.2, 0.2],
            },
            index=["sample1", "sample2", "sample3"],
        )
        summary_data.to_csv(summary_file, sep="\t", index_label="sample")

        result = calculate_distance_matrix(summary_file, threshold=0.5)

        # All off-diagonal elements should be 0 (identical samples)
        for i in range(3):
            for j in range(i + 1, 3):
                assert result.iloc[i, j] == 0.0

    def test_calculate_distance_matrix_all_different_samples(self, tmp_path):
        """Test distance matrix where all samples are completely different."""
        summary_file = tmp_path / "summary.txt"
        summary_data = pd.DataFrame(
            {
                "GCF_001": [1.0, 0.0, 0.0],
                "GCF_002": [0.0, 1.0, 0.0],
                "GCF_003": [0.0, 0.0, 1.0],
            },
            index=["sample1", "sample2", "sample3"],
        )
        summary_data.to_csv(summary_file, sep="\t", index_label="sample")

        result = calculate_distance_matrix(summary_file, threshold=0.5)

        # All off-diagonal elements should be 1.0 (completely different)
        for i in range(3):
            for j in range(i + 1, 3):
                assert result.iloc[i, j] == 1.0

    @patch("prepare.processing.statistics.logger")
    def test_calculate_distance_matrix_logging(self, mock_logger, tmp_path):
        """Test that appropriate logging occurs."""
        summary_file = tmp_path / "summary.txt"
        summary_data = pd.DataFrame(
            {
                "GCF_001": [0.8, 0.2],
            },
            index=["sample1", "sample2"],
        )
        summary_data.to_csv(summary_file, sep="\t", index_label="sample")

        calculate_distance_matrix(summary_file, method="jaccard")

        # Should log completion message
        mock_logger.info.assert_called()
        logged_message = mock_logger.info.call_args[0][0]
        assert "jaccard distances for 2 samples" in logged_message

    @requires_analysis
    def test_perform_hypergeometric_test_boundary_values(self):
        """Test hypergeometric test with boundary values."""
        # Test with minimum valid inputs
        fold_enrichment, p_value = perform_hypergeometric_test(0, 1, 0, 1)
        assert fold_enrichment == 1.0
        assert 0 <= p_value <= 1

        # Test with maximum overlap
        fold_enrichment, p_value = perform_hypergeometric_test(10, 10, 10, 10)
        assert fold_enrichment == 1.0

        # Test with very large population
        # Expected = (1000/1000000) * 100 = 0.1, fold_enrichment = 5/0.1 = 50
        fold_enrichment, p_value = perform_hypergeometric_test(5, 100, 1000, 1000000)
        assert fold_enrichment == 50.0

    @requires_analysis
    def test_perform_hypergeometric_test_negative_inputs(self):
        """Test hypergeometric test with negative inputs."""
        # Test negative sample size
        with pytest.raises(ProcessingError, match="Sample and population sizes must be positive"):
            perform_hypergeometric_test(5, -100, 10, 1000)

        # Test negative population size
        with pytest.raises(ProcessingError, match="Sample and population sizes must be positive"):
            perform_hypergeometric_test(5, 100, 10, -1000)

    @requires_analysis
    @patch("scipy.stats.hypergeom.sf")
    def test_perform_hypergeometric_test_scipy_exception(self, mock_sf):
        """Test handling of scipy exceptions."""
        mock_sf.side_effect = ValueError("Invalid parameters")

        with pytest.raises(ProcessingError, match="Error performing hypergeometric test"):
            perform_hypergeometric_test(5, 100, 10, 1000)

    def test_distance_matrix_memory_efficiency(self, tmp_path):
        """Test distance matrix calculation with moderately large data."""
        # Test with 100 samples to check memory efficiency
        summary_file = tmp_path / "summary.txt"

        # Create synthetic data
        np.random.seed(123)
        n_samples = 100
        n_genomes = 20

        data = {}
        for i in range(n_genomes):
            data[f"GCF_{i:03d}"] = np.random.random(n_samples)

        summary_data = pd.DataFrame(data, index=[f"sample_{i}" for i in range(n_samples)])
        summary_data.to_csv(summary_file, sep="\t", index_label="sample")

        # Should complete without memory issues
        result = calculate_distance_matrix(summary_file, threshold=0.5, method="jaccard")

        assert result.shape == (n_samples, n_samples)

        # Verify properties
        # Convert to numeric array to avoid object dtype issues
        result_numeric = result.astype(float).values
        assert np.allclose(np.diag(result_numeric), 0)  # Diagonal is zero
        assert np.allclose(result_numeric, result_numeric.T)  # Symmetric


class TestNumericPrecision:
    """Test numeric precision and floating point edge cases."""

    @requires_analysis
    def test_enrichment_very_small_numbers(self):
        """Test enrichment with very small floating point numbers."""
        observed = {"gene_A": 1e-10}
        expected = {"gene_A": 1e-15}

        result = calculate_enrichment(observed, expected, min_count=0)

        assert len(result) == 1
        assert result.iloc[0]["fold_enrichment"] == 1e5

    def test_distance_functions_floating_point_precision(self):
        """Test distance functions with floating point precision issues."""
        # Test with numbers that might cause precision issues
        vec1 = np.array([0.1 + 0.2, 0.7, 0.0])  # 0.30000000000000004
        vec2 = np.array([0.3, 0.7, 0.0])

        # Convert to binary based on threshold
        vec1_binary = (vec1 > 0.5).astype(int)
        vec2_binary = (vec2 > 0.5).astype(int)

        distance = _calculate_jaccard_distance(vec1_binary, vec2_binary)
        assert distance == 0.0  # Should be identical after thresholding

    @requires_analysis
    def test_hypergeometric_test_precision_with_large_numbers(self):
        """Test hypergeometric test precision with large numbers."""
        # Test with large numbers that might cause precision issues
        fold_enrichment, p_value = perform_hypergeometric_test(
            success_in_sample=1000, sample_size=10000, success_in_population=5000, population_size=100000
        )

        # Expected = (5000/100000) * 10000 = 500
        # Fold enrichment = 1000/500 = 2.0
        assert abs(fold_enrichment - 2.0) < 1e-10
        assert 0 <= p_value <= 1


if __name__ == "__main__":
    pytest.main([__file__])


class TestCompareGroupMeans:
    """Tests for the t-test / ANOVA choice used when comparing dataset groups."""

    @requires_analysis
    def test_two_groups_use_a_t_test(self):
        from prepare.processing.statistics import compare_group_means

        name, statistic, p_value = compare_group_means([np.array([1.0, 2.0, 3.0]), np.array([7.0, 8.0, 9.0])])
        assert name == "t-test"
        assert statistic < 0
        assert p_value < 0.05

    @requires_analysis
    def test_three_groups_use_anova(self):
        from prepare.processing.statistics import compare_group_means

        groups = [np.array([1.0, 2.0, 3.0]), np.array([1.5, 2.5, 3.5]), np.array([9.0, 10.0, 11.0])]
        name, statistic, p_value = compare_group_means(groups)
        assert name == "ANOVA"
        assert statistic > 0
        assert p_value < 0.05

    def test_missing_scipy_names_the_analysis_extra(self, monkeypatch):
        import sys

        from prepare.core.exceptions import ConfigurationError
        from prepare.processing.statistics import compare_group_means

        monkeypatch.setitem(sys.modules, "scipy", None)
        with pytest.raises(ConfigurationError, match=r"prepare\[analysis\]"):
            compare_group_means([np.array([1.0, 2.0]), np.array([3.0, 4.0])])
