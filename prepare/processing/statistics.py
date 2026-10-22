"""
Statistical analysis for PREPARE.

This module provides functions for statistical analysis of genomic data.
"""

import logging
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Sequence, Tuple, Union

from prepare.core.exceptions import ProcessingError
from prepare.core.optional import require
from prepare.core.utils import get_genome_columns

logger = logging.getLogger(__name__)


def calculate_enrichment(
    observed_counts: Dict[str, int], expected_counts: Dict[str, int], min_count: int = 5
) -> pd.DataFrame:
    """
    Calculate enrichment of observed counts relative to expected.

    Args:
        observed_counts: Dictionary of observed counts
        expected_counts: Dictionary of expected counts
        min_count: Minimum count for inclusion in results

    Returns:
        DataFrame with enrichment statistics

    Raises:
        ProcessingError: If the calculation fails
        ConfigurationError: If scipy is not installed
    """
    stats = require("scipy.stats", "analysis", "An enrichment test")
    try:
        # Prepare counts for common keys
        common_keys = set(observed_counts.keys()) & set(expected_counts.keys())

        if not common_keys:
            logger.warning("No common keys between observed and expected counts")
            return pd.DataFrame()

        results = []

        # Calculate statistics for each key
        for key in common_keys:
            observed = observed_counts.get(key, 0)
            expected = expected_counts.get(key, 0)

            # Skip items with low counts
            if observed < min_count and expected < min_count:
                continue

            # Calculate fold enrichment
            if expected > 0:
                fold_enrichment = observed / expected
            else:
                fold_enrichment = float("inf") if observed > 0 else 1.0

            # Calculate p-value with Fisher's exact test
            # Create contingency table
            # [[observed, total_observed - observed],
            #  [expected, total_expected - expected]]
            total_observed = sum(observed_counts.values())
            total_expected = sum(expected_counts.values())

            contingency_table = np.array(
                [
                    [observed, total_observed - observed],
                    [expected, total_expected - expected],
                ]
            )

            odds_ratio, p_value = stats.fisher_exact(contingency_table)

            # Store results
            results.append(
                {
                    "key": key,
                    "observed": observed,
                    "expected": expected,
                    "fold_enrichment": fold_enrichment,
                    "p_value": p_value,
                }
            )

        # Create DataFrame and sort by enrichment
        if not results:
            logger.warning("No results after filtering")
            return pd.DataFrame()

        df = pd.DataFrame(results)

        # Add multiple testing correction
        if len(df) > 1:
            df["adjusted_p_value"] = stats.false_discovery_control(df["p_value"])
        else:
            df["adjusted_p_value"] = df["p_value"]

        # Sort by fold enrichment
        df = df.sort_values("fold_enrichment", ascending=False)

        return df

    except Exception as e:
        raise ProcessingError(f"Error calculating enrichment: {e}")


def _validate_distance_method(method):
    """
    Validate the distance calculation method.

    Args:
        method: Distance method name

    Raises:
        ProcessingError: If method is invalid
    """
    valid_methods = ("jaccard", "dice", "cosine")
    if method not in valid_methods:
        raise ProcessingError(f"Invalid distance method: {method}. " f"Valid methods: {', '.join(valid_methods)}")


def _create_binary_presence_matrix(summary_df, genome_columns, threshold):
    """
    Create a binary presence/absence matrix.

    Args:
        summary_df: DataFrame with containment data
        genome_columns: List of genome column names
        threshold: Containment threshold

    Returns:
        DataFrame with binary presence/absence data
    """
    presence_df = pd.DataFrame(index=summary_df.index)

    for col in genome_columns:
        presence_df[col] = (summary_df[col] >= threshold).astype(int)

    return presence_df


def _calculate_jaccard_distance(vec1, vec2):
    """
    Calculate Jaccard distance between two vectors.

    Args:
        vec1: First binary vector
        vec2: Second binary vector

    Returns:
        Jaccard distance value
    """
    intersection = np.sum(np.logical_and(vec1, vec2))
    union = np.sum(np.logical_or(vec1, vec2))

    if union > 0:
        return 1.0 - (intersection / union)
    return 1.0


def _calculate_dice_distance(vec1, vec2):
    """
    Calculate Dice distance between two vectors.

    Args:
        vec1: First binary vector
        vec2: Second binary vector

    Returns:
        Dice distance value
    """
    intersection = np.sum(np.logical_and(vec1, vec2))
    sum1 = np.sum(vec1)
    sum2 = np.sum(vec2)

    if (sum1 + sum2) > 0:
        return 1.0 - (2.0 * intersection / (sum1 + sum2))
    return 1.0


def _calculate_cosine_distance(vec1, vec2):
    """
    Calculate Cosine distance between two vectors.

    Args:
        vec1: First binary vector
        vec2: Second binary vector

    Returns:
        Cosine distance value
    """
    dot_product = np.dot(vec1, vec2)
    norm1 = np.linalg.norm(vec1)
    norm2 = np.linalg.norm(vec2)

    if norm1 > 0 and norm2 > 0:
        return 1.0 - (dot_product / (norm1 * norm2))
    return 1.0


def calculate_distance_matrix(
    summary_file: Union[str, Path], threshold: float = 0.1, method: str = "jaccard"
) -> pd.DataFrame:
    """
    Calculate distance matrix between samples based on genome presence.

    Args:
        summary_file: Path to the summary file
        threshold: Minimum containment threshold
        method: Distance method ('jaccard', 'dice', or 'cosine')

    Returns:
        DataFrame with distance matrix

    Raises:
        ProcessingError: If the calculation fails
    """
    try:
        # Validate method
        _validate_distance_method(method)

        # Load summary data
        summary_df = pd.read_csv(summary_file, sep="\t", index_col=0)

        # Get genome columns
        genome_columns = get_genome_columns(summary_df)

        # Create binary presence/absence matrix
        presence_df = _create_binary_presence_matrix(summary_df, genome_columns, threshold)

        # Initialize distance matrix
        sample_count = len(presence_df)
        distance_matrix = pd.DataFrame(index=presence_df.index, columns=presence_df.index)

        # Calculate distances based on method
        distance_functions = {
            "jaccard": _calculate_jaccard_distance,
            "dice": _calculate_dice_distance,
            "cosine": _calculate_cosine_distance,
        }
        distance_function = distance_functions[method]

        # Calculate distances
        for i, sample1 in enumerate(presence_df.index):
            # Self-distance is 0
            distance_matrix.loc[sample1, sample1] = 0.0

            # Calculate distance to other samples
            for sample2 in presence_df.index[i + 1 :]:
                # Get binary vectors
                vec1 = presence_df.loc[sample1].values
                vec2 = presence_df.loc[sample2].values

                # Calculate distance based on method
                distance = distance_function(vec1, vec2)

                # Set distance (symmetric matrix)
                distance_matrix.loc[sample1, sample2] = distance
                distance_matrix.loc[sample2, sample1] = distance

        logger.info(f"Calculated {method} distances for {sample_count} samples")

        return distance_matrix

    except Exception as e:
        if isinstance(e, ProcessingError):
            raise
        raise ProcessingError(f"Error calculating distance matrix: {e}")


def perform_hypergeometric_test(
    success_in_sample: int,
    sample_size: int,
    success_in_population: int,
    population_size: int,
) -> Tuple[float, float]:
    """
    Perform hypergeometric test for enrichment.

    Args:
        success_in_sample: Number of successes in the sample
        sample_size: Size of the sample
        success_in_population: Number of successes in the population
        population_size: Size of the population

    Returns:
        Tuple of (fold_enrichment, p_value)

    Raises:
        ProcessingError: If the calculation fails
        ConfigurationError: If scipy is not installed
    """
    stats = require("scipy.stats", "analysis", "A hypergeometric test")
    try:
        # Validate inputs
        if sample_size <= 0 or population_size <= 0:
            raise ProcessingError("Sample and population sizes must be positive")

        if success_in_sample > sample_size:
            raise ProcessingError("Success count cannot exceed sample size")

        if success_in_population > population_size:
            raise ProcessingError("Success count cannot exceed population size")

        # Calculate expected number of successes
        expected = (success_in_population / population_size) * sample_size

        # Calculate fold enrichment
        if expected > 0:
            fold_enrichment = success_in_sample / expected
        else:
            fold_enrichment = float("inf") if success_in_sample > 0 else 1.0

        # Calculate p-value
        # Probability of observing success_in_sample OR MORE successes
        p_value = stats.hypergeom.sf(
            success_in_sample - 1,  # -1 because sf gives P(X > k)
            population_size,
            success_in_population,
            sample_size,
        )

        return fold_enrichment, p_value

    except Exception as e:
        if isinstance(e, ProcessingError):
            raise
        raise ProcessingError(f"Error performing hypergeometric test: {e}")


def compare_group_means(groups: Sequence[np.ndarray]) -> Tuple[str, float, float]:
    """
    Compare the means of two or more groups: Student's t-test for two groups, one-way ANOVA for more.

    Args:
        groups: One array of observations per group (at least two)

    Returns:
        Tuple of (test name, statistic, p-value), the test name being "t-test" or "ANOVA"

    Raises:
        ConfigurationError: If scipy is not installed
    """
    stats = require("scipy.stats", "analysis", "Comparing dataset groups")
    if len(groups) == 2:
        statistic, p_value = stats.ttest_ind(groups[0], groups[1])
        return "t-test", statistic, p_value
    statistic, p_value = stats.f_oneway(*groups)
    return "ANOVA", statistic, p_value


def holm_adjusted(p_values: Sequence[float]) -> List[float]:
    """Holm step-down adjusted p-values, in the order given, each capped at 1.

    The smallest raw p-value is multiplied by the number of tests, the next by one less, and
    so on; each adjusted value is at least the one before it in that order.
    """
    order = sorted(range(len(p_values)), key=lambda i: p_values[i])
    adjusted = [0.0] * len(p_values)
    running = 0.0
    for rank, index in enumerate(order):
        running = max(running, (len(p_values) - rank) * p_values[index])
        adjusted[index] = min(1.0, running)
    return adjusted
