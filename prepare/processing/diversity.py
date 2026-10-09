"""
Diversity metrics and analysis for PREPARE.

This module provides functions for calculating alpha and beta diversity metrics,
essential for microbiome and metagenomic analysis.
"""

import logging
import numpy as np
import pandas as pd
from typing import Any, Dict, List, Tuple, Union

from prepare.core.exceptions import ProcessingError
from prepare.core.optional import require

logger = logging.getLogger(__name__)


def calculate_alpha_diversity(
    abundance_matrix: Union[pd.DataFrame, np.ndarray],
    metrics: List[str] = ["shannon", "simpson", "chao1", "observed_species"],
) -> pd.DataFrame:
    """
    Calculate alpha diversity metrics for each sample.

    Alpha diversity measures the diversity within individual samples.

    Args:
        abundance_matrix: Sample x Species matrix of abundances
        metrics: List of metrics to calculate

    Returns:
        DataFrame with samples as rows and diversity metrics as columns

    Raises:
        ProcessingError: If calculation fails
    """
    try:
        if isinstance(abundance_matrix, pd.DataFrame):
            data = abundance_matrix.values
            sample_names = abundance_matrix.index
        else:
            data = abundance_matrix
            sample_names = [f"Sample_{i}" for i in range(data.shape[0])]

        # Handle empty data
        if data.size == 0:
            return pd.DataFrame()

        results = {}

        for metric in metrics:
            if metric == "shannon":
                results[metric] = _calculate_shannon_diversity(data)
            elif metric == "simpson":
                results[metric] = _calculate_simpson_diversity(data)
            elif metric == "chao1":
                results[metric] = _calculate_chao1_diversity(data)
            elif metric == "observed_species":
                results[metric] = _calculate_observed_species(data)
            elif metric == "pielou_evenness":
                results[metric] = _calculate_pielou_evenness(data)
            else:
                logger.warning(f"Unknown metric: {metric}")
                continue

        result_df = pd.DataFrame(results, index=sample_names)
        logger.info(f"Calculated alpha diversity for {len(sample_names)} samples")
        return result_df

    except Exception as e:
        raise ProcessingError(f"Failed to calculate alpha diversity: {e}")


def calculate_beta_diversity(
    abundance_matrix: Union[pd.DataFrame, np.ndarray],
    metric: str = "bray_curtis",
    return_dataframe: bool = True,
) -> Union[pd.DataFrame, np.ndarray]:
    """
    Calculate beta diversity (distance/dissimilarity) between samples.

    Beta diversity measures the differences between samples.

    Args:
        abundance_matrix: Sample x Species matrix of abundances
        metric: Distance metric to use
        return_dataframe: Whether to return as DataFrame with sample labels

    Returns:
        Distance matrix as DataFrame or numpy array

    Raises:
        ProcessingError: If calculation fails
    """
    # Bray-Curtis and Jaccard are computed here; every other metric comes from scikit-learn.
    sklearn_metrics = None
    if metric not in ("bray_curtis", "jaccard"):
        sklearn_metrics = require("sklearn.metrics", "analysis", f"The {metric} beta diversity metric")
    try:
        if isinstance(abundance_matrix, pd.DataFrame):
            data = abundance_matrix.values
            sample_names = abundance_matrix.index
        else:
            data = abundance_matrix
            sample_names = [f"Sample_{i}" for i in range(data.shape[0])]

        # Calculate distance matrix
        if metric == "bray_curtis":
            distances = _calculate_bray_curtis_distance(data)
        elif metric == "jaccard":
            distances = _calculate_jaccard_distance(data)
        else:
            assert sklearn_metrics is not None
            distances = sklearn_metrics.pairwise_distances(data, metric=metric)

        if return_dataframe:
            distance_df = pd.DataFrame(distances, index=sample_names, columns=sample_names)
            logger.info(f"Calculated {metric} beta diversity for {len(sample_names)} samples")
            return distance_df
        else:
            return distances

    except Exception as e:
        raise ProcessingError(f"Failed to calculate beta diversity: {e}")


def perform_permanova(
    distance_matrix: Union[pd.DataFrame, np.ndarray],
    metadata: pd.DataFrame,
    formula: str,
    n_permutations: int = 999,
) -> Dict[str, Dict[str, Any]]:
    """
    Perform PERMANOVA (Permutational Multivariate Analysis of Variance).

    Tests whether group differences are statistically significant.

    Args:
        distance_matrix: Beta diversity distance matrix
        metadata: Sample metadata DataFrame
        formula: Formula string (e.g., "treatment + site")
        n_permutations: Number of permutations for p-value calculation

    Returns:
        Dictionary with test statistics and p-values

    Raises:
        ProcessingError: If test fails
    """
    try:
        if isinstance(distance_matrix, pd.DataFrame):
            distances = distance_matrix.values
            sample_names = distance_matrix.index
        else:
            distances = distance_matrix
            sample_names = metadata.index

        # Parse formula to extract variables
        variables = [var.strip() for var in formula.replace("+", " ").split()]

        # Check that all variables exist in metadata
        missing_vars = [var for var in variables if var not in metadata.columns]
        if missing_vars:
            raise ProcessingError(f"Variables not found in metadata: {missing_vars}")

        # Align metadata with distance matrix
        metadata_aligned = metadata.loc[sample_names]

        # Perform PERMANOVA
        results = {}
        for variable in variables:
            groups = metadata_aligned[variable]
            f_stat, p_value = _permanova_test(distances, groups, n_permutations)
            results[variable] = {
                "F_statistic": f_stat,
                "p_value": p_value,
                "significant": bool(p_value < 0.05),
            }

        logger.info(f"PERMANOVA completed for variables: {variables}")
        return results

    except Exception as e:
        raise ProcessingError(f"Failed to perform PERMANOVA: {e}")


# Private helper functions


def _calculate_shannon_diversity(data: np.ndarray) -> np.ndarray:
    """Calculate Shannon diversity index."""
    # Convert to float for calculations
    data_float = data.astype(float)
    row_sums = data_float.sum(axis=1, keepdims=True)

    # Handle zero-sum rows (all zeros)
    zero_rows = (row_sums == 0).flatten()

    # Normalize to relative abundances
    rel_abundance = np.divide(data_float, row_sums, out=np.zeros_like(data_float), where=row_sums != 0)
    rel_abundance = rel_abundance + 1e-12  # Add small constant to avoid log(0)

    # Shannon index: H = -sum(pi * ln(pi))
    shannon = -np.sum(rel_abundance * np.log(rel_abundance), axis=1)

    # Set Shannon to 0 for zero-abundance samples
    shannon[zero_rows] = 0.0

    return shannon


def _calculate_simpson_diversity(data: np.ndarray) -> np.ndarray:
    """Calculate Simpson diversity index (1 - dominance)."""
    # Convert to float for calculations
    data_float = data.astype(float)
    row_sums = data_float.sum(axis=1, keepdims=True)

    # Handle zero-sum rows (all zeros)
    zero_rows = (row_sums == 0).flatten()

    # Normalize to relative abundances
    rel_abundance = np.divide(data_float, row_sums, out=np.zeros_like(data_float), where=row_sums != 0)

    # Simpson index: D = 1 - sum(pi^2)
    simpson = 1 - np.sum(rel_abundance**2, axis=1)

    # Set Simpson to 0 for zero-abundance samples
    simpson[zero_rows] = 0.0

    return simpson


def _calculate_chao1_diversity(data: np.ndarray) -> np.ndarray:
    """Calculate Chao1 richness estimator."""
    chao1 = np.zeros(data.shape[0])

    for i in range(data.shape[0]):
        sample = data[i, :]
        observed = np.sum(sample > 0)
        singletons = np.sum(sample == 1)
        doubletons = np.sum(sample == 2)

        if doubletons > 0:
            chao1[i] = observed + (singletons**2) / (2 * doubletons)
        else:
            chao1[i] = observed + singletons * (singletons - 1) / 2

    return chao1


def _calculate_observed_species(data: np.ndarray) -> np.ndarray:
    """Calculate observed species richness."""
    return np.sum(data > 0, axis=1)


def _calculate_pielou_evenness(data: np.ndarray) -> np.ndarray:
    """Calculate Pielou's evenness index."""
    shannon = _calculate_shannon_diversity(data)
    observed = _calculate_observed_species(data)

    # Pielou's J = H / ln(S)
    evenness = shannon / np.log(observed)
    evenness[observed <= 1] = 0  # Handle cases with <=1 species

    return evenness


def _calculate_bray_curtis_distance(data: np.ndarray) -> np.ndarray:
    """Calculate Bray-Curtis dissimilarity."""
    n_samples = data.shape[0]
    distances = np.zeros((n_samples, n_samples))

    for i in range(n_samples):
        for j in range(i + 1, n_samples):
            numerator = np.sum(np.abs(data[i, :] - data[j, :]))
            denominator = np.sum(data[i, :] + data[j, :])

            if denominator > 0:
                distances[i, j] = numerator / denominator
            else:
                distances[i, j] = 0

            distances[j, i] = distances[i, j]

    return distances


def _calculate_jaccard_distance(data: np.ndarray) -> np.ndarray:
    """Calculate Jaccard dissimilarity (presence/absence)."""
    # Convert to binary
    binary_data = (data > 0).astype(int)

    n_samples = data.shape[0]
    distances = np.zeros((n_samples, n_samples))

    for i in range(n_samples):
        for j in range(i + 1, n_samples):
            intersection = np.sum(binary_data[i, :] & binary_data[j, :])
            union = np.sum(binary_data[i, :] | binary_data[j, :])

            if union > 0:
                distances[i, j] = 1 - (intersection / union)
            else:
                distances[i, j] = 0

            distances[j, i] = distances[i, j]

    return distances


def _permanova_test(distance_matrix: np.ndarray, groups: pd.Series, n_permutations: int) -> Tuple[float, float]:
    """Perform PERMANOVA test."""
    n_samples = distance_matrix.shape[0]

    # Calculate total sum of squares
    upper_triangle_distances = distance_matrix[np.triu_indices(n_samples, k=1)] ** 2
    if len(upper_triangle_distances) > 0:
        grand_mean = np.mean(upper_triangle_distances)
    else:
        grand_mean = 0
    total_ss = np.sum(distance_matrix**2) / (2 * n_samples) - grand_mean

    # Calculate within-group sum of squares
    within_ss = 0
    unique_groups = groups.unique()

    for group in unique_groups:
        group_indices = np.where(groups == group)[0]
        if len(group_indices) > 1:
            group_distances = distance_matrix[np.ix_(group_indices, group_indices)]
            group_upper_triangle = group_distances[np.triu_indices(len(group_indices), k=1)] ** 2
            if len(group_upper_triangle) > 0:
                group_mean = np.mean(group_upper_triangle)
            else:
                group_mean = 0
            within_ss += np.sum(group_distances**2) / (2 * len(group_indices)) - group_mean

    # Calculate between-group sum of squares
    between_ss = total_ss - within_ss

    # Calculate F-statistic
    df_between = len(unique_groups) - 1
    df_within = n_samples - len(unique_groups)

    if df_within > 0 and within_ss > 0:
        f_stat = (between_ss / df_between) / (within_ss / df_within)
    else:
        f_stat = 0

    # Permutation test for p-value
    permuted_f_stats = []
    for _ in range(n_permutations):
        permuted_groups = groups.sample(frac=1).values
        perm_f, _ = _permanova_test(distance_matrix, pd.Series(permuted_groups), 0)
        permuted_f_stats.append(perm_f)

    if len(permuted_f_stats) > 0:
        p_value = np.mean(np.array(permuted_f_stats) >= f_stat)
    else:
        p_value = np.nan  # Return NaN when no permutations were performed

    return f_stat, p_value
