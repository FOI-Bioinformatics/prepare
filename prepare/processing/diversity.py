"""
Diversity metrics and analysis for PREPARE.

This module provides functions for calculating alpha and beta diversity metrics,
essential for microbiome and metagenomic analysis.
"""

import logging
import numpy as np
import pandas as pd
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

from prepare.core.exceptions import ProcessingError, ValidationError
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


def check_metadata_covers_samples(sample_names: Sequence[Any], metadata: pd.DataFrame) -> None:
    """Raise ValidationError naming every sample that has no row in the metadata."""
    present = set(metadata.index)
    missing = [str(name) for name in sample_names if name not in present]
    if missing:
        raise ValidationError(f"Samples missing from the metadata: {', '.join(missing)}")


def perform_permanova(
    distance_matrix: Union[pd.DataFrame, np.ndarray],
    metadata: pd.DataFrame,
    formula: str,
    n_permutations: int = 999,
    seed: Optional[int] = None,
) -> Dict[str, Dict[str, Any]]:
    """
    Perform PERMANOVA (Permutational Multivariate Analysis of Variance).

    Tests whether group differences are statistically significant.

    Args:
        distance_matrix: Beta diversity distance matrix
        metadata: Sample metadata DataFrame
        formula: Formula string (e.g., "treatment + site")
        n_permutations: Number of permutations for p-value calculation
        seed: Seed for the permutation generator; None draws fresh entropy

    Returns:
        Dictionary with test statistics and p-values; each variable is tested on its own

    Raises:
        ValidationError: If a sample in the distance matrix is absent from the metadata
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
        check_metadata_covers_samples(sample_names, metadata)
        # A sample listed more than once is read from its first row, as the plots do.
        unique = metadata[~metadata.index.duplicated(keep="first")]
        metadata_aligned = unique.loc[sample_names]

        # Perform PERMANOVA
        rng = np.random.default_rng(seed)
        results = {}
        for variable in variables:
            groups = metadata_aligned[variable]
            f_stat, p_value = _permanova_test(distances, groups, n_permutations, rng=rng)
            results[variable] = {
                "F_statistic": f_stat,
                "p_value": p_value,
                "significant": bool(p_value < 0.05),
            }

        logger.info(f"PERMANOVA completed for variables: {variables}")
        return results

    except ValidationError:
        raise
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


def _permanova_f(distance_matrix: np.ndarray, labels: np.ndarray) -> float:
    """Pseudo-F of Anderson (2001) for one labelling; nan when it is undefined."""
    n_samples = distance_matrix.shape[0]
    unique_groups = np.unique(labels)
    n_groups = len(unique_groups)
    if n_groups < 2 or n_samples - n_groups < 1:
        return float("nan")

    squared = distance_matrix**2
    upper = np.triu_indices(n_samples, k=1)
    total_ss = squared[upper].sum() / n_samples

    within_ss = 0.0
    for group in unique_groups:
        members = np.flatnonzero(labels == group)
        block = squared[np.ix_(members, members)]
        within_ss += block[np.triu_indices(len(members), k=1)].sum() / len(members)

    between_ss = total_ss - within_ss
    if within_ss <= 0:
        return float("inf") if between_ss > 0 else float("nan")
    return float((between_ss / (n_groups - 1)) / (within_ss / (n_samples - n_groups)))


def _permanova_test(
    distance_matrix: np.ndarray,
    groups: pd.Series,
    n_permutations: int,
    rng: Optional[np.random.Generator] = None,
) -> Tuple[float, float]:
    """One-factor PERMANOVA (Anderson 2001).

    Samples whose group label is missing are removed with their rows and columns.
    SS_T = sum_{i<j} d_ij^2 / N, SS_W = sum_g sum_{i<j in g} d_ij^2 / n_g, SS_A = SS_T - SS_W,
    F = (SS_A / (a - 1)) / (SS_W / (N - a)), and the permutation p-value is
    (1 + #{F_perm >= F}) / (1 + n_permutations). F and p are nan when fewer than two
    groups remain or N - a < 1; p is nan when no permutations are requested.
    """
    if rng is None:
        rng = np.random.default_rng()

    keep = ~pd.isna(np.asarray(groups, dtype=object))
    labels = np.asarray(groups, dtype=object)[keep].astype(str)
    distances = np.asarray(distance_matrix, dtype=float)[np.ix_(keep, keep)]

    f_stat = _permanova_f(distances, labels)
    if np.isnan(f_stat) or n_permutations < 1:
        return f_stat, float("nan")

    n_extreme = 0
    for _ in range(n_permutations):
        if _permanova_f(distances, rng.permutation(labels)) >= f_stat:
            n_extreme += 1
    p_value = (1 + n_extreme) / (1 + n_permutations)

    return f_stat, float(p_value)
