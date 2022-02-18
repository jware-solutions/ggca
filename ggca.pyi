from typing import Optional, Tuple, List


class CorResult:
    gene: str
    gem: str
    cpg_site_id: Optional[str]
    correlation: Optional[float]
    p_value: Optional[float]
    adjusted_p_value: Optional[float]

    def __init__(
		self,
		gene: str,
		gem: str,
		cpg_site_id: Optional[str],
		correlation: Optional[float],
		p_value: Optional[float],
		adjusted_p_value: Optional[float],
	) -> None:
		"""
		Represents a correlation analysis result. Includes Gene, GEM, CpG Site ID (if specified) correlation statistic, p-value and adjusted p-value.

		:param gene: Gene name
		:param gem: Gene Expression Modulator (GEM) name
		:param cpg_site_id: CpG Site ID
		:param correlation: Correlation statistic (Pearson, Spearman or Kendall, as selected)
		:param p_value: P-value
		:param adjusted_p_value: Adjusted p-value (Benjamini-Hochberg, Benjamini-Yekutieli or Bonferroni, as selected)
		"""
		self.gene = gene
		self.gem = gem
		self.cpg_site_id = cpg_site_id
		self.correlation = correlation
		self.p_value = p_value
		self.adjusted_p_value = adjusted_p_value


def correlate(
	gene_file_path: str,
    gem_file_path: str,
    correlation_method: int,
    correlation_threshold: float,
    sort_buf_size: int,
    adjustment_method: int,
    all_vs_all: bool,
    gem_contains_cpg: bool,
    collect_gem_dataset: Optional[bool],
    keep_top_n: Optional[int],
) -> Tuple[List[CorResult], int]:
	"""
	Computes the correlation between both mRNA and GEM files' rows.

	:param gene_file_path: Gene file's path
	:param gem_file_path: Gene Expression Modulator (GEM) file's path
	:param correlation_method: Correlation method to compute (Spearman = 1, Kendall = 2 or Pearson = 3)
	:param correlation_threshold: The threshold to discard all results whose correlation statistic values are below this value
	:param sort_buf_size: Number of elements to sort by block in disk during p-value adjustment process. Greater blocks are faster but consume more memory
	:param adjustment_method: P-value adjustment method (Benjamini-Hochberg = 1, Benjamini-Yekutieli = 2 or Bonferroni = 3)
	:param all_vs_all:  True if all Genes must be evaluated with all GEMs. Otherwise, only matching Genes/GEM will be evaluated (useful for CNA or Methylation analysis)
	:param gem_contains_cpg: Set to True if your GEM data contains CpG Site IDs as the second column to preserve the GEM/CpG Site reference
	:param collect_gem_dataset: True to make the GEM dataset available in memory. This has a HUGE impact in analysis performance. Specify a boolean value to force or use None to allocate in memory automatically when GEM dataset size is small (<= 100MB)
	:param keep_top_n: Specify a number of results to keep or None to return all the resulting combinations
	:return: A tuple with a vec of CorResult and the number of combinations evaluated
	"""
	...


class GGCAError(Exception):
	"""Raises when a general error occurs, such as a read error, file does not exist, among others."""
	...


class GGCADiffSamplesLength(Exception):
	"""Raises when the length of samples in both datasets are different."""
	...


class GGCADiffSamples(Exception):
	"""Raises when Samples in both datasets are different, but they have the same length (maybe they are in different order)."""
	...

class InvalidCorrelationMethod(Exception):
	"""Raises when an invalid correlation method is provided. Only values 1 (Spearman), 2 (Kendall) or 3 (Pearson) are valid."""
	...

class InvalidAdjustmentMethod(Exception):
	"""Raises when an invalid adjustment method is provided. Only values 1 (Benjamini-Hochberg), 2 (Benjamini-Yekutieli) or 3 (Bonferroni) are valid."""
	...