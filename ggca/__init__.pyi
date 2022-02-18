from enum import Enum
from typing import Optional, Tuple, List


class CorrelationMethod(Enum):
	SPEARMAN = 1
	KENDALL = 2
	PEARSON = 3


class AdjustmentMethod(Enum):
	BENJAMINI_HOCHBERG = 1
	BENJAMINI_YEKUTIELI = 2
	BONFERRONI = 3


class CorResult:
    gene: str
    gem: str
    cpg_site_id: Optional[str]
    correlation: Optional[float]
    p_value: Optional[float]
    adjusted_p_value: Optional[float]


def correlate(
	file_1_path: str,
    file_2_path: str,
    correlation_method: CorrelationMethod,
    correlation_threshold: float,
    sort_buf_size: int,
    adjustment_method: AdjustmentMethod,
    all_vs_all: bool,
    gem_contains_cpg: bool,
    collect_gem_dataset: Optional[bool],
    keep_top_n: Optional[int],
) -> Tuple[List[CorResult], int]: ...
