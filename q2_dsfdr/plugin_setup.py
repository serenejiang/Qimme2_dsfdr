import qiime2.plugin
from qiime2.plugin import (SemanticType, Str, Int, Float, Choices,
                          MetadataCategory, Plugin)
from q2_types.feature_table import (
    FeatureTable, Frequency)
from q2_types.sample_data import AlphaDiversity, SampleData
from dsfdr import dsfdr
import pandas as pd


plugin = qiime2.plugin.Plugin(
    name='dsfdr',
    version="0.0.1",
    website='Website for q2-dsfdr',
    package='q2_dsfdr',
    user_support_text=None,
    citation_text=None
)


def permutation_fdr(table : pd.DataFrame,
                    metadata: qiime2.MetadataCategory,
                    statistical_test: str = 'meandiff',
                    transform_function: str = 'log',
                    alpha: float = 0.05,
                    permutations: int=1000) -> pd.Series:

        metadata_series = metadata.to_series()[table.index]
        reject_idx = dsfdr(table.values.T,
			   metadata_series.values,
			   statistical_test,
			   transform_function,
			   alpha, permutations)
        return reject_idx


_statistical_tests = ['meandiff', 'mannwhiteny', 'kruwallis', 'stdmeandiff',
                      'spearman', 'pearson', 'nonzerospearman', 'nonzeropearson']

_transform_functions = ['log', 'rank', 'pa', 'norm']


plugin.methods.register_function(
    function=permutation_fdr,
    inputs={'table': FeatureTable[Frequency]},
    outputs=[('reject', SampleData[AlphaDiversity])],
    parameters={
        'metadata': MetadataCategory,
        'statistical_test': Str % Choices(_statistical_tests),
        'transform_function': Str % Choices(_transform_functions),
        'permutations': Int,
        'alpha': Float,
    },
    name='Discrete FDR',
    description=("Discrete FDR")
)
