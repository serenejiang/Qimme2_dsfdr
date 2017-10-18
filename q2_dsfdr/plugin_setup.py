import qiime2.plugin
from qiime2.plugin import (SemanticType, Str, Int, Float, Choices,
                          MetadataCategory, Plugin)
from q2_types.feature_table import (
    FeatureTable, Frequency)
#from q2_pfdr._pfdr import permutation_fdr
from q2_types.sample_data import AlphaDiversity, SampleData
from dsfdr import dsfdr
import pandas as pd


plugin = qiime2.plugin.Plugin(
    name='dsfdr',
    version="0.0.1",
    website='Website for q2-pfdr',
    package='q2_dsfdr',
    # Information on how to obtain user support should be provided as a free
    # text string via user_support_text. If None is provided, users will
    # be referred to the plugin's website for support.
    user_support_text=None,
    # Information on how the plugin should be cited should be provided as a
    # free text string via citation_text. If None is provided, users
    # will be told to use the plugin's website as a citation.
    citation_text=None
)


def permutation_fdr(table : pd.DataFrame,
                    metadata: qiime2.MetadataCategory,
                    statistical_test: str = 'meandiff',
                    transform_function: str = 'log',
                    alpha: float = 0.05,
                    permutations: int=1000) -> pd.Series:
        # See q2-composition for more details
        # https://github.com/qiime2/q2-composition/blob/master/q2_composition/_ancom.py

        # TODO : Consider renaming the functions to match q2-composition

        metadata_series = metadata.to_series()[table.index]
        # Make sure that metadata and table match up
        reject_idx = dsfdr(table.values.T,
			   metadata_series.values,
			   statistical_test,
			   transform_function,
			   alpha, permutations)
        return reject_idx


_statistical_tests = ['meandiff', 'mannwhiteny', 'kruwallis', 'stdmeandiff',
                      'spearman', 'pearson', 'nonzerospearman', 'nonzeropearson']

_transform_functions = ['log', 'rank', 'pa', 'norm']

# TODO: Pass in a RelativeFrequency, PresenceAbsence, Composition type?
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
    name='Permutation FDR',
    description=("Permutation FDR")
)
