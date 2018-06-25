import qiime2.plugin
from qiime2.plugin import (SemanticType, Str, Int, Float, Choices,
                           MetadataColumn, Metadata, Plugin)
from q2_types.feature_table import (
    FeatureTable, Frequency)
from q2_types.sample_data import AlphaDiversity, SampleData
from dsfdr import dsfdr
import pandas as pd
import numpy as np

_citation = ('Jiang L, Amir A, Morton JT, Heller R, Arias-Castro E, Knight R. 2017. '
             'Discrete False-Discovery Rate Improves Identification of Differentially Abundant Microbes'
             'mSystems.00092-17 '
             'https://doi.org/10.1128/mSystems.00092-17')

_short_description = "Plugin for multiple comparisons in sparse Microbiome Data"

plugin = qiime2.plugin.Plugin(
    name='dsfdr',
    version="0.0.2",
    website='https://github.com/serenejiang/q2_dsfdr',
    package='q2_dsfdr',
    short_description=_short_description,
    description=('This is a QIIME 2 plugin supporting multiple comparisons on sparse microbiome feature tables and metadata'),
    citation_text=_citation
)


def permutation_fdr(table: pd.DataFrame,
                    metadata: qiime2.MetadataColumn,
                    statistical_test: str = 'meandiff',
                    transform_function: str = 'rank',
                    alpha: float = 0.05,
                    permutations: int=1000) -> pd.Series:

        metadata_series = metadata.to_series()[table.index]
        uvals = metadata_series.unique()
        if len(uvals) < 2:
            raise ValueError('Only one value in mapping file data column (%s). Aborting' % uvals[0])
        if len(uvals) > 2:
            raise ValueError('More than two values in mapping file data column (%s). Aborting' % uvals)
        # labels need to be 0/1
        labels = np.zeros(len(metadata_series))
        labels[metadata_series == uvals[1]] = 1

        # allow debug info. q2 takes care of what to show using the --verbose flag
        try:
            dsfdr.logger.setLevel(1)
        except:
            pass

        ret_reject, ret_tstat, ret_pvals = dsfdr.dsfdr(table.values.T,
               labels,
               transform_function,
               statistical_test,
               alpha, permutations)
        return pd.Series(ret_reject, index=table.columns)


_statistical_tests = ['meandiff', 'mannwhiteny', 'kruwallis', 'stdmeandiff',
                      'spearman', 'pearson', 'nonzerospearman', 'nonzeropearson']

_transform_functions = ['rank', 'log', 'pa', 'norm']


plugin.methods.register_function(
    function=permutation_fdr,
    inputs={'table': FeatureTable[Frequency]},
    outputs=[('reject', SampleData[AlphaDiversity])],
    parameters={
        'metadata': Metadata,
        'statistical_test': Str % Choices(_statistical_tests),
        'transform_function': Str % Choices(_transform_functions),
        'permutations': Int,
        'alpha': Float,
    },
    name='Discrete FDR',
    description=("Discrete FDR")
)
