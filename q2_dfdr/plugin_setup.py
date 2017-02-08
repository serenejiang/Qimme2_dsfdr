import qiime.plugin
from qiime.plugin import (SemanticType, Str, Int, Float, Choices,
                          MetadataCategory, Plugin)
from q2_types.feature_table import (
    FeatureTable, Frequency)
import q2_pfdr


plugin = qiime.plugin.Plugin(
    name='pfdr',
    version=q2_pfdr.__version__,
    website='Website for q2-pfdr',
    package='q2_pfdr',
    # Information on how to obtain user support should be provided as a free
    # text string via user_support_text. If None is provided, users will
    # be referred to the plugin's website for support.
    user_support_text=None,
    # Information on how the plugin should be cited should be provided as a
    # free text string via citation_text. If None is provided, users
    # will be told to use the plugin's website as a citation.
    citation_text=None
)

_statistical_tests = q2_pfdr._pfdr.statistical_tests()
_transform_functions = q2_pfdr._pfdr.transform_functions()

# TODO: Pass in a RelativeFrequency, PresenceAbsence, Composition type?
plugin.visualizers.register_function(
    function=q2_pfdr.permutation_fdr,
    inputs={'table': FeatureTable[Frequency]},
    parameters={
        'metadata': MetadataCategory,
        'statistical_test': Str % Choices(_statistical_tests),
        'transform_function': Str % Choices(_transform_functions),
        'permutations': Int
    },
    name='Apply blah ...',
    description=("Apply blah.")
)
