from setuptools import setup, find_packages
import re
import ast

setup(
    name="q2-dsfdr",
    version="0.0.1",
    packages=find_packages(),
    # pandas and q2-dummy-types are only required for the dummy methods and
    # visualizers provided as examples. Remove these dependencies when you're
    # ready to develop your plugin, and add your own dependencies (if there are
    # any).
    author="Serene Jiang",
    author_email="serene1030@gmail.com",
    description="Description of q2-pfdr",
    entry_points={
        "qiime2.plugins":
        ["q2-dsfdr=q2_dsfdr.plugin_setup:plugin"]
    }
)
