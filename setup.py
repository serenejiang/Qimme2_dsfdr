from setuptools import setup, find_packages
import re
import ast

setup(
    name="q2-dsfdr",
    version="0.0.2",
    packages=find_packages(),
    author="Serene Jiang",
    author_email="serene1030@gmail.com",
    description="Description of q2-dsfdr",
    entry_points={
        "qiime2.plugins":
        ["q2-dsfdr=q2_dsfdr.plugin_setup:plugin"]
    }
)
