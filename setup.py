from pathlib import Path
from setuptools import setup

local_path = Path(__file__).resolve().parent

setup(
    install_requires=[
        f"rna_moccasin @ {local_path}/moccasin",
        f"rna_majiq @ {local_path}/majiq",
        f"rna_voila @ {local_path}/voila",
    ]
)