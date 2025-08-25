import setuptools

setuptools.setup(
    name="EpiVirQuant",
    version="0.1",
    author="Jose L. Figueroa III, Sadie M. Hollenack, Richard A. White III",
    description="EpiVirQuant directly counts and measures the size of viral-like particles, bacteria, archaea in epifluorescence microscopy images commonly used in viral ecology, viral estimation and enumeration, and microbial ecology.",
    url="https://github.com/raw-lab/epivirquant",
    scripts=["bin/epivirquant.py"],
    packages=['epivirquant_lib'],
    package_dir={'epivirquant_lib': 'lib'},
    python_requires='>=3.8',
    install_requires=[
        'joblib',
        'configargparse',
        'matplotlib',
        'numpy',
        'scipy',
        'scikit-image',
        'pypher',
        'stardist',
    ]
)