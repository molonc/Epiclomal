from setuptools import setup, find_packages

setup(
    name='Epiclomal',
    packages=find_packages(),
    version='0.0.1',
    description='Epiclomal package, software for clustering and cluster evaluation of sparse DNA methylation data',
    url='https://github.com/shahcompbio/Epiclomal',
    entry_points={'console_scripts': ['epiclomal = epiclomal.epiclomal_run:main', 'evaluate_clustering = epiclomal.evaluate_clustering:main']},
    install_requires=['numba', 'numpy', 'scikit-learn', 'pandas', 'scipy'],
    zip_safe=False,
    )