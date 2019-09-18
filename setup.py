from setuptools import setup, find_packages

setup(
    name='Epiclomal',
    packages=find_packages(where='epiclomal'),
    package_dir={
        '': 'epiclomal',
    },
    version='0.0.1',
    description='Epiclomal package, software for clustering and cluster evaluation of sparse DNA methylation data',
    url='https://github.com/shahcompbio/Epiclomal',
    entry_points={'console_scripts': ['epiclomal = run:main', 'evaluate_clustering = evaluate_clustering:main']},
    )