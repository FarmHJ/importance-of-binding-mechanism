#
# setuptools script
#
from setuptools import setup, find_packages


def get_version():
    """
    Get version number from the modelling module.
    The easiest way would be to just ``import modelling ``, but note that this may  # noqa
    fail if the dependencies have not been installed yet. Instead, we've put
    the version number in a simple version_info module, that we'll import here
    by temporarily adding the oxrse directory to the pythonpath using sys.path.
    """
    import os
    import sys

    sys.path.append(os.path.abspath('modelling'))
    from version_info import VERSION as version
    sys.path.pop()

    return version


def get_readme():
    """
    Load README.md text for use as description.
    """
    with open('README.md') as f:
        return f.read()


setup(
    # Module name (lowercase)
    name='modelling',

    # Version
    version=get_version(),

    description='This is a repository for all the modelling in my DPhil project.',  # noqa

    long_description=get_readme(),

    license='BSD 3-Clause "New" or "Revised" License',

    # author='',

    # author_email='',

    maintainer='',

    maintainer_email='',

    url='https://github.com/FarmHJ/binding-mechanism.git',

    # Packages to include
    packages=find_packages(include=('modelling', 'modelling.*')),

    # List of dependencies
    install_requires=[
        # Dependencies go here!
        'matplotlib',
        'myokit',
        'numpy>=1.8',
        'pandas',
        'pints',
    ],
    extras_require={
        'docs': [
            # Sphinx for doc generation. Version 1.7.3 has a bug:
            'sphinx>=1.5, !=1.7.3',
            # Nice theme for docs
            'sphinx_rtd_theme',
        ],
        'dev': [
            # Flake8 for code style checking
            'flake8>=3',
        ],
    },
)