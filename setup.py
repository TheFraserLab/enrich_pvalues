"""Installation instructions for enrich_pvalues."""
import os
from setuptools import setup

import enrich_pvalues  # For version

VERSION=enrich_pvalues.__version__
GITHUB='https://github.com/MikeDacre/enrich_pvalues'

REQUIREMENTS = ['tabulate', 'tqdm', 'numpy', 'pandas', 'matplotlib', 'seaborn']


def read(fname):
    """Read the contents of a file in this dir."""
    with open(os.path.join(os.path.dirname(__file__), fname)) as fin:
        return fin.read()


# Actual setup instructions
setup(
    name         = 'enrich_pvalues',
    version      = VERSION,
    author       = 'Mike Dacre',
    author_email = 'mike.dacre@gmail.com',
    description  = (
        "Compare one dataset to another at a variety of p-value cutoffs"
    ),
    keywords = (
        "statistics p-values biology molecular-biology console"
    ),
    long_description = read('README.rst'),
    license = 'MIT',

    # URLs
    url = GITHUB,
    download_url='{0}/archive/v{1}.tar.gz'.format(GITHUB, VERSION),

    # Actual packages/modules
    # Packages are directories with __init__.py files
    # Modules are python scripts (minus the .py)
    #  packages=['enrich_pvalues'],
    py_modules=['enrich_pvalues'],

    # Entry points and scripts
    # Entry points are functions that can use sys.argv (e.g. main())
    # Scripts are independent pieces of code intended to be executed as is
    entry_points = {
        'console_scripts': [
            'enrich_pvalues = enrich_pvalues:main',
        ],
    },
    #  scripts = [],

    # Data files
    #  data_files = [],

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 1 - Planning',
        # 'Development Status :: 4 - Beta',
        # 'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: System Administrators',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Topic :: Utilities',
    ],

    # Requirements
    requires=REQUIREMENTS,
    install_requires=REQUIREMENTS
)
