"""Installation instructions for enrich_pvalues."""
import os
from setuptools import setup

import enrich_pvalues  # For version

VERSION=enrich_pvalues.__version__
GITHUB='https://github.com/MikeDacre/enrich_pvalues'

with open('requirements.txt') as fin:
    REQUIREMENTS = [
        i[0] for i in [j.split('>=') for j in fin.read().strip().split('\n')]
    ]


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

    py_modules=['enrich_pvalues'],

    entry_points = {
        'console_scripts': [
            'enrich_pvalues = enrich_pvalues:main',
        ],
    },

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 4 - Beta',
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
