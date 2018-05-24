#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compare one dataset to another at a variety of p-value cutoffs.
"""
import os as _os
import sys as _sys
import bz2 as _bz2
import gzip as _gzip
import argparse as _argparse


__version__ = '0.1'


def open_zipped(infile, mode='r'):
    """Return file handle of file regardless of compressed or not.

    Also returns already opened files unchanged, text mode automatic for
    compatibility with python2.
    """
    # return already open files
    if hasattr(infile, 'write'):
        return infile
    # make text mode automatic
    if len(mode) == 1:
        mode = mode + 't'
    # refuse to handle non-strings that aren't files.
    if not isinstance(infile, str):
        raise ValueError("I cannot open a filename that isn't a string.")
    # treat '-' appropriately
    if infile == '-':
        if 'w' in mode:
            return _sys.stdout
        return _sys.stdin
    # if possible open zipped files
    if infile.endswith('.gz'):
        return _gzip.open(infile, mode)
    if infile.endswith('.bz2'):
        if hasattr(_bz2, 'open'):
            return _bz2.open(infile, mode)
        return _bz2.bz2file(infile, mode)
    # fall back on regular open
    return open(infile, mode)


def main(argv=None):
    """Run as a script."""
    if not argv:
        argv = _sys.argv[1:]

    parser  = _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter
    )

    args = parser.parse_args(argv)

    print('Not implemented yet.')


if __name__ == '__main__' and '__file__' in globals():
    _sys.exit(main())
