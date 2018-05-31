#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compare one dataset to another at a variety of p-value cutoffs.
"""
from __future__ import print_function
import os as _os
import sys as _sys
import bz2 as _bz2
import gzip as _gzip
import json as _json
import argparse as _argparse

import math

import numpy as np
import pandas as pd

from tabulate import tabulate as _tab
from tqdm import tqdm as _tqdm


__version__ = '1.0b2'


###############################################################################
#                          Core Enrichment Algorithm                          #
###############################################################################


def enrich_study(dataset, sig_comp, nsig_comp, simp_fl=False,
                 low_mem=False, conf=None):
    """Compute enrichment of significant data in sig_comp and nsig_comp.

    This is the core algorithm of this script.

    Read in all data from dataset and then take all names that beat a
    significance cutoff (set in conf) and compare to names in sig_comp and
    nsig_comp, computing an enrichment score (percentage in each, sig/nsig).

    Repeats this for every p-value cutoff between conf['max_pval'] (default
    0.05) and conf['min_pval'] (default 1e-15). We check each p-value cutoff
    at intervals of 1e-x and 5e-x for each exponent x between the max and min.

    Params
    ------
    dataset : str
        Path to data table to test, must contain names and p-values. Parse
        instructions in the config
    sig_comp : str or set
        Set of names that are in the significant set, or path to a newline
        separated set (made by split_study)
    nsig_comp : str or set
        Same as sig_comp, but for the non-significant comparison set
    simp_fl : bool, optional
        Treat the study_file as a two column file with no header, where the
        columns are name and pvalue, separated by a tab
    low_mem : bool, optional, not implemented
        Do not load the whole set of names and p-values at the same time,
        parse file over again for every comparison. Slower, but saves memory
        for very larde datasets.
    conf : dict, optional
        A config file to use for parsing.

    Returns
    -------
    [(p-value-cutoff), (enrichment-score)]
    """
    # Build a range of cutoffs to use
    max_exp = abs(math.floor(math.log10(conf['max_pval'])))
    min_exp = abs(math.floor(math.log10(conf['min_pval'])))
    cutoffs = []
    for exp in range(max_exp, min_exp+1):
        # Cumbersome way to handle floating point errors
        cutoffs.append(float('{0:1e}'.format(5*(10**-exp))))
        cutoffs.append(float('{0:1e}'.format(1*(10**-exp))))
    _sys.stderr.write('Testing cutoffs:\n{0}\n'.format(cutoffs))

    # Get the comparison sets
    _sys.stderr.write('Getting comparison data\n')
    scmp = get_set(sig_comp)
    ncmp = get_set(nsig_comp)
    del sig_comp
    del nsig_comp
    _sys.stderr.write(
        'Got {0} sig names and {1} nonsig names.\n'
        .format(len(scmp), len(ncmp))
    )

    # Get the right iterator
    if simp_fl:
        sit = simple_iterator(dataset)
    else:
        sit = study_iterator(
            dataset, conf['test_sep'], conf['test_has_header'],
            conf['test_name_col'], conf['test_pval_col']
        )

    # Keep only the two columns we care about, and only those that beat
    # the max p-value we will test
    data = []
    add_data = data.append
    max_cutoff = conf['max_pval']
    _sys.stderr.write(
        'Reading data file, keeping those less than P {0}\n'.format(max_cutoff)
    )
    for name, pval in _tqdm(sit, unit=' rows', disable=None):
        if pval <= max_cutoff:
            add_data([name, pval])

    # Make our dataframe
    _sys.stderr.write('Converting to DataFrame\n')
    data = pd.DataFrame(data, columns=['name', 'pval'])
    # Make unique, keep most significant
    data = data.sort_values('pval').drop_duplicates(['name', 'pval'])

    # Compute score
    scores = {}
    _sys.stderr.write('Calculating enrichments\n')
    for cutoff in _tqdm(cutoffs, unit=' cuttofs', disable=None):
        test_set = frozenset(data[data.pval <= cutoff].name)
        test_len = len(test_set)
        sigyl  = len(test_set & scmp)
        nsigyl = len(test_set & ncmp)
        sigy   = sigyl/test_len
        nsigy  = nsigyl/test_len
        if nsigy == 0:
            score = np.nan
        else:
            score = sigy/nsigy
        scores['{0:1e}'.format(cutoff)] = {
            'cutoff': cutoff,
            'sig_data': test_len,
            'sig_overlap': sigyl,
            'nonsig_overlap': nsigyl,
            'sig_score': sigy,
            'nonsig_score': nsigy,
            'enrichment_score': score
        }

    # Free memory
    _sys.stderr.write('Done. Clearing memory\n')
    del data, scmp, ncmp

    _sys.stderr.write('Making DataFrame\n')
    scores = pd.DataFrame.from_dict(scores, orient='index')
    scores = scores.sort_values('cutoff', ascending=False)

    print(scores)

    return scores


def plot_scores(scores, outfile=None, figsize=(14,10),
                comp_prefix=None, raw=False):
    """Plot enrichment score and sig count against cutoff.

    Enrichment scores end up on left y axis, total number of significant right
    y axis. x axis is the cutoffs.

    Params
    ------
    scores : DataFrame
        from enrich_study
    outfile : str, optional
        Path to write figure to
    figsize : tuple, optional
        Size of figure
    comp_prefix : str, optional
        The prefix of the comparison data to use for title creation
    raw : bool, optional
        Plot raw counts instead of percentages

    Returns
    -------
    fig, [ax1, ax2]
    """
    from matplotlib import pyplot as plt
    from matplotlib.ticker import FuncFormatter
    import seaborn as sns

    if comp_prefix:
        title = 'Enrichment vs {0} Data'.format(comp_prefix)
    else:
        title = 'Enrichment Scores at Various P-Value Cutoffs'

    scores.reset_index(inplace=True)
    sns.set_style('dark')
    sns.set_palette('deep')

    fig, ax1 = plt.subplots(figsize=figsize)

    # Plot counts
    if raw:
        scores.plot.line(y='sig_data', color='orange', ax=ax1, legend=False)
        ax1.set_ylabel(
            'Number Kept\n(Max {0:,}, Min {1:,})'.format(
                scores.sig_data.max(), scores.sig_data.min()
            ), fontsize=14
        )
    else:
        scores['perc'] = scores.sig_data/scores.sig_data.max()
        scores.plot.line(y='perc', color='orange', ax=ax1, legend=False)
        ax1.set_ylabel(
            'Percentage Kept vs P={0}'.format(scores.cutoff.max()),
            fontsize=14
        )
        ax1.yaxis.set_major_formatter(
            FuncFormatter(lambda y, _: '{0:.0%}'.format(y))
        )

    # Format and label x-axis
    ax1.set_xlabel('p-value cutoff', fontsize=14)
    ax1.set_xticks(range(0, len(scores)+1, 1))
    ax1.set_xticklabels(
        scores.cutoff.apply(
            lambda x: '{0}'.format(float('{0:1e}'.format(x)))
        ), rotation=45
    )

    # Plot enrichment score on opposite x-axis
    ax2 = ax1.twinx()
    scores.plot.line(y='enrichment_score', color='blue', ax=ax2, legend=False)
    ax2.set_ylabel(
        'Enrichment Score\n(sig:non-sig enrichment)',
        fontsize=14
    )
    ax2.set_title(title, fontsize=17)

    # Format and add legend
    fig.tight_layout()
    fig.legend(
        labels=('Percentage Kept', 'Enrichment Score'),
        loc='lower right'
    )

    if outfile:
        fig.savefig(outfile)

    #  return fig, [ax1, ax2]


def get_set(x):
    """Return frozenset from x if x is iterable or newline separated file."""
    if isinstance(x, str):
        with open_zipped(x) as fin:
            return frozenset(fin.read().strip().split('\n'))
    return frozenset(x)


###############################################################################
#                         Splitting Comparison Study                          #
###############################################################################


def split_study(study_file, prefix=None, simp_fl=False, low_mem=False,
                conf=None):
    """Split a sample into a significant file and a non-sig file.

    Params
    ------
    study_file : str
        Path the the file to parse, can be gzipped.
    prefix : str, optional
        A prefix to use for output files, default is study_file name.
    simp_fl : bool, optional
        Treat the study_file as a two column file with no header, where the
        columns are name and pvalue, separated by a tab
    low_mem : bool, optional, not implemented
        Parse file line by line, instead of using pandas. Currently only
        low-mem works
    conf : dict, optional
        A config file to use for parsing.

    Writes
    ------
    <prefix>.sig.txt.gz, <prefix>.non-sig.txt.gz
        Two newline separated of names that are significant or non-significant.
    """
    prefix = prefix if prefix else str(study_file)
    sig_fl = prefix + '.sig.txt.gz'
    nsig_fl = prefix + '.nonsig.txt.gz'

    # Cutoffs
    sig_p = float(conf['comp_sig_pval'])
    non_sig_p = float(conf['comp_nonsig_pval'])
    _sys.stderr.write(
        'Significant set is P <= {0}, non-sig is P >= {1}\n'
        .format(sig_p, non_sig_p)
    )

    if simp_fl:
        sample = simple_iterator(study_file)
    else:
        sample = study_iterator(
            study_file, conf['comp_sep'], conf['comp_has_header'],
            conf['comp_name_col'], conf['comp_pval_col']
        )

    sig_set = set()
    nsig_set = set()
    add_to_sig = sig_set.add
    add_to_nsig = nsig_set.add
    _sys.stderr.write('Splitting dataset\n')
    for name, p_val in _tqdm(sample, unit=' rows', disable=None):
        if p_val <= sig_p:
            add_to_sig(name)
        elif p_val >= non_sig_p:
            add_to_nsig(name)

    _sys.stderr.write('Sorting results and writing\n')
    with open_zipped(sig_fl, 'w') as sigf, open_zipped(nsig_fl, 'w') as nsgf:
        sigf.write('\n'.join(sorted(sig_set)))
        nsgf.write('\n'.join(sorted(nsig_set)))

    _sys.stderr.write(
        'Splitting done, written {} rows to {} and {} rows to {}\n'
        .format(len(sig_set), sig_fl, len(nsig_set), nsig_fl)
    )


###############################################################################
#                                File Iterators                               #
###############################################################################


def study_iterator(infile, sep, has_header, name_col, pval_col):
    """Iterate through infile, yield (name, p-value).

    Params
    ------
    infile : str
        Path to a file to work on.
    sep : str
        Single character to split on
    has_header : bool
        Is first line a header
    name_col : str or int
        Name of col as str if has_header is True, else 0-base column index
    pval_col : str or int
        Name of col as str if has_header is True, else 0-base column index

    Yields
    ------
    name : str
        Name of record
    p-value : float
        P-Value of record
    """
    with open_zipped(infile) as fin:
        if has_header:
            header = fin.readline().strip().split(sep)
            name_idx = header.index(name_col)
            pval_idx = header.index(pval_col)
        else:
            try:
                name_idx = int(name_col)
                pval_idx = int(pval_col)
            except ValueError:
                _sys.stderr.write(
                    'Comp column names must be numbers if no header\n'
                )
                raise
        for line in fin:
            f = line.rstrip().split(sep)
            yield f[name_idx], float(f[pval_idx])


def simple_iterator(infile):
    """Iterate through infile, yield (col 1, col 2)."""
    with open_zipped(infile) as fin:
        for line in fin:
            f = line.rstrip().split('\t')
            assert len(f) == 2
            yield f[0], float(f[1])


###############################################################################
#                                   Config                                    #
###############################################################################


# Contains defaults and help, used to generate a simple key: value dictionary
DEFAULT_CONFIG = {
    'test_sep': {
        'default': '\t', 'help': 'Separator used in the test dataset'
    },
    'comp_sep': {
        'default': '\t', 'help': 'Separator used in the comparison dataset'
    },
    'test_has_header': {
        'default': 1, 'help': '1 if has header, 0 if does not'
    },
    'test_name_col': {
        'default': 'name',
        'help': 'Column name (or number if no header) for names in test data'
    },
    'test_pval_col': {
        'default': 'p_value',
        'help': 'Column name (or number if no header) for pvals in test data'
    },
    'comp_has_header': {
        'default': 1, 'help': '1 if has header, 0 if does not'
    },
    'comp_name_col': {
        'default': 'name',
        'help': 'Column name (or number if no header) for names in comparison data'
    },
    'comp_pval_col': {
        'default': 'p_value',
        'help': 'Column name (or number if no header) for pvals in comparison data'
    },
    'max_pval': {
        'default': 0.05, 'help': 'Max pvalue to test enrichment for'
    },
    'min_pval': {
        'default': 1e-15, 'help': 'Min pvalue to test enrichment for'
    },
    'comp_sig_pval': {
        'default': 1e-4,
        'help': 'pvalue to use as significant cutoff when splitting comparison data'
    },
    'comp_nonsig_pval': {
        'default': 0.98,
        'help': 'pvalue to use as not-significant cutoff when splitting comparison data'
    }
}


def get_default_conf():
    """Return simple dict from DEAFULT_CONFIG."""
    return {k: v['default'] for k, v in DEFAULT_CONFIG.items()}


def conf_help(outstream=_sys.stdout):
    """Print config help to outstream (default STDOUT)."""
    conf = [
        [k, repr(v['default']), v['help']] for k, v in DEFAULT_CONFIG.items()
    ]
    help = _tab(conf, headers=['variable', 'default', 'help'])
    help = 'Config file can contain the following values:\n\n{}\n'.format(
        help
    )
    if outstream:
        outstream.write(help)
    return help


def parse_config_file(conf_file=None):
    """Load a dict from a json file and update it with defaults.

    Params
    ------
    conf_file : str, optional
        Path to a json file. If None, just return default conf.

    Returns
    -------
    config : dict
    """
    conf = get_default_conf()

    if conf_file:
        with open_zipped(conf_file) as fin:
            conf.update(_json.load(fin))

    return conf


###############################################################################
#                                   Helpers                                   #
###############################################################################


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


###########################################################################
#                          Command Line Parsing                           #
###########################################################################


# Mode descriptions
MAIN_DESC = """\
Run the enrichment.

Requires a split comparison dataset, similar to the one produced by the
split_list mode. Inputs are (optionally compressed) tables. The default is to
expect a tab-delimited table with a header line where the column names are
'name' and 'p-value'.  If ``--simple-file`` is passed, the expected input is a
two column file of name\\tp-value, with no header row.

To customize this, pass a json config file, the defaults can be written out by
running dump_config.
"""

SPLIT_DESC = """\
Split an existing study into a highly significant set and a non-significant set.

The expected input is an (optionally compressed) tab delimited file with a
header line and one column labelled 'name' and another labelled 'p-value'. If
``--simple-file`` is passed, the expected input is a two column file of
name\\tp-value, with no header row.

To customize this, pass a json config file, the defaults can be written out by
running dump_config.

The outputs are two newline-separated compressed files where each line is a
name to compare to (case-sensitive). The prefix can be specified with
'--prefix' (default is the same name as the input), the two suffices are
'sig.txt.gz' and 'non-sig.txt.gz', this is non-configurable.

By default, sig gets everything with a p-value smaller than 1e-4, non-sig gets
everything with a p-value greater than 0.99. These values can be configured
in the config.
"""

CONF_DESC = """\
Dump a default config to config file.

Optionally, you can pass an existing config file, and that one will be used
to update the defaults before writing the output conf file.
"""

def core_args(args):
    """Run the enrichment."""
    conf = parse_config_file(args.config_file)
    if args.max_p:
        conf['max_pval'] = args.max_p
    if args.min_p:
        conf['min_pval'] = args.min_p

    sig_comp  = args.prefix + '.sig.txt.gz'
    nsig_comp = args.prefix + '.nonsig.txt.gz'

    bad = []
    for fl in [sig_comp, nsig_comp]:
        if not _os.path.isfile(fl):
            bad.append(fl)
    if bad:
        raise OSError(
            'Need both {0} and {1} to exist, missing {2}'
            .format(sig_comp, nsig_comp, bad)
        )

    scores = enrich_study(
        args.data, sig_comp, nsig_comp, simp_fl=args.simple_file,
        low_mem=False, conf=conf
    )

    if args.output:
        _sys.stderr.write('Writing score table to {0}\n'.format(args.output))
        if args.output.endswith('xls') or args.output.endswith('xlsx'):
            scores.to_excel(
                args.output, index=False, sheet_name='Enrichment Scores'
            )
        elif args.output.endswith('pd') or args.output.endswith('pickle'):
            scores.to_pickle(args.output)
        else:
            with open_zipped(args.output, 'w') as fout:
                scores.to_csv(fout, sep=conf['test_sep'], index=False)

    if args.plot:
        _sys.stderr.write('Plotting scores to {0}\n'.format(args.plot))
        plot_scores(scores, outfile=args.plot, comp_prefix=args.prefix)


def plot_args(args):
    """Run plotting only."""
    conf = parse_config_file(args.config_file)
    _sys.stderr.write('Getting scores\n')
    if args.scores.endswith('xls') or args.scores.endswith('xlsx'):
        scores = pd.read_excel(args.scores, sheet_name='Enrichment Scores')
    elif args.scores.endswith('pd') or args.scores.endswith('pickle'):
        scores = pd.read_pickle(args.scores)
    else:
        with open_zipped(args.scores) as fin:
            scores = pd.read_csv(fin, sep=conf['test_sep'])
    scores['idx'] = scores.cutoff.astype(float)
    scores.set_index('idx', drop=True, inplace=True)
    _sys.stderr.write('Plotting scores to {0}\n'.format(args.plot))
    plot_scores(
        scores, outfile=args.plot, comp_prefix=args.prefix, raw=args.raw
    )


def split_args(args):
    """Run the comp study splitting sub-command."""
    conf = parse_config_file(args.config_file)
    if args.sig_p:
        conf['comp_sig_pval'] = float(args.sig_p)
    if args.non_sig_p:
        conf['comp_nonsig_pval'] = float(args.non_sig_p)

    # Run the splitting algorithm
    split_study(
        args.data_file, prefix=args.prefix,
        simp_fl=args.simple_file, conf=conf
    )


def conf_args(args):
    """Run the config subcommand."""
    conf = parse_config_file(args.config_file)

    _sys.stderr.write('Writing default config to {0}\n\n'.format(args.outfile))
    with open_zipped(args.outfile, 'w') as fout:
        _json.dump(conf, fout, indent=4)

    conf_help(_sys.stdout)
    _sys.stdout.write(
        '\nNone of the above are required, and any extras will be ignored\n\n'
    )


def main(argv=None):
    """Run as a script."""
    if not argv:
        argv = _sys.argv[1:]

    parser  = _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter
    )

    # Shared arguments
    conf_parser = _argparse.ArgumentParser(add_help=False)
    conf_parser.add_argument(
        '-c', dest='config_file',
        help='Optional config file, get by running `dump-config`'
    )
    file_args = _argparse.ArgumentParser(add_help=False)
    file_grp = file_args.add_argument_group('File Overrides')
    file_grp.add_argument(
        '--simple-file', action='store_true',
        help='Treat input file as tab-sep two column file with no header. ' +
        'First column is name, second is p-value.'
    )
    prefix_args = _argparse.ArgumentParser(add_help=False)
    prefix_args.add_argument(
        '-p', '--prefix', default='comp_data',
        help='Optional comparison study prefix'
    )

    # Subparsers
    subparsers = parser.add_subparsers(dest='mode')

    main_mode = subparsers.add_parser(
        'run', description=MAIN_DESC, help='Run the enrichment',
        parents=[conf_parser, file_args, prefix_args],
        formatter_class=_argparse.RawDescriptionHelpFormatter
    )
    files = main_mode.add_argument_group('Optional Output Files')
    files.add_argument(
        '-o', '--output',
        help='Write table to this file in addition to STDOUT, can be Excel fl'
    )
    files.add_argument(
        '--plot', help='Write a plot to this file'
    )
    files2 = main_mode.add_argument_group('Required Inputs')
    files2.add_argument(
        'data', help='Described by config file (dump-config)'
    )
    mconf_override = main_mode.add_argument_group(
        'Config Overrides (Optional, set in config file)'
    )
    mconf_override.add_argument(
        '--max-p', help='Max p-value to consider',
        metavar='maxP', type=float
    )
    mconf_override.add_argument(
        '--min-p', help='Min p-value to consider',
        metavar='maxP', type=float
    )
    main_mode.set_defaults(func=core_args)

    split_mode = subparsers.add_parser(
        'split', description=SPLIT_DESC,
        parents=[conf_parser, file_args, prefix_args],
        help='Split an existing dataset into two files for enrichment',
        formatter_class=_argparse.RawDescriptionHelpFormatter
    )
    split_mode.add_argument(
        'data_file', help='The dataset to parse, must have name and p-value'
    )
    conf_override = split_mode.add_argument_group(
        'Config Overrides (Optional, in config file)'
    )
    conf_override.add_argument(
        '--sig-p', help='P-Value to choose significant records',
        metavar='sigP', type=float
    )
    conf_override.add_argument(
        '--non-sig-p', help='P-Value to chose non-significant records',
        metavar='unsigP', type=float
    )
    split_mode.set_defaults(func=split_args)

    plot_mode = subparsers.add_parser(
        'plot', parents=[conf_parser, prefix_args],
        description='Plot results', help='Plot results'
    )
    plot_mode.add_argument('scores', help='Scores table from run mode')
    plot_mode.add_argument('plot', help='File to write plot to')
    plot_mode.add_argument(
        '--raw', action='store_true',
        help="Plot raw counts instead of percentages"
    )
    plot_mode.set_defaults(func=plot_args)

    conf_mode = subparsers.add_parser(
        'dump-config', description=CONF_DESC, epilog=conf_help(None),
        help='Dump a json config file', parents=[conf_parser],
        formatter_class=_argparse.RawDescriptionHelpFormatter
    )
    conf_mode.add_argument('outfile', help='File to write json config to')
    conf_mode.set_defaults(func=conf_args)


    args = parser.parse_args(argv)
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()
        return 0


if __name__ == '__main__' and '__file__' in globals():
    _sys.exit(main())
