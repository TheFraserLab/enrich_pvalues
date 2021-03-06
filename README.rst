###############
Enrich p-values
###############

A simple script to compare p-values between a test and comparison dataset at a
variety of p-value cutoffs. By plotting the enrichment score at a variety of
cutoffs, it is possible to pick the optimal cutoff for your data.

Version: 1.0-beta2

.. image:: https://github.com/TheFraserLab/enrich_pvalues/raw/master/enrich_score.gif

.. contents:: **Contents**

Algorithm
=========

For each p-value in the interval between ``max_pval`` (default: 0.05) and
``min_pval`` (default: 1e-15), we test at intervals of 1 and 5 for each order of
magnitute, e.g. 0.05, 0.01, 0.005, 0.001, 5e-4, 1e-4, 5e-5, 1e-6, ... 1e-15.

To test, we simply take all identities with a p-value less than the cutoff and
compare them to all identities in the comparison set with p-values below the
``comp_set_pvalue``. We simply ask what percentage or the test set are in the
comparison set. We then do exactly the same with the entire set of identities in
the comparison set that have a p-value greater than 0.98.

The identities are generally going to be gene or SNP names, but they can be
anything (e.g. coordinates) as long as they overlap in the test and comparison
data.

Installation
============

Install via PyPI:

.. code:: shell

    pip install enrich_pvalues

Or install from github:

.. code:: shell

    pip install https://github.com/TheFraserLab/enrich_pvalues/tarball/master

It should work with python 2 or 3, but python 3 is recommended.

Requirements
------------

In ``requirements.txt``, we use numpy, pandas, matplotlib, seaborn, tabulate,
and tqdm.

Usage
=====

After install, run ``enrich_pvalues --help`` to get a full description of all
options. There are four main modes:

- ``dump-config``
- ``split``
- ``run``
- ``plot``

Each has it's own help, so run e.g. ``enrich_pvalues split -h`` to learn how to
split a comparison dataset.

Example Steps
-------------

First, dump a configuration file to describe your data:

.. code:: shell

   enrich_pvalues dump-config enrich_atac.json

This will also print a help table describing each option. You need to describe
your comparison data and your test data, and pick your p-value thresholds.

Next, split your comparison dataset into two tables: significant, and
not-significant:

.. code:: shell

   enrich_pvalues split -c enrich_atac.json --prefix atac /path/to/comp_data.txt.gz

Now, run the enrichment using those two tables and your test data:

.. code:: shell

   enrich_pvalues run -c enrich_atac.json -o atac_scores.xlsx -p atac /path/to/test_data.txt

Note, the second to last argument is the prefix from the second step.

Note: the scores can be excel format, pickled format, or text format, depending
on the suffix. Also, the prefix in this plot step is different, it is used to
title the plot only, and so can be whatever you want.

Finally, plot the data. This can also be done by passing e.g. ``--plot myplot.png``
to the run step, although that has fewer options.

.. code:: shell

   enrich_pvalues plot --prefix caQTL atac_scores.xlsx atac_plot.pdf

The resulting plot will look something like this:

.. image::
    https://github.com/TheFraserLab/enrich_pvalues/raw/master/plot_example.png

To control the name of the comparison dataset, pass ``-p <name>``, this is only
used for title formatting and so does not need to be the same as the prefix used
in earlier steps.

To format the counts as raw numbers instead of a percentage, pass ``--raw``.

Finally, it can be useful to limit the range of cutoffs to zoom the plot into a
region of interest. To do that, pass ``--min-p`` and ``--max-p``. e.g.:

.. code:: shell

    enrich_pvalues plot --min-p 5e-3 --max-p 1e-7 --raw --prefix caQTL atac_scores.xlsx plot_example.png

That command is the one used to create the above example plot.
