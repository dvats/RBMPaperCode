#!/bin/bash

export PKGNAME='bnn_mcmc_examples'
export CONDADIR="$HOME/opt/continuum/miniconda/miniconda3"
export CONDAENV="$CONDADIR/envs/$PKGNAME"
export CONDABIN="$CONDADIR/bin/conda"
export OUTPUTPATH="$HOME/output/bnn_mcmc_examples/mlp/pima/setting2/metropolis_hastings"

qsub \
  -cwd \
  -V \
  -l short \
  -N mlp_pima_mh \
  -o $OUTPUTPATH \
  -e $OUTPUTPATH \
  $CONDABIN run -p $CONDAENV python benchmark_run.py
