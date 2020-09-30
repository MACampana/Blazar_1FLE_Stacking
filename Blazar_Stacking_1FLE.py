#!/usr/bin/env python

import datetime
today_date = datetime.datetime.now().strftime("%b-%d-%Y_%H-%M-%S")

#Imports
import numpy as np
import scipy.interpolate as spint
import pandas as pd

import astropy as ap
import astropy.io.fits as fits

import csky as cy
import histlite as hl
import healpy as hp
import matplotlib.pyplot as plt


plt.ioff()

#ARGS
import argparse

parser = argparse.ArgumentParser(description='Perform stacking analysis on 1FLE Blazars with the 10yrPStracks dataset.')

parser.add_argument('--dat_form', default=None, choices=[None, 'csv', 'txt'], dest='data_save')
parser.add_argument('-g','--gamma', default=2.0, type=float, dest='gamma')
parser.add_argument('-w', '--weight', default='equal', choices=['equal', 'flux'], dest='weighting_scheme')
parser.add_argument('-f', '--flux-pdf', action='store_true', dest='plot_flux_pdf')
parser.add_argument('-n', '--num-trials', default=10000, type=int, dest='num_trials')
parser.add_argument('--seed', default=1, type=int, dest='seed')
parser.add_argument('-t', '--tol', default=.01, type=float, dest='tol')
parser.add_argument('--batches', default=6, type=int, dest='n_batches')
parser.add_argument('--batch-size', default=500, type=int, dest='batch_size')
parser.add_argument('--ns-step', default=5, type=int, dest='n_sig_step')
parser.add_argument('--batch1', default=50, type=int, dest='first_batch_size')
parser.add_argument('--E0', default=100, type=float, dest='ref_E')
parser.add_argument('--funits', default=1e3, type=float, dest='flux_unit')
parser.add_argument('-c', '--do-trials', action='store_true', dest='do_bg_trials')
parser.add_argument('--trfile', default=None, dest='bg_trials_fname')
parser.add_argument('-s', '--sens', action='store_true', dest='do_find_sens')
parser.add_argument('-d', '--disc', action='store_true', dest='do_find_disc')
parser.add_argument('-b', '--bias', action='store_true', dest='do_bias_test')

#Args -> variables
args = parser.parse_args(args=['--trfile', 'flux'])

data_save = args.data_save
gamma = args.gamma
weighting_scheme = args.weighting_scheme
plot_flux_pdf = args.plot_flux_pdf
num_trials = args.num_trials
seed = args.seed
tol = args.tol
n_batches = args.n_batches
batch_size = args.batch_size
n_sig_step = args.n_sig_step
first_batch_size = args.first_batch_size
ref_R = args.ref_E
flux_unit = args.flux_unit
do_bg_trials = args.do_bg_trials
bg_trials_fname = args.bg_trials_fname
do_find_sens = args.do_find_sens
do_find_disc = args.do_find_disc
do_bias_test = args.do_bias_test

print(vars(args))
