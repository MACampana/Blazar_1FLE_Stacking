#!/usr/bin/env python

import datetime
dt = datetime.datetime.now().strftime("%b-%d-%Y_%H-%M-%S")

import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import pickle
import argparse

#parser = argparse.ArgumentParser(description='Plot Sensitivity and Discovery Potential vs Gamma')
#parser.add_argument('-w', '--weight', default='equal', choices=['equal', 'flux'], dest='W')

#args = parser.parse_args()
#W = args.W

eq_snamelist = glob('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/Sens*equal*.pkl')
eq_dnamelist = glob('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/Disc*equal*.pkl')

fx_snamelist = glob('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/Sens*flux*.pkl')
fx_dnamelist = glob('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/Disc*flux*.pkl')

#print(snamelist)

def get_arrays(namelist):
    fxs = []
    gs = []
    for name in namelist:
        f = open(name, "rb")
        dict = pickle.load(f)
        f.close()

        fx = dict['flux_nsig_at1']
        fxs.append(fx)
        g = dict['inj_gamma']
        gs.append(g)
    
    order = np.argsort(gs)
    gs = np.array(gs)[order]
    fxs = np.array(fxs)[order]
    return fxs, gs

eq_sens, eq_s_gammas = get_arrays(eq_snamelist)
eq_disc, eq_d_gammas = get_arrays(eq_dnamelist)
fx_sens, fx_s_gammas = get_arrays(fx_snamelist)
fx_disc, fx_d_gammas = get_arrays(fx_dnamelist)

def plot_sens_v_gamma():
    plt.figure(figsize=(8,6))
    plt.plot(eq_s_gammas, eq_sens, color='b', marker='+', label='Sensitivity, Equal')
    plt.plot(fx_s_gammas, fx_sens, color='g', marker='+', label='Sensitivity, Flux')
    plt.plot(eq_d_gammas, eq_disc, color='b', marker='x', label='DP, Equal', linestyle='None')
    plt.plot(fx_d_gammas, fx_disc, color='g', marker='x', label='DP, Flux', linestyle='None')

    plt.yscale('log')

    plt.ylabel('E^2 dN/dE [TeV /cm^2 /s   @ 1TeV]')
    plt.xlabel('gamma')
    
    #plt.ylim(1e-13, 1e-11)
    
    plt.title('Sensitivities')
    plt.legend()
    plt.grid(b=True, which='both', alpha=.7, linestyle=':')
    plt.savefig('/data/user/mcampana/analysis/Blazar_1FLE/plots/SensDisc_v_gamma_{}'.format(dt))
    plt.close()
    return

#depricated as of 10/23/2020
def plot_sens_v_E():
    xs = np.logspace(-2, 2, 100)
    
    for i in range(len(gammas)):       
        ys = senss[i] * xs**(-gammas[i]) * xs**2
        plt.plot(xs, ys, label='gam={}'.format(gammas[i]), linestyle=':', alpha=1)
    
    plt.vlines(1,1e-14,1e-10,colors=['k'], linestyles=['dashed'], alpha=.5)
    
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('E / 100TeV')
    plt.ylabel('E^2 dN/dE [Tev /cm^2 /s @ 100 TeV]')
    
    plt.ylim(bottom=1e-14)
    plt.xlim(right=2e1)
    
    plt.title('{} Weighting'.format(W.capitalize()))
    plt.savefig('/data/user/mcampana/analysis/Blazar_1FLE/plots/Flux_v_E_{}weighting_{}'.format(W, dt))
    plt.close()
    return

plot_sens_v_gamma()