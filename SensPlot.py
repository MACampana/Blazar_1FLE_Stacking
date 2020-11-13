#!/usr/bin/env python

import datetime
dt = datetime.datetime.now().strftime("%b-%d-%Y_%H-%M-%S")

import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import pickle
#import argparse

def get_arrays_SensGamma(namelist):
    fxs = []
    gs = []
    for name in namelist:
        f = open(name, "rb")
        dict1 = pickle.load(f)
        f.close()

        fx = dict1['flux_nsig_at1']
        fxs.append(fx)
        g = dict1['inj_gamma']
        gs.append(g)

    order = np.argsort(gs)
    gs = np.array(gs)[order]
    fxs = np.array(fxs)[order]
    return fxs, gs

def get_arrays_Diffsens(namelist):
    fxs = []
    Es = []
    for name in namelist:
        f = open(name, "rb")
        dict1 = pickle.load(f)
        f.close()

        fx = dict1['flux_nsig']
        fxs.append(fx)
        E = dict1['flux_E0']
        Es.append(E)
    
    order = np.argsort(Es)
    Es = np.array(Es)[order]
    fxs = np.array(fxs)[order]
    return fxs, Es

#Get file name lists (both weightings, and sens or DP)
eq_snamelist = glob('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/poisson/hemisphere/both/weight/equal/Sens*.pkl')
eq_dnamelist = glob('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/poisson/hemisphere/both/weight/equal/Disc*.pkl')
fx_snamelist = glob('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/poisson/hemisphere/both/weight/flux/Sens*.pkl')
fx_dnamelist = glob('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/poisson/hemisphere/both/weight/flux/Disc*.pkl')
    
eq_sens, eq_s_gammas = get_arrays_SensGamma(eq_snamelist)
eq_disc, eq_d_gammas = get_arrays_SensGamma(eq_dnamelist)
fx_sens, fx_s_gammas = get_arrays_SensGamma(fx_snamelist)
fx_disc, fx_d_gammas = get_arrays_SensGamma(fx_dnamelist)

def plot_sens_v_gamma():
    
    plt.figure(figsize=(8,6))
    plt.plot(eq_s_gammas, eq_sens, color='b', marker='+', label='Sensitivity, Equal')
    plt.plot(fx_s_gammas, fx_sens, color='g', marker='+', label='Sensitivity, Flux')
    plt.plot(eq_d_gammas, eq_disc, color='b', marker='x', label='DP, Equal', linestyle='None')
    plt.plot(fx_d_gammas, fx_disc, color='g', marker='x', label='DP, Flux', linestyle='None')

    plt.yscale('log')

    plt.ylabel('E^2 dN/dE [TeV /cm^2 /s @ 1TeV]')
    plt.xlabel('gamma')
    
    plt.ylim(1e-13, 2e-10)
    
    plt.title('Sensitivities')
    plt.legend()
    plt.grid(b=True, which='major', alpha=.7, linestyle=':')
    plt.savefig('/data/user/mcampana/analysis/Blazar_1FLE/plots/SensDisc_v_gamma_possion_{}'.format(dt))
    plt.close()
    return

#depricated as of 10/23/2020
#def plot_sens_v_E():
#    xs = np.logspace(-2, 2, 100)
#    
#    for i in range(len(gammas)):       
#        ys = senss[i] * xs**(-gammas[i]) * xs**2
#        plt.plot(xs, ys, label='gam={}'.format(gammas[i]), linestyle=':', alpha=1)
#    
#    plt.vlines(1,1e-14,1e-10,colors=['k'], linestyles=['dashed'], alpha=.5)
#    
#    plt.legend()
#    plt.xscale('log')
#    plt.yscale('log')
#    plt.xlabel('E / 100TeV')
#    plt.ylabel('E^2 dN/dE [Tev /cm^2 /s @ 100 TeV]')
#    
#    plt.ylim(bottom=1e-14)
#    plt.xlim(right=2e1)
#    
#    plt.title('{} Weighting'.format(W.capitalize()))
#    plt.savefig('/data/user/mcampana/analysis/Blazar_1FLE/plots/Flux_v_E_{}weighting_{}'.format(W, dt))
#    plt.close()
#    return

#Get file name lists (both weightings, and sens or DP)
N_eq_dsnamelist = glob('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/hemisphere/north/DiffSens*equal*.pkl')
N_fx_dsnamelist = glob('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/hemisphere/north/DiffSens*flux*.pkl')

S_eq_dsnamelist = glob('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/hemisphere/south/DiffSens*equal*.pkl')
S_fx_dsnamelist = glob('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/hemisphere/south/DiffSens*flux*.pkl')

B_eq_dsnamelist = glob('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/hemisphere/both/DiffSens*equal*.pkl')
B_fx_dsnamelist = glob('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/hemisphere/both/DiffSens*flux*.pkl')

N_eq_diffsens, N_eq_s_Es = get_arrays_Diffsens(N_eq_dsnamelist)
N_fx_diffsens, N_fx_s_Es = get_arrays_Diffsens(N_fx_dsnamelist)   

S_eq_diffsens, S_eq_s_Es = get_arrays_Diffsens(S_eq_dsnamelist)
S_fx_diffsens, S_fx_s_Es = get_arrays_Diffsens(S_fx_dsnamelist)   

B_eq_diffsens, B_eq_s_Es = get_arrays_Diffsens(B_eq_dsnamelist)
B_fx_diffsens, B_fx_s_Es = get_arrays_Diffsens(B_fx_dsnamelist)   

def plot_diff_sens():

    plt.figure(figsize=(8,6))
    plt.step(np.append(B_eq_s_Es,10**5), np.append(B_eq_diffsens,B_eq_diffsens[-1]), color='k', label='All-sky, Equal', where='post', linestyle=':')
    plt.step(np.append(B_fx_s_Es,10**5), np.append(B_fx_diffsens,B_fx_diffsens[-1]), color='k', label='All-sky, Flux', where='post')
    
    plt.step(np.append(S_eq_s_Es,10**5), np.append(S_eq_diffsens,S_eq_diffsens[-1]), color='C1', label='South, Equal', where='post', linestyle=':')
    plt.step(np.append(S_fx_s_Es,10**5), np.append(S_fx_diffsens,S_fx_diffsens[-1]), color='C1', label='South, Flux', where='post')
    
    plt.step(np.append(N_eq_s_Es,10**5), np.append(N_eq_diffsens,N_eq_diffsens[-1]), color='C2', label='North, Equal', where='post', linestyle=':')
    plt.step(np.append(N_fx_s_Es,10**5), np.append(N_fx_diffsens,N_fx_diffsens[-1]), color='C2', label='North, Flux', where='post')
    
    plt.hlines(eq_sens[1], 1,10000, linestyles='dashed', alpha=0.7, colors='k', label='All-sky, All-E, Equal')
    plt.hlines(fx_sens[1], 1,10000, linestyles='dashdot', alpha=0.7, colors='k', label='All-sky, All-E, Flux')

    plt.yscale('log')
    plt.xscale('log')
    
    plt.ylabel('Per-Bin E^2 dN/dE [TeV /cm^2 /s]')
    plt.xlabel('Energy [TeV]')
    
    plt.ylim(1e-12, 1e-7)
    plt.xlim(1e-1,1e5)
    
    plt.title('Differential Sensitivities')
    plt.legend(loc='upper right')
    plt.grid(b=True, which='major', alpha=.7, linestyle=':')
    plt.savefig('/data/user/mcampana/analysis/Blazar_1FLE/plots/DiffSens_{}'.format(dt))
    plt.close()
    
    return

plot_sens_v_gamma()
