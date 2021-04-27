#!/usr/bin/env python

import datetime
dt = datetime.datetime.now().strftime("%b-%d-%Y_%H-%M-%S")

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
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
        if 'info' in dict1.keys():
            fx = dict1['info']['flux_nsig']
            fxs.append(fx)
            g = dict1['info']['inj_gamma']
            gs.append(g)
        else:
            fx = dict1['flux_nsig']
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
        if 'info' in dict1.keys():
            fx = dict1['info']['flux_nsig']
            fxs.append(fx)
            E = dict1['info']['flux_E0']
            Es.append(E)
        else:
            fx = dict1['flux_nsig']
            fxs.append(fx)
            E = dict1['flux_E0']
            Es.append(E)
    
    order = np.argsort(Es)
    Es = np.array(Es)[order]
    fxs = np.array(fxs)[order]
    return fxs, Es

#Get file name lists (both weightings, and sens or DP)

eq_snamelist = glob('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/poisson/hemisphere/both/weight/equal/Sens*u1FLEblazars_F*.pkl')
fx_snamelist = glob('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/poisson/hemisphere/both/weight/flux/Sens*u1FLEblazars_F*.pkl')

eq_sens, eq_s_gammas = get_arrays_SensGamma(eq_snamelist)
fx_sens, fx_s_gammas = get_arrays_SensGamma(fx_snamelist)
    
eq_dnamelist = glob('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/poisson/hemisphere/both/weight/equal/Disc*u1FLEblazars_F*.pkl')
fx_dnamelist = glob('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/poisson/hemisphere/both/weight/flux/Disc*u1FLEblazars_F*.pkl')

eq_disc, eq_d_gammas = get_arrays_SensGamma(eq_dnamelist)
fx_disc, fx_d_gammas = get_arrays_SensGamma(fx_dnamelist)


#Get file name lists (both weightings, and sens or DP)
N_eq_dsnamelist = glob('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/poisson/hemisphere/north/weight/equal/DiffSens*u1FLE*Apr-22*.pkl')
N_fx_dsnamelist = glob('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/poisson/hemisphere/north/weight/flux/DiffSens*u1FLE*Apr-22*.pkl')
S_eq_dsnamelist = glob('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/poisson/hemisphere/south/weight/equal/DiffSens*u1FLE*Apr-22*.pkl')
S_fx_dsnamelist = glob('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/poisson/hemisphere/south/weight/flux/DiffSens*u1FLE*Apr-22*.pkl')
B_eq_dsnamelist = glob('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/poisson/hemisphere/both/weight/equal/DiffSens*u1FLE*Apr-22*.pkl')
B_fx_dsnamelist = glob('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/poisson/hemisphere/both/weight/flux/DiffSens*u1FLE*Apr-22*.pkl')

N_eq_diffsens, N_eq_s_Es = get_arrays_Diffsens(N_eq_dsnamelist)
N_fx_diffsens, N_fx_s_Es = get_arrays_Diffsens(N_fx_dsnamelist)   
S_eq_diffsens, S_eq_s_Es = get_arrays_Diffsens(S_eq_dsnamelist)
S_fx_diffsens, S_fx_s_Es = get_arrays_Diffsens(S_fx_dsnamelist)   
B_eq_diffsens, B_eq_s_Es = get_arrays_Diffsens(B_eq_dsnamelist)
B_fx_diffsens, B_fx_s_Es = get_arrays_Diffsens(B_fx_dsnamelist)  


def plot_sens_v_gamma():
    
    figax1 = plt.figure(figsize=(8,6))

    figax1.plot(eq_s_gammas, eq_sens, color='orangered', marker='x', label='Sens, Equal', linestyle='--')
    figax1.plot(fx_s_gammas, fx_sens, color='orangered', marker='x', label='Sens, Flux')

    figax1.plot(eq_d_gammas, eq_disc, color='lightsalmon', marker='x', label='DP, Equal', linestyle='--')
    figax1.plot(fx_d_gammas, fx_disc, color='lightsalmon', marker='x', label='DP, Flux')
    
    figax1.set_yscale('log')
    figax1.set_ylabel('E^2 dN/dE [TeV /cm^2 /s @ 100TeV]')
    figax1.set_ylim(1e-15, 1e-13)
    figax1.legend()
    figax1.set_title('Sensitivities: 1FLE Blazar Stacking')
    figax1.grid(b=True, which='major', alpha=.7, linestyle=':')
    figax1.tick_params(axis='x', which='both', labelbottom=False)
    figax1.set_xlabel('gamma')    
    
    plt.savefig('/data/user/mcampana/analysis/Blazar_1FLE/plots/u1FLE_SensDisc_v_gamma_{}'.format(dt))
    plt.close()
    return

def plot_diff_sens(): #

    plt.figure(figsize=(10,6))
    plt.step(np.clip(np.append(B_eq_s_Es[:-1], 10**4.5), 14.8, 14797.3), np.append(B_eq_diffsens[:-1], B_eq_diffsens[-2]), color='k', label='All-sky, Equal', where='post', linestyle=':')
    plt.step(np.clip(np.append(B_fx_s_Es[:-1], 10**4.5), 10.5, 12188.9), np.append(B_fx_diffsens[:-1], B_fx_diffsens[-2]), color='k', label='All-sky, Flux', where='post')
    
    plt.step(np.clip(np.append(S_eq_s_Es[:-1], 10**4.5), 14.8, 14797.3), np.append(S_eq_diffsens[:-1], S_eq_diffsens[-2]), color='C1', label='South, Equal', where='post', linestyle=':')
    plt.step(np.clip(np.append(S_fx_s_Es[:-1], 10**4.5), 10.5, 12188.9), np.append(S_fx_diffsens[:-1], S_fx_diffsens[-2]), color='C1', label='South, Flux', where='post')
    
    plt.step(np.clip(np.append(N_eq_s_Es[:-1], 10**4.5), 14.8, 14797.3), np.append(N_eq_diffsens[:-1], N_eq_diffsens[-2]), color='C2', label='North, Equal', where='post', linestyle=':')
    plt.step(np.clip(np.append(N_fx_s_Es[:-1], 10**4.5), 10.5, 12188.9), np.append(N_fx_diffsens[:-1], N_fx_diffsens[-2]), color='C2', label='North, Flux', where='post')
    
    plt.hlines(eq_sens[1], 14.8, 14797.3, linestyles='dashed', alpha=0.7, colors='k', label='All-sky, All-E, Equal')
    plt.hlines(fx_sens[1], 10.5, 12188.9, linestyles='dashdot', alpha=0.7, colors='k', label='All-sky, All-E, Flux')

    plt.yscale('log')
    plt.xscale('log')
    
    plt.ylabel('Per-Bin E^2 dN/dE [TeV /cm^2 /s]')
    plt.xlabel('Energy [TeV]')
    
    plt.ylim(1e-12, 1e-9)
    plt.xlim(9e0,2e4)
    
    plt.title('Differential Sensitivities: 1FLE Stacking (gamma=2.0)')
    plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
    plt.grid(b=True, which='major', alpha=.7, linestyle=':')
    plt.tight_layout()
    plt.savefig('/data/user/mcampana/analysis/Blazar_1FLE/plots/u1FLE_DiffSens_{}'.format(dt))
    plt.close()
    
    return

#plot_sens_v_gamma()
plot_diff_sens()