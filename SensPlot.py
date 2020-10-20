#!/usr/bin/env python

import datetime
dt = datetime.datetime.now().strftime("%b-%d-%Y_%H-%M-%S")

import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import pickle
import argparse

parser = argparse.ArgumentParser(description='Plot Sensitivity and Discovery Potential vs Gamma')
parser.add_argument('-w', '--weight', default='equal', choices=['equal', 'flux'], dest='W')

args = parser.parse_args()
W = args.W

snamelist = glob('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/Sens*{}*.pkl'.format(W))
dnamelist = glob('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/Disc*{}*.pkl'.format(W))

#print(snamelist)

gammas = []
senss = []
gammad = []
discs = []
for sname in snamelist:
    f = open(sname, "rb")
    s_dict = pickle.load(f)
    f.close()

    sens = s_dict['flux_nsig']
    senss.append(sens)
    gamma = s_dict['inj_gamma']
    gammas.append(gamma)

for dname in dnamelist:
    f = open(dname, "rb")
    d_dict = pickle.load(f)
    f.close()

    disc = d_dict['flux_nsig']
    discs.append(disc)
    gamma = d_dict['inj_gamma']
    gammad.append(gamma)

orders = np.argsort(gammas)
gammas = np.array(gammas)[orders]
senss = np.array(senss)[orders]

orderd = np.argsort(gammad)
gammad = np.array(gammad)[orderd]
discs = np.array(discs)[orderd]

print('gamma (s): ', gammas)
print('sens: ', senss)
print('----------------')
print('gamma (d): ', gammad)
print('disc: ', discs)

def plot_sens_v_gamma():
    plt.figure(figsize=(8,6))
    plt.plot(gammas, senss, color='k', marker='+', label='Sensitivity')
    plt.plot(gammad, discs, color='r', marker='+', label='DP @ gamma=2')

    plt.yscale('log')

    plt.ylabel('E^2 dN/dE [TeV /cm^2 /s   @ 100TeV]')
    plt.xlabel('gamma')
    
    plt.ylim(1e-13, 1e-11)
    
    plt.title('{} Weighting'.format(W.capitalize()))
    plt.legend()
    plt.grid(b=True, which='both', alpha=.7, linestyle=':')
    plt.savefig('/data/user/mcampana/analysis/Blazar_1FLE/plots/SensDisc_v_gamma_{}weighting_{}'.format(W, dt))
    plt.close()
    return

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

plot_sens_v_E()