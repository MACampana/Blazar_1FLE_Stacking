#!/usr/bin/env python

# See https://icecube.wisc.edu/~mrichman/docs/csky/about.html
import sys
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

import pickle

plt.ioff()

#ARGS
import argparse

parser = argparse.ArgumentParser(description='Perform stacking analysis on 1FLE Blazars with the 10yrPStracks dataset.')

parser.add_argument('--dat-form', default=None, choices=[None, 'csv', 'txt'], dest='data_save', help='csv or txt for saving sources. None for no save.')
parser.add_argument('-g','--gamma', default=[2.0], type=float, dest='gamma', nargs='+', help='Spectral index for injections.')
parser.add_argument('-w', '--weight', default=['equal'], choices=['equal', 'flux'], dest='weighting_scheme', nargs='+', help='Weighting scheme, either flux or equal.')
parser.add_argument('-f', '--flux-pdf', action='store_true', dest='plot_flux_pdf', help='Use this option to plot the flux pdf.')
parser.add_argument('-n', '--num-trials', default=1000, type=int, dest='num_trials', help='Number of background trials to compute.')
parser.add_argument('--cpus', default=1, type=int, dest='cpus', help='Number of CPUS for parallel trial running.')
parser.add_argument('--seed', default=[1], type=int, dest='seed', nargs='+', help='Seed.')
parser.add_argument('-t', '--tol', default=.01, type=float, dest='tol', help='Tolerance for sensitivity and discovery potential.')
parser.add_argument('--batches', default=6, type=int, dest='n_batches', help='Number of batches to compute for sens and disc.')
parser.add_argument('--batch-size', default=500, type=int, dest='batch_size', help='Batch size for computing sens and disc.')
parser.add_argument('--ns-step', default=50, type=int, dest='n_sig_step', help='Step for initial injection check for sens and disc.')
parser.add_argument('--batch1', default=200, type=int, dest='first_batch_size', help='First batch size.')
parser.add_argument('--E0', default=100, type=float, dest='ref_E', help='Reference energy for converting ns to flux, in units of flux_unit.')
parser.add_argument('--funits', default=1e3, type=float, dest='flux_unit', help='Energy units for flux, wrt GeV.')
parser.add_argument('-c', '--dotr', action='store_true', dest='do_trials', help='Use this option to perform trials.')
parser.add_argument('-l', '--lotr', action='store_true', dest='load_trials', help='Use this option to load trials for calculations.')
parser.add_argument('-s', '--sensdisc', action='store_true', dest='do_find_sensdisc', help='Use this option to do sensitivity and discovery calculation.')
parser.add_argument('-b', '--bias', action='store_true', dest='do_bias_test', help='Use this option to perform and plot bias tests, for given weights and gammas.')
parser.add_argument('-i', '--chi2', action='store_true', dest='plot_TSchi2', help='Use this option to plot the TS dist and chi2 fit of loaded trials.')

#Args -> variables
args = parser.parse_args()

data_save = args.data_save
gamma = args.gamma
if type(gamma)==int: gamma=[gamma]
weighting_scheme = args.weighting_scheme
plot_flux_pdf = args.plot_flux_pdf
num_trials = args.num_trials
seed = args.seed
tol = args.tol
n_batches = args.n_batches
batch_size = args.batch_size
n_sig_step = args.n_sig_step
first_batch_size = args.first_batch_size
ref_E = args.ref_E
flux_unit = args.flux_unit
do_trials = args.do_trials
load_trials = args.load_trials
do_find_sensdisc = args.do_find_sensdisc
do_bias_test = args.do_bias_test
plot_TSchi2 = args.plot_TSchi2
cpus = args.cpus

print('Loading Fits...')
#Open Catalog fits file
fits_1fle = fits.open('/data/user/mcampana/analysis/Blazar_1FLE/1fle.fits')
#Get Column Names
cols = fits_1fle[1].data.columns.names
#Extract Blazars
b_mask = (fits_1fle[1].data['CLASS1'] == 'bll') | (fits_1fle[1].data['CLASS1'] == 'fsrq')
blaz_1fle = fits_1fle[1].data[b_mask]

#29 BL LACs -> seems like 2 were changed to FSRQ from paper
#102 FSRQs -> 2 from BL LAC and 2 from Unclassified from paper
#131 sources total

#Make into dictionary
data = {}
for n in cols:
    data[n] = blaz_1fle.field(n)
#Close Fits file
fits_1fle.close()
print('Fits closed...')

print('Considering saving data...')
#Save source list with all data as CSV and Text Files
if data_save == 'csv':
    pd.DataFrame.from_dict(data).to_csv('/data/user/mcampana/analysis/Blazar_1FLE/1FLE_Blazars_Data.csv')
elif data_save == 'txt':
    f = open("/data/user/mcampana/analysis/Blazar_1FLE/1FLE_Blazars_Data.txt","w")
    f.write(str(data))
    f.close()
elif data_save == None:
    pass

#Weighting schemes
num_blaz = len(data['Name'])
#print(num_blaz)

print('Considering plotting flux pdf...')
#Plot the Flux weight PDF
if plot_flux_pdf:
    plt.figure(figsize=(6,4))
    plt.hist(data['EF30-100'], bins='auto', histtype='step', linewidth=2)
    plt.xlabel("Energy Flux from 30-100 MeV [erg /cm^2 /s]")
    plt.ylabel("Counts")
    plt.title("Energy Flux PDF for 1FLE Blazars")
    #plt.semilogy()
    plt.savefig('/data/user/mcampana/analysis/Blazar_1FLE/plots/EnergyFlux_30-100MeV_PDF_{}.png'.format(today_date))
    plt.close()
    sys.exit("Plotted Flux PDF, exiting program...")
    

print('Getting data selection...')
#Data selection
selection = cy.selections.PSDataSpecs.ps_10yr

print('Getting analysis arrays...')
#Create analysis with data selection (saved in working dir)
#ana = cy.get_analysis(cy.selections.repo, selection)
#ana.save('./')
ana_dir = cy.utils.ensure_dir('/data/user/mcampana/analysis/Blazar_1FLE/subanalyses')
#ana_repo = cy.selections.Repository(ana_dir)
ana = cy.get_analysis(cy.selections.repo, selection, dir=ana_dir)
#ana = cy.get_analysis(cy.selections.repo, cy.selections.PSDataSpecs.ps_2011, dir='/data/user/mcampana/analysis/Blazar_1FLE')

#Trial save and load directories
trials_dir = cy.utils.ensure_dir('/data/user/mcampana/analysis/Blazar_1FLE/trials')
sig_dir = cy.utils.ensure_dir('{}/sig'.format(trials_dir))
bg_dir = cy.utils.ensure_dir('{}/bg'.format(trials_dir))

#=========================================================
#       FUNCTION DEFINITIONS 
#=========================================================
#Do BG Trials
def bg_trials():
    
    trials = tr.get_many_fits(num_trials, seed=s)
    #Save Trials to numpy array file
    dir_ = cy.utils.ensure_dir('{}/weight/{}/gamma/{}'.format(bg_dir,w,g))
    save_trials_fname = '{}/BG_{}trials_{}weight_gamma{}_10yrPStracks_1FLEblazars_seed{}_{}.npy'.format(dir_, num_trials, w, g, s, today_date)
    np.save(save_trials_fname, trials.as_array)
    print("BG Trials saved ->", save_trials_fname)
    
    return

#Do Signal Trials
def sig_trials():
    
    trials = tr.get_many_fits(num_trials/2, n_sig, seed=s)
    #Save Trials to numpy array file    
    dir_ = cy.utils.ensure_dir('{}/weight/{}/gamma/{}/nsig/{}'.format(sig_dir,w,g,n_sig))
    save_trials_fname = '{}/SIG_{}trials_{}weight_gamma{}_nsig{}_10yrPStracks_1FLEblazars_seed{}_{}.npy'.format(dir_, int(num_trials/2), w, g, n_sig, s, today_date)
    np.save(save_trials_fname, trials.as_array)
    print("Signal Trials saved ->", save_trials_fname)
      
    return

#Plot TS and chi2 distribution

def TSchi2():
    fig, ax = plt.subplots()

    b = cy.bk.get_best(bg_chi2, 'weight', w, 'gamma', g)
    h = b.get_hist(bins=30)
    hl.plot1d(ax, h, crosses=True, label='{} bg trials'.format(b.n_total))

    x = h.centers[0]
    norm = h.integrate().values
    ax.semilogy(x, norm * b.pdf(x), lw=1, ls='--',
                label='Chi2 fit: {} dof, eta={}'.format(np.round(b.ndof,3), np.round(b.eta,3)))

    ax.set_xlabel('TS')
    ax.set_ylabel('number of trials')
    ax.legend()
    ax.set_title('{} Weighted, Gamma={}'.format(w.capitalize(),g))
    ax.text(10, 5e2, 'gamma={}'.format(g), ha='right', va='center')
    plt.tight_layout()
    plt.savefig('/data/user/mcampana/analysis/Blazar_1FLE/plots/Chi2TS_{}trials_{}weighting_gamma{}_10yrPStracks_1FLEblazars_{}.png'.format(b.n_total, w, g, today_date))
    plt.close()
    
    return


def multi_TSchi2():
    nrow, ncol = 3, 3
    fig, aaxs = plt.subplots(nrow, ncol, figsize=(12,12), sharex=True, sharey=True, dpi=200)
    axs = np.ravel(aaxs)
    # keep track of which ax's we already used
    used_axs = []
    for (i, gam) in enumerate(gamma):
        ax = axs[i]
        # plot histogram
        b = cy.bk.get_best(bg_chi2, 'weight', w, 'gamma', gam)
        h = b.get_hist(bins=30, range=(0, 10))
        hl.plot1d(ax, h, crosses=True, label='{} bg trials'.format(b.n_total))
        # plot chi2 fit to nonzero values
        norm = h.integrate().values
        ts = np.linspace(.1, h.range[0][1], 100)
        ax.plot(ts, norm * b.pdf(ts), label='Chi2 fit: {} dof, eta={}'.format(np.round(b.ndof,2), np.round(b.eta,3)))
        # set limits and label dec
        ax.semilogy(nonposy='clip')
        ax.set_ylim(.3, 3e3)
        ax.text(10, 5e2, 'gamma={}'.format(gam), ha='right', va='center')
        ax.legend(loc='best')
        used_axs.append(ax)
    # hide unused ax's
    for ax in axs:
        if ax not in used_axs:
            ax.set_visible(False)
    # add x and y labels
    for ax in aaxs[-1]:
        if ax in used_axs:
            ax.set_xlabel('TS')
    for ax in aaxs[:,0]:
        ax.set_ylabel('trials per bin')
    plt.tight_layout()
    plt.savefig('/data/user/mcampana/analysis/Blazar_1FLE/plots/Chi2TS_{}weighting_10yrPStracks_1FLEblazars_{}.png'.format(w, today_date))
    plt.close()
    
    return


def ndarray_to_Chi2TSD(tri):
    return cy.dists.Chi2TSD(cy.utils.Arrays(tri))

def get_n_sig(beta=0.9, nsigma=None):
    # get signal trials, background distribution, and trial runner
    sig_trials = cy.bk.get_best(sig, 'weight', w, 'gamma', g, 'nsig')
    b = cy.bk.get_best(bg_chi2, 'weight', w, 'gamma', g)
    # determine ts threshold
    if nsigma is not None:
        ts = b.isf_nsigma(nsigma)
        kind = 'Disc'
    else:
        ts = b.median()
        kind = 'Sens'
    # include background trials in calculation
    trials = {0: b.trials}
    trials.update(sig_trials)
    # get number of signal events
    # (arguments prevent additional trials from being run)
    result = tr.find_n_sig(ts, beta, trials=trials, max_batch_size=0, logging=True, n_bootstrap=1) #,tol=tol)
#    result = tr.find_n_sig(
#        ts, beta,
#        n_sig_step=n_sig_step,
#        first_batch_size=first_batch_size,
#        batch_size=batch_size,
#        # 10% tolerance -- let's get an estimate quick!
#        tol=tol,
#        # number of signal signal strengths (default 6, i'm really tryina rush here)
#        n_batches=n_batches
#    )

    #Get flux, and add parameters to dictionary for saving
    flux_nsig = tr.to_E2dNdE(result, E0=ref_E, unit=flux_unit)   # TeV/cm2/s  @  100TeV by default
    dnde_nsig = tr.to_dNdE(result, E0=ref_E, unit=flux_unit)     # 1/TeV/cm2/s  @  100TeV by default
    flux_nsig_at1 = tr.to_E2dNdE(result, E0=1, unit=flux_unit)   # TeV/cm2/s  @  1 TeV 
    
    result['info']['flux_nsig'] = flux_nsig
    result['info']['flux_nsig_at1'] = flux_nsig_at1
    result['info']['dnde_nsig'] = dnde_nsig
    result['info']['inj_gamma'] = g
    result['info']['flux_E0'] = ref_E
    result['info']['flux_Eunit'] = flux_unit

    save_nsig_fname = '/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/{}Info_{}weight_gamma{}_10yrPStracks_1FLEblazars_{}.pkl'.format(kind, w, g, today_date)
    f = open(save_nsig_fname, "wb")
    pickle.dump(result['info'],f)
    f.close()
    print("{} info saved ->".format(kind), save_nsig_fname)
    
    # return flux
    return flux_nsig

def plot_sensdisc():
    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.plot([2.0], fluxes_disc, linestyle='', marker='+', markersize=5, color='r', label='Discovery Potential (5sigma)')
    ax.plot(gamma, fluxes_sens, linestyle='-', color='k', label='Sensitivity')
    ax.set_xlabel('gamma')
    ax.set_ylabel('E^2 dN/dE  [TeV / cm^2 / s  @ 100 TeV]')
    ax.legend()
    ax.set_title('{} weighting'.format(w))
    ax.grid()
    plt.tight_layout()
    plt.savefig('/data/user/mcampana/analysis/Blazar_1FLE/plots/SensDisc_{}weighting_10yrPStracks_1FLEblazars_{}.png'.format(w, today_date))
    plt.close()
    
    return

#Doing Bias tests:
def bias_test():
    #Making new trials...
    #n_sigs = np.r_[:101:10]
    #trials = [tr.get_many_fits(100, n_sig=n_sig, logging=False, seed=n_sig) for n_sig in n_sigs]
    
    #Using loaded trials
    n_sigs = list(cy.bk.get_best(sig, 'weight', w, 'gamma', g, 'nsig').keys())
    n_sigs.sort()
    trials = [cy.bk.get_best(sig, 'weight', w, 'gamma', g, 'nsig', n_sig) for n_sig in n_sigs]
    b = cy.bk.get_best(bg_chi2, 'weight', w, 'gamma', g)
    trials.insert(0, b.trials)
    n_sigs.insert(0,0)
        
    for (n_sig, t) in zip(n_sigs, trials):
        t['ntrue'] = np.repeat(n_sig, len(t))

    allt = cy.utils.Arrays.concatenate(trials)

    fig, axs = plt.subplots(1, 2, figsize=(6,3))

    dns = np.mean(np.diff(n_sigs))
    ns_bins = np.r_[n_sigs - 0.5*dns, n_sigs[-1] + 0.5*dns]
    expect_kw = dict(color='C0', ls='--', lw=1, zorder=-10)
    expect_gamma = tr.sig_injs[0].flux[0].gamma
    
    ax = axs[0]
    h = hl.hist((allt.ntrue, allt.ns), bins=(ns_bins, 100))
    hl.plot1d(ax, h.contain_project(1),errorbands=True, drawstyle='default')

    lim = ns_bins[[0, -1]]
    ax.set_xlim(ax.set_ylim(lim))
    ax.plot(lim, lim, **expect_kw)
    ax.set_aspect('equal')


    ax = axs[1]
    h = hl.hist((allt.ntrue, allt.gamma), bins=(ns_bins, 100))
    hl.plot1d(ax, h.contain_project(1),errorbands=True, drawstyle='default')
    ax.axhline(expect_gamma, **expect_kw)
    ax.set_xlim(axs[0].get_xlim())

    for ax in axs:
        ax.set_xlabel('n_inj')
        ax.grid()
    axs[0].set_ylabel('n_s')
    axs[1].set_ylabel('gamma')

    plt.tight_layout()
    plt.savefig('/data/user/mcampana/analysis/Blazar_1FLE/plots/BiasTest_{}weighting_gamma{}_10yrPStracks_1FLEblazars_{}.png'.format(w, g, today_date))
    plt.close()
    
    return
#=========================================================
#=========================================================

if do_trials:
    
    for w in weighting_scheme:
        if w == 'equal':
            src_weights = np.ones(num_blaz) / num_blaz
        elif w == 'flux':
            flux_data = np.copy(data['EF30-100'])
            src_weights = flux_data / np.sum(flux_data)
            
        srcs = cy.utils.Sources(ra=data['RAdeg'], dec=data['DEdeg'], deg=True, weight=src_weights)
        
        for g in gamma:
            #Get Trials Runner
            tr = cy.get_trial_runner(src=srcs, ana=ana, flux=cy.hyp.PowerLawFlux(g), mp_cpus=cpus)
            
            for s in seed:
                #print('Doing {} BG trials for {} weight, gamma = {}, seed = {} ...'.format(num_trials,w,g,s))
                #bg_trials()  
                
                if w == 'equal':
                    if g == 1.75: 
                        n_sigs = np.r_[2:22.1:2] 
                    elif g == 2.0: 
                        n_sigs = np.r_[2:42.1:4, 60:100.1:4] 
                    elif g == 2.25: 
                        n_sigs = np.r_[10:50.1:4] 
                    elif g == 2.5: 
                        n_sigs = np.r_[40:80.1:4] 
                    elif g == 2.75: 
                        n_sigs = np.r_[70:110.1:4] 
                    elif g == 3.0:
                        n_sigs = np.r_[110:150.1:4] 
                    elif g == 3.25:
                        n_sigs = np.r_[160:200.1:4] 
                        
                elif w == 'flux':
                    if g == 1.75: 
                        n_sigs = np.r_[2:22.1:2] 
                    elif g == 2.0: 
                        n_sigs = np.r_[4:24.1:2, 40:80.1:4] 
                    elif g == 2.25: 
                        n_sigs = np.r_[2:42.1:4] 
                    elif g == 2.5: 
                        n_sigs = np.r_[20:60.1:4] 
                    elif g == 2.75: 
                        n_sigs = np.r_[40:80.1:4] 
                    elif g == 3.0:
                        n_sigs = np.r_[70:110.1:4]
                    elif g == 3.25:
                        n_sigs = np.r_[100:140.1:4] 
                    
                for n_sig in n_sigs:
                    print('Doing {} Signal trials for {} weight, gamma = {}, seed = {}, n_sig = {} ...'.format(num_trials/2,w,g,s,n_sig))
                    sig_trials()
                    
elif load_trials:
    
    #print('Loading all BG Trials ...')
    #bg_chi2 = cy.bk.get_all(
    #    # disk location
    #    bg_dir,
    #    # filename pattern
    #    'BG*.npy',
    #    # how to combine items within each directory
    #    merge=np.concatenate,
    #    # what to do with items after merge
    #    post_convert=ndarray_to_Chi2TSD)
    
    #print('Loading all SIG Trials ...')
    #sig = cy.bk.get_all(
    #    # disk location
    #    sig_dir,
    #    # filename pattern
    #    'SIG*.npy',
    #    # how to combine items within each directory
    #    merge=np.concatenate,
    #    # what to do with items after merge
    #    post_convert=cy.utils.Arrays)
            
    if do_find_sensdisc:
        for w in weighting_scheme:
            if w == 'equal':
                src_weights = np.ones(num_blaz) / num_blaz
            elif w == 'flux':
                flux_data = np.copy(data['EF30-100'])
                src_weights = flux_data / np.sum(flux_data)
            
            srcs = cy.utils.Sources(ra=data['RAdeg'], dec=data['DEdeg'], deg=True, weight=src_weights)
            
            #fluxes_sens = []
            #fluxes_disc = []
            
            for g in gamma:
                
                print('Calculating Sensitivity for: {} Weighting and Gamma={} ...'.format(w,g))
                tr = cy.get_trial_runner(src=srcs, ana=ana, flux=cy.hyp.PowerLawFlux(g), mp_cpus=cpus)
                
                f_sens = get_n_sig(beta=0.9, nsigma=None)
                print('')
                print('Sensitivity Flux: ', f_sens)
                
                if g==2.0:
                    
                    f_disc = get_n_sig(beta=0.5, nsigma=5)
                    print('')
                    print('Discovery Potential Flux: ', f_disc)
                    #fluxes_disc.append(f_disc)

                #fluxes_sens.append(f_sens)
                print('')
                
            #plot_sensdisc()
    
    if plot_TSchi2:
        for w in weighting_scheme:
            g = 2.0
            print('Plotting TS distribution with Chi2 fit for BG trials with {} weighting and gamma = {} ...'.format(w,g))
            TSchi2()
        
    if do_bias_test:

        for w in weighting_scheme:
            if w == 'equal':
                src_weights = np.ones(num_blaz) / num_blaz
            elif w == 'flux':
                flux_data = np.copy(data['EF30-100'])
                src_weights = flux_data / np.sum(flux_data)

            srcs = cy.utils.Sources(ra=data['RAdeg'], dec=data['DEdeg'], deg=True, weight=src_weights)          

            for g in gamma:
                print('Plotting bias tests for {} weighting and gamma={}'.format(w,g))
                tr = cy.get_trial_runner(src=srcs, ana=ana, flux=cy.hyp.PowerLawFlux(g), mp_cpus=cpus, update_bg=True)
                bias_test()


else:
    raise ValueError('No trials given or loaded in arguments...')


print('===========================')
print('===== Script Finished =====')
