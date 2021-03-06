#!/usr/bin/env python

# See https://icecube.wisc.edu/~mrichman/docs/csky/about.html

#=========================================================
#       IMPORTS
#=========================================================
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

import argparse

#=========================================================
#       ARGUMENTS
#=========================================================

parser = argparse.ArgumentParser(description='Perform stacking analysis on 1FLE Blazars.')

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
parser.add_argument('-s', '--sens', action='store_true', dest='do_find_sens', help='Use this option to do sensitivity calculation.')
parser.add_argument('-d', '--disc', action='store_true', dest='do_find_disc', help='Use this option to do discovery potential calculation.')
parser.add_argument('-b', '--bias', action='store_true', dest='do_bias_test', help='Use this option to perform and plot bias tests, for given weights and gammas.')
parser.add_argument('-i', '--chi2', action='store_true', dest='plot_TSchi2', help='Use this option to plot the TS dist and chi2 fit of loaded trials.')
parser.add_argument('--diff-sens', action='store_true', dest='diff_sens', help='Use this option to calculate differential sensitvities.')
parser.add_argument('--hemi', default='both', choices=['both', 'north', 'south'], type=str, dest='hemisphere', help='ONE of both, north, OR south. Hemisphere for sources.')
parser.add_argument('--no-poisson', action='store_false', dest='poisson', help='Use this option to NOT inject with poisson distribution during signal trial making.')
parser.add_argument('--dataset', default='ps-v3p2', choices=['ps-v3p2', 'ps-v4'], type=str, dest='dataset', help='Dataset.')
parser.add_argument('--ecut', action='store_true', dest='ecut', help='Use this option to set limits on MC energies (then set the limits with --elo and --ehi).')
parser.add_argument('--elo', default=0.0, type=float, dest='elo', help='Low Energy cutoff for MC.')
parser.add_argument('--ehi', default=np.inf, type=float, dest='ehi', help='High Energy cutoff for MC.')
parser.add_argument('--no-bg', action='store_false', dest='do_bg', help='Use this option to NOT perform background trials.')

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
do_find_sens = args.do_find_sens
do_find_disc = args.do_find_disc
do_bias_test = args.do_bias_test
plot_TSchi2 = args.plot_TSchi2
cpus = args.cpus
diff_sens = args.diff_sens
hemisphere = args.hemisphere
poisson = args.poisson
dSetName = args.dataset
E_cut = args.ecut
ethresh_high = args.ehi
ethresh_low = args.elo
do_bg_trials = args.do_bg

#=========================================================
#       CATALOG/DATA EXTRACTION
#=========================================================

#FOR LOADING FITS FILE (NOT currently in use)
if False:
    print('Loading Fits...')
    #Open Catalog fits file
    fits_1fle = fits.open('/data/user/mcampana/analysis/Blazar_1FLE/1fle.fits')
    #Get Column Names
    cols = fits_1fle[1].data.columns.names
    #Extract Blazars (all, north, south)
    if hemisphere == 'both':
        b_mask = (fits_1fle[1].data['CLASS1'] == 'bll') | (fits_1fle[1].data['CLASS1'] == 'fsrq')
    elif hemisphere == 'north':
        b_mask = ((fits_1fle[1].data['CLASS1'] == 'bll') | (fits_1fle[1].data['CLASS1'] == 'fsrq')) & (fits_1fle[1].data['DEdeg'] > -5.0)
    elif hemisphere == 'south':
        b_mask = ((fits_1fle[1].data['CLASS1'] == 'bll') | (fits_1fle[1].data['CLASS1'] == 'fsrq')) & (fits_1fle[1].data['DEdeg'] <= -5.0)

    blaz_1fle = fits_1fle[1].data[b_mask]

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
        pd.DataFrame.from_dict(data).to_csv('/data/user/mcampana/analysis/Blazar_1FLE/1FLE_Blazars_Data_{}.csv'.format(hemisphere))
    elif data_save == 'txt':
        f = open("/data/user/mcampana/analysis/Blazar_1FLE/1FLE_Blazars_Data_{}.txt".format(hemisphere),"w")
        f.write(str(data))
        f.close()
    elif data_save == None:
        pass

#FOR LOADING PICKLED DICTIONARY (currently used)
if True:
    if hemisphere == 'both':
        data = pickle.load( open('/data/user/mcampana/analysis/Blazar_1FLE/u1FLE_Blazars_Catalog.pkl', 'rb') )
    elif hemisphere == 'north':
        data = pickle.load( open('/data/user/mcampana/analysis/Blazar_1FLE/u1FLE_North_Blazars_Catalog.pkl', 'rb') )
    elif hemisphere == 'south':
        data = pickle.load( open('/data/user/mcampana/analysis/Blazar_1FLE/u1FLE_South_Blazars_Catalog.pkl', 'rb') )

#For Weighting
num_blaz = len(data['Name'])
print("Stacking {} 1FLE Blazars in the {} hemisphere...".format(num_blaz, hemisphere))

print('Considering plotting flux pdf...')
#Plot the Flux weight PDF
if plot_flux_pdf:
    plt.figure(figsize=(6,4))
    plt.hist(data['EF30-100'], bins='auto', histtype='step', linewidth=2)
    plt.xlabel("Energy Flux from 30-100 MeV [erg /cm^2 /s]")
    plt.ylabel("Counts")
    plt.title("Energy Flux PDF for u1FLE Blazars")
    #plt.semilogy()
    plt.savefig('/data/user/mcampana/analysis/Blazar_1FLE/plots/EnergyFlux_30-100MeV_PDF_{}_{}.png'.format(hemisphere,today_date))
    plt.close()
    sys.exit("Plotted Flux PDF, exiting program...")

#=========================================================
#       CSKY SETUP
#=========================================================

if poisson:
    pDirName = 'poisson'
else:
    pDirName = 'nopoisson'
    
#********************Data selection***********************
# Currently using PS Tracks v3p2, because v4 is not fully ready

print('Getting data selection...')

if dSetName == 'ps-v3p2': #Previously calling 10yrPStracks
    selection = cy.selections.PSDataSpecs.ps_10yr
    
    if not E_cut:
        trials_dir = cy.utils.ensure_dir('/data/user/mcampana/analysis/Blazar_1FLE/trials_u1FLE')
    elif E_cut:
        
        # --------Setting Energy Range for truncated MC energies--------
        print('Truncating MC energy range to {}-{} GeV'.format(ethresh_low, ethresh_high))
        
        trials_dir = cy.utils.ensure_dir('/data/user/mcampana/analysis/Blazar_1FLE/trials_Erange_{}-{}'.format(ethresh_low, ethresh_high))
        
    version = 'version-003-p02'
    
elif dSetName == 'ps-v4':
    selection = cy.selections.PSDataSpecs.ps_v4
    trials_dir = cy.utils.ensure_dir('/data/user/mcampana/analysis/Blazar_1FLE/trials_psv4')
    version = 'version-004-p00'

print('Getting analysis arrays...')
#Create analysis with data selection (saved in working dir)
#ana = cy.get_analysis(cy.selections.repo, selection)
#ana.save('./')
ana_dir = cy.utils.ensure_dir('/data/user/mcampana/analysis/Blazar_1FLE/subanalyses/{}'.format(dSetName))
ana = cy.get_analysis(cy.selections.repo, version, selection, dir=ana_dir)

#Trial save and load directories
sig_dir = cy.utils.ensure_dir('{}/sig'.format(trials_dir))
if not E_cut:
    bg_dir = cy.utils.ensure_dir('{}/bg'.format(trials_dir))
elif E_cut:
    bg_dir = cy.utils.ensure_dir('/data/user/mcampana/analysis/Blazar_1FLE/trials_u1FLE/bg')

#=========================================================
#       FUNCTION DEFINITIONS 
#=========================================================

#Do BG Trials
def bg_trials():
    
    trials = tr.get_many_fits(num_trials, seed=s)
    #Save Trials to numpy array file
    dir_ = cy.utils.ensure_dir('{}/hemisphere/{}/weight/{}'.format(bg_dir,h,w))
    save_trials_fname = '{}/BG_{}trials_u1FLEblazars_seed{}_{}.npy'.format(dir_, num_trials, s, today_date)
    np.save(save_trials_fname, trials.as_array)
    print("BG Trials saved ->", save_trials_fname)
    
    return

#Do Signal Trials
def sig_trials():
    
    trials = tr.get_many_fits(num_trials/2, n_sig, seed=s, poisson=poisson)
    #Save Trials to numpy array file    
    dir_ = cy.utils.ensure_dir('{}/{}/hemisphere/{}/weight/{}/gamma/{}/nsig/{}'.format(sig_dir,pDirName,h,w,g,n_sig))
    save_trials_fname = '{}/SIG_{}trials_u1FLEblazars_seed{}_{}.npy'.format(dir_, int(num_trials/2), s, today_date)
    np.save(save_trials_fname, trials.as_array)
    print("Signal Trials saved ->", save_trials_fname)
      
    return

#Plot TS and chi2 distribution

def TSchi2():
    fig, ax = plt.subplots()

    b = cy.bk.get_best(bg_chi2, 'hemisphere', h, 'weight', w, 'gamma', g)
    hst = b.get_hist(bins=30)
    hl.plot1d(ax, hst, crosses=True, label='{} bg trials'.format(b.n_total))

    x = hst.centers[0]
    norm = hst.integrate().values
    ax.semilogy(x, norm * b.pdf(x), lw=1, ls='--',
                label='Chi2 fit: {} dof, eta={}'.format(np.round(b.ndof,3), np.round(b.eta,3)))

    ax.set_xlabel('TS')
    ax.set_ylabel('number of trials')
    ax.legend()
    ax.set_title('{} Weighted, {} Hemisphere(s)'.format(w.capitalize(), h.capitalize()))

    plt.tight_layout()
    plt.savefig('/data/user/mcampana/analysis/Blazar_1FLE/plots/Chi2TS_{}trials_{}weighting_{}hemi_{}_{}_u1FLEblazars_{}.png'.format(b.n_total, w, h, pDirName, dSetName, today_date))
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
    plt.savefig('/data/user/mcampana/analysis/Blazar_1FLE/plots/Chi2TS_{}weighting_{}hemi_{}_{}_u1FLEblazars_{}.png'.format(w,h, pDirName, dSetName, today_date))
    plt.close()
    
    return


def ndarray_to_Chi2TSD(tri):
    return cy.dists.Chi2TSD(cy.utils.Arrays(tri))

def get_n_sig(beta=0.9, nsigma=None):
    # get signal trials, background distribution, and trial runner
    sig_trials = cy.bk.get_best(sig, pDirName, 'hemisphere', h, 'weight', w, 'gamma', g, 'nsig')
    b = cy.bk.get_best(bg_chi2, 'hemisphere', h, 'weight', w)
    # determine ts threshold
    if nsigma is not None:
        ts = b.isf_nsigma(nsigma)
        kind = 'Disc{}'.format(nsigma)
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
#        tol=tol,
#        # number of signal signal strengths
#        n_batches=n_batches
#    )

    #Get flux, and add parameters to dictionary for saving
    flux_nsig = tr.to_E2dNdE(result, E0=ref_E, unit=flux_unit)   # TeV/cm2/s  @  100TeV by default when not diff_sens
    
    if not diff_sens:
        flux_nsig_at1 = tr.to_E2dNdE(result, E0=1, unit=flux_unit)   # TeV/cm2/s  @  1 TeV 
        result['info']['flux_nsig_at1'] = flux_nsig_at1
        dnde_nsig = tr.to_dNdE(result, E0=ref_E, unit=flux_unit)     # 1/TeV/cm2/s  @  100TeV by default
        result['info']['dnde_nsig'] = dnde_nsig
        
    result['info']['flux_nsig'] = flux_nsig
    result['info']['inj_gamma'] = g
    result['info']['flux_E0'] = ref_E
    result['info']['flux_Eunit'] = flux_unit
    
    sensSaveDir = cy.utils.ensure_dir('/data/user/mcampana/analysis/Blazar_1FLE/sens_disc/{}/hemisphere/{}/weight/{}'.format(pDirName, h, w))    
    if diff_sens:
        save_nsig_fname = '{}/Diff{}Info_gamma{}_Emin{}_Emax{}_{}_u1FLEblazars_{}.pkl'.format(sensSaveDir, kind, g, np.round(E_min,1), np.round(E_max,1), dSetName, today_date)
    elif E_cut:
        save_nsig_fname = '{}/{}Info_gamma{}_Erange_{}-{}_{}_u1FLEblazars_{}.pkl'.format(sensSaveDir, kind, g, ethresh_low, ethresh_high, dSetName, today_date) 
    else:
        save_nsig_fname = '{}/{}Info_gamma{}_{}_u1FLEblazars_{}.pkl'.format(sensSaveDir, kind, g, dSetName, today_date) 
    
    f = open(save_nsig_fname, "wb")
    pickle.dump(result,f)
    f.close()
    print("{} dict saved ->".format(kind), save_nsig_fname)
    
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
    plt.savefig('/data/user/mcampana/analysis/Blazar_1FLE/plots/SensDisc_{}weighting_{}hemi_{}_{}_u1FLEblazars_{}.png'.format(w,h,pDirName, dSetName, today_date))
    plt.close()
    
    return

#Doing Bias tests:
def bias_test():
    #Making new trials...
    n_sigs = np.r_[:201:10]
    trials = [tr.get_many_fits(1000, n_sig=n_sig, logging=False, seed=n_sig, poisson=False) for n_sig in n_sigs]
        
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

    plt.title('Gamma={}, {} Weighted'.format(g,w.capitalize()))
    plt.tight_layout()
    plt.savefig('/data/user/mcampana/analysis/Blazar_1FLE/plots/BiasTest_{}weighting_gamma{}_{}hemi_{}_{}_u1FLEblazars_{}.png'.format(w, g, hemisphere, pDirName, dSetName, today_date))
    plt.close()
    
    return
#=========================================================
#      DOING THE WORK
#=========================================================

h = hemisphere

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
            
            #Set signal flux (for injections)
            #   Power law from gamma, with or without some cutoff (for determining sensitivity changes and energy range)
            if not E_cut:
                flux = cy.hyp.PowerLawFlux(g)
            elif E_cut:
                flux = cy.hyp.PowerLawFlux(g, energy_range=(ethresh_low, ethresh_high))
                
            tr = cy.get_trial_runner(src=srcs, ana=ana, flux=flux, mp_cpus=cpus)
            
            for s in seed:
                print('Doing {} BG trials for {} data, {} weight, gamma = {}, seed = {} ...'.format(num_trials,dSetName,w,g,s))
                
                if do_bg_trials:
                    if g == 2.0:
                        bg_trials()
                
                if w == 'equal':
                    if g == 1.75: 
                        n_sigs = np.r_[2:22.1:2, 30:70.1:4] 
                    elif g == 2.0: 
                        n_sigs = np.r_[2:42.1:4, 60:100.1:4] 
                    elif g == 2.25: 
                        n_sigs = np.r_[10:50.1:4, 120:160.1:4]
                    elif g == 2.5: 
                        n_sigs = np.r_[40:80.1:4, 200:240.1:4] 
                    elif g == 2.75: 
                        n_sigs = np.r_[70:110.1:4, 330:370.1:4]
                    elif g == 3.0:
                        n_sigs = np.r_[110:150.1:4, 520:560.1:4] 
                    elif g == 3.25:
                        n_sigs = np.r_[160:200.1:4, 700:740.1:4] 
                        
                elif w == 'flux':
                    if g == 1.75: 
                        n_sigs = np.r_[2:22.1:2, 20:60.1:4] 
                    elif g == 2.0: 
                        n_sigs = np.r_[4:24.1:2, 40:80.1:4] 
                    elif g == 2.25: 
                        n_sigs = np.r_[2:42.1:4, 80:120.1:4] 
                    elif g == 2.5: 
                        n_sigs = np.r_[20:60.1:4, 140:180.1:4] 
                    elif g == 2.75: 
                        n_sigs = np.r_[40:80.1:4, 230:270.1:4] 
                    elif g == 3.0:
                        n_sigs = np.r_[70:110.1:4, 340:380.1:4]
                    elif g == 3.25:
                        n_sigs = np.r_[100:140.1:4, 440:480.1:4] 
                        
                    
                for n_sig in n_sigs:
                    print('Doing {} Signal trials for {} data, {} hemisphere, {} weight, gamma = {}, seed = {}, n_sig = {} ...'.format(num_trials/2,dSetName,h,w,g,s,n_sig))
                    sig_trials()
                    
elif load_trials:
    
    print('Loading all BG Trials ...')
    bg_chi2 = cy.bk.get_all(
        # disk location
        bg_dir,
        # filename pattern
        'BG*.npy',
        # how to combine items within each directory
        merge=np.concatenate,
        # what to do with items after merge
        post_convert=ndarray_to_Chi2TSD)
    
    print('Loading all SIG Trials ...')
    sig = cy.bk.get_all(
        # disk location
        sig_dir,
        # filename pattern
        'SIG*.npy',
        # how to combine items within each directory
        merge=np.concatenate,
        # what to do with items after merge
        post_convert=cy.utils.Arrays)
            
    if do_find_sens or do_find_disc:
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
                
                if diff_sens:
                    
                    if w == 'equal':
                        #Energy cutoffs for 90% sens range determined separately (in GeV to the nearest 100 GeV)
                        E_low_cut = 14800.0
                        E_high_cut = 14797300.0
                    elif w == 'flux':
                        #Energy cutoffs for 90% sens range determined separately (in GeV to the nearest 100 GeV)
                        E_low_cut = 10500.0
                        E_high_cut = 12188900.0
                        
                    #logE_range = np.logspace(np.log10(E_low_cut), np.log10(E_high_cut), 21) #n bins in log (E/ 1 GeV)
                    logE_range = np.logspace(4, 7.5, 15)    #E range rounded for plotting conventions
                    for i in range(len(logE_range[:-1])):
                        E_min = logE_range[i]
                        E_max = logE_range[i+1]
                        
                        ref_E = E_min * 10**-3         #in TeV
                        flux = cy.hyp.PowerLawFlux(gamma=g, energy_range=(E_min, E_max))

                        print('Calculating Differential Sensitivity for: Emin={}, {} data, {} Hermisphere, {} Weighting, and Gamma={} ...'.format(E_min,dSetName,h,w,g))
                        tr = cy.get_trial_runner(src=srcs, ana=ana, flux=flux, mp_cpus=cpus, update_bg=True)

                        f_sens = get_n_sig(beta=0.9, nsigma=None)
                        print('')
                        print('Sensitivity Flux: ', f_sens)
                    
                else:
                    
                    if not E_cut:
                        flux = cy.hyp.PowerLawFlux(g)
                    elif E_cut:
                        flux = cy.hyp.PowerLawFlux(g, energy_range=(ethresh_low, ethresh_high))
                
                    print('Calculating Sensitivity for: {} data, {} Hemisphere, {} Weighting, and Gamma={} ...'.format(dSetName,h,w,g))
                    tr = cy.get_trial_runner(src=srcs, ana=ana, flux=flux, mp_cpus=cpus, update_bg=True)
                    
                    if do_find_sens:

                        f_sens = get_n_sig(beta=0.9, nsigma=None)
                        print('')
                        print('Sensitivity Flux: ', f_sens)

                    if do_find_disc:

                        f_disc = get_n_sig(beta=0.5, nsigma=5)
                        print('')
                        print('Discovery Potential Flux: ', f_disc)
                        
                    print('')
    
    if plot_TSchi2:
        for w in weighting_scheme:
            g = 2.0
            print('Plotting TS distribution with Chi2 fit for {} data, {} hemisphere BG trials with {} weighting ...'.format(dSetName,h,w))
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
                print('Plotting bias tests for {} data, {}hemisphere, {} weighting, and gamma={}'.format(dSetName,h,w,g))
                tr = cy.get_trial_runner(src=srcs, ana=ana, flux=cy.hyp.PowerLawFlux(g), mp_cpus=cpus, update_bg=True)
                bias_test()


else:
    raise ValueError('No trials given or loaded in arguments...')


print('===========================')
print('===== Script Finished =====')
