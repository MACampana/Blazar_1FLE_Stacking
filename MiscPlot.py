#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import csky as cy

bg_dir = '/data/user/mcampana/analysis/Blazar_1FLE/trials/bg'

print('Compiling BG Trials...')
bg_trials = cy.bk.get_all(
        # disk location
        bg_dir,
        # filename pattern
        'BG*.npy',
        # how to combine items within each directory
        merge=np.concatenate,
        # what to do with items after merge
        post_convert=cy.utils.Arrays)

g = 2.0

bg_eq = cy.bk.get_best(bg_trials, 'weight', 'equal', 'gamma', g)
bg_fx = cy.bk.get_best(bg_trials, 'weight', 'flux', 'gamma', g)
    
plt.figure(figsize=(8,6))
h, bins, m = plt.hist(bg_eq['ns'], bins='auto', histtype='step', linewidth=2, label='Equal')
plt.hist(bg_fx['ns'], bins=bins, histtype='step', linewidth=2, label='Flux')
plt.semilogy()
plt.title('BG Trial n_s Distribution')
plt.ylabel('Number of BG Trials')
plt.xlabel('Fitted n_s')
plt.legend()
plt.savefig('/data/user/mcampana/analysis/Blazar_1FLE/plots/BG_ns_hist.png')
    
plt.clf()
    
plt.figure(figsize=(8,6))
h, bins, m = plt.hist(bg_eq[bg_eq['ns']>0]['gamma'], bins='auto', histtype='step', linewidth=2, label='Equal')
plt.hist(bg_fx[bg_fx['ns']>0]['gamma'], bins=bins, histtype='step', linewidth=2, label='Flux')
plt.semilogy()
plt.title('BG Trial Spectral Index Distribution')
plt.ylabel('Number of BG Trials')
plt.xlabel('Fitted $\gamma$')
plt.legend()
plt.savefig('/data/user/mcampana/analysis/Blazar_1FLE/plots/BG_gamma_hist.png')
