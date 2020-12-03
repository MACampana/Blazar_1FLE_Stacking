#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
#import csky as cy
import pandas as pd

df_1fle = pd.read_csv('1FLE_Blazars_Data.csv')

pFlux = df_1fle['F30-100']
eFlux = df_1fle['EF30-100']

fig, ax = plt.subplots(facecolor='w', figsize=(6,6))
ax.scatter(pFlux/np.sum(pFlux), eFlux/np.sum(eFlux), marker='+', color='b')
ax.set_title('1FLE Blazars Flux Weight Comparison')
ax.set_ylabel('Normalized Energy Flux')
ax.set_xlabel('Normalized Photon Flux')

lim = 0, .1
ax.plot(lim, lim, 'k--', zorder=-10, alpha=.3)
ax.set(xlim=lim, ylim=lim, aspect='equal')

plt.savefig('./plots/1FLE_EvPFlux.png')

print(np.max(np.abs((pFlux/np.sum(pFlux)) - (eFlux/np.sum(eFlux)))))
print(np.min(np.abs((pFlux/np.sum(pFlux)) - (eFlux/np.sum(eFlux)))))