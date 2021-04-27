#!/usr/bin/env python

import numpy as np
import astropy as ap
import astropy.io.fits as fits
import healpy as hp
import matplotlib.pyplot as plt
import pandas as pd
import pickle as pkl

data = pkl.load( open('./u1FLE_Blazars_Catalog.pkl', 'rb') )
RAdeg = np.array(data['RAdeg'])
DEdeg = np.array(data['DEdeg'])
Flux = np.array(data['EF30-100'])

coord = 'C'
rot = (-180,0,0)
size = Flux * 1e12 / np.min(Flux * 1e12) * 5

plt.subplots()

hp.visufunc.mollview(np.zeros(12) + np.inf, cbar=False, coord=coord, rot=rot, margins=(1,1,1,1))
hp.graticule(30, color='k', alpha=0.7)
plt.title('1FLE Blazars, 30-100 MeV Flux Scaling', pad=25)

hp.projplot(np.arange(0,360), np.zeros(360), coord='G', color='r', lonlat=True, alpha=.5)
hp.projscatter(266.42, -28.99, lonlat=True, coord=coord, marker='+', color='k', label='Galactic Center')

hp.projscatter(RAdeg, DEdeg, lonlat=True, coord=coord, s=size, alpha=0.8, c='b')

hp.projtext(0,0,'90', size='large', ha='center', va='bottom', coord=coord)
hp.projtext(np.pi/2,2.0*np.pi - .001,'360', size='large', ha='right', coord=coord)
hp.projtext(np.pi/2,0,'0', size='large', ha='left', coord=coord)
hp.projtext(np.pi,0,'-90', size='large', ha='center', va='top', coord=coord)

plt.savefig('./plots/u1FLE_FluxScaled_SkyMap.png', bbox_inches='tight')
plt.close()