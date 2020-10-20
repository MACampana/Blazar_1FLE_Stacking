#!/usr/bin/env python

#Comparing Sources in a few different FermiLAT Catalogs

import astropy.io.fits as fits
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib_venn as mpv

#Open fits files
fle1 = fits.open('./1fle.fits')
fhl3 = fits.open('./3fhl.fit')
fgl3 = fits.open('./3fgl.fit')

#fle1.info()
#fhl3.info()
#fgl3.info()

#Get catalogs
catList_1fle = fle1[1]
catList_3fhl = fhl3[1]
catList_3fgl = fgl3[1]

#print(catList_1fle.columns[:])
#print(catList_3fhl.columns[:])
#print(catList_3fgl.columns[:])

#Get important data columns into arrays
data_1fle = catList_1fle.data
name_1fle = np.array(data_1fle.field('Name'))
name3fgl_1fle = np.array(data_1fle.field('Assoc3FGL'))
class_1fle = np.array(data_1fle.field('Class1'))
flux30_100_1fle = np.array(data_1fle.field('EF30-100'))


data_3fhl = catList_3fhl.data
name_3fhl = np.array(data_3fhl.field('Source_Name'))
name3fgl_3fhl = np.array(data_3fhl.field('ASSOC_GAM'))
class_3fhl = np.array(data_3fhl.field('CLASS'))
fluxTotal_3fhl = np.array(data_3fhl.field('Energy_Flux')).byteswap().newbyteorder()


data_3fgl = catList_3fgl.data
name_3fgl = np.array(data_3fgl.field('Source_Name'))
class_3fgl = np.array(data_3fgl.field('CLASS1'))

#replace non 3fgl associations with separate values for 1fle and 3fhl
#This is to avoid merging dataframes on NaN values

for i in range(len(name3fgl_3fhl)):
    if '3FGL' not in name3fgl_3fhl[i]:
        name3fgl_3fhl[i] = None
        #print(name3fgl_3fhl[i])

for i in range(len(name3fgl_1fle)):
    if '3FGL' not in name3fgl_1fle[i]:
        name3fgl_1fle[i] = False
        #print(name3fgl_3fhl[i])

#2LAC(AGN), 3FHL(north only)

#Make dataframes with relevant data columns

df_3fgl = pd.DataFrame(name_3fgl)
df_3fgl.columns = ['Name_3FGL']
df_3fgl['Class_3FGL'] = class_3fgl


df_3fhl = pd.DataFrame(name_3fhl)
df_3fhl.columns = ['Name_3FHL']
df_3fhl['Class_3FHL'] = class_3fhl
df_3fhl['Name_3FGL'] = name3fgl_3fhl
df_3fhl['IntFlux_10Gto1T'] = fluxTotal_3fhl


df_1fle = pd.DataFrame(name_1fle)
df_1fle.columns = ['Name_1FLE']
df_1fle['Class_1FLE'] = class_1fle
df_1fle['Name_3FGL'] = name3fgl_1fle
df_1fle['MeVFlux_1FLE'] = flux30_100_1fle


#Merge 1fle and 3fhl on shared 3fgl names

df_1fle_3fhl = pd.merge(df_3fhl, df_1fle, on='Name_3FGL', how='inner')
#display(df_1fle_3fhl)
#with pd.option_context('display.max_rows', None, 'display.max_columns', None):
#    display(df_1fle_3fhl.sort_values('Class_1FLE'))


overlap_class_list = np.array(df_1fle_3fhl['Class_1FLE'])
unique, counts = np.unique(overlap_class_list, return_counts=True)
overlap_class_dict = dict(zip(unique, counts))
print('The following are the number of each class that is included in the overlapped sources.')
print(overlap_class_dict)

print('The 3FHL Catalog seems to have two double counted sources--ones that share a 3FGL association.')
print('These sources are not double counted in the 1FLE catalog, but will be in this comparison.')
print(df_1fle_3fhl.loc[(df_1fle_3fhl.loc[:,'Name_3FGL'] == '3FGL J1801.3-2326e') | (df_1fle_3fhl.loc[:,'Name_3FGL'] == '3FGL J1633.0-4746e')])

#Get and plot integral fluxes of the overlapped sources
#And compare to hist of full 3fhl source fluxes

flux_overlap132 = df_1fle_3fhl['IntFlux_10Gto1T'].to_numpy()
MeVFlux_overlap = df_1fle_3fhl['MeVFlux_1FLE'].to_numpy()
print(np.min(flux_overlap132), np.max(flux_overlap132))
print(np.min(fluxTotal_3fhl), np.max(fluxTotal_3fhl))

fig, ax1 = plt.subplots(facecolor='w', figsize=(8,6))

logbins = np.logspace(np.log10(1.0e-13),np.log10(1.0e-9), 25)
ax1.hist(fluxTotal_3fhl, bins=logbins, alpha=1, color='r', label='3FHL', histtype='step', linestyle='-')
ax1.hist(flux_overlap132, bins=logbins, alpha=1, color='r', label='Shared', histtype='step', linestyle='--')

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylabel('Counts')
ax1.set_xlabel('Integrated Flux (10GeV-1TeV) [erg/cm^2/s]', color='r')
ax1.tick_params(axis='x', labelcolor='r')

#ax2 = ax1.twiny()

#logbins = np.logspace(np.log10(6.0e-8),np.log10(1.0e-4), 25)
#ax2.hist(flux30_100_1fle, bins=logbins, alpha=1, color='b', label='1FLE', histtype='step', linestyle='-')

#ax2.set_xscale('log')
#ax2.set_yscale('log')
#ax2.set_ylabel('Counts')
#ax2.set_xlabel('Integrated Photon Flux (30-100 MeV) [1/cm^2/s]', color='b')
#ax2.tick_params(axis='x', labelcolor='b')

plt.legend(loc='best')
plt.savefig('./plots/HE_FluxPDF_3FHL_wOverlap.png')
plt.close()

fig, ax = plt.subplots(facecolor='w', figsize=(6,6))
ax.scatter(MeVFlux_overlap/np.sum(flux30_100_1fle), flux_overlap132/np.sum(fluxTotal_3fhl), marker='+', color='purple')
ax.set_xlabel('1FLE Normalized Flux (30-100 MeV)', color='b')
ax.set_ylabel('3FHL Normalized Flux (10 GeV-100 TeV)', color='r')
#plt.xlim(-0.001,.2)
#plt.ylim(-0.001,.2)

lim = .00001, 1
ax.plot(lim, lim, 'k--', zorder=-10)
ax.set(xlim=lim, ylim=lim, aspect='equal')

ax.set_yscale('log')
ax.set_xscale('log')

ax.set_title('Normalized Fluxes of 1FLE/3FHL Shared Sources')
plt.savefig('./plots/NormFlux_OverlapComparison_3FHL.png')
plt.close()

plt.figure(facecolor='w', figsize=(6,6))
plt.title('3FGL Associated Sources')
mpv.venn2([set(name3fgl_1fle),set(name3fgl_3fhl)], set_labels=('1FLE', '3FHL'), set_colors=('b','r'))
plt.savefig('./plots/VennOverlap_3FHL.png')
plt.close()

fle1.close()
fhl3.close()