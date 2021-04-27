#!/usr/bin/env python
#Comparing Sources in a few different FermiLAT Catalogs

#=====IMPORTS=====
import astropy.io.fits as fits
import astropy as ap
import pandas as pd
import pickle 
import numpy as np
import matplotlib_venn as mpv
from icecube import astro
import argparse
import matplotlib.pyplot as plt
plt.ioff()

#=====Arguments=====
parser = argparse.ArgumentParser(description='Compare 1FLE Catalog to 3FHL and 2LAC.')

parser.add_argument('-s', action='store_true', dest='save_cat', help='Use this option to save the constructed updated 1FLE catalog as pickled dict (Saves to PWD).')
parser.add_argument('--3FHL', action='store_true', dest='do_3FHL', help='Use this option to compare to 3FHL.')
parser.add_argument('--2LAC', action='store_true', dest='do_2LAC', help='Use this option to compare to 2LAC.')
parser.add_argument('--store', default='./plots', type=str, dest='direc', help='Directory to save plots.')

args = parser.parse_args()

save_cat = args.save_cat
do_3FHL = args.do_3FHL
do_2LAC = args.do_2LAC
direc = args.direc

#=====Get 1FLE data and 3FGL (For associating)=====
#Open fits files
fle1 = fits.open('./1fle.fits')
fgl3 = fits.open('./3fgl.fit')

#Get catalogs
catList_1fle = fle1[1]
catList_3fgl = fgl3[1]

#Get important data columns into arrays
data_1fle = catList_1fle.data
name_1fle = np.array(data_1fle.field('Name'))
name3fgl_1fle = np.array(data_1fle.field('Assoc3FGL'))
class_1fle = np.array(data_1fle.field('Class1'))
flux30_100_1fle = np.array(data_1fle.field('EF30-100'))
unc_flux30_100_1fle = np.array(data_1fle.field('e_EF30-100'))
ph_flux30_100_1fle = np.array(data_1fle.field('F30-100'))
RA_1fle = np.array(data_1fle.field('RAdeg'))
DEC_1fle = np.array(data_1fle.field('DEdeg'))
z_1fle = np.array(data_1fle.field('z'))

data_3fgl = catList_3fgl.data
name_3fgl = np.array(data_3fgl.field('Source_Name'))
class_3fgl = np.array(data_3fgl.field('CLASS1'))
name2fgl_3fgl = np.array(data_3fgl.field('2FGL_Name'))
ph_flux30_100_3fgl = np.array(data_3fgl.field('Flux30_100')).byteswap().newbyteorder()
RA_3fgl = np.array(data_3fgl.field('RAJ2000')).byteswap().newbyteorder()
DEC_3fgl = np.array(data_3fgl.field('DEJ2000')).byteswap().newbyteorder()

#Make dataframes with relevant data columns
df_3fgl = pd.DataFrame(name_3fgl)
df_3fgl.columns = ['Name_3FGL']
df_3fgl['Class_3FGL'] = class_3fgl
df_3fgl['MeVPhoFlux_3FGL'] = ph_flux30_100_3fgl
df_3fgl['RA_3FGL'] = RA_3fgl
df_3fgl['DEC_3FGL'] = DEC_3fgl
fgl3.close()

#replace non 3fgl associations with separate values for 1fle
#This is to avoid merging dataframes on NaN values
for i in range(len(name3fgl_1fle)):
    if '3FGL' not in name3fgl_1fle[i]:
        name3fgl_1fle[i] = False
        #print(name3fgl_3fhl[i])
        
df_1fle = pd.DataFrame(name_1fle)
df_1fle.columns = ['Name_1FLE']
df_1fle['Class_1FLE'] = class_1fle
df_1fle['Name_3FGL'] = name3fgl_1fle
df_1fle['MeVFlux_1FLE'] = flux30_100_1fle
df_1fle['Unc_MeVFlux_1FLE'] = unc_flux30_100_1fle
df_1fle['MeVPhoFlux_1FLE'] = ph_flux30_100_1fle
df_1fle['RA_1FLE'] = RA_1fle
df_1fle['DEC_1FLE'] = DEC_1fle
df_1fle['z_1FLE'] = z_1fle
fle1.close()


#=====Merge 1FLE with 4FGL to update Classifications=====
#Open fits
fgl4 = fits.open('./4FGL_DR2.fit') 

#Put needed data into arrays
name_4fgl = np.array(fgl4[1].data.field('Source_Name'))
ph_flux50_100_4fgl = np.array(fgl4[1].data.field('Flux_Band')[:,0]).byteswap().newbyteorder()
flux50_100_4fgl = np.array(fgl4[1].data.field('nuFnu_Band')[:,0]).byteswap().newbyteorder()
RA_4fgl = np.array(fgl4[1].data.field('RAJ2000')).byteswap().newbyteorder()
DEC_4fgl = np.array(fgl4[1].data.field('DEJ2000')).byteswap().newbyteorder()
name3fgl_4fgl = np.array(fgl4[1].data.field('ASSOC_FGL'))
class_4fgl = np.array(fgl4[1].data.field('CLASS1'))

#replace non 3fgl associations with separate values for 4fgl
#This is to avoid merging dataframes on NaN values

for i in range(len(name3fgl_4fgl)):
    if '3FGL' not in name3fgl_4fgl[i]:
        name3fgl_4fgl[i] = None
        
#Dataframe it and close file
df_4fgl = pd.DataFrame(name_4fgl)
df_4fgl.columns = ['Name_4FGL']
df_4fgl['Class_4FGL'] = class_4fgl
df_4fgl['Name_3FGL'] = name3fgl_4fgl
df_4fgl['MeVFlux_4FGL'] = flux50_100_4fgl
df_4fgl['MeVPhoFlux_4FGL'] = ph_flux50_100_4fgl
df_4fgl['RA_4FGL'] = RA_4fgl
df_4fgl['DEC_4FGL'] = DEC_4fgl
fgl4.close()

df_1fle_4fgl = pd.merge(df_4fgl, df_1fle, on='Name_3FGL', how='inner')

#This data frame will be used for class comparisons. Additionally, this data frame could be 
#saved and used for comparisons on its own if desired.

if save_cat:
    #Select blazars only
    DF = df_1fle_4fgl
    blazMask = ((DF['Class_4FGL'].str.lower().str.strip() == 'fsrq') | (DF['Class_4FGL'].str.lower().str.strip() == 'bll') | (DF['Class_1FLE'].str.lower().str.strip() == 'fsrq') | (DF['Class_1FLE'].str.lower().str.strip() == 'bll'))
    blaz = DF[blazMask]
    
    #Select Columns of Interest (Dataframe)
    cols = ['Name_1FLE', 'Name_4FGL', 'Class_4FGL', 'RA_1FLE', 'DEC_1FLE', 'MeVFlux_1FLE', 'Unc_MeVFlux_1FLE', 'z_1FLE']
    u1FLE = blaz[cols].copy()
    u1FLE.rename(columns={'Name_1FLE': 'Name', 'Class_4FGL': 'CLASS1', 'RA_1FLE': 'RAdeg', 'DEC_1FLE': 'DEdeg', 'MeVFlux_1FLE': 'EF30-100', 'Unc_MeVFlux_1FLE': 'e_EF30-100'}, inplace=True)
    u1FLE['CLASS1'] = u1FLE['CLASS1'].apply(lambda s: s.lower().strip())
    u1FLE.reset_index(drop=True, inplace=True)
    
    #Make into dictionary
    u1FLE_d = u1FLE.to_dict(orient='list')
    
    #Also make dictionaries for South and North sources
    nmask = (u1FLE['DEdeg'] >= -5.0)
    smask = (u1FLE['DEdeg'] < -5.0)

    u1FLE_s = u1FLE[smask]
    u1FLE_n = u1FLE[nmask]

    u1FLE_s_d = u1FLE_s.to_dict(orient='list')
    u1FLE_n_d = u1FLE_n.to_dict(orient='list')
    
    #Save as pickled files
    pickle.dump(u1FLE_d, open( './u1FLE_Blazars_Catalog.pkl' , 'wb' ) )
    pickle.dump(u1FLE_s_d, open( './u1FLE_South_Blazars_Catalog.pkl' , 'wb' ) )
    pickle.dump(u1FLE_n_d, open( './u1FLE_North_Blazars_Catalog.pkl' , 'wb' ) )
    

#=====3FHL Comparison=====
if do_3FHL:
    print('Doing 3FHL Comparison...')

    #Open fits and get catalog
    fhl3 = fits.open('./3fhl.fit')
    catList_3fhl = fhl3[1]

    #Put data into arrays
    data_3fhl = catList_3fhl.data
    name_3fhl = np.array(data_3fhl.field('Source_Name'))
    name3fgl_3fhl = np.array(data_3fhl.field('ASSOC_GAM'))
    class_3fhl = np.array(data_3fhl.field('CLASS'))
    fluxTotal_3fhl = np.array(data_3fhl.field('Energy_Flux')).byteswap().newbyteorder() #erg/cm^2/s

    ph_fluxTotal_3fhl = np.array(data_3fhl.field('Flux')).byteswap().newbyteorder() #/cm^2/s
    pivE_3fhl = np.array(data_3fhl.field('Pivot_Energy')) #GeV
    specI_3fhl = np.array(data_3fhl.field('Spectral_Index')) 
    fluxD_3fhl = np.array(data_3fhl.field('Flux_Density')) #/cm^2/GeV/s
    specType_3fhl = np.array(data_3fhl.field('SpectrumType')) 

    #replace non 3fgl associations with separate values for 3fhl
    #This is to avoid merging dataframes on NaN values

    for i in range(len(name3fgl_3fhl)):
        if '3FGL' not in name3fgl_3fhl[i]:
            name3fgl_3fhl[i] = None
            #print(name3fgl_3fhl[i])

    #Dataframe and close fits
    df_3fhl = pd.DataFrame(name_3fhl)
    df_3fhl.columns = ['Name_3FHL']
    df_3fhl['Class_3FHL'] = class_3fhl
    df_3fhl['Name_3FGL'] = name3fgl_3fhl
    df_3fhl['IntFlux_10Gto1T'] = fluxTotal_3fhl
    df_3fhl['IntPhoFlux_10Gto1T'] = ph_fluxTotal_3fhl
    fhl3.close()

    #Merge 1fle and 3fhl on shared 3fgl names
    df_1fle_3fhl = pd.merge(df_3fhl, df_1fle_4fgl, on='Name_3FGL', how='inner')

    #Get counts of overlapped classes
    overlap_class_list = np.array(df_1fle_3fhl['Class_4FGL'].apply(lambda s: s.lower().strip()))
    unique, counts = np.unique(overlap_class_list, return_counts=True)
    overlap_class_dict = dict(zip(unique, counts))
    print('The following are the number of each class that is included in the overlapped sources as per 4FGL.')
    print(overlap_class_dict)

    #Get and plot integral fluxes of the overlapped sources
    #And compare to hist of full 3fhl source fluxes
    flux_overlap132 = df_1fle_3fhl['IntFlux_10Gto1T'].to_numpy()
    MeVFlux_overlap = df_1fle_3fhl['MeVFlux_1FLE'].to_numpy()

    #-----------------------------------

    #HE Histogram
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
    plt.title('3FHL Source Fluxes with Shared 1FLE Sources')
    #plt.show()
    plt.savefig('{}/u1FLE_3FHL_HEFlux_Overlap.png'.format(direc))
    plt.clf()

    #-----------------------------------

    #Flux weight scatter
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
    plt.savefig('{}/u1FLE_3FHL_FluxWeight_Comparison.png'.format(direc))
    plt.clf()
    
    #-----------------------------------

    #Venn Diagram
    plt.figure(facecolor='w', figsize=(6,6))
    plt.title('3FGL Associated Sources')
    mpv.venn2([set(df_1fle['Name_3FGL']),set(df_3fhl['Name_3FGL'])], set_labels=('1FLE', '3FHL'), set_colors=('b','r'))
    plt.savefig('{}/u1FLE_3FHL_Venn.png'.format(direc))
    plt.clf()
    
    print('Done 3FHL Comparison...')
    
    
#=====2LAC Comparison=====
if do_2LAC:
    print('Doing 2LAC Comparison...')

    #Get 2FGL catalog that Thorsten used for 2LAC E Fluxes
    thorst = fits.open('./fermi_2year_ps_catalogue_official.fits')
    catList_thorst = thorst[1]
    data_thorst = catList_thorst.data
    name_thorst = np.array(data_thorst.field('Source_Name'))
    EFlux_thorst = np.array(data_thorst.field('Energy_Flux100')).byteswap().newbyteorder()

    #Make into DF and close fits
    df_2fgl = pd.DataFrame(name_thorst)
    df_2fgl.columns = ['Name_2FGL']
    df_2fgl['GeVFlux_2FGL'] = EFlux_thorst
    thorst.close()

    #Add 2FGL association to 3FGL DF
    df_3fgl['Name_2FGL'] = name2fgl_3fgl

    #Read 2LAC into Dataframe
    import csv
    df_2lac = pd.read_table('./2lac.txt', sep='|')

    #Rename Columns
    df_2lac.columns = ["id","Name_2FGL","RA","DEC","THETA","Class_2LAC", "Redshift", "SED", "GeVPhoFlux_2LAC", "Spectral_Ind"]

    #Get rid of whitespace so naming conventions match
    df_2lac['Name_2FGL'] = df_2lac['Name_2FGL'].apply(lambda x: x.replace(' ',''))
    df_2fgl['Name_2FGL'] = df_2fgl['Name_2FGL'].apply(lambda x: x.replace(' ',''))
    df_3fgl['Name_2FGL'] = df_3fgl['Name_2FGL'].apply(lambda x: x.replace(' ',''))

    #Keeping only the columns of current interest
    df_2lac = df_2lac[['Name_2FGL', 'Class_2LAC']]

    #First, merge 1fle (w/ 4fgl classes) and 3fgl on shared 3fgl name
    df_1fle_3fgl = pd.merge(df_3fgl, df_1fle_4fgl, on='Name_3FGL', how='inner')

    #Now merge 2LAC and 2FGL (thorst)
    df_2lac_2fgl = pd.merge(df_2lac, df_2fgl, on='Name_2FGL', how='inner')

    #Then, merge the two merged DFs
    df_1fle_3fgl_2lac_2fgl = pd.merge(df_1fle_3fgl, df_2lac_2fgl, on='Name_2FGL', how='inner')

    #Keep only columns of current interest
    df_1fle_3fgl_2lac_2fgl = df_1fle_3fgl_2lac_2fgl[['Name_3FGL', 'Name_2FGL', 'Name_1FLE', 'Class_4FGL', 'GeVFlux_2FGL', 'MeVFlux_1FLE']]

    #Get overlapped classes
    overlap_class_list = np.array(df_1fle_3fgl_2lac_2fgl['Class_4FGL'].apply(lambda s: s.lower().strip()))
    unique, counts = np.unique(overlap_class_list, return_counts=True)
    overlap_class_dict = dict(zip(unique, counts))
    print('The following are the number of each class that is included in the overlapped sources.')
    print(overlap_class_dict)

    #-----------------------------------

    #2LAC Flux Hist
    fig, ax1 = plt.subplots(facecolor='w', figsize=(8,6))

    logbins = np.logspace(np.log10(1.0e-12),np.log10(1.0e-8), 25)
    ax1.hist(df_2fgl['GeVFlux_2FGL'], bins=logbins, alpha=1, color='g', label='2LAC', histtype='step', linestyle='-')
    ax1.hist(df_1fle_3fgl_2lac_2fgl['GeVFlux_2FGL'], bins=logbins, alpha=1, color='g', label='Shared', histtype='step', linestyle='--')

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylabel('Counts')
    ax1.set_xlabel('Energy Flux (.1-100GeV) [erg/cm^2/s]', color='g')
    ax1.tick_params(axis='x', labelcolor='g')

    plt.legend()
    plt.savefig('{}/u1FLE_2LAC_HEFlux_Overlap.png'.format(direc))
    plt.clf()

    #-----------------------------------

    #Flux weight scatter
    fig, ax = plt.subplots(facecolor='w', figsize=(6,6))
    ax.scatter(df_1fle_3fgl_2lac_2fgl['MeVFlux_1FLE']/np.sum(flux30_100_1fle), df_1fle_3fgl_2lac_2fgl['GeVFlux_2FGL']/np.sum(df_2fgl['GeVFlux_2FGL']), marker='+', color='cadetblue')
    ax.set_xlabel('1FLE Normalized Energy Flux (30-100 MeV)', color='b')
    ax.set_ylabel('2LAC Normalized Energy Flux (.1-100 GeV)', color='g')
    #plt.xlim(-0.001,.2)
    #plt.ylim(-0.001,.2)

    lim = .00004, .04
    ax.plot(lim, lim, 'k--', zorder=-10)
    ax.set(xlim=lim, ylim=lim, aspect='equal')

    ax.set_yscale('log')
    ax.set_xscale('log')

    ax.set_title('Normalized Energy Fluxes of 1FLE/2LAC Shared Sources')
    plt.savefig('{}/u1FLE_2LAC_FluxWeight_Comparison.png'.format(direc))
    plt.clf()

    #-----------------------------------

    #Overlap Venn Diagram
    plt.figure(facecolor='w', figsize=(6,6))
    plt.title('2FGL Associated Sources')
    mpv.venn2([set(df_1fle_3fgl['Name_2FGL']),set(df_2lac['Name_2FGL'])], set_labels=('1FLE', '2LAC'), set_colors=('b','g'))
    plt.savefig('{}/u1FLE_2LAC_Venn.png'.format(direc))

    print('Done 2LAC Comparison...')