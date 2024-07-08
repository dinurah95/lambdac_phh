#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import root_pandas
from root_pandas import read_root
import ROOT
import matplotlib.pyplot as plt
import numpy as np
#from matplotlib.patches import HatchPatternFill
get_ipython().run_line_magic('matplotlib', 'inline')


# In[2]:


# Loading sample

lambda_c = "/group/belle2/users2022/dhettiar/ntuples/MC/d/k3pi_complete.root"


# In[3]:


# coloum defined in the root file to load

c_v_1 =['M']
c_v_2 =['isSignal','CMS_cosTheta','CMS_p','pi_1_p','pi_2_p','pi_3_p']
c_v_3 =['K_binaryID'] 


# In[4]:


Loading_cut = '1.80 < M & M < 1.92'


# In[5]:


columns= (c_v_1 + 
         c_v_2 +
         c_v_3 
         )


# In[6]:


# data frame defined

df_lam = root_pandas.read_root(lambda_c, key='D0tree', columns=columns, where=Loading_cut) 


# In[7]:


def mass_plot_all_Stack(var, cuts, nbins, myrange, myXlabel, mylog):
    
    figM = plt.figure(figsize=(15,8))
    
    # define a numpy array from the var column in the dataframe
    npbkg = df_lam.query(cuts+ ' & isSignal!=1')[var].to_numpy()
    nptrue = df_lam.query(cuts+' & isSignal==1')[var].to_numpy()
    
    # create a matplotlib histogram
    #plt.hist([nptrue,npbkg], bins=nbins, range=myrange,label=["Signal","Bkg"],histtype='step',log=False, stacked=True)
    plt.hist([npbkg, nptrue], 
             bins=nbins, range=myrange,
             label=["background","signal"], 
             histtype='stepfilled', 
             log=mylog, stacked=True)
    
    # set plot features
    plt.xlim(myrange)
    plt.xlabel(myXlabel, fontsize=16)
    epb = ((max(myrange) - min(myrange))*1000)/ nbins  #events per bins
    plt.ylabel(r'Events / (%.3f MeV/$c^2$)' % epb, fontsize=16)
    
 
    sig_len = len(nptrue)
    sig_mean = np.mean(nptrue)
    sig_std = np.std(nptrue)

    bkg_len = len(npbkg)
    bkg_mean = np.mean(npbkg)
    bkg_std = np.std(npbkg)
    
    print('Bkg_Entries:', bkg_len)
    print('Bkg_Mean:', '%.3f' %bkg_mean)
    print('Bkg_Std:', '%.3f' %bkg_std)
    
    print('--------')
    
    print('Sig_Entries:', sig_len)
    print('Sig_Mean:', r'%.3f' %sig_mean)
    print('Sig_Std:', r'%.3f' %sig_std)
    
    print('--------')
    
    plt.legend(fontsize=16)
    plt.show()


# In[8]:


def mass_plot_all_noStack(var, cuts, nbins, myrange, myXlabel, mylog):
    
    figM = plt.figure(figsize=(15,8))
    
    # define a numpy array from the var column in the dataframe
    npbkg = df_lam.query(cuts+ ' & isSignal!=1')[var].to_numpy()
    nptrue = df_lam.query(cuts+' & isSignal==1')[var].to_numpy()
    
    # create a matplotlib histogram
    plt.hist(npbkg, bins=nbins, range=myrange,label=["background",],histtype='stepfilled', log=mylog,  alpha=0.8)
    plt.hist(nptrue, bins=nbins, range=myrange,label=["signal"],histtype='stepfilled',log=mylog, alpha=0.5 )
    
    # set plot features
    plt.xlim(myrange)
    plt.xlabel(myXlabel, fontsize=16)
    epb = ((max(myrange) - min(myrange))*1000)/ nbins  #events per bins
    plt.ylabel(r'Events / (%.3f MeV/$c^2$)' % epb, fontsize=16)
   
    
    sig_len = len(nptrue)
    sig_mean = np.mean(nptrue)
    sig_std = np.std(nptrue)

    bkg_len = len(npbkg)
    bkg_mean = np.mean(npbkg)
    bkg_std = np.std(npbkg)
    
    print('Bkg_Entries:', bkg_len)
    print('Bkg_Mean:', '%.3f' %bkg_mean)
    print('Bkg_Std:', '%.3f' %bkg_std)
    
    print('--------')
    
    print('Sig_Entries:', sig_len)
    print('Sig_Mean:', r'%.3f' %sig_mean)
    print('Sig_Std:', r'%.3f' %sig_std)
    
    print('--------')
   

    plt.legend(fontsize=14)
    plt.show()


# In[9]:


def Individual_signal_mass_plot(var, cuts, nbins, myrange, myXlabel, mylog):
    
    figM = plt.figure(figsize=(15,8))
    
    # define a numpy array from the var column in the dataframe
    nptrue = df_lam.query(cuts+' & isSignal==1')[var].to_numpy()
    
    # create a matplotlib histogram
    #plt.hist([nptrue,npbkg], bins=nbins, range=myrange,label=["Signal","Bkg"],histtype='step',log=False, stacked=True)
    plt.hist(nptrue, bins=nbins, range=myrange,label=["Signal"],histtype='stepfilled',log=mylog, color='darkorange')
    
    # set plot features
    plt.xlim(myrange)
    plt.xlabel(myXlabel, fontsize=16)
    epb = ((max(myrange) - min(myrange))*1000)/ nbins  #events per bins
    plt.ylabel(r'Events / (%.3f MeV/$c^2$)' % epb, fontsize=16)
    
 
    sig_len = len(nptrue)
    sig_mean = np.mean(nptrue)
    sig_std = np.std(nptrue)

    
    print('Sig_Entries:', sig_len)
    print('Sig_Mean:', r'%.3f' %sig_mean)
    print('Sig_Std:', r'%.3f' %sig_std)
    
    print('--------')
    
    plt.legend(fontsize=14)
    plt.show()


# In[10]:


def Individual_bkg_mass_plot(var, cuts, nbins, myrange, myXlabel, mylog):
    
    figM = plt.figure(figsize=(15,8))
    
    # define a numpy array from the var column in the dataframe
    npbkg = df_lam.query(cuts+ ' & Lambdac_isSignal!=1')[var].to_numpy()
        
    # create a matplotlib histogram
    plt.hist(npbkg, bins=nbins, range=myrange,label=["Signal!=1",],histtype='stepfilled', log=mylog, color='royalblue')
    
    
    # set plot features
    plt.xlim(myrange)
    plt.xlabel(myXlabel, fontsize=16)
    epb = ((max(myrange) - min(myrange))*1000)/ nbins  #events per bins
    plt.ylabel(r'Events / (%.3f MeV/$c^2$)' % epb, fontsize=16)
   
    
    bkg_len = len(npbkg)
    bkg_mean = np.mean(npbkg)
    bkg_std = np.std(npbkg)
    
    print('Bkg_Entries:', bkg_len)
    print('Bkg_Mean:', '%.3f' %bkg_mean)
    print('Bkg_Std:', '%.3f' %bkg_std)
    
    print('--------')
    
 

    plt.legend(fontsize=14)
    plt.show()


# In[11]:


def get_purity_and_sigeff(var, cuts):
    # Create numpy arrays for event types
    nptrue = df_lam.query(cuts + '& 1.846 < M < 1.882'+'  & isSignal==1')[var].to_numpy()
    npbkg = df_lam.query(cuts + ' & 1.846 < M < 1.882'+' & isSignal!=1')[var].to_numpy()

    sig_events = len(nptrue)
    bkg_events = len(npbkg)
    tot_events = sig_events + bkg_events

    # Calculate purity
    purity = sig_events / tot_events * 100

    # Create numpy arrays for signal before/after cut
    nptrue_before = df_lam.query('isSignal == 1')[var].to_numpy()
    sig_before = len(nptrue_before)
    nptrue_after = df_lam.query(cuts + ' & isSignal == 1')[var].to_numpy()
    sig_after = len(nptrue_after)

    # Calculate signal efficiency
    sigeff = sig_after / sig_before * 100

    return purity, sigeff


# In[12]:


mass_plot_all_Stack(var='M',
           cuts='1.80 < M < 1.92',           
           nbins=100, myrange=(1.80,1.92),
           myXlabel=r'$D^0$ Mass ($D^0\to K^- \pi^+ \pi^+ \pi^- $) [GeV/$c^2$]', mylog=False)


# In[13]:


get_purity_and_sigeff('M', '1.80 < M < 1.92')


# In[14]:


mass_plot_all_Stack(var='M',
           cuts='1.80 < M < 1.92 \
                   & K_binaryID > 0.7',           
           nbins=100, myrange=(1.80,1.92),
           myXlabel=r'$D^0$ Mass ($D^0\to K^- \pi^+ \pi^+ \pi^- $) [GeV/$c^2$]', mylog=False)


# In[15]:


get_purity_and_sigeff('M', '1.80 < M < 1.92                   & K_binaryID > 0.7')


# In[16]:


mass_plot_all_Stack(var='M',
           cuts='1.80 < M < 1.92 \
                   & K_binaryID > 0.7 \
                   & pi_1_p > 0.30 \
                   & pi_2_p > 0.30 \
                   & pi_3_p > 0.30',           
           nbins=100, myrange=(1.80,1.92),
           myXlabel=r'$D^0$ Mass ($D^0\to K^- \pi^+ \pi^+ \pi^- $) [GeV/$c^2$]', mylog=False)


# In[17]:


get_purity_and_sigeff('M', '1.80 < M < 1.92                   & K_binaryID > 0.7                    & CMS_p > 2.5                    & pi_1_p > 0.30                    & pi_2_p > 0.30                    & pi_3_p > 0.30')


# In[18]:


mass_plot_all_noStack(var='K_binaryID',
           cuts='1.80 < M < 1.92',           
           nbins=100, myrange=(0.2,1.0),
           myXlabel=r'K binaryID ($D^0\to K^- \pi^+ \pi^+ \pi^- $)', mylog=True)


# In[ ]:




