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

lambda_c = "/group/belle2/users2022/dhettiar/ntuples/MC/lcp/pkk_dis.root"


# In[3]:


# coloum defined in the root file to load

c_v_1 =['Lambdac_isSignal','Lambdac_M']
c_v_2 =['p_trinaryID',"Lambdac_p_cms",'Lambdac_flightDistance','p_p','K1_p','K2_p']
c_v_3 =['K1_binaryID','K2_binaryID']


# In[4]:


Loading_cut = '2.24 < Lambdac_M & Lambdac_M < 2.34'


# In[5]:


columns= (c_v_1 + 
         c_v_2+  
         c_v_3  
         )


# In[6]:


# data frame defined

df_lam = root_pandas.read_root(lambda_c, key='lcp_pkk', columns=columns, where=Loading_cut) 


# In[7]:


def mass_plot_all_Stack(var, cuts, nbins, myrange, myXlabel, mylog):
    
    figM = plt.figure(figsize=(15,8))
    
    # define a numpy array from the var column in the dataframe
    npbkg = df_lam.query(cuts+ ' & Lambdac_isSignal!=1')[var].to_numpy()
    nptrue = df_lam.query(cuts+' & Lambdac_isSignal==1')[var].to_numpy()
    
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
    npbkg = df_lam.query(cuts+ ' & Lambdac_isSignal!=1')[var].to_numpy()
    nptrue = df_lam.query(cuts+' & Lambdac_isSignal==1')[var].to_numpy()
    
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
   

    plt.legend(fontsize=16)
    plt.show()


# In[9]:


def Individual_signal_mass_plot(var, cuts, nbins, myrange, myXlabel, mylog):
    
    figM = plt.figure(figsize=(15,8))
    
    # define a numpy array from the var column in the dataframe
    nptrue = df_lam.query(cuts+' & Lambdac_isSignal==1')[var].to_numpy()
    
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


def mass_plot(var, cuts, nbins, myrange, myXlabel, mylog):
    
    figM = plt.figure(figsize=(15,8))
    
    # define a numpy array from the var column in the dataframe
    nptotal = df_lam.query(cuts)[var].to_numpy()
    
    # create a matplotlib histogram
    #plt.hist([nptrue,npbkg], bins=nbins, range=myrange,label=["Signal","Bkg"],histtype='step',log=False, stacked=True)
    plt.hist(nptotal, bins=nbins, range=myrange,label='Generic',histtype='stepfilled',log=mylog)
    
    # set plot features
    plt.xlim(myrange)
    plt.xlabel(myXlabel, fontsize=16)
    epb = ((max(myrange) - min(myrange))*1000)/ nbins  #events per bins
    plt.ylabel(r'Events / (%.3f MeV/$c^2$)' % epb, fontsize=16)
    
 
    sig_len = len(nptotal)
    sig_mean = np.mean(nptotal)
    sig_std = np.std(nptotal)

    
    print('Entries:', sig_len)
    print('Mean:', r'%.3f' %sig_mean)
    print('Std:', r'%.3f' %sig_std)
    
    print('--------')
    
    plt.legend(fontsize=14)
    plt.show()


# In[12]:


def Abs_signal_mass_plot(var, cuts, nbins, myrange, myXlabel, mylog):
    
    figM = plt.figure(figsize=(15,8))
    
    # define a numpy array from the var column in the dataframe
    nptrue = df_lam.query(cuts+' & Lambdac_isSignal==1')
    
    # Create separate arrays for positive and negative values
    nptrue_positive = nptrue[nptrue[var] > 0][var].to_numpy()
    nptrue_negative = np.abs(nptrue[nptrue[var] < 0][var].to_numpy())
    
    
    
    # create a matplotlib histogram
    #plt.hist([nptrue,npbkg], bins=nbins, range=myrange,label=["Signal","Bkg"],histtype='step',log=False, stacked=True)
    plt.hist(nptrue_positive, bins=nbins, range=myrange,label=["Signal_Positive"],histtype='step',log=mylog, color='darkorange', linewidth=2.0)
    plt.hist(nptrue_negative, bins=nbins, range=myrange,label=["Signal_Negative"],histtype='step',log=mylog, color='steelblue', linewidth=2.0)
    
    
    
    # set plot features
    plt.xlim(myrange)
    plt.xlabel(myXlabel, fontsize=16)
    epb = ((max(myrange) - min(myrange))*1000)/ nbins  #events per bins
    plt.ylabel(r'Events / (%.3f MeV/$c^2$)' % epb, fontsize=16)
    
 
    pos_sig_len = len(nptrue_positive)
    pos_sig_mean = np.mean(nptrue_positive)
    pos_sig_std = np.std(nptrue_positive)
    
    neg_sig_len = len(nptrue_negative)
    neg_sig_mean = np.mean(nptrue_negative)
    neg_sig_std = np.std(nptrue_negative)

    
    print('Pos_Sig_Entries:', pos_sig_len)
    print('Pos_Sig_Mean:', r'%.3f' %pos_sig_mean)
    print('Pos_Sig_Std:', r'%.3f' %pos_sig_std)
    
    print('--------')
    
    print('Neg_Sig_Entries:', neg_sig_len)
    print('Neg_Sig_Mean:', r'%.3f' %neg_sig_mean)
    print('Neg_Sig_Std:', r'%.3f' %neg_sig_std)
    
    plt.legend(fontsize=14)
    plt.show()


# In[13]:


def Abs_bkg_mass_plot(var, cuts, nbins, myrange, myXlabel, mylog):
    
    figM = plt.figure(figsize=(15,8))
    
    # define a numpy array from the var column in the dataframe
    npbkg = df_lam.query(cuts+ ' & Lambdac_isSignal!=1')
    
    # Create separate arrays for positive and negative values
    npbkg_positive = npbkg[npbkg[var] > 0][var].to_numpy()
    npbkg_negative = np.abs(npbkg[npbkg[var] < 0][var].to_numpy())
        
    # create a matplotlib histogram
    plt.hist(npbkg_positive, bins=nbins, range=myrange,label=["Bkg_Signal!=1",],histtype='step', log=mylog, color='darkorange', linewidth=2.0)
    plt.hist(npbkg_negative, bins=nbins, range=myrange,label=["Bkg_Signal!=1",],histtype='step', log=mylog, color='steelblue', linewidth=2.0)
    
    # set plot features
    plt.xlim(myrange)
    plt.xlabel(myXlabel, fontsize=16)
    epb = ((max(myrange) - min(myrange))*1000)/ nbins  #events per bins
    plt.ylabel(r'Events / (%.3f MeV/$c^2$)' % epb, fontsize=16)
   
    
    pos_bkg_len = len(npbkg_positive)
    pos_bkg_mean = np.mean(npbkg_positive)
    pos_bkg_std = np.std(npbkg_positive)
    
    neg_bkg_len = len(npbkg_negative)
    neg_bkg_mean = np.mean(npbkg_negative)
    neg_bkg_std = np.std(npbkg_negative)

    
    print('Pos_Bkg_Entries:', pos_bkg_len)
    print('Pos_Bkg_Mean:', r'%.3f' %pos_bkg_mean)
    print('Pos_Bkg_Std:', r'%.3f' %pos_bkg_std)
    
    print('--------')
    
    print('Neg_Bkg_Entries:', neg_bkg_len)
    print('Neg_Bkg_Mean:', r'%.3f' %neg_bkg_mean)
    print('Neg_Bkg_Std:', r'%.3f' %neg_bkg_std)
    
 

    plt.legend(fontsize=14)
    plt.show()


# In[14]:


def get_purity_and_sigeff(var, cuts):
    # Create numpy arrays for event types
    nptrue = df_lam.query(cuts + '& 2.278 < Lambdac_M < 2.294'+'  & Lambdac_isSignal==1')[var].to_numpy()
    npbkg = df_lam.query(cuts + ' &2.278 < Lambdac_M < 2.294'+' & Lambdac_isSignal!=1')[var].to_numpy()

    sig_events = len(nptrue)
    bkg_events = len(npbkg)
    tot_events = sig_events + bkg_events

    # Calculate purity
    purity = sig_events / tot_events * 100

    # Create numpy arrays for signal before/after cut
    nptrue_before = df_lam.query('Lambdac_isSignal == 1')[var].to_numpy()
    sig_before = len(nptrue_before)
    nptrue_after = df_lam.query(cuts + ' & Lambdac_isSignal == 1')[var].to_numpy()
    sig_after = len(nptrue_after)

    # Calculate signal efficiency
    sigeff = sig_after / sig_before * 100

    return purity, sigeff


# In[15]:


mass_plot_all_Stack(var='Lambdac_M',
           cuts='2.240 < Lambdac_M < 2.340',           
           nbins=100, myrange=(2.24,2.34),
           myXlabel=r'$\Lambda_c$ Mass ($\Lambda_c^+\to p^+ K^- K^+$) [GeV/$c^2$]', mylog=False)


# In[16]:


get_purity_and_sigeff('Lambdac_M', '2.240 < Lambdac_M < 2.340')


# In[17]:


mass_plot_all_Stack(var='Lambdac_M',
           cuts='2.240 < Lambdac_M < 2.340 \
                    & 0.7 < K1_binaryID \
                    & 0.7 < K2_binaryID\
                    & 0 < Lambdac_flightDistance  \
                    & 0.9 < p_trinaryID \
                    & 2.5 < Lambdac_p_cms\
                    & 4.0 > K1_p\
                    & 4.0 > K2_p\
                    & 5.0 > p_p',           
           nbins=100, myrange=(2.24,2.34),
           myXlabel=r'$\Lambda_c$ Mass ($\Lambda_c^+\to p^+ K^- K^+$) [GeV/$c^2$]', mylog=False)


# In[18]:


get_purity_and_sigeff('Lambdac_M', '2.240 < Lambdac_M < 2.340                    & 0.7 < K1_binaryID                     & 0.7 < K2_binaryID                    & 0 < Lambdac_flightDistance                      & 0.9 < p_trinaryID                     & 2.5 < Lambdac_p_cms                    & 4.0 > K1_p                    & 4.0 > K2_p                    & 5.0 > p_p')


# In[19]:


mass_plot_all_noStack(var='Lambdac_p_cms',
           cuts='2.240 < Lambdac_M < 2.340',           
           nbins=100, myrange=(2,5),
           myXlabel=r'$\Lambda_c$ momentum_cms ($\Lambda_c^+\to p^+ K^- K^+$) [GeV/c]', mylog=True)


# In[20]:


get_purity_and_sigeff('Lambdac_M', '2.240 < Lambdac_M < 2.340                    & 2.5 < Lambdac_p_cms')


# In[21]:


mass_plot_all_noStack(var='p_trinaryID',
           cuts='2.240 < Lambdac_M < 2.340',           
           nbins=100, myrange=(0.2,1),
           myXlabel=r'proton trinaryID ($\Lambda_c^+\to p^+ K^- K^+$)', mylog=True)


# In[22]:


get_purity_and_sigeff('Lambdac_M', '2.240 < Lambdac_M < 2.340                    & 2.5 < Lambdac_p_cms                    & 0.9 < p_trinaryID')


# In[23]:


mass_plot_all_noStack(var='K1_binaryID',
           cuts='2.240 < Lambdac_M < 2.340',           
           nbins=100, myrange=(0,1),
           myXlabel=r'kaon binaryID ($\Lambda_c^+\to p^+ K^- K^+$)', mylog=True)


# In[24]:


get_purity_and_sigeff('Lambdac_M', '2.240 < Lambdac_M < 2.340                    & 2.5 < Lambdac_p_cms                    & 0.9 < p_trinaryID                    & 0.7 < K1_binaryID                    & 0.7 < K2_binaryID')


# In[25]:


mass_plot_all_noStack(var='Lambdac_flightDistance',
           cuts='2.240 < Lambdac_M < 2.340',           
           nbins=100, myrange=(-0.2,0.2),
           myXlabel=r'$\Lambda_c$ flightDistance ($\Lambda_c^+\to p^+ K^- K^+$) [cm]', mylog=True)


# In[26]:


get_purity_and_sigeff('Lambdac_M', '2.240 < Lambdac_M < 2.340                    & 0.7 < K1_binaryID                    & 0.7 < K2_binaryID                    & 0.9 < p_trinaryID                    & 2.5 < Lambdac_p_cms                    & 0 < Lambdac_flightDistance')


# In[27]:


(/////////////)

