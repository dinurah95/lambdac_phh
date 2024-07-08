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


# Loading all lambda_c sample
lambda_c = "/group/belle2/users2022/dhettiar/ntuples/MC/lcp/ppipi_dis.root"


# In[3]:


# coloum defined in the root file to load

c_v_1 =['Lambdac_isSignal','Lambdac_M']
c_v_2 =['Lambdac_p_cms','p_trinaryID','Lambdac_significanceOfDistance','Lambdac_flightDistance']
c_v_3 =['p_p', 'pi1_p','pi2_p'] 


# In[4]:


Loading_cut = ('2.24 < Lambdac_M & Lambdac_M < 2.34')


# In[5]:


columns= (c_v_1 + 
         c_v_2+  
         c_v_3   
         )


# In[6]:


df_lam = read_root(lambda_c, key='lcp_ppipi', columns=columns, where=Loading_cut) 


# In[7]:


def plotFom(var, _min, _max, nbins, _greaterthan, label, cuts):
    
    #define cuts that will always be applied, mostly to reduce processing time
    
    basecuts = cuts 
 
    #use the size of the variable's numpy arrays to define a quantity that represents the total number 
    #of signal events and background events, respectively
    
    bkgtot = df_lam.query(basecuts+' & Lambdac_isSignal!=1')[var].to_numpy().size
    sigtot = df_lam.query(basecuts+' & Lambdac_isSignal==1')[var].to_numpy().size
    
    #define empty arrays which will be appended in the for loop
    
    testcuts = []
    globalsig = []
    globalbkg = []
    fom = []
    
    interval = (_max - _min)/nbins #ddefine an interval based on the number of bins and the range of our variable
    
    
    
    for bin in range(0, nbins):
        testcut = interval*bin+(_min) #define a test cut to apply for each bin
        testcuts.append(testcut) #append that test cut to our array
        globalcuts = ""
        if(_greaterthan == True): #if the figure of merit is being used to find the best ">" cut,
            globalcuts = basecuts + " & " + var + " > " + str(testcut) #define our overall cuts as base cuts
            #plus the variable greater than the test cut (e.g apply test cut for each bin)
        else:
            globalcuts = basecuts + " & " + var + " < " + str(testcut) #similarly for < cuts
            
        globalsig.append(df_lam.query(globalcuts+' & Lambdac_isSignal==1')[var].to_numpy().size)
        globalbkg.append(df_lam.query(globalcuts+' & Lambdac_isSignal!=1')[var].to_numpy().size) #check signal/bg
        #events after the test cut is applied
        
        #append the figure of merit array with the figure of merit calculation for this bin, which is
        #(number of sig events)/Sqrt((number of sig events + number of bkg events))
        fom.append(globalsig[bin]/np.sqrt(globalsig[bin]+globalbkg[bin]))
     
    
    fig, ax = plt.subplots()
    # Twin the x-axis twice to make independent y-axes.
    axes = [ax, ax.twinx(), ax.twinx()]
    # Make some space on the right side for the extra y-axis.
    fig.subplots_adjust(right=0.75)
    # Move the last y-axis spine over to the right by 20% of the width of the axes
    axes[-1].spines['right'].set_position(('axes', 1.2))
    # To make the border of the right-most axis visible, we need to turn the frame
    # on. This hides the other plots, however, so we need to turn its fill off.
    axes[-1].set_frame_on(True)
    axes[-1].patch.set_visible(False)
    
    sigeff = []
    purity = []
    for bin in range(0, (nbins - 1)):
        sigeff.append(globalsig[bin]/sigtot) #calculate the efficiency of each cut (how much signal did we lose?)
        purity.append(globalsig[bin]/(globalbkg[bin]+globalsig[bin])) #how much of our sample is signal?
    sigeff.append(sigeff[nbins - 2])
    purity.append(purity[nbins - 2])
    # And finally we get to plot things...
    axes[0].plot(testcuts, fom, color='Red')
    axes[0].set_ylabel('Figure of merit', color='Red')
    axes[1].plot(testcuts, sigeff, color='Blue')
    axes[1].set_ylabel('Signal efficiency', color='Blue')
    
    axes[2].plot(testcuts, purity, color='Green')
    axes[2].set_ylabel('Purity', color='Green')
    axes[0].set_xlabel(label)
    ax.grid()
    # only one line may be specified; full height
    #plt.axvline(x=-0.7, color='black', lw=1, ls='--')
    plt.show()


# In[8]:


plotFom('Lambdac_p_cms', 
        2.0, 
        5.0, 
        100, 
        True, 
        r'$\Lambda_c$ momentum_cms ($\Lambda_c^+ \to p^+ \pi^+ \pi^-$)', 
        cuts ='2.240 < Lambdac_M < 2.340')


# In[9]:


plotFom('p_trinaryID', 
        0.2, 
        1.0, 
        100, 
        True, 
        r'proton trinaryID ($\Lambda_c^+ \to p^+ \pi^+ \pi^-$)', 
        cuts ='2.240 < Lambdac_M < 2.340\
                & 2.5 < Lambdac_p_cms')


# In[10]:


plotFom('Lambdac_flightDistance', 
        -0.03, 
        0.03, 
        100, 
        True, 
        r'$\Lambda_c$ flightDistance ($\Lambda_c^+ \to p^+ \pi^+ \pi^-$)', 
        cuts ='2.240 < Lambdac_M < 2.340\
                & 0.9 < p_trinaryID\
                & 2.5 < Lambdac_p_cms')


# In[14]:


plotFom('Lambdac_significanceOfDistance', 
        0, 
        1.0, 
        100, 
        True, 
        r'$\Lambda_c$ significanceOfDistance ($\Lambda_c^+ \to p^+ \pi^+ \pi^-$)', 
        cuts ='2.240 < Lambdac_M < 2.340\
                & 0.9 < p_trinaryID\
                & 2.5 < Lambdac_p_cms\
                & 0 < Lambdac_flightDistance')


# In[12]:


plotFom('pi1_p', 
        0, 
        0.5, 
        100, 
        True, 
        r'$\pi$ momentum ($\Lambda_c^+ \to p^+ \pi^+ \pi^-$)', 
        cuts ='2.240 < Lambdac_M < 2.340\
                & 0.9 < p_trinaryID\
                & 2.5 < Lambdac_p_cms\
                & 0 < Lambdac_flightDistance\
                & 0.25 < Lambdac_significanceOfDistance')


# In[13]:


plotFom('p_p', 
        0, 
        1.5, 
        100, 
        True, 
        r'proton momentum ($\Lambda_c^+ \to p^+ \pi^+ \pi^-$)', 
        cuts ='2.240 < Lambdac_M < 2.340\
                & 0.9 < p_trinaryID\
                & 2.5 < Lambdac_p_cms\
                & 0 < Lambdac_flightDistance\
                & 0.30 < pi1_p\
                & 0.30 < pi2_p')


# In[ ]:




