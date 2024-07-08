#!/usr/bin/env python
# coding: utf-8

# In[1]:


import zfit
from zfit import z
import pandas as pd
import root_pandas
from root_pandas import read_root
import ROOT
import matplotlib.pyplot as plt
import tensorflow as tf
#ZFIT_DISABLE_TF_WARNINGS=1
import numpy as np
import mplhep
from hepstats.splot import compute_sweights


# In[2]:


import ROOT as r
r.gROOT.LoadMacro('/group/belle2/users2022/dhettiar/belle2style/Belle2Style.C') 
r.SetBelle2Style()

# Make nice looking plots
plt.rcParams.update({
          'font.size': 20,
          'figure.figsize': (12, 8),
})


# In[3]:


lccol = ["Lambdac_M","reWeights","p_charge","K_p","pi_p","K_cosTheta","pi_cosTheta"]
dcol = ["M","K_charge","pi_2_p","K_p","CMS_cosTheta","K_cosTheta","pi_2_cosTheta","isSignal"]

lcinfile = '/group/belle2/users2022/dhettiar/ntuples/MC/lcp/weights/pkpi_pkk_weighted.root'
dinfile = '/group/belle2/users2022/dhettiar/ntuples/MC/d/k3pi_pkk_625.root'


# In[4]:


lcdf = root_pandas.read_root(lcinfile, key='lcp_pkpi', columns=lccol, where='K_p < 4.0 && pi_p <4.0')


# In[5]:


ddf = root_pandas.read_root(dinfile, key='D0tree', columns=dcol, where='K_p < 4.0 && pi_2_p < 4.0')


# In[6]:


# Standard plot settings
lupper = 2.34
llower = 2.24
dupper = 1.92
dlower = 1.80               


# In[7]:


# Standard cuts for lambda_c and D^0
lc_cut = str(llower)+' < Lambdac_M < '+str(lupper)
d_cut = str(dlower)+ ' < M < ' +str(dupper) 


# In[8]:


# Simple plotting function
def plotVar(mydf, var, cuts, nbins=100, myrange=(llower,lupper), mylabel="", log=False):

    if mylabel=="":
        mylabel=var

    ax = plt.subplot()

    # define a numpy array from the given column in the dataframe
    npdata = mydf.query(cuts)[var].to_numpy()

    # create histograms
    ydata, bin_edges = np.histogram(npdata, bins=nbins, range=myrange)
    bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])
    ax.errorbar(bin_centers, ydata, yerr = np.sqrt(ydata), fmt='k.', label="Data")
    
    # set plot features
    plt.ylim(0,None)
    if log==True:
        plt.ylim(0.1,None)
        plt.yscale("log")
    plt.xlim(myrange)
    plt.xlabel(mylabel)
    plt.legend(loc=0)
    plt.show()


# In[9]:


# Method to overlay fitted pdfs and data sample
def plot_model(model, mydata, nevents, nbins=100, myrange=(llower,lupper), mylabel="", plot_data=True):

    fig = plt.figure()
    fig.subplots_adjust(hspace=0.4)

    from matplotlib import gridspec
    gs = gridspec.GridSpec(2, 1, height_ratios=[5, 1]) 
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    lower, upper = mydata.data_range.limit1d

    x = tf.linspace(lower, upper, num=nbins)  # np.linspace also works
    for mod, frac in zip(model.pdfs, model.params.values()):
        y = mod.pdf(x) * nevents / nbins * mydata.data_range.area()
        y *= frac
        ax0.plot(x, y)
    data_plot = zfit.run(z.unstack_x(mydata))  # we could also use the `to_pandas` method
    y = model.pdf(x) * nevents / nbins * mydata.data_range.area()
    ax0.plot(x, y)

    counts, bin_edges = np.histogram(data_plot, bins=nbins, range=myrange)
    bin_centers = (bin_edges[:-1] + bin_edges[1:])/2.
    yerrs = counts**0.5
    ax0.errorbar(bin_centers, counts, yerr=yerrs, fmt='k.', linestyle='')
    ax0.set_xlabel(mylabel)
    ax0.set_xlim(myrange)
    ax0.set_ylim(0,None)

    X = np.ma.masked_equal(yerrs,0)
    ypull = (y-counts)/X
    ax1.errorbar(bin_centers, ypull, yerr=ypull*[0], fmt='k.', linestyle='')
    ax1.set_xlim(myrange)
    ax1.set_ylim(-5,5)


# In[10]:


plotVar(lcdf,"Lambdac_M",lc_cut,mylabel=r'M($\Lambda_c^+\rightarrow p K^- \pi^+$)')


# In[11]:


plotVar(ddf,"M",d_cut,myrange=(1.80,1.92),mylabel=r'M($D^0\rightarrow K^- \pi^+ \pi^- \pi^+$)')


# In[12]:


# Stack plot of d_M
np_d_M_sw = ddf.query(d_cut)['M'].to_numpy()
np_d_K_p_sw = ddf.query(d_cut)['K_p'].to_numpy()
np_d_pi_p_sw = ddf.query(d_cut)['pi_2_p'].to_numpy()
np_d_K_cosTheta_sw = ddf.query(d_cut)['K_cosTheta'].to_numpy()
np_d_pi_cosTheta_sw = ddf.query(d_cut)['pi_2_cosTheta'].to_numpy()


# In[13]:


# Stack plot of lc_M
np_lc_M_sw = lcdf.query(lc_cut)['Lambdac_M'].to_numpy()
np_lc_K_p_sw = lcdf.query(lc_cut)['K_p'].to_numpy()
np_lc_pi_p_sw = lcdf.query(lc_cut)['pi_p'].to_numpy()
np_lc_K_cosTheta_sw = lcdf.query(lc_cut)['K_cosTheta'].to_numpy()
np_lc_pi_cosTheta_sw = lcdf.query(lc_cut)['pi_cosTheta'].to_numpy()


# # ZFit

# In[14]:


# Define default parameters for zfit
lcrange = (llower,lupper)
lcobs = zfit.Space('Lambdac_M', lcrange)


# In[15]:


# Lambdac fit parameters
mu = zfit.Parameter("mu", 2.289, 2.27, 2.3)
s1 = zfit.param.Parameter("s1", 0.005, 0.0001, 0.01)
s2 = zfit.param.Parameter("s2", 0.002, 0.0001, 0.01)
fg1 = zfit.param.Parameter("fg1", 0.20, 0., 1.)
a1 = zfit.Parameter("a1", 0.01, -1e6, 1e6)
a2 = zfit.Parameter("a2", 0.01, -1e6, 1e6)

# define PDFs
gaus1 = zfit.pdf.Gauss(obs=lcobs, mu=mu, sigma=s1)
gaus2 = zfit.pdf.Gauss(obs=lcobs, mu=mu, sigma=s2)
gaus = zfit.pdf.SumPDF(pdfs=[gaus1,gaus2], fracs=[fg1])

poly = zfit.pdf.Chebyshev(obs=lcobs, coeffs=[a1, a2])
sig_yield = zfit.Parameter('sig_yield', 1000000, 0, 1e7, step_size=1)
bkg_yield = zfit.Parameter('bkg_yield', 4000000, 0, 1e8, step_size=1)

sig_ext = gaus.create_extended(sig_yield)
bkg_ext = poly.create_extended(bkg_yield)

pdf_ext = zfit.pdf.SumPDF(pdfs=[sig_ext,bkg_ext])


# In[16]:


# D ----> Kpi Define default parameters for zfit
drange = (dlower,dupper)
dobs = zfit.Space('M', drange)


# In[17]:


# D----> Kpi fit parameters
dmu = zfit.Parameter("dmu", 1.864, 1.85, 1.88)
ds1 = zfit.param.Parameter("ds1", 0.005, 0.0001, 0.02)
ds2 = zfit.param.Parameter("ds2", 0.003, 0.0001, 0.02)
dfg1 = zfit.param.Parameter("dfg1", 0.20, 0., 1.)
da1 = zfit.Parameter("da1", 0.01, -1e6, 1e6)
da2 = zfit.Parameter("da2", 0.01, -1e6, 1e6)

# define PDFs
dgaus1 = zfit.pdf.Gauss(obs=dobs, mu=dmu, sigma=ds1)
dgaus2 = zfit.pdf.Gauss(obs=dobs, mu=dmu, sigma=ds2)
dgaus = zfit.pdf.SumPDF(pdfs=[dgaus1,dgaus2], fracs=[dfg1])

dpoly = zfit.pdf.Chebyshev(obs=dobs, coeffs=[da1, da2])
dsig_yield = zfit.Parameter('dsig_yield', 3000000, 0, 1e8, step_size=1)
dbkg_yield = zfit.Parameter('dbkg_yield', 30000000, 0, 1e9, step_size=1)

dsig_ext = dgaus.create_extended(dsig_yield)
dbkg_ext = dpoly.create_extended(dbkg_yield)

dpdf_ext = zfit.pdf.SumPDF(pdfs=[dsig_ext,dbkg_ext])


# # Apply weights for $\Lambda_c$

# In[18]:


# find the momentum and cos(theta) bin for each event
pkpidf = lcdf.query(lc_cut)


# In[19]:


pkpidf


# In[20]:


# this is a little dangerous since we replace keys with indices that depend on the order 
def getWeights(df):
    np_arr = df.to_numpy()
    return np_arr[:,0], np_arr[:,1]


# In[21]:


pkpi_M, pkpi_weights = getWeights(pkpidf)


# ## $\Lambda_c$ Integrated fit

# In[22]:


# Fill an array with the data to be fit
data_np = lcdf.query(lc_cut).Lambdac_M.to_numpy()
data = zfit.data.Data.from_numpy(obs=lcobs, array=data_np, weights=pkpi_weights)
data.set_data_range(lcrange) 
data_size = data.n_events


# In[23]:


# Define loss function
nll_simultaneous = zfit.loss.ExtendedUnbinnedNLL(model=pdf_ext,data=data)
minimizer = zfit.minimize.Minuit()
result = minimizer.minimize(nll_simultaneous)  
print(result)


# In[24]:


# Define loss function
nll_simultaneous = zfit.loss.ExtendedUnbinnedNLL(model=pdf_ext,data=data)
minimizer = zfit.minimize.Minuit()
result = minimizer.minimize(nll_simultaneous)  
print(result)


# In[25]:


plot_model(model=pdf_ext, mydata=data, myrange=lcrange, nevents=int(data_size), mylabel=r'M($\Lambda_{c}$)  [GeV/$c^{2}$]')


# In[26]:


lcweights = compute_sweights(pdf_ext, data)


# ## $D^0$ Integrated fit

# In[27]:


# Fill an array with the data to be fit
data_np = ddf.query(d_cut).M.to_numpy()
data = zfit.data.Data.from_numpy(obs=dobs, array=data_np)
data.set_data_range(drange) 
data_size = data.n_events


# In[28]:


# Define loss function
nll_simultaneous = zfit.loss.ExtendedUnbinnedNLL(model=dpdf_ext,data=data)
minimizer = zfit.minimize.Minuit()
result = minimizer.minimize(nll_simultaneous)  
print(result)


# In[29]:


plot_model(model=dpdf_ext, mydata=data, myrange=(1.80,1.92), nevents=int(data_size), mylabel=r'M($D^{0}$)  [GeV/$c^{2}$]')


# In[30]:


dweights = compute_sweights(dpdf_ext, data)


# # sPlot

# In[31]:


nbins = 10
nrange = (0,4.0)
ncrange= (-1,1)


# In[32]:


dscale = lcweights[sig_yield].sum()/dweights[dsig_yield].sum()
print(dscale)


# In[33]:


plt.figure(figsize=(24, 8))

# 2D Plot for lambdac
plt.subplot(1, 2, 1)
plt.hist2d(np_lc_pi_p_sw, np_lc_K_p_sw, bins=nbins, range=[nrange,nrange], weights=lcweights[sig_yield], alpha=1, density=False)
plt.colorbar()  # Show color scale
plt.xlabel(r"$\pi$ momentum")
plt.ylabel(r"$K$ momentum")
plt.title(r"$\Lambda_c^+\rightarrow p K^- \pi^+$")

# 2D Plot for D
plt.subplot(1, 2, 2)
plt.hist2d(np_d_pi_p_sw, np_d_K_p_sw, bins=nbins, range=[nrange,nrange], weights=dweights[dsig_yield], alpha=1, density=False)
plt.colorbar()  # Show color scale
plt.xlabel(r"$\pi$ momentum")
plt.ylabel(r"$K$ momentum")
plt.title(r"$D^0\rightarrow K^- \pi^- \pi^+ \pi^+$")

plt.tight_layout()
plt.show()


# In[34]:


plt.figure(figsize=(24, 8))

# 2D Plot for lambdac
plt.subplot(1, 2, 1)
plt.hist2d(np_lc_pi_cosTheta_sw, np_lc_K_cosTheta_sw, bins=nbins, range=[ncrange,ncrange], weights=lcweights[sig_yield], alpha=1, density=False)
plt.colorbar()  # Show color scale
plt.xlabel(r"$\pi$ cosTheta")
plt.ylabel(r"$K$ cosTheta")
plt.title(r"$\Lambda_c^+\rightarrow p K^- \pi^+$")

# 2D Plot for D
plt.subplot(1, 2, 2)
plt.hist2d(np_d_pi_cosTheta_sw, np_d_K_cosTheta_sw, bins=nbins, range=[ncrange,ncrange], weights=dweights[dsig_yield], alpha=1, density=False)
plt.colorbar()  # Show color scale
plt.xlabel(r"$\pi$ cosTheta")
plt.ylabel(r"$K$ cosTheta")
plt.title(r"$D^0\rightarrow K^- \pi^- \pi^+ \pi^+$")

plt.tight_layout()
plt.show()


# # Reweighting

# In[35]:


dscale = lcweights[sig_yield].sum()/dweights[dsig_yield].sum()
print(dscale)


# In[36]:


# Weighting
lc_counts_per_bin_p = np.histogram2d(np_lc_pi_p_sw, np_lc_K_p_sw, bins=nbins, range=[nrange,nrange], weights=lcweights[sig_yield], density=False)
d_counts_per_bin_p = np.histogram2d(np_d_pi_p_sw, np_d_K_p_sw, bins=nbins, range=[nrange,nrange], weights=[dscale]*dweights[dsig_yield], density=False)

# Calculate the event ratio for each bin
lc_counts_per_bin_p = lc_counts_per_bin_p[0]  
d_counts_per_bin_p = d_counts_per_bin_p[0]
event_ratio_p = lc_counts_per_bin_p / d_counts_per_bin_p

# Calculate the bin edges
bin_edges_p = np.linspace(nrange[0], nrange[1], nbins + 1)

# Digitize for each dimension
pi_bins_p = np.digitize(np_d_pi_p_sw, bin_edges_p) - 1  
k_bins_p = np.digitize(np_d_K_p_sw, bin_edges_p) - 1 

ddf_pweighted = ddf.query(d_cut)

d0_pweights = []
for event in range(len(ddf_pweighted)):
    pi_bin_index_p = pi_bins_p[event]
    k_bin_index_p = k_bins_p[event]
    
    # Ensure the bin indices are within the valid range
    if pi_bin_index_p >= nbins or k_bin_index_p >= nbins: 
        weight = 1
    else:
        # Calculate weight from the 2D event_ratio for both pi_p and K_p
        pweight = event_ratio_p[pi_bin_index_p, k_bin_index_p]  # Correctly indexing the 2D event_ratio
    
    d0_pweights.append(pweight)


# In[37]:


d0_pweights_sw = d0_pweights*dweights[dsig_yield]


# In[38]:


# Weighting
lc_counts_per_bin_cosTheta = np.histogram2d(np_lc_pi_cosTheta_sw, np_lc_K_cosTheta_sw, bins=nbins, range=[ncrange,ncrange], weights=lcweights[sig_yield], density=False)
d_counts_per_bin_cosTheta = np.histogram2d(np_d_pi_cosTheta_sw, np_d_K_cosTheta_sw, bins=nbins, range=[ncrange,ncrange], weights=[dscale]*d0_pweights_sw, density=False)

# Calculate the event ratio for each bin
lc_counts_per_bin_cosTheta = lc_counts_per_bin_cosTheta[0]  
d_counts_per_bin_cosTheta = d_counts_per_bin_cosTheta[0]
event_ratio_cosTheta = lc_counts_per_bin_cosTheta / d_counts_per_bin_cosTheta

# Calculate the bin edges
bin_edges_cosTheta = np.linspace(ncrange[0], ncrange[1], nbins + 1)

# Digitize for each dimension
pi_bins_cosTheta = np.digitize(np_d_pi_cosTheta_sw, bin_edges_cosTheta) - 1  
k_bins_cosTheta = np.digitize(np_d_K_cosTheta_sw, bin_edges_cosTheta) - 1 

ddf_cweighted = ddf.query(d_cut)

d0_cweights = []
for event in range(len(ddf_cweighted)):
    pi_bin_index_cosTheta = pi_bins_cosTheta[event]
    k_bin_index_cosTheta = k_bins_cosTheta[event]
    
    # Ensure the bin indices are within the valid range
    if pi_bin_index_cosTheta >= nbins or k_bin_index_cosTheta >= nbins: 
        cweight = 1
    else:
        # Calculate weight from the 2D event_ratio for both pi_p and K_p
        cweight = event_ratio_cosTheta[pi_bin_index_cosTheta, k_bin_index_cosTheta]  # Correctly indexing the 2D event_ratio
    
    d0_cweights.append(cweight)


# In[39]:


reWeights = [a*b for a,b in zip(d0_pweights,d0_cweights)]


# In[40]:


ddf_weighted = ddf.query(d_cut)
ddf_weighted['reWeights'] = reWeights


# In[41]:


# Plot adjusted D0 plot
plt.figure(figsize=(24, 8))

# Plot for Lambdac for reference
plt.subplot(1, 2, 1)
plt.hist2d(np_lc_pi_p_sw, np_lc_K_p_sw, bins=nbins, range=[nrange, nrange], weights=lcweights[sig_yield], alpha=1, density=False)
plt.colorbar()  # Show color scale
plt.xlabel(r"$\pi$ momentum")
plt.ylabel(r"$K$ momentum")
plt.title(r"$\Lambda_c^+\rightarrow p K^- \pi^+$")

# Adjusted plot for D0 using reweighted events
plt.subplot(1, 2, 2)
plt.hist2d(np_d_pi_p_sw, np_d_K_p_sw, bins=nbins, range=[nrange, nrange], weights=reWeights*dweights[dsig_yield], alpha=1, density=False)  # Using pi_p_reweights, Not sure how to use weights seperately
plt.colorbar()  # Show color scale
plt.xlabel(r"$\pi$ momentum")
plt.ylabel(r"$K$ momentum")
plt.title(r"Weighted $D^0\rightarrow K^- \pi^- \pi^+ \pi^+$")

plt.tight_layout()
plt.show()


# In[42]:


# Plot adjusted D0 plot
plt.figure(figsize=(24, 8))

# Plot for Lambdac for reference
plt.subplot(1, 2, 1)
plt.hist2d(np_lc_pi_cosTheta_sw, np_lc_K_cosTheta_sw, bins=nbins, range=[ncrange, ncrange], weights=lcweights[sig_yield], alpha=1, density=False)
plt.colorbar()  # Show color scale
plt.xlabel(r"$\pi$ cosTheta")
plt.ylabel(r"$K$ cosTheta")
plt.title(r"$\Lambda_c^+\rightarrow p K^- \pi^+$")

# Adjusted plot for D0 using reweighted events
plt.subplot(1, 2, 2)
plt.hist2d(np_d_pi_cosTheta_sw, np_d_K_cosTheta_sw, bins=nbins, range=[ncrange, ncrange], weights=reWeights*dweights[dsig_yield], alpha=1, density=False)  # Using pi_p_reweights, Not sure how to use weights seperately
plt.colorbar()  # Show color scale
plt.xlabel(r"$\pi$ cosTheta")
plt.ylabel(r"$K$ cosTheta")
plt.title(r"Weighted $D^0\rightarrow K^- \pi^- \pi^+ \pi^+$")

plt.tight_layout()
plt.show()


# In[43]:


# Specify the directory to write the root file
path = "/group/belle2/users2022/dhettiar/ntuples/MC/d/weights/k3pi_pkk_weighted.root"


# In[44]:


root_pandas.to_root(ddf_weighted, path, key='D0tree', mode='w')


# # Plots for Belle2Note

# In[45]:


plt.figure(figsize=(40, 10))

# 2D Plot for lambdac
plt.subplot(1, 3, 1)
plt.hist2d(np_lc_pi_p_sw, np_lc_K_p_sw, bins=nbins, range=[nrange,nrange], weights=lcweights[sig_yield], alpha=1, density=False)
plt.colorbar()  # Show color scale
plt.xlabel(r"$\pi$ momentum")
plt.ylabel(r"$K$ momentum")
plt.title(r"$\Lambda_c^+\rightarrow p K^- \pi^+$")

# 2D Plot for D
plt.subplot(1, 3, 2)
plt.hist2d(np_d_pi_p_sw, np_d_K_p_sw, bins=nbins, range=[nrange,nrange], weights=dweights[dsig_yield], alpha=1, density=False)
plt.colorbar()  # Show color scale
plt.xlabel(r"$\pi$ momentum")
plt.ylabel(r"$K$ momentum")
plt.title(r"$D^0\rightarrow K^- \pi^- \pi^+ \pi^+$")

# Adjusted plot for D0 using reweighted events
plt.subplot(1, 3, 3)
plt.hist2d(np_d_pi_p_sw, np_d_K_p_sw, bins=nbins, range=[nrange, nrange], weights=reWeights*dweights[dsig_yield], alpha=1, density=False)  # Using pi_p_reweights, Not sure how to use weights seperately
plt.colorbar()  # Show color scale
plt.xlabel(r"$\pi$ momentum")
plt.ylabel(r"$K$ momentum")
plt.title(r"Weighted $D^0\rightarrow K^- \pi^- \pi^+ \pi^+$")

plt.tight_layout()
plt.show()


# In[46]:


plt.figure(figsize=(40, 10))

# 2D Plot for lambdac
plt.subplot(1, 3, 1)
plt.hist2d(np_lc_pi_cosTheta_sw, np_lc_K_cosTheta_sw, bins=nbins, range=[ncrange,ncrange], weights=lcweights[sig_yield], alpha=1, density=False)
plt.colorbar()  # Show color scale
plt.xlabel(r"$\pi$ cosTheta")
plt.ylabel(r"$K$ cosTheta")
plt.title(r"$\Lambda_c^+\rightarrow p K^- \pi^+$")

# 2D Plot for D
plt.subplot(1, 3, 2)
plt.hist2d(np_d_pi_cosTheta_sw, np_d_K_cosTheta_sw, bins=nbins, range=[ncrange,ncrange], weights=dweights[dsig_yield], alpha=1, density=False)
plt.colorbar()  # Show color scale
plt.xlabel(r"$\pi$ cosTheta")
plt.ylabel(r"$K$ cosTheta")
plt.title(r"$D^0\rightarrow K^- \pi^- \pi^+ \pi^+$")

# Adjusted plot for D0 using reweighted events
plt.subplot(1, 3, 3)
plt.hist2d(np_d_pi_cosTheta_sw, np_d_K_cosTheta_sw, bins=nbins, range=[ncrange, ncrange], weights=reWeights*dweights[dsig_yield], alpha=1, density=False)  # Using pi_p_reweights, Not sure how to use weights seperately
plt.colorbar()  # Show color scale
plt.xlabel(r"$\pi$ cosTheta")
plt.ylabel(r"$K$ cosTheta")
plt.title(r"Weighted $D^0\rightarrow K^- \pi^- \pi^+ \pi^+$")


plt.tight_layout()
plt.show()

