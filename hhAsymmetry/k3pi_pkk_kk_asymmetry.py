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


dcol = ["M","K_charge","CMS_cosTheta","isSignal","pi_1_p","pi_3_p","pi_1_cosTheta","pi_3_cosTheta"]
dinfile = '/group/belle2/users2022/dhettiar/ntuples/MC/d/k3pi_pkk_625.root'


# In[4]:


ddf = root_pandas.read_root(dinfile, key='D0tree', columns=dcol, where='pi_1_p < 4.0 && pi_3_p <4.0')


# In[5]:


# Standard plot settings
dupper = 1.92
dlower = 1.80             


# In[6]:


# Standard cuts for lambda_c and D^0
d_cut = str(dlower)+ ' < M < ' +str(dupper)


# In[7]:


# Simple plotting function
def plotVar(mydf, var, cuts, nbins=100, myrange=(dlower,dupper), mylabel="", log=False):

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


# In[8]:


# Method to overlay fitted pdfs and data sample
def plot_model(model, mydata, nevents, nbins=100, myrange=(dlower,dupper), mylabel="", plot_data=True):

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


# In[9]:


plotVar(ddf,"M",d_cut,mylabel=r'M($D^0\rightarrow K^- \pi^+ \pi^- \pi^+$)')


# In[10]:


# Stack plot of lc_M
np_d_M_sw = ddf.query(d_cut)['M'].to_numpy()
np_d_pi_1_p_sw = ddf.query(d_cut)['pi_1_p'].to_numpy()
np_d_pi_3_p_sw = ddf.query(d_cut)['pi_3_p'].to_numpy()
np_d_pi_1_cosTheta_sw = ddf.query(d_cut)['pi_1_cosTheta'].to_numpy()
np_d_pi_3_cosTheta_sw = ddf.query(d_cut)['pi_3_cosTheta'].to_numpy()


# # ZFit

# In[11]:


# D ----> Kpi Define default parameters for zfit
drange = (dlower,dupper)
dobs = zfit.Space('M', drange)


# In[12]:


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
dsig_yield = zfit.Parameter('dsig_yield', 3000000, 0, 1e7, step_size=1)
dbkg_yield = zfit.Parameter('dbkg_yield', 30000000, 0, 1e8, step_size=1)

dsig_ext = dgaus.create_extended(dsig_yield)
dbkg_ext = dpoly.create_extended(dbkg_yield)

dpdf_ext = zfit.pdf.SumPDF(pdfs=[dsig_ext,dbkg_ext])


# ## $D^0$ Integrated fit

# In[13]:


# Fill an array with the data to be fit
data_np = ddf.query(d_cut).M.to_numpy()
data = zfit.data.Data.from_numpy(obs=dobs, array=data_np)
data.set_data_range(drange) 
data_size = data.n_events


# In[14]:


# Define loss function
nll_simultaneous = zfit.loss.ExtendedUnbinnedNLL(model=dpdf_ext,data=data)
minimizer = zfit.minimize.Minuit()
result = minimizer.minimize(nll_simultaneous)  
print(result)


# In[15]:


plot_model(model=dpdf_ext, mydata=data, myrange=(1.80,1.92), nevents=int(data_size), mylabel=r'M($D^{0}$)  [GeV/$c^{2}$]')


# In[16]:


dweights = compute_sweights(dpdf_ext, data)


# # sPlot

# In[17]:


nbins = 25
nrange = (0,4.0)
cosTheta_range= (-1,1)


# In[18]:


fig, axs = plt.subplots(1, 2, figsize=(30, 10))
#nbins = 10

axs[0].hist(np_d_pi_3_p_sw, bins=nbins, range=nrange, weights=dweights[dsig_yield], label=r"$D^0 \ (\pi^+ \ momentum)$", alpha=.5, density=False)
axs[0].hist(np_d_pi_1_p_sw, bins=nbins, range=nrange, weights=dweights[dsig_yield], label=r"$D^0 \ (\pi^- \ momentum)$", alpha=.5, density=False)
axs[0].set_xlabel(r"$D^0 \ (\pi \ momentum)$")
axs[0].legend(fontsize=20)

axs[1].hist(np_d_pi_3_cosTheta_sw, bins=nbins, range=cosTheta_range, weights=dweights[dsig_yield], label=r"$D^0 \ (\pi^+ \ cosTheta)$", alpha=.5, density=False)
axs[1].hist(np_d_pi_1_cosTheta_sw, bins=nbins, range=cosTheta_range, weights=dweights[dsig_yield], label=r"$D^0 \ (\pi^- \ cosTheta)$", alpha=.5, density=False)
axs[1].set_xlabel(r"$D^0 \ (\pi \ cosTheta)$")
axs[1].legend(fontsize=20)


# In[19]:


# Compute ratio for each bin with uncertainity
hist_d_pi_3_p, bin_edges = np.histogram(np_d_pi_3_p_sw, bins=nbins, range=nrange, weights=dweights[dsig_yield], density=False)
hist_d_pi_1_p, _ = np.histogram(np_d_pi_1_p_sw, bins=nbins, range=nrange, weights=dweights[dsig_yield], density=False)

pi_3_pi_1_p_ratio = (hist_d_pi_3_p/hist_d_pi_1_p)

hist_d_pi_3_p_err = np.sqrt(hist_d_pi_3_p)
hist_d_pi_1_p_err = np.sqrt(hist_d_pi_1_p)

pi_3_pi_1_p_ratio_err = pi_3_pi_1_p_ratio * np.sqrt((hist_d_pi_3_p_err / hist_d_pi_3_p)**2 + (hist_d_pi_1_p_err / hist_d_pi_1_p)**2)

bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

# Plotting with uncertainties
plt.errorbar(bin_centers, pi_3_pi_1_p_ratio, yerr=pi_3_pi_1_p_ratio_err, fmt='o')
plt.xlabel(r"$\pi \ momentum (D^0 \rightarrow K^- \pi^- \pi^+ \pi^+)$")
plt.ylabel('Ratio $\pi^+ / \pi^-$')
plt.legend(fontsize=12)
plt.show()


# In[20]:


ratio_with_uncertainty = [f"{value:.3f} +/- {error:.3f}" for value, error in zip(pi_3_pi_1_p_ratio, pi_3_pi_1_p_ratio_err)]

# Print or use the formatted values
for i, ratio in enumerate(ratio_with_uncertainty):
    print(f"bin {i+1}: {ratio}")


# In[21]:


# Compute ratio for each bin with uncertainity
hist_d_pi_3_cosTheta, bin_edges = np.histogram(np_d_pi_3_cosTheta_sw, bins=nbins, range=cosTheta_range, weights=dweights[dsig_yield], density=False)
hist_d_pi_1_cosTheta, _ = np.histogram(np_d_pi_1_cosTheta_sw, bins=nbins, range=cosTheta_range, weights=dweights[dsig_yield], density=False)

pi_3_pi_1_cosTheta_ratio = (hist_d_pi_3_cosTheta/hist_d_pi_1_cosTheta)

hist_d_pi_3_cosTheta_err = np.sqrt(hist_d_pi_3_cosTheta)
hist_d_pi_1_cosTheta_err = np.sqrt(hist_d_pi_1_cosTheta)

pi_3_pi_1_cosTheta_ratio_err = pi_3_pi_1_cosTheta_ratio * np.sqrt((hist_d_pi_3_cosTheta_err / hist_d_pi_3_cosTheta)**2 + (hist_d_pi_1_cosTheta_err / hist_d_pi_1_cosTheta)**2)

bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

# Plotting with uncertainties
plt.errorbar(bin_centers, pi_3_pi_1_cosTheta_ratio, yerr=pi_3_pi_1_cosTheta_ratio_err, fmt='o')
plt.xlabel(r"$\pi \ cosTheta (D^0 \rightarrow K^- \pi^- \pi^+ \pi^+)$")
plt.ylabel('Ratio $\pi^+ / \pi^-$')
plt.legend(fontsize=12)
plt.show()


# In[22]:


ratio_with_uncertainty = [f"{value:.3f} +/- {error:.3f}" for value, error in zip(pi_3_pi_1_cosTheta_ratio, pi_3_pi_1_cosTheta_ratio_err)]

# Print or use the formatted values
for i, ratio in enumerate(ratio_with_uncertainty):
    print(f"bin {i+1}: {ratio}")


# In[ ]:





# In[ ]:





# In[ ]:




