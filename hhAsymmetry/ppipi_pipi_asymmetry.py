#!/usr/bin/env python
# coding: utf-8

# In[30]:


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
from uncertainties import unumpy as un


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


lccol = ["Lambdac_M","Lambdac_isSignal","p_charge","pi1_p","pi2_p","pi1_cosTheta","pi2_cosTheta"]
lcinfile = '/group/belle2/users2022/dhettiar/ntuples/MC/lcp/ppipi_s.root'


# In[4]:


lcdf = root_pandas.read_root(lcinfile, key='lcp_ppipi', columns=lccol, where='pi1_p < 4.0 && pi2_p <4.0')


# In[5]:


# Standard plot settings
lupper = 2.34
llower = 2.24             


# In[6]:


# Standard cuts for lambda_c and D^0
lc_cut = str(llower)+' < Lambdac_M < '+str(lupper)


# In[7]:


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


# In[8]:


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


# In[9]:


plotVar(lcdf,"Lambdac_M",lc_cut,mylabel=r'M($\Lambda_c^+\rightarrow p \pi^- \pi^+$)')


# In[10]:


# Stack plot of lc_M
np_lc_M_sw = lcdf.query(lc_cut)['Lambdac_M'].to_numpy()
np_lc_pi1_p_sw = lcdf.query(lc_cut)['pi1_p'].to_numpy()
np_lc_pi2_p_sw = lcdf.query(lc_cut)['pi2_p'].to_numpy()
np_lc_pi1_cosTheta_sw = lcdf.query(lc_cut)['pi1_cosTheta'].to_numpy()
np_lc_pi2_cosTheta_sw = lcdf.query(lc_cut)['pi2_cosTheta'].to_numpy()


# # ZFit

# In[11]:


# Define default parameters for zfit
lcrange = (llower,lupper)
lcobs = zfit.Space('Lambdac_M', lcrange)
issignal = 'Lambdac_isSignal==1'

# Get the signal and background for reference
signal = lcdf.query(lc_cut + ' and Lambdac_isSignal==1').Lambdac_M.to_numpy()
bkg = lcdf.query(lc_cut + ' and Lambdac_isSignal!=1').Lambdac_M.to_numpy()
print('signal: ' + str(len(signal)))
print('bkg: ' + str(len(bkg)))


# In[12]:


# Lambdac fit parameters
mu = zfit.Parameter("mu", 2.289, 2.27, 2.3)
s1 = zfit.param.Parameter("s1", 0.020, 0.0001, 0.04)
s2 = zfit.param.Parameter("s2", 0.005, 0.0001, 0.04)
fg1 = zfit.param.Parameter("fg1", 0.20, 0., 1.)
a1 = zfit.Parameter("a1", 0.01, -1e6, 1e6)

# define PDFs
gaus1 = zfit.pdf.Gauss(obs=lcobs, mu=mu, sigma=s1)
gaus2 = zfit.pdf.Gauss(obs=lcobs, mu=mu, sigma=s2)
gaus = zfit.pdf.SumPDF(pdfs=[gaus1,gaus2], fracs=[fg1])

poly = zfit.pdf.Chebyshev(obs=lcobs, coeffs=[a1])
sig_yield = zfit.Parameter('sig_yield', 500000, 0, 1e7, step_size=1)
bkg_yield = zfit.Parameter('bkg_yield', 500000, 0, 1e9, step_size=1)

sig_ext = gaus.create_extended(sig_yield)
bkg_ext = poly.create_extended(bkg_yield)

pdf_ext = zfit.pdf.SumPDF(pdfs=[sig_ext,bkg_ext])


# ## $\Lambda_c$ signal fit

# In[13]:


# Fill an array with the data to be fit
data_np = lcdf.query(lc_cut+' and '+issignal).Lambdac_M.to_numpy()
data = zfit.data.Data.from_numpy(obs=lcobs, array=data_np)
data.set_data_range(lcrange) 
data_size = data.n_events


# In[14]:


# Fix background to zero
a1.floating = False
bkg_yield.set_value(0)
bkg_yield.floating = False

# Define loss function
nll_simultaneous = zfit.loss.ExtendedUnbinnedNLL(model=pdf_ext,data=data)
minimizer = zfit.minimize.Minuit()
result = minimizer.minimize(nll_simultaneous)  
print(result)


# In[15]:


plot_model(model=pdf_ext, mydata=data, nevents=int(data_size), mylabel=r'M($\Lambda_c^+\rightarrow p \pi^- \pi^+$)  [GeV/$c^{2}$]')


# In[16]:


# Save the signal parameters
lcsigmean=result.params[mu].get('value')
lcsigsigma1=result.params[s1].get('value')
lcsigsigma2=result.params[s2].get('value')
lcsigfg1=result.params[fg1].get('value')

print('lcsigmean =',r'%.5f' %lcsigmean)
print('lcsigsigma1 =',r'%.5f' %lcsigsigma1)
print('lcsigsigma2 =',r'%.5f' %lcsigsigma2)
print('lcsigfg1 =',r'%.5f' %lcsigfg1)


# ## $\Lambda_c$ Integrated fit

# In[17]:


# Fill an array with the data to be fit
data_np = lcdf.query(lc_cut).Lambdac_M.to_numpy()
data = zfit.data.Data.from_numpy(obs=lcobs, array=data_np)
data.set_data_range(lcrange) 
data_size = data.n_events


# In[18]:


# Float background
a1.floating = True
bkg_yield.set_value(5000000)
bkg_yield.floating = True

# Fix signal shape
s1.floating = False
s1.set_value(lcsigsigma1)
s2.floating = False
s2.set_value(lcsigsigma2)
fg1.floating = False
fg1.set_value(lcsigfg1)

# Define loss function
nll_simultaneous = zfit.loss.ExtendedUnbinnedNLL(model=pdf_ext,data=data)
minimizer = zfit.minimize.Minuit()
result = minimizer.minimize(nll_simultaneous)  
print(result)


# In[19]:


plot_model(model=pdf_ext, mydata=data, myrange=lcrange, nevents=int(data_size), mylabel=r'M($\Lambda_c^+\rightarrow p \pi^- \pi^+$)  [GeV/$c^{2}$]')


# In[20]:


lcweights = compute_sweights(pdf_ext, data)


# # sPlot

# In[21]:


nbins = 25
nrange = (0,4.0)
cosTheta_range= (-1,1)


# In[22]:


# pkk matter channel kaon momentum distributions
fig, axs = plt.subplots(1, 2, figsize=(30, 10))
#nbins = 10

axs[0].hist(np_lc_pi1_p_sw, bins=nbins, range=nrange, weights=lcweights[sig_yield], label=r"$\Lambda_c \ (\pi^+ \ momentum)$", alpha=.5, density=False)
axs[0].hist(np_lc_pi2_p_sw, bins=nbins, range=nrange, weights=lcweights[sig_yield], label=r"$\Lambda_c \ (\pi^- \ momentum)$", alpha=.5, density=False)
axs[0].set_xlabel(r"$\Lambda_c \ (\pi \ momentum)$")
axs[0].legend(fontsize=20)

axs[1].hist(np_lc_pi1_cosTheta_sw, bins=nbins, range=cosTheta_range, weights=lcweights[sig_yield], label=r"$\Lambda_c \ (\pi^+ \ cosTheta)$", alpha=.5, density=False)
axs[1].hist(np_lc_pi2_cosTheta_sw, bins=nbins, range=cosTheta_range, weights=lcweights[sig_yield],label=r"$\Lambda_c \ (\pi^- \ cosTheta)$", alpha=.5, density=False)
axs[1].set_xlabel(r"$\Lambda_c \ (\pi \ cosTheta)$")
axs[1].legend(fontsize=20)


# # Difference

# In[48]:


# Compute ratio for each bin with uncertainity
hist_lc_pi1_p, bin_edges = np.histogram(np_lc_pi1_p_sw, bins=nbins, range=nrange, weights=lcweights[sig_yield], density=False)
hist_lc_pi2_p, _ = np.histogram(np_lc_pi2_p_sw, bins=nbins, range=nrange, weights=lcweights[sig_yield], density=False)

pi1_pi2_p_ratio = (hist_lc_pi1_p/hist_lc_pi2_p)

hist_lc_pi1_p_err = np.sqrt(hist_lc_pi1_p)
hist_lc_pi2_p_err = np.sqrt(hist_lc_pi2_p)

pi1_pi2_p_ratio_err = pi1_pi2_p_ratio * np.sqrt((hist_lc_pi1_p_err / hist_lc_pi1_p)**2 + (hist_lc_pi2_p_err / hist_lc_pi2_p)**2)

bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

# Plotting with uncertainties
plt.errorbar(bin_centers, pi1_pi2_p_ratio, yerr=pi1_pi2_p_ratio_err, fmt='o')
plt.xlabel(r"$\pi \ momentum (\Lambda_c^+\rightarrow p \pi^- \pi^+)$")
plt.ylabel('Ratio $\pi^+ / \pi^-$')
plt.legend(fontsize=12)
plt.show()


# In[44]:


ratio_with_uncertainty = [f"{value:.3f} +/- {error:.3f}" for value, error in zip(pi1_pi2_p_ratio, pi1_pi2_p_ratio_err)]

# Print or use the formatted values
for i, ratio in enumerate(ratio_with_uncertainty):
    print(f"bin {i+1}: {ratio}")


# In[50]:


# Compute ratio for each bin with uncertainity
hist_lc_pi1_cosTheta, bin_edges = np.histogram(np_lc_pi1_cosTheta_sw, bins=nbins, range=cosTheta_range, weights=lcweights[sig_yield], density=False)
hist_lc_pi2_cosTheta, _ = np.histogram(np_lc_pi2_cosTheta_sw, bins=nbins, range=cosTheta_range, weights=lcweights[sig_yield], density=False)

pi1_pi2_cosTheta_ratio = (hist_lc_pi1_cosTheta/hist_lc_pi2_cosTheta)

hist_lc_pi1_cosTheta_err = np.sqrt(hist_lc_pi1_cosTheta)
hist_lc_pi2_cosTheta_err = np.sqrt(hist_lc_pi2_cosTheta)

pi1_pi2_cosTheta_ratio_err = pi1_pi2_cosTheta_ratio * np.sqrt((hist_lc_pi1_cosTheta_err / hist_lc_pi1_cosTheta)**2 + (hist_lc_pi2_cosTheta_err / hist_lc_pi2_cosTheta)**2)

bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

# Plotting with uncertainties
plt.errorbar(bin_centers, pi1_pi2_cosTheta_ratio, yerr=pi1_pi2_cosTheta_ratio_err, fmt='o')
plt.xlabel(r"$\pi \ cosTheta (\Lambda_c^+\rightarrow p \pi^- \pi^+)$")
plt.ylabel('Ratio $(\pi^+ / \pi^-)$')
plt.legend(fontsize=12)
plt.show()


# In[51]:


ratio_with_uncertainty = [f"{value:.3f} +/- {error:.3f}" for value, error in zip(pi1_pi2_cosTheta_ratio, pi1_pi2_cosTheta_ratio_err)]

# Print or use the formatted values
for i, ratio in enumerate(ratio_with_uncertainty):
    print(f"bin {i+1}: {ratio}")


# In[ ]:




