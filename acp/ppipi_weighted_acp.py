# %%
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
import uncertainties as un
from uncertainties import unumpy


# %%
import ROOT as r
r.gROOT.LoadMacro('/belle2work/BelleII/belle2style/Belle2Style.C') 
r.SetBelle2Style()

# Make nice looking plots
plt.rcParams.update({
          'font.size': 20,
          'figure.figsize': (12, 8),
})


# %%
# Define columns to read into memory
col = ["Lambdac_M","Lambdac_cosTheta_cms","p_charge"]
ctrcol = ["reWeights"]
mccol =["Lambdac_isSignal"]

siginfile = '/b2diska/dinura/rootfiles/ntuples/MC/lcp/ppipi_s.root'
ctrinfile = '/b2diska/dinura/rootfiles/ntuples/MC/lcp/weights/pkpi_ppipi_weighted.root'


# %%
sigdf = root_pandas.read_root(siginfile, key='lcp_ppipi', columns=col+mccol, where='Lambdac_M > 2.24 && Lambdac_M < 2.34')


# %%
ctrdf = root_pandas.read_root(ctrinfile, key='lcp_pkpi', columns=col+ctrcol+mccol, where= 'Lambdac_M > 2.24 && Lambdac_M < 2.34')


# %%
# Standard plot settings
lupper = 2.34
llower = 2.24


# %%
# Standard cuts for lambda_c and D^0
lc_cut = str(llower)+' < Lambdac_M < '+str(lupper)


# %%
binedges=[-1,0,1]
bin1 = str(binedges[0])+"<=Lambdac_cosTheta_cms<"+str(binedges[1])
bin2 = str(binedges[1])+"<=Lambdac_cosTheta_cms<"+str(binedges[2])
print("bin range")
print(bin1)
print(bin2)


# %%
# Method to overlay fitted pdfs and data sample
def plot_model(model, mydata, nevents, nbins=100, myrange=(llower,lupper), mylabel="", plot_data=True, save_dir=None, save_name="plot"):

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
    
    # Save the plot if save_dir is provided
    if save_dir:
        save_path = f"{save_dir}{save_name}.png"
        plt.savefig(save_path)
        print(f"Plot saved at {save_path}")




# %%
# # Apply weights for A$_{CP}$ fit
# find the momentum and cos(theta) bin for each event
pkpidf = ctrdf.query(lc_cut)
pkpidf_bin1 = ctrdf.query(lc_cut+" and "+bin1+' and p_charge==1')
pkpidf_bin2 = ctrdf.query(lc_cut+" and "+bin2+' and p_charge==1')
apkpidf_bin1 = ctrdf.query(lc_cut+" and "+bin1+' and p_charge==-1')
apkpidf_bin2 = ctrdf.query(lc_cut+" and "+bin2+' and p_charge==-1')


# %%
# this is a little dangerous since we replace keys with indices that depend on the order 
def getWeights(df):
    np_arr = df.to_numpy()
    return np_arr[:,0], np_arr[:,3]


# %%
pkpi_M, pkpi_weights = getWeights(pkpidf)
pkpi_bin1_M, pkpi_bin1_weights = getWeights(pkpidf_bin1)
pkpi_bin2_M, pkpi_bin2_weights = getWeights(pkpidf_bin2)
apkpi_bin1_M, apkpi_bin1_weights = getWeights(apkpidf_bin1)
apkpi_bin2_M, apkpi_bin2_weights = getWeights(apkpidf_bin2)

# %%
# # ZFit for $\Lambda_c^+$

# Define default parameters for zfit
lcrange = (llower,lupper)
lcobs = zfit.Space('Lambdac_M', lcrange)
issignal = 'Lambdac_isSignal==1'

# Get the signal and background for reference
signal = sigdf.query(lc_cut + ' and Lambdac_isSignal==1').Lambdac_M.to_numpy()
bkg = sigdf.query(lc_cut + ' and Lambdac_isSignal!=1').Lambdac_M.to_numpy()
print('signal: ' + str(len(signal)))
print('bkg: ' + str(len(bkg)))


# %%
# Signal fit parameters
mu = zfit.Parameter("mu", 2.289, 2.27, 2.3)
s1 = zfit.param.Parameter("s1", 0.020, 0.0001, 0.04)
s2 = zfit.param.Parameter("s2", 0.005, 0.0001, 0.04)
fg1 = zfit.param.Parameter("fg1", 0.20, 0., 1.)
a1 = zfit.Parameter("a1", 0.01, -1e6, 1e6)
a2 = zfit.Parameter("a2", 0.01, -1e6, 1e6)

# define PDFs
gaus1 = zfit.pdf.Gauss(obs=lcobs, mu=mu, sigma=s1)
gaus2 = zfit.pdf.Gauss(obs=lcobs, mu=mu, sigma=s2)
gaus = zfit.pdf.SumPDF(pdfs=[gaus1,gaus2], fracs=[fg1])

poly = zfit.pdf.Chebyshev(obs=lcobs, coeffs=[a1, a2])
sig_yield = zfit.Parameter('sig_yield', 100000, 0, 1e9, step_size=1)
bkg_yield = zfit.Parameter('bkg_yield', 500000, 0, 1e9, step_size=1)

sig_ext = gaus.create_extended(sig_yield)
bkg_ext = poly.create_extended(bkg_yield)

pdf_ext = zfit.pdf.SumPDF(pdfs=[sig_ext,bkg_ext])


# # Signal Fit for $\Lambda_c^+\$

# %%
# Fill an array with the data to be fit
data_np = sigdf.query(lc_cut+' and '+issignal).Lambdac_M.to_numpy()
data = zfit.data.Data.from_numpy(obs=lcobs, array=data_np)
data.set_data_range(lcrange) 
data_size = data.n_events


# %%
# Fix background to zero
a1.floating = False
a2.floating = False
bkg_yield.set_value(0)
bkg_yield.floating = False

# Define loss function
nll_simultaneous = zfit.loss.ExtendedUnbinnedNLL(model=pdf_ext,data=data)
minimizer = zfit.minimize.Minuit()
result = minimizer.minimize(nll_simultaneous)  
print("signal fit - signal channel")
print(result)


# %%
plot_model(model=pdf_ext, mydata=data, nevents=int(data_size), mylabel=r'$\Lambda_c$ Mass ($\Lambda_c^+\rightarrow p \pi^- \pi^+$)  [GeV/$c^{2}$]', plot_data=True, save_dir="/b2diska/dinura/scripts/acp/py/ppipi/plots/weighted/", save_name="signal_fit_signalchannel")


# %%
# Save the signal parameters
lcsigmean=result.params[mu].get('value')
lcsigsigma1=result.params[s1].get('value')
lcsigsigma2=result.params[s2].get('value')
lcsigfg1=result.params[fg1].get('value')

print("signal parameters")
print('lcsigmean =',r'%.5f' %lcsigmean)
print('lcsigsigma1 =',r'%.5f' %lcsigsigma1)
print('lcsigsigma2 =',r'%.5f' %lcsigsigma2)
print('lcsigfg1 =',r'%.5f' %lcsigfg1)

# %%
# # Intergrated Fits 
# ## $\Lambda_c^+ \rightarrow pK^-\pi^+$ - control channel integrated fit

# Fill an array with the data to be fit
data_np = ctrdf.query(lc_cut).Lambdac_M.to_numpy()
data = zfit.data.Data.from_numpy(obs=lcobs, array=data_np, weights=pkpi_weights)
data.set_data_range(lcrange)
data_size = data.n_events

# %%
# Float background
a1.floating = True
a2.floating = True

sig_yield.set_value(2000000)
sig_yield.floating = True
bkg_yield.set_value(5000000)
bkg_yield.floating = True

# Fix signal shape
s1.floating = True
s1.set_value(lcsigsigma1)
s2.floating = True
s2.set_value(lcsigsigma2)
fg1.floating = True
fg1.set_value(lcsigfg1)

# Define loss function
nll_simultaneous = zfit.loss.ExtendedUnbinnedNLL(model=pdf_ext,data=data)
minimizer = zfit.minimize.Minuit()
result = minimizer.minimize(nll_simultaneous)  
print("intergrated fit - control channel")
print(result)


# %%
plot_model(model=pdf_ext, mydata=data, myrange=lcrange, nevents=int(data_size), mylabel=r'$\Lambda_c$ Mass ($\Lambda_c^+\rightarrow p K^- \pi^+$)  [GeV/$c^{2}$]', plot_data=True, save_dir="/b2diska/dinura/scripts/acp/py/ppipi/plots/weighted/", save_name="intergrated_fit_controlchannel")


# %%
# Save the control parameters
lcctrmean=mu
lcctrsigma1=s1
lcctrsigma2=s2
lcctrfg1=fg1

mu.floating = False
s1.floating = False
s2.floating = False
fg1.floating = False

print("control parameters")
print('lcctrmean =',r'%.5f' %lcctrmean)
print('mu =',r'%.5f' %lcctrmean)
print('lcctrsigma1 =',r'%.5f' %lcctrsigma1)
print('lcctrsigma2 =',r'%.5f' %lcctrsigma2)
print('lcctrfg1 =',r'%.5f' %lcctrfg1)

# %%
# # Fits in cos($\theta$) bins for $\Lambda_c^+\rightarrow p \pi^- \pi^+$ & $\Lambda_c^+\rightarrow p K^- \pi^+$

nsigtrue=[]
nasigtrue=[]
nsig=[]
nasig=[]

nsigtruee=[]
nasigtruee=[]
nsige=[]
nasige=[]


# %%
nctrtrue=[]
nactrtrue=[]
nctr=[]
nactr=[]

nctrtruee=[]
nactrtruee=[]
nctre=[]
nactre=[]


# %%
# Common shift/scale factors
smu = zfit.param.Parameter("smu", 0., -0.01, 0.01)
ss1 = zfit.param.Parameter("ss1", 1.0, 0., 5.0)
ss2 = zfit.param.Parameter("ss2", 1.0, 0., 5.0)

# %%
class LcSigGauss1(zfit.pdf.ZPDF):
    _N_OBS = 1  # dimension, can be omitted
    _PARAMS = ['mean', 'std']  # the name of the parameters

    def _unnormalized_pdf(self, x):
        x = z.unstack_x(x)  # returns a list with the columns: do x, y, z = z.unstack_x(x) for 3D
        mean = self.params['mean']
        std = self.params['std']
        return 1/(std*lcsigsigma1*np.sqrt(2*np.pi))*z.exp(- ((x - (mean+lcsigmean)) / (std*lcsigsigma1)) ** 2)


# %%
class LcSigGauss2(zfit.pdf.ZPDF):
    _N_OBS = 1  # dimension, can be omitted
    _PARAMS = ['mean', 'std']  # the name of the parameters

    def _unnormalized_pdf(self, x):
        x = z.unstack_x(x)  # returns a list with the columns: do x, y, z = z.unstack_x(x) for 3D
        mean = self.params['mean']
        std = self.params['std']
        return 1/(std*lcsigsigma2*np.sqrt(2*np.pi))*z.exp(- ((x - (mean+lcsigmean)) / (std*lcsigsigma2)) ** 2)


# %%
class LcCtrGauss1(zfit.pdf.ZPDF):
    _N_OBS = 1  # dimension, can be omitted
    _PARAMS = ['mean', 'std']  # the name of the parameters

    def _unnormalized_pdf(self, x):
        x = z.unstack_x(x)  # returns a list with the columns: do x, y, z = z.unstack_x(x) for 3D
        mean = self.params['mean']
        std = self.params['std']
        return 1/(std*lcctrsigma1*np.sqrt(2*np.pi))*z.exp(- ((x - (mean+lcctrmean)) / (std*lcctrsigma1)) ** 2)


# %%
class LcCtrGauss2(zfit.pdf.ZPDF):
    _N_OBS = 1  # dimension, can be omitted
    _PARAMS = ['mean', 'std']  # the name of the parameters

    def _unnormalized_pdf(self, x):
        x = z.unstack_x(x)  # returns a list with the columns: do x, y, z = z.unstack_x(x) for 3D
        mean = self.params['mean']
        std = self.params['std']
        return 1/(std*lcctrsigma2*np.sqrt(2*np.pi))*z.exp(- ((x - (mean+lcctrmean)) / (std*lcctrsigma2)) ** 2)


# %%
# signal fit parameters
sig_fg1 = zfit.Parameter("sig_fg1", 0.20, 0., 1.)
asig_fg1 = zfit.Parameter("asig_fg1", 0.20, 0., 1.)

sig_a1 = zfit.Parameter("sig_a1", -0.03, -1e6, 1e6)
asig_a1 = zfit.Parameter("asig_a1", -0.03, -1e6, 1e6)

sig_a2 = zfit.Parameter("sig_a2", 0.01, -1e6, 1e6)
asig_a2 = zfit.Parameter("asig_a2", 0.01, -1e6, 1e6)

# define PDFs
sig_gaus1 = LcSigGauss1(obs=lcobs, mean=smu, std=ss1)
sig_gaus2 = LcSigGauss2(obs=lcobs, mean=smu, std=ss2)
sig_gaus = zfit.pdf.SumPDF(pdfs=[sig_gaus1,sig_gaus2], fracs=[sig_fg1])

sig_poly = zfit.pdf.Chebyshev(obs=lcobs, coeffs=[sig_a1, sig_a2])

sig_sig_yield = zfit.Parameter('sig_sig_yield', 85000, 0, 0.2e6, step_size=1)
sig_bkg_yield = zfit.Parameter('sig_bkg_yield', 2000000, 0, 0.5e7, step_size=1)

sig_sig_ext = sig_gaus.create_extended(sig_sig_yield)
sig_bkg_ext = sig_poly.create_extended(sig_bkg_yield)

sig_pdf_ext = zfit.pdf.SumPDF(pdfs=[sig_sig_ext,sig_bkg_ext])

asig_gaus1 = LcSigGauss1(obs=lcobs, mean=smu, std=ss1)
asig_gaus2 = LcSigGauss2(obs=lcobs, mean=smu, std=ss2)
asig_gaus = zfit.pdf.SumPDF(pdfs=[asig_gaus1,asig_gaus2], fracs=[asig_fg1])

asig_poly = zfit.pdf.Chebyshev(obs=lcobs, coeffs=[asig_a1, asig_a2])

asig_sig_yield = zfit.Parameter('asig_sig_yield', 85000, 0, 0.2e6, step_size=1)
asig_bkg_yield = zfit.Parameter('asig_bkg_yield', 2000000, 0, 0.5e7, step_size=1)

asig_sig_ext = asig_gaus.create_extended(asig_sig_yield)
asig_bkg_ext = asig_poly.create_extended(asig_bkg_yield)

asig_pdf_ext = zfit.pdf.SumPDF(pdfs=[asig_sig_ext,asig_bkg_ext])


# %%
# Fit parameters
ctr_fg1 = zfit.Parameter("ctr_fg1", 0.02, 0., 1.)
actr_fg1 = zfit.Parameter("actr_fg1", 0.02, 0., 1.)

ctr_a1 = zfit.Parameter("ctr_a1", 0.01, -1e6, 1e6)
actr_a1 = zfit.Parameter("actr_a1", 0.01, -1e6, 1e6)

ctr_a2 = zfit.Parameter("ctr_a2", 0.01, -1e6, 1e6)
actr_a2 = zfit.Parameter("actr_a2", 0.01, -1e6, 1e6)

# define PDFs
ctr_gaus1 = LcCtrGauss1(obs=lcobs, mean=smu, std=ss1)
ctr_gaus2 = LcCtrGauss2(obs=lcobs, mean=smu, std=ss2)
ctr_gaus = zfit.pdf.SumPDF(pdfs=[ctr_gaus1,ctr_gaus2], fracs=[ctr_fg1])

ctr_poly = zfit.pdf.Chebyshev(obs=lcobs, coeffs=[ctr_a1, ctr_a2])

ctr_sig_yield = zfit.Parameter('ctr_sig_yield', 500000, 0, 1e6, step_size=1)
ctr_bkg_yield = zfit.Parameter('ctr_bkg_yield', 800000, 0, 0.2e7, step_size=1)

ctr_sig_ext = ctr_gaus.create_extended(ctr_sig_yield)
ctr_bkg_ext = ctr_poly.create_extended(ctr_bkg_yield)

ctr_pdf_ext = zfit.pdf.SumPDF(pdfs=[ctr_sig_ext,ctr_bkg_ext])

actr_gaus1 = LcCtrGauss1(obs=lcobs, mean=smu, std=ss1)
actr_gaus2 = LcCtrGauss2(obs=lcobs, mean=smu, std=ss2)
actr_gaus = zfit.pdf.SumPDF(pdfs=[actr_gaus1,actr_gaus2], fracs=[actr_fg1])

actr_poly = zfit.pdf.Chebyshev(obs=lcobs, coeffs=[actr_a1, actr_a2])

actr_sig_yield = zfit.Parameter('actr_sig_yield', 500000, 0, 1e6, step_size=1)
actr_bkg_yield = zfit.Parameter('actr_bkg_yield', 800000, 0, 0.2e7, step_size=1)

actr_sig_ext = actr_gaus.create_extended(actr_sig_yield)
actr_bkg_ext = actr_poly.create_extended(actr_bkg_yield)

actr_pdf_ext = zfit.pdf.SumPDF(pdfs=[actr_sig_ext,actr_bkg_ext])


# %%
# Fix values from signal fit
sig_fg1.set_value(lcsigfg1)
sig_fg1.floating = False
asig_fg1.set_value(lcsigfg1)
asig_fg1.floating = False

# %%
# ## bin 1

print(bin1)


# %%
# Get the signal and background for reference
print("truthed-matched values bin1")

sig_signal = sigdf.query(lc_cut+" and "+bin1+' and Lambdac_isSignal==1 and p_charge==1').Lambdac_M.to_numpy()
sig_bkg = sigdf.query(lc_cut+" and "+bin1+' and Lambdac_isSignal!=1 and p_charge==1').Lambdac_M.to_numpy()
print('sig signal bin1: ' + str(len(sig_signal)))
print('sig bkg bin1: ' + str(len(sig_bkg)))
nsigtrue.append(len(sig_signal))

asig_signal = sigdf.query(lc_cut+" and "+bin1+' and Lambdac_isSignal==1 and p_charge==-1').Lambdac_M.to_numpy()
asig_bkg = sigdf.query(lc_cut+" and "+bin1+' and Lambdac_isSignal!=1 and p_charge==-1').Lambdac_M.to_numpy()
print('asig signal bin1: ' + str(len(asig_signal)))
print('asig bkg bin1: ' + str(len(asig_bkg)))
nasigtrue.append(len(asig_signal))

# Get the signal and background for reference
ctr_signal = ctrdf.query(lc_cut+" and "+bin1+' and Lambdac_isSignal==1 and p_charge==1').Lambdac_M.to_numpy()
ctr_bkg = ctrdf.query(lc_cut+" and "+bin1+' and Lambdac_isSignal!=1 and p_charge==1').Lambdac_M.to_numpy()
print('ctr signal bin1 - not weighted: ' + str(len(ctr_signal)))
print('ctr bkg bin1 - not weighted: ' + str(len(ctr_bkg)))
nctrtrue.append(len(ctr_signal))

actr_signal = ctrdf.query(lc_cut+" and "+bin1+' and Lambdac_isSignal==1 and p_charge==-1').Lambdac_M.to_numpy()
actr_bkg = ctrdf.query(lc_cut+" and "+bin1+' and Lambdac_isSignal!=1 and p_charge==-1').Lambdac_M.to_numpy()
print('actr signal bin1 - not weighted: ' + str(len(actr_signal)))
print('actr bkg bin1 - not weighted: ' + str(len(actr_bkg)))
nactrtrue.append(len(actr_signal))


# %%
# Fill an array with the data to be fit
sig_data_np = sigdf.query(lc_cut+" and "+bin1+' and p_charge==1').Lambdac_M.to_numpy()
sig_data = zfit.data.Data.from_numpy(obs=lcobs, array=sig_data_np)
sig_data.set_data_range(lcrange)
sig_data_size = sig_data.n_events

asig_data_np = sigdf.query(lc_cut+" and "+bin1+' and p_charge==-1').Lambdac_M.to_numpy()
asig_data = zfit.data.Data.from_numpy(obs=lcobs, array=asig_data_np)
asig_data.set_data_range(lcrange)
asig_data_size = asig_data.n_events

#ctr_data_np = ctrdf.query(lc_cut+" and "+bin1+' and p_charge==1').Lambdac_M.to_numpy()
ctr_data_np = pkpi_bin1_M #weighted mass
ctr_data = zfit.data.Data.from_numpy(obs=lcobs, array=ctr_data_np, weights=pkpi_bin1_weights)
ctr_data.set_data_range(lcrange)
ctr_data_size = ctr_data.n_events

#actr_data_np = ctrdf.query(lc_cut+" and "+bin1+' and p_charge==-1').Lambdac_M.to_numpy()
actr_data_np = apkpi_bin1_M #weighted mass
actr_data = zfit.data.Data.from_numpy(obs=lcobs, array=actr_data_np, weights=apkpi_bin1_weights)
actr_data.set_data_range(lcrange)
actr_data_size = actr_data.n_events

print("bin1 data sizes")
print(int(sig_data_size))
print(int(asig_data_size))
print(int(ctr_data_size))
print(int(actr_data_size))


# %%
# Simultaenous loss
nll_simultaneous = zfit.loss.ExtendedUnbinnedNLL(model=[sig_pdf_ext,asig_pdf_ext,ctr_pdf_ext,actr_pdf_ext], data=[sig_data,asig_data,ctr_data,actr_data])
minimizer = zfit.minimize.Minuit()
result = minimizer.minimize(nll_simultaneous)

# %%
# Simultaenous loss
nll_simultaneous = zfit.loss.ExtendedUnbinnedNLL(model=[sig_pdf_ext,asig_pdf_ext,ctr_pdf_ext,actr_pdf_ext], data=[sig_data,asig_data,ctr_data,actr_data])
minimizer = zfit.minimize.Minuit()
result = minimizer.minimize(nll_simultaneous)

# %%
# Simultaenous loss
nll_simultaneous = zfit.loss.ExtendedUnbinnedNLL(model=[sig_pdf_ext,asig_pdf_ext,ctr_pdf_ext,actr_pdf_ext], data=[sig_data,asig_data,ctr_data,actr_data])
minimizer = zfit.minimize.Minuit()
result = minimizer.minimize(nll_simultaneous)
print(result)


# %%
result.hesse()

errors_sig, new_result_sig = result.errors(params=sig_pdf_ext.get_params())
print("signal component completed")
nsig.append(result.params[sig_sig_yield].get('value'))
nsige.append(result.params[sig_sig_yield].get('hesse').get('error'))

# %%
errors_ctr, new_result_ctr = result.errors(params=ctr_pdf_ext.get_params())
print("control component completed")
nctr.append(result.params[ctr_sig_yield].get('value'))
nctre.append(result.params[ctr_sig_yield].get('hesse').get('error'))

errors_asig, new_result_asig = result.errors(params=asig_pdf_ext.get_params())
print("antimatter signal component completed")
nasig.append(result.params[asig_sig_yield].get('value'))
nasige.append(result.params[asig_sig_yield].get('hesse').get('error'))

# %%
errors_actr, new_result_actr = result.errors(params=actr_pdf_ext.get_params())
print("antimatter control component completed")
nactr.append(result.params[actr_sig_yield].get('value'))
nactre.append(result.params[actr_sig_yield].get('hesse').get('error'))

print(result)

# %%
plot_model(model=sig_pdf_ext, mydata=sig_data, nevents=int(sig_data_size), mylabel=r'$\Lambda_c$ Mass ($\Lambda_{c}^{+}\rightarrow p \pi^- \pi^+$)  [GeV/$c^{2}$]', plot_data=True, save_dir="/b2diska/dinura/scripts/acp/py/ppipi/plots/weighted/", save_name="signalchannel_matter_bin1")


# %%
plot_model(model=asig_pdf_ext, mydata=asig_data, nevents=int(asig_data_size), mylabel=r'$\Lambda_c$ Mass ($\bar{\Lambda}_{c}^{-}\rightarrow \barp \pi^- \pi^+$)  [GeV/$c^{2}$]', plot_data=True, save_dir="/b2diska/dinura/scripts/acp/py/ppipi/plots/weighted/", save_name="signalchannel_anitmatter_bin1")


# %%
plot_model(model=ctr_pdf_ext, mydata=ctr_data, nevents=int(ctr_data_size), mylabel=r'$\Lambda_c$ Mass ($\Lambda_{c}^{+}\rightarrow p K^- \pi^+$)  [GeV/$c^{2}$]', plot_data=True, save_dir="/b2diska/dinura/scripts/acp/py/ppipi/plots/weighted/", save_name="controlchannel_matter_bin1")


# %%
plot_model(model=actr_pdf_ext, mydata=actr_data, nevents=int(actr_data_size), mylabel=r'$\Lambda_c$ Mass ($\bar{\Lambda}_{c}^{-}\rightarrow \barp K^+ \pi^-$)  [GeV/$c^{2}$]', plot_data=True, save_dir="/b2diska/dinura/scripts/acp/py/ppipi/plots/weighted/", save_name="controlchannel_antimatter_bin1")

# %%
# ## bin 2

print(bin2)


# %%
# Get the signal and background for reference
print("truthed-matched values bin2")

sig_signal = sigdf.query(lc_cut+" and "+bin2+' and Lambdac_isSignal==1 and p_charge==1').Lambdac_M.to_numpy()
sig_bkg = sigdf.query(lc_cut+" and "+bin2+' and Lambdac_isSignal!=1 and p_charge==1').Lambdac_M.to_numpy()
print('sig signal bin2: ' + str(len(sig_signal)))
print('sig bkg bin2: ' + str(len(sig_bkg)))
nsigtrue.append(len(sig_signal))

asig_signal = sigdf.query(lc_cut+" and "+bin2+' and Lambdac_isSignal==1 and p_charge==-1').Lambdac_M.to_numpy()
asig_bkg = sigdf.query(lc_cut+" and "+bin2+' and Lambdac_isSignal!=1 and p_charge==-1').Lambdac_M.to_numpy()
print('asig signal bin2: ' + str(len(asig_signal)))
print('asig bkg bin2: ' + str(len(asig_bkg)))
nasigtrue.append(len(asig_signal))

# Get the signal and background for reference
ctr_signal = ctrdf.query(lc_cut+" and "+bin2+' and Lambdac_isSignal==1 and p_charge==1').Lambdac_M.to_numpy()
ctr_bkg = ctrdf.query(lc_cut+" and "+bin2+' and Lambdac_isSignal!=1 and p_charge==1').Lambdac_M.to_numpy()
print('ctr signal bin2 - not weighted: ' + str(len(ctr_signal)))
print('ctr bkg bin2 - not weighted: ' + str(len(ctr_bkg)))
nctrtrue.append(len(ctr_signal))

actr_signal = ctrdf.query(lc_cut+" and "+bin2+' and Lambdac_isSignal==1 and p_charge==-1').Lambdac_M.to_numpy()
actr_bkg = ctrdf.query(lc_cut+" and "+bin2+' and Lambdac_isSignal!=1 and p_charge==-1').Lambdac_M.to_numpy()
print('actr signal bin2 - not weighted: ' + str(len(actr_signal)))
print('actr bkg bin2 - not weighted: ' + str(len(actr_bkg)))
nactrtrue.append(len(actr_signal))


# %%
# Fill an array with the data to be fit
sig_data_np = sigdf.query(lc_cut+" and "+bin2+' and p_charge==1').Lambdac_M.to_numpy()
sig_data = zfit.data.Data.from_numpy(obs=lcobs, array=sig_data_np)
sig_data.set_data_range(lcrange)
sig_data_size = sig_data.n_events

asig_data_np = sigdf.query(lc_cut+" and "+bin2+' and p_charge==-1').Lambdac_M.to_numpy()
asig_data = zfit.data.Data.from_numpy(obs=lcobs, array=asig_data_np)
asig_data.set_data_range(lcrange)
asig_data_size = asig_data.n_events

#ctr_data_np = ctrdf.query(lc_cut+" and "+bin2+' and p_charge==1').Lambdac_M.to_numpy()
ctr_data_np = pkpi_bin2_M #weighted mass
ctr_data = zfit.data.Data.from_numpy(obs=lcobs, array=ctr_data_np, weights=pkpi_bin2_weights)
ctr_data.set_data_range(lcrange)
ctr_data_size = ctr_data.n_events

#actr_data_np = ctrdf.query(lc_cut+" and "+bin2+' and p_charge==-1').Lambdac_M.to_numpy()
actr_data_np = apkpi_bin2_M
actr_data = zfit.data.Data.from_numpy(obs=lcobs, array=actr_data_np, weights=apkpi_bin2_weights)
actr_data.set_data_range(lcrange)
actr_data_size = actr_data.n_events

print("bin2 data sizes")
print(int(sig_data_size))
print(int(asig_data_size))
print(int(ctr_data_size))
print(int(actr_data_size))


# %%
# Simultaenous loss
nll_simultaneous = zfit.loss.ExtendedUnbinnedNLL(model=[sig_pdf_ext,asig_pdf_ext,ctr_pdf_ext,actr_pdf_ext], data=[sig_data,asig_data,ctr_data,actr_data])
minimizer = zfit.minimize.Minuit()
result = minimizer.minimize(nll_simultaneous)

# %%
# Simultaenous loss
nll_simultaneous = zfit.loss.ExtendedUnbinnedNLL(model=[sig_pdf_ext,asig_pdf_ext,ctr_pdf_ext,actr_pdf_ext], data=[sig_data,asig_data,ctr_data,actr_data])
minimizer = zfit.minimize.Minuit()
result = minimizer.minimize(nll_simultaneous)

# %%
# Simultaenous loss
nll_simultaneous = zfit.loss.ExtendedUnbinnedNLL(model=[sig_pdf_ext,asig_pdf_ext,ctr_pdf_ext,actr_pdf_ext], data=[sig_data,asig_data,ctr_data,actr_data])
minimizer = zfit.minimize.Minuit()
result = minimizer.minimize(nll_simultaneous)
print(result)

# %%
result.hesse()

# %%
errors_sig, new_result_sig = result.errors(params=sig_pdf_ext.get_params())
print("signal component completed")
nsig.append(result.params[sig_sig_yield].get('value'))
nsige.append(result.params[sig_sig_yield].get('hesse').get('error'))

# %%
errors_ctr, new_result_ctr = result.errors(params=ctr_pdf_ext.get_params())
print("control component completed")
nctr.append(result.params[ctr_sig_yield].get('value'))
nctre.append(result.params[ctr_sig_yield].get('hesse').get('error'))

# %%
errors_asig, new_result_asig = result.errors(params=asig_pdf_ext.get_params())
print("antimatter signal component completed")
nasig.append(result.params[asig_sig_yield].get('value'))
nasige.append(result.params[asig_sig_yield].get('hesse').get('error'))

# %%
errors_actr, new_result_actr = result.errors(params=actr_pdf_ext.get_params())
print("antimatter control component completed")
nactr.append(result.params[actr_sig_yield].get('value'))
nactre.append(result.params[actr_sig_yield].get('hesse').get('error'))


# %%
print(result)

# %%
plot_model(model=sig_pdf_ext, mydata=sig_data, nevents=int(sig_data_size), mylabel=r'$\Lambda_c$ Mass ($\Lambda_{c}^{+}\rightarrow p \pi^- \pi^+$)  [GeV/$c^{2}$]', plot_data=True, save_dir="/b2diska/dinura/scripts/acp/py/ppipi/plots/weighted/", save_name="signalchannel_matter_bin2")


# %%
plot_model(model=asig_pdf_ext, mydata=asig_data, nevents=int(asig_data_size), mylabel=r'$\Lambda_c$ Mass ($\bar\Lambda_{c}^{-}\rightarrow \barp \pi^- \pi^+$)  [GeV/$c^{2}$]', plot_data=True, save_dir="/b2diska/dinura/scripts/acp/py/ppipi/plots/weighted/", save_name="signalchannel_antimatter_bin2")


# %%
plot_model(model=ctr_pdf_ext, mydata=ctr_data, nevents=int(ctr_data_size), mylabel=r'$\Lambda_c$ Mass ($\Lambda_{c}^{+}\rightarrow p K^- \pi^+$)  [GeV/$c^{2}$]', plot_data=True, save_dir="/b2diska/dinura/scripts/acp/py/ppipi/plots/weighted/", save_name="controlchannel_matter_bin2")


# %%
plot_model(model=actr_pdf_ext, mydata=actr_data, nevents=int(actr_data_size), mylabel=r'$\Lambda_c$ Mass ($\bar\Lambda_{c}^{-}\rightarrow \barp K^+ \pi^-$)  [GeV/$c^{2}$]', plot_data=True, save_dir="/b2diska/dinura/scripts/acp/py/ppipi/plots/weighted/", save_name="controlchannel_antimatter_bin2")


# %%
sigar = unumpy.uarray(nsig,nsige)
asigar = unumpy.uarray(nasig,nasige)
ctrar = unumpy.uarray(nctr,nctre)
actrar = unumpy.uarray(nactr,nactre)

# %%
nsigtruee = np.sqrt(nsigtrue)
nasigtruee = np.sqrt(nasigtrue)
nctrtruee = np.sqrt(nctrtrue)
nactrtruee = np.sqrt(nactrtrue)

# %%
trsigar = unumpy.uarray(nsigtrue,nsigtruee)
trasigar = unumpy.uarray(nasigtrue,nasigtruee)
trctrar = unumpy.uarray(nctrtrue,nctrtruee)
tractrar = unumpy.uarray(nactrtrue,nactrtruee)

# %%
# # Apply weights for true A$_{CP}$

nctrweighted=[]
nactrweighted=[]
nctrweighted_un=[]
nactrweighted_un=[]


# %%
# Get the weighted true Acp value
ctr_signal = ctrdf.query(lc_cut+" and "+bin1+' and Lambdac_isSignal==1 and p_charge==1')['reWeights'].to_numpy()
ctr_signal_2 = ctr_signal*ctr_signal
ctr_bkg = ctrdf.query(lc_cut+" and "+bin1+' and Lambdac_isSignal!=1 and p_charge==1')['reWeights'].to_numpy()
print('weighted ctr signal bin1: ' + str(np.sum(ctr_signal)))
print('weighted ctr bkg bin1: ' + str(np.sum(ctr_bkg)))
nctrweighted.append(np.sum(ctr_signal))
nctrweighted_un.append(np.sum(ctr_signal_2))

actr_signal = ctrdf.query(lc_cut+" and "+bin1+' and Lambdac_isSignal==1 and p_charge==-1')['reWeights'].to_numpy()
actr_signal_2 = actr_signal*actr_signal
actr_bkg = ctrdf.query(lc_cut+" and "+bin1+' and Lambdac_isSignal!=1 and p_charge==-1')['reWeights'].to_numpy()
print('weighted actr signal bin1: ' + str(np.sum(actr_signal)))
print('weighted actr bkg bin1: ' + str(np.sum(actr_bkg)))
nactrweighted.append(np.sum(actr_signal))
nactrweighted_un.append(np.sum(actr_signal_2))


# %%
ctr_signal = ctrdf.query(lc_cut+" and "+bin2+' and Lambdac_isSignal==1 and p_charge==1')['reWeights'].to_numpy()
ctr_signal_2 = ctr_signal*ctr_signal
ctr_bkg = ctrdf.query(lc_cut+" and "+bin2+' and Lambdac_isSignal!=1 and p_charge==1')['reWeights'].to_numpy()
print('weighted ctr signal bin2: ' + str(np.sum(ctr_signal)))
print('weighted ctr bkg bin2: ' + str(np.sum(ctr_bkg)))
nctrweighted.append(np.sum(ctr_signal))
nctrweighted_un.append(np.sum(ctr_signal_2))

actr_signal = ctrdf.query(lc_cut+" and "+bin2+' and Lambdac_isSignal==1 and p_charge==-1')['reWeights'].to_numpy()
actr_signal_2 = actr_signal*actr_signal
actr_bkg = ctrdf.query(lc_cut+" and "+bin2+' and Lambdac_isSignal!=1 and p_charge==-1')['reWeights'].to_numpy()
print('weighted actr signal bin2: ' + str(np.sum(actr_signal)))
print('weighted actr bkg bin2: ' + str(np.sum(actr_bkg)))
nactrweighted.append(np.sum(actr_signal))
nactrweighted_un.append(np.sum(actr_signal_2))


# %%
nctrweightede = np.sqrt(nctrweighted_un)
nactrweightede = np.sqrt(nactrweighted_un)
    
trctrwar = unumpy.uarray(nctrweighted,nctrweightede)
tractrwar = unumpy.uarray(nactrweighted,nactrweightede)


# %%
# # A$_{CP}$ Determination

def getAsym(sig,asig,ctr,actr):
    a1 = ((sig[0]-asig[0])/(sig[0]+asig[0])+(sig[1]-asig[1])/(sig[1]+asig[1]))/2
          
    a2 = ((ctr[0]-actr[0])/(ctr[0]+actr[0])+(ctr[1]-actr[1])/(ctr[1]+actr[1]))/2
     
    acp = a1-a2     
    
    return a1, a2, acp


# %%
print("true acp")
print(getAsym(trsigar,trasigar,trctrar,tractrar)) #True Acp Values


# %%
print("true weighted acp")
print(getAsym(trsigar,trasigar,trctrwar,tractrwar)) #True Weighted Acp


# %%
print("fit acp")
print(getAsym(sigar,asigar,ctrar,actrar)) #Weighted Fit Acp Values

