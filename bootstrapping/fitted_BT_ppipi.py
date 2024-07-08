# %%
import zfit
from zfit import z
import root_pandas
import matplotlib.pyplot as plt
import tensorflow as tf
#ZFIT_DISABLE_TF_WARNINGS=1
import numpy as np
import mplhep
import uncertainties as un
from uncertainties import unumpy
from hepstats.splot import compute_sweights

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
col = ["Lambdac_M","Lambdac_cosTheta_cms","p_charge","p_cosTheta","p_p"]
mccol =["Lambdac_isSignal"]

siginfile = '/b2diska/dinura/rootfiles/ntuples/MC/lcp/ppipi_s.root'
ctrinfile = '/b2diska/dinura/rootfiles/ntuples/MC/lcp/pkpi_ppipi_s.root'

# %%
#bootstrap sample size in %
smplsize = 0.25
# how many bootstrap samples
b = 200

# %%
# Standard plot settings
lupper = 2.34
llower = 2.24

# %%
# Standard cuts for lambda_c
lc_cut = str(llower)+' < Lambdac_M < '+str(lupper)


# %%
binedges=[-1,0,1]
bin1 = str(binedges[0])+"<=Lambdac_cosTheta_cms<"+str(binedges[1])
bin2 = str(binedges[1])+"<=Lambdac_cosTheta_cms<"+str(binedges[2])
print(bin1)
print(bin2)

# %%
sigdf_main = root_pandas.read_root(siginfile, key='lcp_ppipi', columns=col+mccol, where='Lambdac_M > 2.24 && Lambdac_M < 2.34')


# %%
ctrdf_main = root_pandas.read_root(ctrinfile, key='lcp_pkpi', columns=col+mccol, where= 'Lambdac_M > 2.24 && Lambdac_M < 2.34')

# %%
# Define default parameters for zfit
lcrange = (llower,lupper)
lcobs = zfit.Space('Lambdac_M', lcrange)
issignal = 'Lambdac_isSignal==1'


# %%
# Lambdac fit parameters
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
sig_yield = zfit.Parameter('sig_yield', 500000, 0, 1e9, step_size=1)
bkg_yield = zfit.Parameter('bkg_yield', 20000000, 0, 1e9, step_size=1)

sig_ext = gaus.create_extended(sig_yield)
bkg_ext = poly.create_extended(bkg_yield)

pdf_ext = zfit.pdf.SumPDF(pdfs=[sig_ext,bkg_ext])


# %%
# Common shift/scale factors
smu = zfit.param.Parameter("smu", 0.0, -0.01, 0.01)
ss1 = zfit.param.Parameter("ss1", 1.0, 0., 2.5)
ss2 = zfit.param.Parameter("ss2", 1.0, 0., 2.5)


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

sig_a1 = zfit.Parameter("sig_a1", 0.01, -1e6, 1e6)
asig_a1 = zfit.Parameter("asig_a1", 0.01, -1e6, 1e6)

sig_a2 = zfit.Parameter("sig_a2", 0.01, -1e6, 1e6)
asig_a2 = zfit.Parameter("asig_a2", 0.01, -1e6, 1e6)

# define PDFs
sig_gaus1 = LcSigGauss1(obs=lcobs, mean=smu, std=ss1)
sig_gaus2 = LcSigGauss2(obs=lcobs, mean=smu, std=ss2)
sig_gaus = zfit.pdf.SumPDF(pdfs=[sig_gaus1,sig_gaus2], fracs=[sig_fg1])

sig_poly = zfit.pdf.Chebyshev(obs=lcobs, coeffs=[sig_a1,sig_a2])

sig_sig_yield = zfit.Parameter('sig_sig_yield', 20000, 0, 1e5, step_size=1)
sig_bkg_yield = zfit.Parameter('sig_bkg_yield', 600000, 0, 2e6, step_size=1)

sig_sig_ext = sig_gaus.create_extended(sig_sig_yield)
sig_bkg_ext = sig_poly.create_extended(sig_bkg_yield)

sig_pdf_ext = zfit.pdf.SumPDF(pdfs=[sig_sig_ext,sig_bkg_ext])

asig_gaus1 = LcSigGauss1(obs=lcobs, mean=smu, std=ss1)
asig_gaus2 = LcSigGauss2(obs=lcobs, mean=smu, std=ss2)
asig_gaus = zfit.pdf.SumPDF(pdfs=[asig_gaus1,asig_gaus2], fracs=[asig_fg1])

asig_poly = zfit.pdf.Chebyshev(obs=lcobs, coeffs=[asig_a1,asig_a2])

asig_sig_yield = zfit.Parameter('asig_sig_yield', 20000, 0, 1e5, step_size=1)
asig_bkg_yield = zfit.Parameter('asig_bkg_yield', 600000, 0, 2e6, step_size=1)

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

ctr_poly = zfit.pdf.Chebyshev(obs=lcobs, coeffs=[ctr_a1,ctr_a2])

ctr_sig_yield = zfit.Parameter('ctr_sig_yield', 100000, 0, 0.5e6, step_size=1)
ctr_bkg_yield = zfit.Parameter('ctr_bkg_yield', 400000, 0, 1e6, step_size=1)

ctr_sig_ext = ctr_gaus.create_extended(ctr_sig_yield)
ctr_bkg_ext = ctr_poly.create_extended(ctr_bkg_yield)

ctr_pdf_ext = zfit.pdf.SumPDF(pdfs=[ctr_sig_ext,ctr_bkg_ext])

actr_gaus1 = LcCtrGauss1(obs=lcobs, mean=smu, std=ss1)
actr_gaus2 = LcCtrGauss2(obs=lcobs, mean=smu, std=ss2)
actr_gaus = zfit.pdf.SumPDF(pdfs=[actr_gaus1,actr_gaus2], fracs=[actr_fg1])

actr_poly = zfit.pdf.Chebyshev(obs=lcobs, coeffs=[actr_a1,actr_a2])

actr_sig_yield = zfit.Parameter('actr_sig_yield', 100000, 0, 0.5e6, step_size=1)
actr_bkg_yield = zfit.Parameter('actr_bkg_yield', 400000, 0, 1e6, step_size=1)

actr_sig_ext = actr_gaus.create_extended(actr_sig_yield)
actr_bkg_ext = actr_poly.create_extended(actr_bkg_yield)

actr_pdf_ext = zfit.pdf.SumPDF(pdfs=[actr_sig_ext,actr_bkg_ext])


# %%
# Sample data frames generate==>
print("<==sample data frames are generated==>")

# Print lengths of the original DataFrames
print("length of sigdf_main:", len(sigdf_main))
print("length of ctrdf_Main:", len(ctrdf_main))

# Sample 25% of events (using integer division for floor)
sample_size1 = int(smplsize * len(sigdf_main))
sample_size2 = int(smplsize * len(ctrdf_main))

# %%
# Sample 5 times and create corresponding DataFrames
print("<==Text File Defined==>")

with open("ppipi_fitted_a1.txt", "w") as file_1, open("ppipi_fitted_a2.txt", "w") as file_2:

    for i in range(b):
        # Sample dataframes randomly with replacement
        sigdf = sigdf_main.sample(sample_size1, replace=False)
        ctrdf = ctrdf_main.sample(sample_size2, replace=False)


        # Print information about the sampled DataFrames
        # for i in range(b):
        print(f"\n sigdf_{i} (sampled):")
        print(sigdf.head())  # Print the first few rows

        print(f"\n ctrdf_{i} (sampled):")
        print(ctrdf.head())  # Print the first few rows

        print(f"\n length of sigdf_{i}:", len(sigdf))
        print(f"\n length of ctrdf_{i}:", len(ctrdf))

        # %%
        # Stack plot of signal channel
        np_ppipi_M_sw = sigdf.query(lc_cut)['Lambdac_M'].to_numpy()    
        np_ppipi_p_sw = sigdf.query(lc_cut)['p_p'].to_numpy()       
        np_ppipi_cosTheta_sw = sigdf.query(lc_cut)['p_cosTheta'].to_numpy()     

        # %%
        # Stack plot of control channel
        np_pkpi_M_sw = ctrdf.query(lc_cut)['Lambdac_M'].to_numpy()
        np_pkpi_p_sw = ctrdf.query(lc_cut)['p_p'].to_numpy()      
        np_pkpi_cosTheta_sw = ctrdf.query(lc_cut)['p_cosTheta'].to_numpy()  

        # %%
        # Get the signal and background for reference
        signal = sigdf.query(lc_cut + ' and Lambdac_isSignal==1').Lambdac_M.to_numpy()
        bkg = sigdf.query(lc_cut + ' and Lambdac_isSignal!=1').Lambdac_M.to_numpy()
        print('signal: ' + str(len(signal)))
        print('bkg: ' + str(len(bkg)))

        # %%
        # Fill an array with the signal channel data to be fit - Signal Channel Signal Fit 
        data_np = sigdf.query(lc_cut+' and '+issignal).Lambdac_M.to_numpy()
        data = zfit.data.Data.from_numpy(obs=lcobs, array=data_np)
        data.set_data_range(lcrange) 
        data_size = data.n_events

        # Fix background to zero
        a1.floating = False
        a2.floating = False
        bkg_yield.set_value(0)
        bkg_yield.floating = False

        mu.floating = True
        s1.floating = True
        s2.floating = True
        fg1.floating = True

        # %%
        # Define loss function
        nll_simultaneous = zfit.loss.ExtendedUnbinnedNLL(model=pdf_ext,data=data)
        minimizer = zfit.minimize.Minuit()
        result = minimizer.minimize(nll_simultaneous)  
        print(result)

        # %%
        # Save the signal parameters
        lcsigmean=result.params[mu].get('value')
        lcsigsigma1=result.params[s1].get('value')
        lcsigsigma2=result.params[s2].get('value')
        lcsigfg1=result.params[fg1].get('value')

        # %%
        # Fill an array with the data to be fit - Integrated Fit Signal Channel
        data_np = sigdf.query(lc_cut).Lambdac_M.to_numpy()
        data = zfit.data.Data.from_numpy(obs=lcobs, array=data_np)
        data.set_data_range(lcrange) 
        data_size = data.n_events

        # %%
        # Float background
        a1.floating = True
        bkg_yield.set_value(2500000)
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

        # %%
        ppipiweights = compute_sweights(pdf_ext, data)

        # %%
        # Fill an array with the data to be fit
        data_np = ctrdf.query(lc_cut).Lambdac_M.to_numpy()
        data = zfit.data.Data.from_numpy(obs=lcobs, array=data_np)
        data.set_data_range(lcrange) 
        data_size = data.n_events

        # %%
        # Float background
        a1.floating = True
        a2.floating = True
        
        sig_yield.set_value(500000)
        sig_yield.floating = True
        bkg_yield.set_value(1000000)
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
        print(result)

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

        # %%
        pkpiweights = compute_sweights(pdf_ext, data)

        # %%
        #reweighting
        nbins = 25
        p_range = (0, 5.0)
        cosTheta_range = (-1,1)

        # %%
        pkpiscale = ppipiweights[sig_yield].sum()/pkpiweights[sig_yield].sum()

        # %%
        # Reweighting for proton momentum
        ppipi_counts_per_bin_p = np.histogram(np_ppipi_p_sw, bins=nbins, weights=ppipiweights[sig_yield], range=p_range)
        pkpi_counts_per_bin_p = np.histogram(np_pkpi_p_sw, bins=nbins, weights=[pkpiscale]*pkpiweights[sig_yield], range=p_range)

        # Calculate the event ratio for each bin
        ppipi_counts_per_bin_p = ppipi_counts_per_bin_p[0]  # Extract the counts from the tuple
        pkpi_counts_per_bin_p = pkpi_counts_per_bin_p[0]
        event_ratio_p = ppipi_counts_per_bin_p / pkpi_counts_per_bin_p

        # Calculate the bin edges
        bin_edges_p = np.linspace(p_range[0], p_range[1], nbins + 1)

        # Assign each event to a bin.
        event_bins_p = np.digitize(np_pkpi_p_sw, bin_edges_p)
        event_bins_p -= 1

        ctrdf_pweighted = ctrdf.query(lc_cut)

        # Create an empty list to store weighted events
        pkpi_pweights = []

        for event in range(len(ctrdf_pweighted)):
            event_bin_index_p = event_bins_p[event]

            # Ensure the bin indices are within the valid range
            if event_bin_index_p >= nbins: 
                weight = 1

            else:
                # Calculate weight from the 2D event_ratio for both pi_p and K_p
                pweight = event_ratio_p[event_bin_index_p] 

            pkpi_pweights.append(pweight)    

        # %%
        pkpi_pweights_sw = pkpi_pweights*pkpiweights[sig_yield]

        # %%
        # Reweighting for the proton costheta
        ppipi_counts_per_bin_cosTheta = np.histogram(np_ppipi_cosTheta_sw, bins=nbins, weights=ppipiweights[sig_yield], range=cosTheta_range)
        pkpi_counts_per_bin_cosTheta = np.histogram(np_pkpi_cosTheta_sw, bins=nbins, weights=[pkpiscale]*pkpi_pweights_sw, range=cosTheta_range)

        # Calculate the event ratio for each bin
        ppipi_counts_per_bin_cosTheta = ppipi_counts_per_bin_cosTheta[0]  # Extract the counts from the tuple
        pkpi_counts_per_bin_cosTheta = pkpi_counts_per_bin_cosTheta[0]
        event_ratio_cosTheta = ppipi_counts_per_bin_cosTheta / pkpi_counts_per_bin_cosTheta

        # Calculate the bin edges
        bin_edges_cosTheta = np.linspace(cosTheta_range[0], cosTheta_range[1], nbins + 1)

        # I want to scale up and compair with "np_xic_sig_p_sw"
        event_bins_cosTheta = np.digitize(np_pkpi_cosTheta_sw, bin_edges_cosTheta)
        event_bins_cosTheta -= 1

        ctrdf_cweighted = ctrdf.query(lc_cut)

        # Create an empty list to store weighted events
        pkpi_cweights = []

        for event in range(len(ctrdf_cweighted)):
            event_bin_index_cosTheta = event_bins_cosTheta[event]

            # Ensure the bin indices are within the valid range
            if event_bin_index_cosTheta >= nbins: 
                weight = 1

            else:
                # Calculate weight from the 2D event_ratio for both pi_p and K_p
                cweight = event_ratio_cosTheta[event_bin_index_cosTheta] 

            pkpi_cweights.append(cweight)    

        # %%
        reWeights = [a*b for a,b in zip(pkpi_pweights,pkpi_cweights)]

        # %%
        ctrdf['reWeights'] = reWeights

        # %%
        reWeights_sw = reWeights*pkpiweights[sig_yield]

        # %%
        # find the momentum and cos(theta) bin for each event
        pkpidf = ctrdf.query(lc_cut)
        pkpidf_bin1 = ctrdf.query(lc_cut+" and "+bin1+' and p_charge==1')
        pkpidf_bin2 = ctrdf.query(lc_cut+" and "+bin2+' and p_charge==1')
        apkpidf_bin1 = ctrdf.query(lc_cut+" and "+bin1+' and p_charge==-1')
        apkpidf_bin2 = ctrdf.query(lc_cut+" and "+bin2+' and p_charge==-1')

        # %%
        pkpidf_bin1

        # %%
        def getWeights(df):
            np_arr = df.to_numpy()
            return np_arr[:,0], np_arr[:,6]

        # %%
        pkpi_bin1_M, pkpi_bin1_weights = getWeights(pkpidf_bin1)
        pkpi_bin2_M, pkpi_bin2_weights = getWeights(pkpidf_bin2)
        apkpi_bin1_M, apkpi_bin1_weights = getWeights(apkpidf_bin1)
        apkpi_bin2_M, apkpi_bin2_weights = getWeights(apkpidf_bin2)

        # %%
        nsig=[]
        nasig=[]

        nsige=[]
        nasige=[]

        # %%
        nctr=[]
        nactr=[]

        nctre=[]
        nactre=[]

        # %%
        # Fix values from signal fit
        sig_fg1.set_value(lcsigfg1)
        sig_fg1.floating = False
        asig_fg1.set_value(lcsigfg1)
        asig_fg1.floating = False

        # %%
        print(bin1)

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
        nsig.append(result.params[sig_sig_yield].get('value'))
        nasig.append(result.params[asig_sig_yield].get('value'))
        nctr.append(result.params[ctr_sig_yield].get('value'))
        nactr.append(result.params[actr_sig_yield].get('value'))

        # %%
        print(bin2)

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
        nsig.append(result.params[sig_sig_yield].get('value'))
        nasig.append(result.params[asig_sig_yield].get('value'))
        nctr.append(result.params[ctr_sig_yield].get('value'))
        nactr.append(result.params[actr_sig_yield].get('value'))

        # %%
        # Initialize empty arrays for a1, a2, acp
        fitted_a1Noerr_values = []
        fitted_a2Noerr_values = []

        # %%
        def Fitted_a1_getAsymNoErr(sig1,sig2):
            a1_er_array = np.divide(np.subtract(sig1,sig2),np.add(sig1,sig2))
            a1_er = (np.sum(a1_er_array))/2

            fitted_a1Noerr_values.append(a1_er)
            print('Fitted_sig_a1Noerr=', format(a1_er, ".6f"))
            return a1_er

        # %%
        def Fitted_a2_getAsymNoErr(ctr1,ctr2):
            a2_er_array = np.divide(np.subtract(ctr1,ctr2),np.add(ctr1,ctr2))
            a2_er = (np.sum(a2_er_array))/2

            fitted_a2Noerr_values.append(a2_er)
            print('Fitted_ctr_a2Noerr=', format(a2_er, ".6f"))
            return a2_er

        # %%
        print("Fitted_a1_values_NoErr")
        Fitted_a1_NoErr = Fitted_a1_getAsymNoErr(nsig,nasig)
        print(Fitted_a1_NoErr)

        # %%
        print("Fitted_a2_values_NoErr")
        Fitted_a2_NoErr = Fitted_a2_getAsymNoErr(nctr,nactr)
        print(Fitted_a2_NoErr)

        # %%
        # Write value to the text file
        file_1.write(str(Fitted_a1_NoErr) + "\n")
        file_2.write(str(Fitted_a2_NoErr) + "\n")

        # %%
        print("End_____loop", i)

# %%
print("******LOOP____END******")
print("******Final Acp values From Each Loop ******")

print("Fitted_a1_values : ",fitted_a1Noerr_values)
print("Fitted_a2_values : ",fitted_a2Noerr_values)

print("******Thank You******")
print("******Have a blessed day******")
print("Special Note: In the output Root file, Real num of Entries = Entries/num of samples(b)")
