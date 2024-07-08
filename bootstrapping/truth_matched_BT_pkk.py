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
def True_a1_getAsymNoErr(sig1,sig2):
    a1_er_array = np.divide(np.subtract(sig1,sig2),np.add(sig1,sig2))
    a1_er = (np.sum(a1_er_array))/2
    
    true_a1Noerr_values.append(a1_er)
    
    print('true_sig_a1Noerr=', format(a1_er, ".6f"))
    return a1_er


# %%
def True_a2_getAsymNoErr(ctr1,ctr2):
    a2_er_array = np.divide(np.subtract(ctr1,ctr2),np.add(ctr1,ctr2))
    a2_er = (np.sum(a2_er_array))/2

    true_a2Noerr_values.append(a2_er)

    print('true_ctr_a2Noerr=', format(a2_er, ".6f"))
    return a2_er


# %%
# Define columns to read into memory
col = ["Lambdac_M","Lambdac_cosTheta_cms","p_charge","p_p","p_cosTheta"]
mccol =["Lambdac_isSignal"]

siginfile = '/b2diska/dinura/rootfiles/ntuples/MC/lcp/pkk_s.root'
ctrinfile = '/b2diska/dinura/rootfiles/ntuples/MC/lcp/pkpi_pkk_s.root'

# %%
#bootstrap sample size in %
smplsize = 0.25
# how many bootstrap samples
b = 500

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
# Initialize empty arrays for a1, a2, acp
true_a1Noerr_values = []
true_a2Noerr_values = []

# %%
sigdf_main = root_pandas.read_root(siginfile, key='lcp_pkk', columns=col+mccol, where='Lambdac_M > 2.24 && Lambdac_M < 2.34')


# %%
ctrdf_main = root_pandas.read_root(ctrinfile, key='lcp_pkpi', columns=col+mccol, where= 'Lambdac_M > 2.24 && Lambdac_M < 2.34')

# %%
# Define default parameters for zfit
lcrange = (llower,lupper)
lcobs = zfit.Space('Lambdac_M', lcrange)
issignal = 'Lambdac_isSignal==1'

# Get the signal and background for reference
signal = sigdf_main.query(lc_cut + ' and Lambdac_isSignal==1').Lambdac_M.to_numpy()
bkg = sigdf_main.query(lc_cut + ' and Lambdac_isSignal!=1').Lambdac_M.to_numpy()
print('signal: ' + str(len(signal)))
print('bkg: ' + str(len(bkg)))


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

with open("pkk_true_a1.txt", "w") as file_1, open("pkk_true_a2.txt", "w") as file_2:

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
        np_pkk_M_sw = sigdf.query(lc_cut)['Lambdac_M'].to_numpy()    
        np_pkk_p_sw = sigdf.query(lc_cut)['p_p'].to_numpy()       
        np_pkk_cosTheta_sw = sigdf.query(lc_cut)['p_cosTheta'].to_numpy()    

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
        bkg_yield.set_value(500000)
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

        # %%
        pkkweights = compute_sweights(pdf_ext, data)

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
        
        sig_yield.set_value(700000)
        sig_yield.floating = True
        bkg_yield.set_value(2000000)
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
        pkpiscale = pkkweights[sig_yield].sum()/pkpiweights[sig_yield].sum()

        # %%
        # Reweighting for proton momentum
        pkk_counts_per_bin_p = np.histogram(np_pkk_p_sw, bins=nbins, weights=pkkweights[sig_yield], range=p_range)
        pkpi_counts_per_bin_p = np.histogram(np_pkpi_p_sw, bins=nbins, weights=[pkpiscale]*pkpiweights[sig_yield], range=p_range)

        # Calculate the event ratio for each bin
        pkk_counts_per_bin_p = pkk_counts_per_bin_p[0]  # Extract the counts from the tuple
        pkpi_counts_per_bin_p = pkpi_counts_per_bin_p[0]
        event_ratio_p = pkk_counts_per_bin_p / pkpi_counts_per_bin_p

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
        pkk_counts_per_bin_cosTheta = np.histogram(np_pkk_cosTheta_sw, bins=nbins, weights=pkkweights[sig_yield], range=cosTheta_range)
        pkpi_counts_per_bin_cosTheta = np.histogram(np_pkpi_cosTheta_sw, bins=nbins, weights=[pkpiscale]*pkpi_pweights_sw, range=cosTheta_range)

        # Calculate the event ratio for each bin
        pkk_counts_per_bin_cosTheta = pkk_counts_per_bin_cosTheta[0]  # Extract the counts from the tuple
        pkpi_counts_per_bin_cosTheta = pkpi_counts_per_bin_cosTheta[0]
        event_ratio_cosTheta = pkk_counts_per_bin_cosTheta / pkpi_counts_per_bin_cosTheta

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
        # Apply Weights
        # find the momentum and cos(theta) bin for each event
        pkpidf_bin1 = ctrdf.query(lc_cut+" and "+bin1+' and p_charge==1')
        pkpidf_bin2 = ctrdf.query(lc_cut+" and "+bin2+' and p_charge==1')
        apkpidf_bin1 = ctrdf.query(lc_cut+" and "+bin1+' and p_charge==-1')
        apkpidf_bin2 = ctrdf.query(lc_cut+" and "+bin2+' and p_charge==-1')

        # %%
        pkpidf_bin1

        # %%
        # this is a little dangerous since we replace keys with indices that depend on the order 
        def getWeights(df):
            np_arr = df.to_numpy()
            return np_arr[:,0], np_arr[:,6]

        # %%
        pkpi_bin1_M, pkpi_bin1_weights = getWeights(pkpidf_bin1)
        pkpi_bin2_M, pkpi_bin2_weights = getWeights(pkpidf_bin2)
        apkpi_bin1_M, apkpi_bin1_weights = getWeights(apkpidf_bin1)
        apkpi_bin2_M, apkpi_bin2_weights = getWeights(apkpidf_bin2)

        # %%
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
        print(bin1)

        # %%
        # Get the signal and background for reference
        sig_signal = sigdf.query(lc_cut+" and "+bin1+' and Lambdac_isSignal==1 and p_charge==1').Lambdac_M.to_numpy()
        print('sig signal: ' + str(len(sig_signal)))
        nsigtrue.append(len(sig_signal))

        asig_signal = sigdf.query(lc_cut+" and "+bin1+' and Lambdac_isSignal==1 and p_charge==-1').Lambdac_M.to_numpy()
        print('asig signal: ' + str(len(asig_signal)))
        nasigtrue.append(len(asig_signal))

        # Get the signal and background for reference
        ctr_signal = ctrdf.query(lc_cut+" and "+bin1+' and Lambdac_isSignal==1 and p_charge==1')['reWeights'].to_numpy()
        print('ctr signal: ' + str(np.sum(ctr_signal)))
        nctrtrue.append(np.sum(ctr_signal))

        actr_signal = ctrdf.query(lc_cut+" and "+bin1+' and Lambdac_isSignal==1 and p_charge==-1')['reWeights'].to_numpy()
        print('actr signal: ' + str(np.sum(actr_signal)))
        nactrtrue.append(np.sum(actr_signal))

        # %%
        print(bin2)

        # %%
        # Get the signal and background for reference
        sig_signal = sigdf.query(lc_cut+" and "+bin2+' and Lambdac_isSignal==1 and p_charge==1').Lambdac_M.to_numpy()
        print('sig signal: ' + str(len(sig_signal)))
        nsigtrue.append(len(sig_signal))

        asig_signal = sigdf.query(lc_cut+" and "+bin2+' and Lambdac_isSignal==1 and p_charge==-1').Lambdac_M.to_numpy()
        print('asig signal: ' + str(len(asig_signal)))
        nasigtrue.append(len(asig_signal))

        # Get the signal and background for reference
        ctr_signal = ctrdf.query(lc_cut+" and "+bin2+' and Lambdac_isSignal==1 and p_charge==1')['reWeights'].to_numpy()
        print('ctr signal: ' + str(np.sum(ctr_signal)))
        nctrtrue.append(np.sum(ctr_signal))

        actr_signal = ctrdf.query(lc_cut+" and "+bin2+' and Lambdac_isSignal==1 and p_charge==-1')['reWeights'].to_numpy()
        print('actr signal: ' + str(np.sum(actr_signal)))
        nactrtrue.append(np.sum(actr_signal))

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
        print("true_a1_values_NoErr")
        True_a1_NoErr = True_a1_getAsymNoErr(nsigtrue,nasigtrue)
        print(True_a1_NoErr)

        # %%
        print("true_a2_values_NoErr")
        True_a2_NoErr = True_a2_getAsymNoErr(nctrtrue,nactrtrue)
        print(True_a2_NoErr)

        # %%
        # Write value to the text file
        file_1.write(str(True_a1_NoErr) + "\n")
        file_2.write(str(True_a2_NoErr) + "\n")

        # %%
        print("End_____loop", i)

# %%
print("******LOOP____END******")
print("******Final Acp values From Each Loop ******")
print("True_a1_values : ",true_a1Noerr_values)
print("True_a2_values : ",true_a2Noerr_values)

print("******Thank You******")
print("******Have a blessed day******")
print("Special Note: In the output Root file, Real num of Entries = Entries/num of samples(b)")

# %%
