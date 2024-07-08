# %%
import zfit
from zfit import z
import root_pandas
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
dcol = ["M","reWeights","CMS_cosTheta","K_charge"]
dmccol = ["isSignal"]

dinfile = '/b2diska/dinura/rootfiles/ntuples/MC/dk3pi/weights/k3pi_ppipi_weighted.root'

ddf = root_pandas.read_root(dinfile, key='D0tree', columns=dcol+dmccol, where='1.80 < M && M < 1.92')


# %%
dupper = 1.92
dlower = 1.80


# %%
d_cut = str(dlower)+ ' < M < ' +str(dupper)


# %%
binedges=[-1,0,1]
dbin1 = str(binedges[0])+"<=CMS_cosTheta<"+str(binedges[1])
dbin2 = str(binedges[1])+"<=CMS_cosTheta<"+str(binedges[2])
print(dbin1)
print(dbin2)


# %%


# D -----> Kpi Method to overlay fitted pdfs and data sample
def Dplot_model(model, mydata, nevents, nbins=100, myrange=(dlower,dupper), mylabel="", plot_data=True, save_dir=None, save_name="plot"):

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
# # Apply weights for $D^0$

# find the momentum and cos(theta) bin for each event
ddf = ddf.query(d_cut)
ddf_bin1 = ddf.query(d_cut+" and "+dbin1+' and K_charge==-1')
ddf_bin2 = ddf.query(d_cut+" and "+dbin2+' and K_charge==-1')
addf_bin1 = ddf.query(d_cut+" and "+dbin1+' and K_charge==1')
addf_bin2 = ddf.query(d_cut+" and "+dbin2+' and K_charge==1')


# %%
# this is a little dangerous since we replace keys with indices that depend on the order 
def getWeights(df):
    np_arr = df.to_numpy()
    return np_arr[:,0], np_arr[:,1]


# %%
d_M, d_weights = getWeights(ddf)
d_bin1_M, d_bin1_weights = getWeights(ddf_bin1)
d_bin2_M, d_bin2_weights = getWeights(ddf_bin2)
ad_bin1_M, ad_bin1_weights = getWeights(addf_bin1)
ad_bin2_M, ad_bin2_weights = getWeights(addf_bin2)

# %%
# # ZFit for $D^0$

# D ----> Kpi Define default parameters for zfit
drange = (dlower,dupper)
dobs = zfit.Space('M', drange)
dissignal = 'isSignal==1'

# Get the signal and background for reference
dsignal = ddf.query(d_cut + ' and isSignal==1').M.to_numpy()
dbkg = ddf.query(d_cut + ' and isSignal!=1').M.to_numpy()
print('signal: ' + str(len(dsignal)))
print('bkg: ' + str(len(dbkg)))


# %%
# D----> Kpi Signal fit parameters
dmu = zfit.Parameter("dmu", 1.864, 1.85, 1.88)
ds1 = zfit.param.Parameter("ds1", 0.005, 0.0001, 0.015)
ds2 = zfit.param.Parameter("ds2", 0.003, 0.0001, 0.015)
dfg1 = zfit.param.Parameter("dfg1", 0.2, 0., 1.)
da1 = zfit.param.Parameter("da1", 0.01, -1e6, 1e6)
da2 = zfit.param.Parameter("da2", 0.01, -1e6, 1e6)

# define PDFs
dgaus1 = zfit.pdf.Gauss(obs=dobs, mu=dmu, sigma=ds1)
dgaus2 = zfit.pdf.Gauss(obs=dobs, mu=dmu, sigma=ds2)
dgaus = zfit.pdf.SumPDF(pdfs=[dgaus1,dgaus2], fracs=[dfg1])

dpoly = zfit.pdf.Chebyshev(obs=dobs, coeffs=[da1, da2])

dsig_yield = zfit.Parameter('dsig_yield', 6000000, 0, 1e8, step_size=1)
dbkg_yield = zfit.Parameter('dbkg_yield', 25000000, 0, 1e9, step_size=1)

dsig_ext = dgaus.create_extended(dsig_yield)
dbkg_ext = dpoly.create_extended(dbkg_yield)

dpdf_ext = zfit.pdf.SumPDF(pdfs=[dsig_ext,dbkg_ext])

# %%
# # Intergrated Fit

# Fill an array with the data to be fit
data_np = ddf.query(d_cut).M.to_numpy()
data = zfit.data.Data.from_numpy(obs=dobs, array=data_np, weights = d_weights)
data.set_data_range(drange) 
data_size = data.n_events


# %%
# Define loss function
nll_simultaneous = zfit.loss.ExtendedUnbinnedNLL(model=dpdf_ext,data=data)
minimizer = zfit.minimize.Minuit()
result = minimizer.minimize(nll_simultaneous)  
print(result)


# %%
Dplot_model(model=dpdf_ext, mydata=data, myrange=(1.80,1.92), nevents=int(data_size), mylabel=r'$D^0$ Mass ($D^{0} \to K^- \pi^+ \pi^+ \pi^-$)  [GeV/$c^{2}$]', plot_data=True, save_dir="/b2diska/dinura/scripts/acp/py/dk3pi/plots/weighted/", save_name="dk3pi_ppipi_intergrated")


# %%
# # Fits in cos($\theta$) bins for $D^0\rightarrow K^- \pi^- \pi^+ \pi^+$ 

ndtrue=[]
nadtrue=[]
nd=[]
nad=[]

ndtruee=[]
nadtruee=[]
nde=[]
nade=[]

ndweighted=[]
nadweighted=[]
ndweighted_un=[]
nadweighted_un=[]


# %%
# signal fit parameters
d_fg1 = zfit.Parameter("d_fg1", 0.20, 0., 1.)
d_fg1.floating = True
ad_fg1 = zfit.Parameter("ad_fg1", 0.20, 0., 1.)
ad_fg1.floating = True

d_a1 = zfit.Parameter("d_a1", 0.01, -1e6, 1e6)
ad_a1 = zfit.Parameter("ad_a1", 0.01, -1e6, 1e6)

d_a2 = zfit.Parameter("d_a2", 0.01, -1e6, 1e6)
ad_a2 = zfit.Parameter("ad_a2", 0.01, -1e6, 1e6)

# define PDFs
d_gaus1 = zfit.pdf.Gauss(obs=dobs, mu=dmu, sigma=ds1)
d_gaus2 = zfit.pdf.Gauss(obs=dobs, mu=dmu, sigma=ds2)
d_gaus = zfit.pdf.SumPDF(pdfs=[d_gaus1,d_gaus2], fracs=[d_fg1])

d_poly = zfit.pdf.Chebyshev(obs=dobs, coeffs=[d_a1, d_a2])

d_sig_yield = zfit.Parameter('d_sig_yield', 1500000, 0, 0.5e7, step_size=1)
d_bkg_yield = zfit.Parameter('d_bkg_yield', 4000000, 0, 0.5e7, step_size=1)

d_sig_ext = d_gaus.create_extended(d_sig_yield)
d_bkg_ext = d_poly.create_extended(d_bkg_yield)

d_pdf_ext = zfit.pdf.SumPDF(pdfs=[d_sig_ext,d_bkg_ext])

ad_gaus1 = zfit.pdf.Gauss(obs=dobs, mu=dmu, sigma=ds1)
ad_gaus2 = zfit.pdf.Gauss(obs=dobs, mu=dmu, sigma=ds2)
ad_gaus = zfit.pdf.SumPDF(pdfs=[ad_gaus1,ad_gaus2], fracs=[ad_fg1])

ad_poly = zfit.pdf.Chebyshev(obs=dobs, coeffs=[ad_a1, ad_a2])

ad_sig_yield = zfit.Parameter('ad_sig_yield', 1500000, 0, 0.5e7, step_size=1)
ad_bkg_yield = zfit.Parameter('ad_bkg_yield', 4000000, 0, 0.5e7, step_size=1)

ad_sig_ext = ad_gaus.create_extended(ad_sig_yield)
ad_bkg_ext = ad_poly.create_extended(ad_bkg_yield)

ad_pdf_ext = zfit.pdf.SumPDF(pdfs=[ad_sig_ext,ad_bkg_ext])

# %%
# ## bin 1

print(dbin1)


# %%
# Get the signal and background for reference
d_signal = ddf.query(d_cut+" and "+dbin1+' and isSignal==1 and K_charge==-1').M.to_numpy()
d_bkg = ddf.query(d_cut+" and "+dbin1+' and isSignal!=1 and K_charge==-1').M.to_numpy()
print('d signal: ' + str(len(d_signal)))
print('d bkg: ' + str(len(d_bkg)))
ndtrue.append(len(d_signal))

ad_signal = ddf.query(d_cut+" and "+dbin1+' and isSignal==1 and K_charge==1').M.to_numpy()
ad_bkg = ddf.query(d_cut+" and "+dbin1+' and isSignal!=1 and K_charge==1').M.to_numpy()
print('ad signal: ' + str(len(ad_signal)))
print('ad bkg: ' + str(len(ad_bkg)))
nadtrue.append(len(ad_signal))


# %%
# Fill an array with the data to be fit
d_data_np = ddf.query(d_cut+" and "+dbin1+' and K_charge==-1').M.to_numpy()
d_data = zfit.data.Data.from_numpy(obs=dobs, array=d_data_np, weights= d_bin1_weights) #applyng weights 
d_data.set_data_range(drange)
d_data_size = d_data.n_events

ad_data_np = ddf.query(d_cut+" and "+dbin1+' and K_charge==1').M.to_numpy()
ad_data = zfit.data.Data.from_numpy(obs=dobs, array=ad_data_np, weights=ad_bin1_weights) #applying weights
ad_data.set_data_range(drange)
ad_data_size = ad_data.n_events

print(int(d_data_size))
print(int(ad_data_size))


# %%
# Simultaenous loss
nll_simultaneous = zfit.loss.ExtendedUnbinnedNLL(model=[d_pdf_ext,ad_pdf_ext], data=[d_data,ad_data])
minimizer = zfit.minimize.Minuit()
result = minimizer.minimize(nll_simultaneous)


# %%
# Simultaenous loss
nll_simultaneous = zfit.loss.ExtendedUnbinnedNLL(model=[d_pdf_ext,ad_pdf_ext], data=[d_data,ad_data])
minimizer = zfit.minimize.Minuit()
result = minimizer.minimize(nll_simultaneous)
print(result)

# %%
result.hesse()


# %%
errors_d, new_result_d = result.errors(params=d_pdf_ext.get_params())
print("bin1 matter completed")

# %%
errors_ad, new_result_ad = result.errors(params=ad_pdf_ext.get_params())
print("bin1 antimatter completed")

# %%
print(result)


# %%
Dplot_model(model=d_pdf_ext, mydata=d_data, nevents=int(d_data_size), mylabel=r'$D^0$ Mass ($D^0\rightarrow K^- \pi^+ \pi^+ \pi^-$)  [GeV/$c^{2}$]', plot_data=True, save_dir="/b2diska/dinura/scripts/acp/py/dk3pi/plots/weighted/", save_name="dk3pi_ppipi_matter_bin1")


# %%
Dplot_model(model=ad_pdf_ext, mydata=ad_data, nevents=int(ad_data_size), mylabel=r'$D^0$ Mass ($\barD^0\rightarrow K^+ \pi^- \pi^- \pi^+$)  [GeV/$c^{2}$]', plot_data=True, save_dir="/b2diska/dinura/scripts/acp/py/dk3pi/plots/weighted/", save_name="dk3pi_ppipi_antimatter_bin1")


# %%
nd.append(result.params[d_sig_yield].get('value'))
nde.append(result.params[d_sig_yield].get('hesse').get('error'))
nad.append(result.params[ad_sig_yield].get('value'))
nade.append(result.params[ad_sig_yield].get('hesse').get('error'))

# %%
# ## bin 2

print(dbin2)


# %%
# Get the signal and background for reference
d_signal = ddf.query(d_cut+" and "+dbin2+' and isSignal==1 and K_charge==-1').M.to_numpy()
d_bkg = ddf.query(d_cut+" and "+dbin2+' and isSignal!=1 and K_charge==-1').M.to_numpy()
print('d signal: ' + str(len(d_signal)))
print('d bkg: ' + str(len(d_bkg)))
ndtrue.append(len(d_signal))

ad_signal = ddf.query(d_cut+" and "+dbin2+' and isSignal==1 and K_charge==1').M.to_numpy()
ad_bkg = ddf.query(d_cut+" and "+dbin2+' and isSignal!=1 and K_charge==1').M.to_numpy()
print('ad signal: ' + str(len(ad_signal)))
print('ad bkg: ' + str(len(ad_bkg)))
nadtrue.append(len(ad_signal))


# %%
# Fill an array with the data to be fit
d_data_np = ddf.query(d_cut+" and "+dbin2+' and K_charge==-1').M.to_numpy()
d_data = zfit.data.Data.from_numpy(obs=dobs, array=d_data_np, weights= d_bin2_weights) #applying weights
d_data.set_data_range(drange)
d_data_size = d_data.n_events

ad_data_np = ddf.query(d_cut+" and "+dbin2+' and K_charge==1').M.to_numpy()
ad_data = zfit.data.Data.from_numpy(obs=dobs, array=ad_data_np, weights= ad_bin2_weights) #applying weights
ad_data.set_data_range(drange)
ad_data_size = ad_data.n_events

print(int(d_data_size))
print(int(ad_data_size))


# %%
# Simultaenous loss
nll_simultaneous = zfit.loss.ExtendedUnbinnedNLL(model=[d_pdf_ext,ad_pdf_ext], data=[d_data,ad_data])
minimizer = zfit.minimize.Minuit()
result = minimizer.minimize(nll_simultaneous)


# %%
# Simultaenous loss
nll_simultaneous = zfit.loss.ExtendedUnbinnedNLL(model=[d_pdf_ext,ad_pdf_ext], data=[d_data,ad_data])
minimizer = zfit.minimize.Minuit()
result = minimizer.minimize(nll_simultaneous)
print(result)

# %%
result.hesse()


# %%
errors_d, new_result_d = result.errors(params=d_pdf_ext.get_params())
print("bin2 matter completed")

# %%
errors_ad, new_result_ad = result.errors(params=ad_pdf_ext.get_params())
print("bin2 antimatter completed")

# %%
print(result)


# %%
Dplot_model(model=d_pdf_ext, mydata=d_data, nevents=int(d_data_size), mylabel=r'$D^0$ Mass ($D^0\rightarrow K^- \pi^+ \pi^+ \pi^-$)  [GeV/$c^{2}$]', plot_data=True, save_dir="/b2diska/dinura/scripts/acp/py/dk3pi/plots/weighted/", save_name="dk3pi_ppipi_matter_bin2")


# %%
Dplot_model(model=ad_pdf_ext, mydata=ad_data, nevents=int(ad_data_size), mylabel=r'$D^0$ Mass ($\barD^0\rightarrow K^+ \pi^- \pi^- \pi^+$)  [GeV/$c^{2}$]', plot_data=True, save_dir="/b2diska/dinura/scripts/acp/py/dk3pi/plots/weighted/", save_name="dk3pi_ppipi_antimatter_bin2")


# %%
nd.append(result.params[d_sig_yield].get('value'))
nde.append(result.params[d_sig_yield].get('hesse').get('error'))
nad.append(result.params[ad_sig_yield].get('value'))
nade.append(result.params[ad_sig_yield].get('hesse').get('error'))


# %%
dar = unumpy.uarray(nd,nde)
adar = unumpy.uarray(nad,nade)


# %%
ndtruee = np.sqrt(ndtrue)
nadtruee = np.sqrt(nadtrue)


# %%
trdar = unumpy.uarray(ndtrue,ndtruee)
tradar = unumpy.uarray(nadtrue,nadtruee)


# %%
# Get the weighted true Acp value
d_signal = (ddf.query(d_cut + " and " + dbin1 + ' and isSignal==1 and K_charge==-1')['reWeights']).to_numpy()
d_signal_2 = d_signal*d_signal
d_bkg = (ddf.query(d_cut + " and " + dbin1 + ' and isSignal!=1 and K_charge==-1')['reWeights']).to_numpy()
ndweighted.append(np.sum(d_signal))
ndweighted_un.append(np.sum(d_signal_2))
print('weighted d signal bin1: ' + str(np.sum(d_signal)))
print('weighted d bkg bin1: ' + str(np.sum(d_bkg)))


ad_signal = (ddf.query(d_cut + " and " + dbin1 + ' and isSignal==1 and K_charge==1')['reWeights']).to_numpy()
ad_signal_2 = ad_signal*ad_signal
ad_bkg = (ddf.query(d_cut + " and " + dbin1 + ' and isSignal!=1 and K_charge==1')['reWeights']).to_numpy()
nadweighted.append(np.sum(ad_signal))
nadweighted_un.append(np.sum(ad_signal_2))
print('weighted ad signal bin1: ' + str(np.sum(ad_signal)))
print('weighted ad bkg bin1: ' + str(np.sum(ad_bkg)))


# %%
# Get the weighted true Acp value
d_signal = (ddf.query(d_cut + " and " + dbin2 + ' and isSignal==1 and K_charge==-1')['reWeights']).to_numpy()
d_signal_2 = d_signal*d_signal
d_bkg = (ddf.query(d_cut + " and " + dbin2 + ' and isSignal!=1 and K_charge==-1')['reWeights']).to_numpy()
ndweighted.append(np.sum(d_signal))
ndweighted_un.append(np.sum(d_signal_2))
print('weighted d signal bin2: ' + str(np.sum(d_signal)))
print('weighted d bkg bin12 ' + str(np.sum(d_bkg)))


ad_signal = (ddf.query(d_cut + " and " + dbin2 + ' and isSignal==1 and K_charge==1')['reWeights']).to_numpy()
ad_signal_2 = ad_signal*ad_signal
ad_bkg = (ddf.query(d_cut + " and " + dbin2 + ' and isSignal!=1 and K_charge==1')['reWeights']).to_numpy()
nadweighted.append(np.sum(ad_signal))
nadweighted_un.append(np.sum(ad_signal_2))
print('weighted ad signal bin2: ' + str(np.sum(ad_signal)))
print('weighted ad bkg bin2: ' + str(np.sum(ad_bkg)))


# %%
ndweightede = np.sqrt(ndweighted_un)
nadweightede = np.sqrt(nadweighted_un)
    
trdwar = unumpy.uarray(ndweighted,ndweightede)
tradwar = unumpy.uarray(nadweighted,nadweightede)

print("More decimal places")

def getAsym(d,ad):

    a3 = ((d[0]-ad[0])/(d[0]+ad[0])+(d[1]-ad[1])/(d[1]+ad[1]))/2

    return '{:.10f}'.format(a3)


# %%
print("True Acp")
print(getAsym(trdar,tradar)) #True Araw


# %%
print("Weighted True Acp")
print(getAsym(trdwar,tradwar)) #Weighted True Araw


# %%
print("Fitted Acp")
print(getAsym(dar,adar)) #Weighted Araw

