import basf2 as b2
import modularAnalysis as ma
import variables as va
import variables.collections as vc
import variables.utils as vu
import stdCharged as stdc
from stdPi0s import stdPi0s
import vertex as vx
import stdPhotons
import sys

b2.conditions.append_globaltag(ma.getAnalysisGlobaltag())

mypath = b2.create_path()

# load input ROOT file
ma.inputMdst(environmentType='default', filename='/belle2work/janaka/CharmAnalysis/Recons_test_delete/sub00/mdst_000001_prod00027826_task251967000001.root', path=mypath)

# apply loose track quality selection
trackQuality = 'thetaInCDCAcceptance and nCDCHits>0'
goodTrack = trackQuality
ma.fillParticleList('p+:good', 'protonID > 0.2 and ' + goodTrack, True, path=mypath)
ma.fillParticleList('pi+:good', goodTrack, True, path=mypath)
ma.fillParticleList('K+:good', goodTrack, True, path=mypath)

# reconstruct pi0s
stdPi0s(listtype='eff40_May2020', path=mypath, beamBackgroundMVAWeight='MC15rd', fakePhotonMVAWeight='MC15rd', biasCorrectionTable='')

# build Sigma+ candidates
ma.reconstructDecay(decayString='Sigma+:ppi0loose -> p+:good pi0:eff40_May2020', cut='1.159 < M < 1.219', path=mypath)

# save the Sigma+ mass before the vertex fit
ma.variablesToExtraInfo("Sigma+:ppi0loose", variables={'M': 'M_before_fit'}, path=mypath)
va.variables.addAlias("sigma_M_BeforeFit", "extraInfo(M_before_fit)")


# build candidates and apply vertex fit with Sigma+ mass constraint
xic_cuts = '2.35 < M < 2.6 and  useCMSFrame(p)>=2.0'
ma.reconstructDecay(decayString='Xi_c+:sigpipi -> Sigma+:ppi0loose pi+:good pi-:good', cut=xic_cuts, path=mypath)
ma.reconstructDecay(decayString='Xi_c+:sigkk -> Sigma+:ppi0loose K+:good K-:good', cut=xic_cuts, path=mypath)
vx.treeFit('Xi_c+:sigpipi', conf_level=0.001, ipConstraint=True, updateAllDaughters=True, massConstraint=[3222], path=mypath)
vx.treeFit('Xi_c+:sigkk', conf_level=0.001, ipConstraint=True, updateAllDaughters=True, massConstraint=[3222], path=mypath)



lc_cuts = '2.15 < M < 2.4 and  useCMSFrame(p)>=2.0'
ma.reconstructDecay(decayString='Lambda_c+:sigpipi -> Sigma+:ppi0loose pi+:good pi-:good', cut=lc_cuts, path=mypath)
ma.reconstructDecay(decayString='Lambda_c+:sigkk -> Sigma+:ppi0loose K+:good K-:good', cut=lc_cuts, path=mypath)
ma.reconstructDecay(decayString='Lambda_c+:pkpi -> p+:good K-:good pi+:good', cut=lc_cuts, path=mypath)
ma.reconstructDecay(decayString='Lambda_c+:ppipi -> p+:good pi+:good pi-:good', cut=lc_cuts, path=mypath)
ma.reconstructDecay(decayString='Lambda_c+:pkk -> p+:good K+:good K-:good', cut=lc_cuts, path=mypath)
vx.treeFit('Lambda_c+:sigpipi', conf_level=0.001, ipConstraint=True, updateAllDaughters=True, massConstraint=[3222], path=mypath)
vx.treeFit('Lambda_c+:sigkk', conf_level=0.001, ipConstraint=True, updateAllDaughters=True, massConstraint=[3222], path=mypath)
vx.treeFit('Lambda_c+:pkpi', conf_level=0.001, ipConstraint=True, updateAllDaughters=True, path=mypath)
vx.treeFit('Lambda_c+:ppipi', conf_level=0.001, ipConstraint=True, updateAllDaughters=True, path=mypath)
vx.treeFit('Lambda_c+:pkk', conf_level=0.001, ipConstraint=True, updateAllDaughters=True, path=mypath)

# rank candidates
ma.rankByHighest(particleList='Xi_c+:sigpipi',variable='chiProb',outputVariable='vert_chi2_rank',allowMultiRank=True, path=mypath)
ma.rankByHighest(particleList='Xi_c+:sigkk',variable='chiProb',outputVariable='vert_chi2_rank',allowMultiRank=True, path=mypath)
ma.rankByHighest(particleList='Lambda_c+:sigpipi',variable='chiProb',outputVariable='vert_chi2_rank',allowMultiRank=True, path=mypath)
ma.rankByHighest(particleList='Lambda_c+:sigkk',variable='chiProb',outputVariable='vert_chi2_rank',allowMultiRank=True, path=mypath)
ma.rankByHighest(particleList='Lambda_c+:pkpi',variable='chiProb',outputVariable='vert_chi2_rank',allowMultiRank=True, path=mypath)
ma.rankByHighest(particleList='Lambda_c+:ppipi',variable='chiProb',outputVariable='vert_chi2_rank',allowMultiRank=True, path=mypath)
ma.rankByHighest(particleList='Lambda_c+:pkk',variable='chiProb',outputVariable='vert_chi2_rank',allowMultiRank=True, path=mypath)
va.variables.addAlias("vertex_chi2rank","extraInfo(vert_chi2_rank)")

va.variables.addAlias('vtxNDF','extraInfo(ndf)')
va.variables.addAlias('vtxChi2','extraInfo(chiSquared)')

# do MC matching
ma.matchMCTruth(list_name='Xi_c+:sigpipi', path=mypath)
ma.matchMCTruth(list_name='Xi_c+:sigkk', path=mypath)
ma.matchMCTruth(list_name='Lambda_c+:sigpipi', path=mypath)
ma.matchMCTruth(list_name='Lambda_c+:sigkk', path=mypath)
ma.matchMCTruth(list_name='Lambda_c+:pkpi', path=mypath)
ma.matchMCTruth(list_name='Lambda_c+:ppipi', path=mypath)
ma.matchMCTruth(list_name='Lambda_c+:pkk', path=mypath)

# define useful variables
va.variables.addAlias('p_cms', 'useCMSFrame(p)')
va.variables.addAlias('cosTheta_cms', 'useCMSFrame(cosTheta)')
#va.variables.addAlias('labAngleBetweeengammas', 'useLabFrame(daughterAngle(0,1))')

# photon variables
va.variables.addAlias('gamma0_cpsd', 'daughter(0,clusterPulseShapeDiscriminationMVA)')
va.variables.addAlias('gamma0_bbs', 'daughter(0,beamBackgroundSuppression)')
va.variables.addAlias('gamma0_fps', 'daughter(0,fakePhotonSuppression)')
va.variables.addAlias('gamma0_mcErrors', 'daughter(0,mcErrors)')
va.variables.addAlias('gamma0_genmother_0', 'daughter(0,genMotherPDG(0))')
va.variables.addAlias('gamma0_e9e21', 'daughter(0,clusterE9E21)')
va.variables.addAlias('gamma0_mcPDG', 'daughter(0,mcPDG)')
#va.variables.addAlias('gamma0_mcPrimary', 'daughter(0,mcPrimary)')

va.variables.addAlias('gamma1_cpsd', 'daughter(1,clusterPulseShapeDiscriminationMVA)')
va.variables.addAlias('gamma1_bbs', 'daughter(1,beamBackgroundSuppression)')
va.variables.addAlias('gamma1_fps', 'daughter(1,fakePhotonSuppression)')
va.variables.addAlias('gamma1_mcErrors', 'daughter(1,mcErrors)')
va.variables.addAlias('gamma1_genmother_0', 'daughter(1,genMotherPDG(0))')
va.variables.addAlias('gamma1_e9e21', 'daughter(1,clusterE9E21)')
va.variables.addAlias('gamma1_mcPDG', 'daughter(1,mcPDG)')
#va.variables.addAlias('gamma1_mcPrimary', 'daughter(1,mcPrimary)')

gamma_var = ['gamma0_cpsd', 'gamma0_bbs', 'gamma0_fps', 'gamma0_e9e21', 
             'gamma1_cpsd', 'gamma1_bbs', 'gamma1_fps', 'gamma1_e9e21']

gamma_var_mc = ['gamma0_mcErrors', 'gamma0_genmother_0', 'gamma0_mcPDG',  
                'gamma1_mcErrors', 'gamma1_genmother_0', 'gamma1_mcPDG']

xicvars = ['M', 'cosTheta', 'isSignal', 'vertex_chi2rank', 'vtxChi2','charge','chiProb','PDG']

# event level variables
#eventWiseVariables = ['nTracks','beamE','beamPx','beamPy','beamPz']
#eventWiseVariables += ['IPX','IPY','IPZ']
#eventWiseVariables += ['eventRandom']
eventWiseVariables = ['eventRandom']

# common variables
commonVariables = vc.kinematics
commonVariables += vc.mc_variables + vc.mc_truth

J_commonVariables = vc.kinematics

# track variables
trackVariables = vc.track + vc.track_hits + vc.pid + ['protonID','pionID','kaonID'] 
trackVariables =  ['protonID','pionID','kaonID'] 
trackVariables += ['charge']
trackVariables += ['cosTheta']
J_trackVariables = vc.track_hits
J_trackVariables += ['dr','dz','protonID','pionID','kaonID']

# composite variables
compositeVariables = vc.inv_mass
compositeVariables += ['distance']
compositeVariables += ['prodVertexX', 'prodVertexY', 'prodVertexZ', 'significanceOfDistance']
compositeVariables += ['cosAngleBetweenMomentumAndVertexVector']
compositeVariables += ['mcFlightTime', 'mcFlightDistance','flightDistance']

# cluster variables                                                                                                                                                                               
clusterVariables = vc.cluster

# additional variables
#va.variables.addAlias('px_cms','useCMSFrame(px)')
#va.variables.addAlias('py_cms','useCMSFrame(py)')
#va.variables.addAlias('pz_cms','useCMSFrame(pz)')
#cms_variables = ['px_cms','py_cms','pz_cms','p_cms','cosTheta_cms']
cms_variables = ['p_cms','cosTheta_cms']

# define the variables we want to include for the final state particles

xicp_sigpipi_vars = vu.create_aliases_for_selected(list_of_variables=J_commonVariables,
                                           decay_string='^Xi_c+:sigpipi -> [^Sigma+:ppi0loose -> ^p+ ^pi0] ^pi- ^pi+') + \
    vu.create_aliases_for_selected(list_of_variables=J_trackVariables,
                                           decay_string='^Xi_c+:sigpipi -> [^Sigma+:ppi0loose -> ^p+ pi0] ^pi- ^pi+') + \
    vu.create_aliases_for_selected(list_of_variables=gamma_var + cms_variables,
                                           decay_string='Xi_c+:sigpipi -> [Sigma+:ppi0loose -> p+ ^pi0] pi- pi+') + \
    vu.create_aliases_for_selected(list_of_variables=['sigma_M_BeforeFit'],
                                           decay_string='Xi_c+:sigpipi -> [^Sigma+:ppi0loose -> p+ pi0] pi- pi+') + \
    vu.create_aliases_for_selected(list_of_variables=cms_variables + xicvars,
                                           decay_string='^Xi_c+:sigpipi -> [^Sigma+:ppi0loose -> ^p+ ^pi0] ^pi- ^pi+') +\
    vu.create_aliases_for_selected(list_of_variables=compositeVariables,
                                           decay_string='^Xi_c+:sigpipi -> [^Sigma+:ppi0loose -> ^p+ pi0] pi- pi+')


xicp_sigkk_vars = vu.create_aliases_for_selected(list_of_variables=J_commonVariables,
                                           decay_string='^Xi_c+:sigkk -> [^Sigma+:ppi0loose -> ^p+ ^pi0] ^K- ^K+') + \
   vu.create_aliases_for_selected(list_of_variables=J_trackVariables,
                                           decay_string='^Xi_c+:sigkk -> [^Sigma+:ppi0loose -> ^p+ pi0] ^K- ^K+') + \
    vu.create_aliases_for_selected(list_of_variables=gamma_var + cms_variables,
                                           decay_string='Xi_c+:sigkk -> [Sigma+:ppi0loose -> p+ ^pi0] K- K+') + \
    vu.create_aliases_for_selected(list_of_variables=['sigma_M_BeforeFit'],
                                           decay_string='Xi_c+:sigkk -> [^Sigma+:ppi0loose -> p+ pi0] K- K+') + \
    vu.create_aliases_for_selected(list_of_variables=cms_variables + xicvars,
                                           decay_string='^Xi_c+:sigkk -> [^Sigma+:ppi0loose -> ^p+ ^pi0] ^K- ^K+') + \
    vu.create_aliases_for_selected(list_of_variables=compositeVariables,
                                           decay_string='^Xi_c+:sigkk -> [^Sigma+:ppi0loose -> ^p+ pi0] K- K+')


lcp_sigpipi_vars = vu.create_aliases_for_selected(list_of_variables=J_commonVariables,
                                           decay_string='^Lambda_c+:sigpipi -> [^Sigma+:ppi0loose -> ^p+ ^pi0] ^pi- ^pi+') + \
    vu.create_aliases_for_selected(list_of_variables=J_trackVariables,
                                           decay_string='^Lambda_c+:sigpipi -> [^Sigma+:ppi0loose -> ^p+ pi0] ^pi- ^pi+') + \
    vu.create_aliases_for_selected(list_of_variables=gamma_var + cms_variables,
                                           decay_string='Lambda_c+:sigpipi -> [Sigma+:ppi0loose -> p+ ^pi0] pi- pi+') + \
    vu.create_aliases_for_selected(list_of_variables=['sigma_M_BeforeFit'],
                                           decay_string='Lambda_c+:sigpipi -> [^Sigma+:ppi0loose -> p+ pi0] pi- pi+') + \
    vu.create_aliases_for_selected(list_of_variables=cms_variables + xicvars,
                                           decay_string='^Lambda_c+:sigpipi -> [^Sigma+:ppi0loose -> ^p+ ^pi0] ^pi- ^pi+') +\
    vu.create_aliases_for_selected(list_of_variables=compositeVariables,
                                           decay_string='^Lambda_c+:sigpipi -> [^Sigma+:ppi0loose -> ^p+ pi0] pi- pi+')

lcp_sigkk_vars = vu.create_aliases_for_selected(list_of_variables=J_commonVariables,
                                           decay_string='^Lambda_c+:sigkk -> [^Sigma+:ppi0loose -> ^p+ ^pi0] ^K- ^K+') + \
    vu.create_aliases_for_selected(list_of_variables=J_trackVariables,
                                           decay_string='^Lambda_c+:sigkk -> [^Sigma+:ppi0loose -> ^p+ pi0] ^K- ^K+') + \
    vu.create_aliases_for_selected(list_of_variables=gamma_var + cms_variables,
                                           decay_string='Lambda_c+:sigkk -> [Sigma+:ppi0loose -> p+ ^pi0] K- K+') + \
    vu.create_aliases_for_selected(list_of_variables=['sigma_M_BeforeFit'],
                                           decay_string='Lambda_c+:sigkk -> [^Sigma+:ppi0loose -> p+ pi0] K- K+') + \
    vu.create_aliases_for_selected(list_of_variables=cms_variables + xicvars,
                                           decay_string='^Lambda_c+:sigkk -> [^Sigma+:ppi0loose -> ^p+ ^pi0] ^K- ^K+') + \
    vu.create_aliases_for_selected(list_of_variables=compositeVariables,
                                           decay_string='^Lambda_c+:sigkk -> [^Sigma+:ppi0loose -> ^p+ pi0] K- K+') 

lcp_pkpi_vars = vu.create_aliases_for_selected(list_of_variables=J_commonVariables,
                                           decay_string='^Lambda_c+:pkpi -> ^p+ ^K- ^pi+',
                                           prefix=['Lambdac','p','K','pi']) + \
    vu.create_aliases_for_selected(list_of_variables=J_trackVariables,
                                           decay_string='^Lambda_c+:pkpi -> ^p+ ^K- ^pi+',
                                           prefix=['Lambdac','p','K','pi']) + \
    vu.create_aliases_for_selected(list_of_variables=compositeVariables + cms_variables + xicvars,
                                           decay_string='^Lambda_c+:pkpi -> ^p+ ^K- ^pi+',
                                           prefix=['Lambdac','p','K','pi'])


lcp_ppipi_vars = vu.create_aliases_for_selected(list_of_variables=J_commonVariables,
                                           decay_string='^Lambda_c+:ppipi -> ^p+ ^pi- ^pi+',
                                           prefix=['Lambdac','p','pi1','pi2']) + \
    vu.create_aliases_for_selected(list_of_variables=J_trackVariables,
                                           decay_string='^Lambda_c+:ppipi -> ^p+ ^pi- ^pi+',
                                           prefix=['Lambdac','p','pi1','pi2']) + \
    vu.create_aliases_for_selected(list_of_variables=compositeVariables + cms_variables + xicvars,
                                           decay_string='^Lambda_c+:ppipi -> ^p+ ^pi- ^pi+',
                                           prefix=['Lambdac','p','pi1','pi2'])

lcp_pkk_vars = vu.create_aliases_for_selected(list_of_variables=J_commonVariables,
                                           decay_string='^Lambda_c+:pkk -> ^p+ ^K- ^K+',
                                           prefix=['Lambdac','p','K1','K2']) + \
    vu.create_aliases_for_selected(list_of_variables=J_trackVariables,
                                           decay_string='^Lambda_c+:pkk -> ^p+ ^K- ^K+',
                                           prefix=['Lambdac','p','K1','K2']) + \
    vu.create_aliases_for_selected(list_of_variables=compositeVariables + cms_variables + xicvars,
                                           decay_string='^Lambda_c+:pkpi -> ^p+ ^K- ^K+',
                                           prefix=['Lambdac','p','K1','K2'])


# write out the ntuple file
ma.variablesToNtuple(decayString='Xi_c+:sigpipi', variables=xicp_sigpipi_vars+eventWiseVariables, filename="ntuple.root",
                     treename='xicp_sigpipi', path=mypath)
ma.variablesToNtuple(decayString='Xi_c+:sigkk', variables=xicp_sigkk_vars+eventWiseVariables, filename="ntuple.root",
                     treename='xicp_sigkk', path=mypath)


ma.variablesToNtuple(decayString='Lambda_c+:sigpipi', variables=lcp_sigpipi_vars+eventWiseVariables, filename="ntuple.root",
                     treename='lcp_sigpipi', path=mypath)
ma.variablesToNtuple(decayString='Lambda_c+:sigkk', variables=lcp_sigkk_vars+eventWiseVariables, filename="ntuple.root",
                     treename='lcp_sigkk', path=mypath)
ma.variablesToNtuple(decayString='Lambda_c+:pkpi', variables=lcp_pkpi_vars+eventWiseVariables, filename="ntuple.root",
                     treename='lcp_pkpi', path=mypath)
ma.variablesToNtuple(decayString='Lambda_c+:ppipi', variables=lcp_ppipi_vars+eventWiseVariables, filename="ntuple.root",
                     treename='lcp_ppipi', path=mypath)
ma.variablesToNtuple(decayString='Lambda_c+:pkk', variables=lcp_pkk_vars+eventWiseVariables, filename="ntuple.root",
                     treename='lcp_pkk', path=mypath)

# Process the events
b2.process(mypath)

# print out the summary
print(b2.statistics)



