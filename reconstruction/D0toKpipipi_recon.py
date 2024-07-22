import basf2 as b2
import modularAnalysis as ma
import variables.collections as vc
import variables.utils as vu
import stdCharged as stdc
import variables as va
from stdV0s import stdKshorts
from stdPi0s import stdPi0s
from stdPhotons import stdPhotons
from datetime import *
from ROOT import Belle2
import charmFlavorTagger as cft
import vertex as vx
import os
import math
import sys

my_path = b2.create_path()
b2.set_log_level(b2.LogLevel.ERROR)

#append analysis global tag where the CFT payload is stored
b2.conditions.prepend_globaltag(ma.getAnalysisGlobaltag())

input_file = []
output_file = "k3pi_ntuple_MC.root"
ma.inputMdstList(environmentType='default', filelist=input_file, path=my_path)

ma.fillParticleList(
    "pi+:D0", "thetaInCDCAcceptance and nCDCHits>0 and dr < 1 and abs(dz) < 3", path=my_path
)
ma.fillParticleList(
    "K+:D0", "thetaInCDCAcceptance and nCDCHits>0 and dr < 1 and abs(dz) < 3", path=my_path
)

kinvars = ['E', 'p', 'cosTheta','theta','pt','px','py','pz']
cmsvars1 = vu.create_aliases(kinvars, "useCMSFrame({variable})", "CMS")
cmsvars = ['beamE','beamPx','beamPy','beamPz']
swapmass = vu.create_aliases('M',"useAlternativeDaughterHypothesis(M, 0:pi-, 1:K+)","swap")
va.variables.addAlias('kbinaryID', 'formula(kaonID/(pionID+kaonID))')
va.variables.addAlias('pibinaryID', 'formula(pionID/(pionID+kaonID))')
va.variables.addAlias('flightDistanceErrRatio', 'formula(flightDistance/flightDistanceErr)')

ma.cutAndCopyList('K+:D0_cut', 'K+:D0', 'kbinaryID > 0.2', path=my_path)
ma.cutAndCopyList('pi+:D0_cut', 'pi+:D0', 'pionIDNN > 0.1', path=my_path)

ma.reconstructDecay(
    decayString="D0:Kpipipi -> K-:D0_cut pi-:D0_cut pi+:D0_cut pi+:D0_cut",
    cut="1.70 < M < 2.00 and useCMSFrame(p)>=2.0",
    path=my_path,
)

vx.treeFit('D0:Kpipipi', conf_level=0.001, updateAllDaughters=True, ipConstraint=True, path=my_path)
ma.matchMCTruth(list_name="D0:Kpipipi", path=my_path)

va.variables.addAlias("genGrandmotherPDG", "genMotherPDG(1)")
vertexvars = vc.flight_info
vertexvars += ["M","eventRandom","chiProb"]
vertexvars += ["dr","dz","distance","significanceOfDistance"]
vertexvars += ["flightDistanceErr","flightDistanceErrRatio","flightTime","flightTimeErr"]
vertexvars += ["cosAngleBetweenMomentumAndVertexVector"]

truthvars = vc.mc_truth + vc.mc_variables + ["genGrandmotherPDG"]
vertextruthvars = vc.mc_flight_info
vertexvars += truthvars + vertextruthvars

trackvars = ['dr','dz','d0', 'z0', 'phi0', 'omega', 'charge', 'ndf', 'pValue','nPXDHits', 'nSVDHits', 'nCDCHits', 'M']
pidvars = ['pionID', 'kaonID', 'pionIDNN', 'kaonIDNN','kbinaryID','pibinaryID']

d0vars = vu.create_aliases_for_selected(
    list_of_variables = vertexvars + kinvars + cmsvars + cmsvars1 +swapmass, 
    decay_string="^D0 -> K- pi- pi+ pi+",
)

fsvars = vu.create_aliases_for_selected(
    list_of_variables = truthvars + kinvars + trackvars + cmsvars1 + pidvars, 
    decay_string="D0 -> ^K- ^pi- ^pi+ ^pi+ ",
)

ma.variablesToNtuple(
    "D0:Kpipipi",
    d0vars + fsvars,
    filename=output_file,
    treename="D0tree",
    path=my_path,
)

b2.process(my_path)
