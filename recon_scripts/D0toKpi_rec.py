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
output_file = "out_CFT.root"
ma.inputMdstList(environmentType='default', filelist=input_file, path=my_path)

ma.fillParticleList(
    "pi+:D0", "thetaInCDCAcceptance and dr < 1 and abs(dz) < 3", path=my_path
)
ma.fillParticleList(
    "K+:D0", "thetaInCDCAcceptance and dr < 1 and abs(dz) < 3", path=my_path
)

kinvars = ['E', 'p', 'cosTheta','pt']
boostvars = ['beamE','beamPx','beamPy','beamPz']
swapmass = vu.create_aliases('M',"useAlternativeDaughterHypothesis(M, 0:pi-, 1:K+)","swap")

ma.reconstructDecay(
    decayString="D0:Kpi -> K-:D0 pi+:D0",
    cut="1.66 < M < 2.06",
    path=my_path,
)

ma.matchMCTruth(list_name="D0:Kpi", path=my_path)

#build the rest of the event associated to the D0
#ma.buildRestOfEvent(target_list_name='D0:Kpi',
#                    path=my_path)

#Charm Flavor Tagging Function
#cft.charmFlavorTagger(
#    'D0:Kpi',
#    path=my_path)

vx.treeFit('D0:Kpi', conf_level=0.001, updateAllDaughters=True, ipConstraint=True, path=my_path)
ma.cutAndCopyList('D0:K+pi-', 'D0:Kpi', '1.75 < M < 1.95', writeOut=False, path=my_path)

va.variables.addAlias("genGrandmotherPDG", "genMotherPDG(1)")
vertexvars = vc.flight_info
vertexvars += ["M", "chiProb"]

truthvars = vc.mc_truth + vc.mc_variables + ["genGrandmotherPDG", "M"]
vertextruthvars = vc.mc_flight_info
vertexvars += truthvars + vertextruthvars

trackvars = ['d0', 'z0', 'phi0', 'omega', 'charge', 'ndf', 'pValue','nPXDHits', 'nSVDHits', 'nCDCHits']
pidvars = ['pionID', 'kaonID', 'pionID_noSVD', 'kaonID_noSVD', 'pionIDNN', 'kaonIDNN']

#CFT = ["CFT_qr"]

d0vars = vu.create_aliases_for_selected(
    list_of_variables = vertexvars + kinvars + boostvars + swapmass, 
#+ CFT,
    decay_string="^D0 -> K- pi+",
)

fsvars = vu.create_aliases_for_selected(
    list_of_variables = truthvars + kinvars + trackvars + pidvars, 
    decay_string="D0 -> ^K- ^pi+",
)

ma.variablesToNtuple(
    "D0:K+pi-",
    d0vars + fsvars,
    filename=output_file,
    treename="D0tree",
    path=my_path,
)

b2.process(my_path)
