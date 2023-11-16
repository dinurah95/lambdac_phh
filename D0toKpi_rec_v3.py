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

# append analysis global tag where the CFT payload is stored                                                                                                                                                       
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

kin_variables = ['E', 'p', 'cosTheta','px', 'py', 'pz', 'pt']
cms_kinematics = vu.create_aliases(kin_variables, "useCMSFrame({variable})", "CMS")

ma.reconstructDecay(
    decayString="D0:Kpi -> K-:D0 pi+:D0",
    cut="1.66 < M < 2.06 and CMS_p>2",
    path=my_path,
)

ma.matchMCTruth(list_name="D0:Kpi", path=my_path)

#build the rest of the event associated to the D0
ma.buildRestOfEvent(target_list_name='D0:Kpi',
                    path=my_path)


#Charm Flavor Tagging Function
cft.charmFlavorTagger(
    'D0:Kpi',
    path=my_path)

vx.treeFit('D0:Kpi', conf_level=0.001, updateAllDaughters=True, ipConstraint=True, path=my_path)

va.variables.addAlias("genGrandmotherPDG", "genMotherPDG(1)")


track =["dr", "dx","dy","dz"]
other_sig = ["M","genMotherPDG","genGrandmotherPDG","charge"]
particleIDs = ['pionID', 'kaonID', 'protonID','electronID', 'muonID', 'deuteronID', 'pidIsMostLikely()']
trk_variables = ['d0', 'z0', 'phi0', 'omega', 'ndf', 'pValue','nPXDHits', 'nSVDHits', 'nCDCHits','lastCDCLayer', 'firstCDCLayer','nVXDHits']    
NN_variables = ['kaonIDNN','pionIDNN']
flight_variables = ['mcFlightTime', 'mcFlightDistance','flightDistance','flightTime',"significanceOfDistance","distance"]
CFT = ["CFT_qr","CFT_prob"]

vars = vu.create_aliases_for_selected(list_of_variables=vc.mc_truth, decay_string="^D0 -> ^K- ^pi+") + \
	vu.create_aliases_for_selected(list_of_variables=track, decay_string="^D0 -> ^K- ^pi+") + \
	vu.create_aliases_for_selected(list_of_variables=other_sig, decay_string="^D0 -> ^K- ^pi+") + \
	vu.create_aliases_for_selected(list_of_variables=particleIDs, decay_string="D0 -> ^K- ^pi+") + \
	vu.create_aliases_for_selected(list_of_variables=kin_variables, decay_string="^D0 -> ^K- ^pi+") + \
	vu.create_aliases_for_selected(list_of_variables=trk_variables, decay_string="^D0 -> ^K- ^pi+") + \
	vu.create_aliases_for_selected(list_of_variables=NN_variables, decay_string="D0 -> ^K- ^pi+") + \
	vu.create_aliases_for_selected(list_of_variables=flight_variables, decay_string="^D0 -> ^K- ^pi+") + \
	vu.create_aliases_for_selected(list_of_variables=CFT, decay_string="^D0 -> ^K- ^pi+") + \
	vu.create_aliases_for_selected(list_of_variables=cms_kinematics, decay_string="^D0 -> ^K- ^pi+") 

   

ma.variablesToNtuple(
    "D0:Kpi",
    vars,
    filename=output_file,
    treename="D0tree",
    path=my_path,
)


b2.process(my_path)
