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
import vertex as vx
import os
import math
import sys

my_path = b2.create_path()
b2.set_log_level(b2.LogLevel.ERROR)

input_file = []
output_file = "out_CFT.root"
ma.inputMdstList(input_file, path=my_path)


ma.fillParticleList(
    "pi+:D0", "thetaInCDCAcceptance and dr < 1 and abs(dz) < 3 and pionID>0.9", path=my_path
)
ma.fillParticleList(
    "K+:D0", "thetaInCDCAcceptance and dr < 1 and abs(dz) < 3 and kaonID>0.9", path=my_path
)


ma.reconstructDecay(
    decayString="D0:full -> K-:D0 pi+:D0",
    cut="1.78 < InvM < 1.92",
    path=my_path,
)

ma.cutAndCopyList('D0:Kpi', 'D0:full', 'eventRandom < 0.1', path=my_path)

ma.matchMCTruth(list_name="D0:Kpi", path=my_path)

vx.treeFit('D0:Kpi',0.001,path=my_path)

va.variables.addAlias("genGrandmotherPDG", "genMotherPDG(1)")
va.variables.addAlias("binaryPID_K_pi", "binaryPID(321,211)")
va.variables.addAlias("binaryPID_mu_pi", "binaryPID(13,211)")

kinematics = ["pt", "p", "E","cosTheta","theta","phi"]
cms_kinematics = vu.create_aliases(kinematics, "useCMSFrame({variable})", "CMS")

other_sig = ["M","genMotherPDG","genGrandmotherPDG","charge","protonID","electronID","kaonID","pionID","muonID"]

roe_vars = []

ma.buildRestOfEvent('D0:Kpi', path=my_path)

cleanMask = ('cleanMask', 'dr < 1 and abs(dz) < 3','')
ma.appendROEMasks('D0:Kpi', [cleanMask], path = my_path)



va.variables.addAlias("nROE_Kaons","nROE_Charged(cleanMask,321)")
va.variables.addAlias("nROE_Muons","nROE_Charged(cleanMask,13)")
va.variables.addAlias("nROE_Tracks","nROE_Tracks(cleanMask)")
va.variables.addAlias("ROE_M","roeM(cleanMask)")
va.variables.addAlias("ROE_E","roeE(cleanMask)")
va.variables.addAlias("ROE_p","roeP(cleanMask)")
va.variables.addAlias("ROE_pt","roePt(cleanMask)")
va.variables.addAlias("ROE_ptheta","roePTheta(cleanMask)")
va.variables.addAlias("ROE_DeltaE","roeDeltae(cleanMask)")

d_vars = kinematics+cms_kinematics+vc.mc_truth+["M","genMotherPDG","charge"]
vars = vu.create_aliases_for_selected(
    list_of_variables=vc.mc_truth
    + kinematics
    + cms_kinematics
    + other_sig,
    decay_string="D0 -> ^K- ^pi+",
)

ROE_event_vars = ['nROE_Kaons','nROE_Muons','nROE_Tracks','ROE_M','ROE_E','ROE_p','ROE_pt','ROE_DeltaE','ROE_ptheta']


ma.variablesToNtuple(
    "D0:Kpi",
    d_vars+vars+ROE_event_vars,
    filename=output_file,
    treename="D0tree",
    path=my_path,
)


b2.process(my_path)
