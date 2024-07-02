from TIMBER.Analyzer import *
from TIMBER.Tools.Common import *

import ROOT
from ROOT import TFile
import sys, os
import gc

gc.disable()

from TIMBER.Tools.RestFramesHandler import load_restframes

# From https://gist.github.com/pieterdavid/a560e65658386d70a1720cb5afe4d3e9#file-df-py  example
import correctionlib
correctionlib.register_pyroot_binding()

sys.path.append('../../')
sys.path.append('../../../')

num_threads = 1
#file_name = 'root://cmsxrootd.fnal.gov//store/mc/RunIISummer20UL18NanoAODv9/TprimeTprime_M-1500_TuneCP5_13TeV-madgraph-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/447AD74F-034B-FA42-AD05-CD476A98C43D.root'
file_name = 'ourtestfile.root'

# Import the C++
CompileCpp('TIMBER/Framework/include/common.h') # Compile (via gInterpreter) commonly used c++ code
CompileCpp('TIMBER/Framework/Tprime1lep/cleanjet.cc') # Compile Our vlq c++ code
CompileCpp('TIMBER/Framework/Tprime1lep/utilities.cc') # Compile Our vlq c++ code
CompileCpp('TIMBER/Framework/Tprime1lep/lumiMask.cc')
CompileCpp('TIMBER/Framework/Tprime1lep/selfDerived_corrs.cc')
CompileCpp('TIMBER/Framework/Tprime1lep/corrlib_funcs.cc') 
ROOT.gInterpreter.ProcessLine('#include "TString.h"')

handler_name = 'Tprime_handler.cc'
class_name = 'Tprime_RestFrames_Container'

# Enable using 4 threads
ROOT.ROOT.EnableImplicitMT(num_threads)

# load rest frames handler
load_restframes(num_threads, handler_name, class_name, 'rfc')

# Create analyzer instance
a = analyzer(file_name)

print('==========================INITIALIZED ANALYZER========================')

# ------------------ Command Line Arguments ------------------
campaign = sys.argv[1]

# set variable for year depending on campaign input, preserves some lines calling year
#will this work? I think so... selfDerived_corr (where year is defined) is compiled above this (?) 
if (campaign == "Summer22"): year = "2022";

elif (campaign == "Summer22EE"): year = "2022";

#elif (campaign == "Prompt"): year = "2022";

elif (campaign == "Summer23"): year = "2023";

elif (campaign == "Summer23BPix"): year = "2023";

#selfDerived_corr.cc --> replace year block with 2022/2023 corrections?
else: print(f'ERROR: Can\'t parse the RUN 3 campaign to assign correctionLib json file. Expected Summer22, Summer22EE, Prompt, Summer23, or Summer23BPix. Got: {campaign}\n')

# ------------------ Important Variables ------------------

#TODO isMC? isVV? isSig? etc.
isMC = True
'''if "/mc/" in file_name:
  isMC = True
else: 
  isMC = False'''
#   isMC = !(sampleName.Contains("Single") || sampleName.Contains("Data18") || sampleName.Contains("EGamma"));

isData = False
#   if(inputFile.find("Single") != std::string::npos || inputFile.find("EGamma") != std::string::npos) isData = true;

jesvar = "Nominal"
if not isData:
  jesvar; # will have it run through each of the options
  #"Nominal","JECup","JECdn","JERup","JERdn"

debug = False

ROOT.gInterpreter.Declare("""string year = "' + year + '"; bool isMC = \""""+str(isMC)+"""\"; bool debug = \""""+str(debug)+"""\"; string jesvar = "' + jesvar + '"; """)

# ------------------ Golden JSON Data ------------------
# change the jsonfile path to somewhere they have it in TIMBER
jsonfile = "../TIMBER/data/LumiJSON/"
if (year == "2022"): jsonfile = jsonfile #+ "(...)_Collisions22_*JSON.txt"
elif (year == "2023"): jsonfile = jsonfile #+ "(...)_Collisions23_*JSON.txt"
else: print(f'ERROR: Can\'t parse the year to assign a golden json file. Try changing the \"campaign\" initial argument.')

ROOT.gInterpreter.Declare("""
const auto myLumiMask = lumiMask::fromJSON(\"""" + jsonfile + """\");
//  std::cout << "Testing the JSON! Known good run/lumi returns: " << myLumiMask.accept(315257, 10) << ", and known bad run returns: " << myLumiMask.accept(315257, 90) << std::endl;
""")

# ------------------ Self-derived corrections ------------------

#TODO more things here

    # Lepton scale factors not in correctionLib
ROOT.gInterpreter.ProcessLine('initialize(campaign);')
#cami - initialize with campaign becaus ethe corrections are unique for each one
#change in selfderived -- not necessary bc it's calling the function -- setup corr?

#can do things like this inside:   include <iostream>
#using namespace std; 
#std::cout << elid_pts.size() << elid_pts.at(2); //TODO remove
#std::cout << elecidsfs.size() << elecidsfs.at(0).size() << elecidsfs.at(0).at(0) << endl;
#""")

# ------------------ correctionsLib corrections ------------------
mutrig = "TkMu50";
#missing what changes to do for yr, jecyr, jeryr, and jecver
if (campaign == "Summer22"): #CAMI: check where else yrstr is used, change to "era_num"
# deepjetL = "0.0508"; era_num = "2022"; yr = "16"; jecyr = "UL16APV"; jeryr = "Summer20UL16APV_JRV3"; jecver = "V7"
elif (campaign == "Summer22EE"): 
# deepjetL = "0.0480"; era_num = "2022"; yr = "16"; jecyr = "UL16"; jeryr = "Summer20UL16_JRV3"; jecver = "V7"
#elif (campaign == "Prompt"): #only used in jetvetomaps.json
# deepjetL = "0.0480"; era_num = "2022"; yr = "16"; jecyr = "UL16"; jeryr = "Summer20UL16_JRV3"; jecver = "V7"
elif (campaign == "Summer23"): 
# mutrig = "OldMu100_or_TkMu100"; deepjetL = "0.0532"; era_num = "2023"; yr = "17"; jecyr = "UL17"; jeryr = "Summer19UL17_JRV2"; jecver = "V5"
elif (campaign == "Summer23BPix"): 
# mutrig = "OldMu100_or_TkMu100"; deepjetL = "0.0490"; era_num = "2023"; yr = "18"; jecyr = "UL18"; jeryr = "Summer19UL18_JRV2"; jecver = "V5"
else: print(f'ERROR: Can\'t parse the year to assign correctionLib json files. Expected 2022 or 2023. Try changing the \"campaign\" initial argument.')


#change era_num to year cami
ROOT.gInterpreter.Declare("""
string year = \""""+era_num+"""\"; 
string yr = \""""+yr+"""\"; 
string jecyr = \""""+jecyr+"""\"; 
string jeryr = \""""+jeryr+"""\"; 
string jecver = \""""+jecver+"""\"; 
string mutrig = \""""+mutrig+"""\";
float deepjetL = """+deepjetL+""";
""")


ROOT.gInterpreter.Declare("""
auto csetPU = correction::CorrectionSet::from_file("jsonpog-integration/POG/LUM/"+era_num+"_"+campaign+"/puWeights.json");
auto electroncorrset = correction::CorrectionSet::from_file("jsonpog-integration/POG/EGM/"+era_num+"_"+campaign+"/electron.json");
auto muoncorrset = correction::CorrectionSet::from_file("jsonpog-integration/POG/MUO/"+era_num+"_"+campaign+"/muon_Z.json");
auto jetvetocorrset = correction::CorrectionSet::from_file("jsonpog-integration/POG/JME/"+era_num+"_"+campaign+"/jetvetomaps.json");

# auto metcorrset = correction::CorrectionSet::from_file("jsonpog-integration/POG/JME/"+era_num+"_"+campaign+"/met.json");
# didn't find a met.json correction. Only jetvetomap, jet_jerc and fatJet_jerc

#cami change the if else to only change the dd=string so we run the correction once.
if (campaign == "Summer22"){ 
   auto corrPU = csetPU->at("Collisions2022_355100_357900_eraBCD_GoldenJson");
   

   };
elif (campaign == "Summer22EE"){ 

   auto corrPU = csetPU->at("Collisions2022_359022_362760_eraEFG_GoldenJson");

   };
elif (campaign == "Summer23"){ #only used in jetvetomaps.json

   auto corrPU = csetPU->at("Collisions2023_366403_369802_eraBC_GoldenJson");

   };
elif (campaign == "Summer23BPix"){ 

   auto corrPU = csetPU->at("Collisions2023_369803_370790_eraD_GoldenJson");

   };
elif (campaign == "Prompt"){

   };

else { False };

Electron-ID-SFauto electroncorr = electroncorrset->at("Electron-ID-SF"); #for 2022 it's "v2" and for 2023 it's "v3"
 
auto muoncorr = muoncorrset->at("NUM_TrackerMuons_DEN_genTracks"); ##no name matches closely? cami -- comment out (?)

## in recofunc, change the muon part of the func to return value one instead of figuring our value from corr

auto muonidcorr = muoncorrset->at("NUM_MediumID_DEN_TrackerMuons");

auto muonhltcorr = muoncorrset->at("NUM_Mu50_or_"+mutrig+"_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose"); 

auto jetvetocorr = jetvetocorrset->at("Summer19UL"+yr+"_V1");

auto metptcorr = metcorrset->at("pt_metphicorr_pfmet_mc");

auto metphicorr = metcorrset->at("phi_metphicorr_pfmet_mc");
//if (!isMC) {
//  metptcorr = metcorrset->at("pt_metphicorr_pfmet_data");
//  metphicorr = metcorrset->at("phi_metphicorr_pfmet_data"); };

auto ak4corrset = correction::CorrectionSet::from_file("jsonpog-integration/POG/JME/"+yrstr+"_UL/jet_jerc.json"); 

auto ak8corrset = correction::CorrectionSet::from_file("jsonpog-integration/POG/JME/"+yrstr+"_UL/fatJet_jerc.json"); 

auto ak4corr = ak4corrset->compound().at("Summer19"+jecyr+"_"+jecver+"_MC_L1L2L3Res_AK4PFchs");

auto ak4corrL1 = ak4corrset->at("Summer19"+jecyr+"_"+jecver+"_MC_L1FastJet_AK4PFchs"); 
//if(!isMC){ ak4corr = ak4corrset->compound().at("Summer19"+jecyr+"_Run"+jecera+"_"+jecver+"_DATA_L1L2L3Res_AK4PFchs"); };

auto ak4corrUnc = ak4corrset->at("Summer19"+jecyr+"_"+jecver+"_MC_Total_AK4PFchs"); 

auto ak4ptres = ak4corrset->at(jeryr+"_MC_PtResolution_AK4PFchs"); 

auto ak4jer = ak4corrset->at(jeryr+"_MC_ScaleFactor_AK4PFchs"); 

auto ak8corr = ak8corrset->compound().at("Summer19"+jecyr+"_"+jecver+"_MC_L1L2L3Res_AK8PFPuppi"); 
//if(!isMC){ ak8corr = ak8corrset->compound().at("Summer19"+jecyr+"_Run"+jecera+"_"+jecver+"_DATA_L1L2L3Res_AK8PFPuppi"); };

auto ak8corrUnc = ak8corrset->at("Summer19"+jecyr+"_"+jecver+"_MC_Total_AK8PFPuppi"); 
""") #TODO if statement doesn't work in this .Declare() ???

#from muonhltcorr => std::cout << "\t loaded muon trig" << std::endl; // REDO ME (Do we need to change something?)

# ------------------ MET Cuts ------------------
metCuts = CutGroup('METCuts')
metCuts.Add('MET Filters', 'Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && Flag_goodVertices == 1 && Flag_HBHENoiseFilter == 1 && Flag_HBHENoiseIsoFilter == 1 && Flag_eeBadScFilter == 1 && Flag_globalSuperTightHalo2016Filter == 1 && Flag_BadPFMuonFilter == 1 && Flag_ecalBadCalibFilter == 1')
metCuts.Add('Pass MET > 50', 'MET_pt > 50')
metCuts.Add('Event has jets',        'nJet > 0 && nFatJet > 0') # need jets 

# ------------------ Golden JSON (Data) || GEN Info (MC) ------------------
gjsonVars = VarGroup('GoldenJsonVars')
gjsonCuts = CutGroup('GoldenJsonCuts')
if isMC is False: # apply golden json to data
  gjsonVars.Add("passesJSON", "goldenjson(myLumiMask, run, luminosityBlock)")
  gjsonCuts.Add("Data passes Golden JSON", "passesJSON == 1") 
else:
  gjsonVars.Add("PileupWeights", "pufunc(corrPU, Pileup_nTrueInt)")

# ------------------ LEPTON Definitions ------------------
lVars = VarGroup('LeptonVars')

'''
lVars.Add("TightMu", "abs(Muon_eta) < 2.4 && Muon_tightId == true && Muon_miniIsoId >= 3 && Muon_pt > 50")
lVars.Add("TightEl", "abs(Electron_eta) < 2.5 && Electron_mvaFall17V2noIso_WP90 == true && Electron_miniPFRelIso_all < 0.1 && Electron_pt > 50")
lVars.Add("VetoMu", "abs(Muon_eta) < 2.4 && Muon_looseId == true && Muon_miniIsoId >= 1 && Muon_pt > 10 && TightMu == false")
lVars.Add("VetoEl", "abs(Electron_eta) < 2.5 && Electron_mvaFall17V2noIso_WPL == true && Electron_miniPFRelIso_all < 0.4 && Electron_pt > 10 && TightEl == false")
lVars.Add("nVetoLep", "Sum(VetoMu)+Sum(VetoEl)")
lVars.Add("nTightMu", "Sum(TightMu)")
lVars.Add("nTightEl", "Sum(TightEl)")
lVars.Add("TMuon_pt", "Muon_pt[TightMu == true]")
lVars.Add("TMuon_eta", "Muon_eta[TightMu == true]")
lVars.Add("TMuon_phi", "Muon_phi[TightMu == true]")
lVars.Add("TMuon_mass", "Muon_mass[TightMu == true]")
lVars.Add("TElectron_pt", "Electron_pt[TightEl == true]")
lVars.Add("TElectron_eta", "Electron_eta[TightEl == true]")
lVars.Add("TElectron_phi", "Electron_phi[TightEl == true]")
lVars.Add("TElectron_mass", "Electron_mass[TightEl == true]")
lVars.Add("TMuon_P4", "fVectorConstructor(TMuon_pt,TMuon_eta,TMuon_phi,TMuon_mass)")
lVars.Add("TElectron_P4", "fVectorConstructor(TElectron_pt,TElectron_eta,TElectron_phi,TElectron_mass)")
lVars.Add("TMuon_jetIdx", "Muon_jetIdx[TightMu == true]")
lVars.Add("TElectron_jetIdx", "Electron_jetIdx[TightEl == true]")

'''
if year == "2018": elHEMcut = " && (Electron_eta > -1.479 || (Electron_phi < -1.57 || Electron_phi > -0.87))"
ROOT.gInterpreter.Declare('string elHEMcut = "'+elHEMcut+'"; ')

lVars.Add("Electron_cutBasedIdNoIso_tight", "Electron_cutBasedIdNoIso_tight(nElectron, Electron_vidNestedWPBitmap)")
lVars.Add("TPassMu", "abs(Muon_eta)<2.4 && Muon_mediumId==1 && Muon_miniIsoId>=3 && abs(Muon_dz) < 0.5 && Muon_dxy < 0.2")
#lVars.Add("TPassEl", "Form(\"(abs(Electron_eta)<1.442 || (abs(Electron_eta)>1.566 && abs(Electron_eta)<2.5)) && Electron_cutBasedIdNoIso_tight==1 && Electron_miniPFRelIso_all<0.1%s\",elHEMcut.c_str())")
lVars.Add("TPassEl", "(abs(Electron_eta)<1.442 || (abs(Electron_eta)>1.566 && abs(Electron_eta)<2.5)) && Electron_cutBasedIdNoIso_tight==1 && Electron_miniPFRelIso_all<0.1" + elHEMcut)
lVars.Add("VetoMu", "TPassMu && (Muon_pt>25)")
lVars.Add("VetoEl", "TPassEl && (Electron_pt>25)")
lVars.Add("SignalIsoMu", "TPassMu && (Muon_pt>=55)")
lVars.Add("SignalIsoEl", "TPassEl && (Electron_pt>=55)")
lVars.Add("nVetoLep", "(int) (Sum(VetoMu)+Sum(VetoEl))")
lVars.Add("SMuon_pt", "Muon_pt[SignalIsoMu == true]")
lVars.Add("SMuon_eta", "Muon_eta[SignalIsoMu == true]")
lVars.Add("SMuon_phi", "Muon_phi[SignalIsoMu == true]")
lVars.Add("SMuon_mass", "Muon_mass[SignalIsoMu == true]")
lVars.Add("SElectron_pt", "Electron_pt[SignalIsoEl == true]")
lVars.Add("SElectron_eta", "Electron_eta[SignalIsoEl == true]")
lVars.Add("SElectron_phi", "Electron_phi[SignalIsoEl == true]")
lVars.Add("SElectron_mass", "Electron_mass[SignalIsoEl == true]")
lVars.Add("Muon_P4", "fVectorConstructor(Muon_pt,Muon_eta,Muon_phi,Muon_mass)")
lVars.Add("SMuon_P4", "fVectorConstructor(SMuon_pt,SMuon_eta,SMuon_phi,SMuon_mass)")
lVars.Add("SElectron_P4", "fVectorConstructor(SElectron_pt,SElectron_eta,SElectron_phi,SElectron_mass)")
lVars.Add("SMuon_jetIdx", "Muon_jetIdx[SignalIsoMu == true]")
lVars.Add("SElectron_jetIdx", "Electron_jetIdx[SignalIsoEl]")
lVars.Add("nSignalIsoMu", "(int) Sum(SignalIsoMu)")
lVars.Add("nSignalIsoEl", "(int) Sum(SignalIsoEl)")
lVars.Add("VetoIsoMu", "(VetoMu == true && Muon_pt < 55)")
lVars.Add("VetoIsoEl", "(VetoEl == true && Electron_pt < 55)")
lVars.Add("nVetoIsoLep", "(int) (Sum(VetoIsoMu)+Sum(VetoIsoEl))")

# ------------------ LEPTON SELECTION ------------------

# ||||||||||||||||||||| TODO ||||||||||||||||||||||||
### WE need to find out which triggers exist in run 3
# ___________________________________________________

tkmutrig = " || HLT_OldMu100 || HLT_TkMu100"
eltrig = "HLT_Ele115_CaloIdVT_GsfTrkIdT || HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 || HLT_Photon200"
if(year == "2017" and era == "B"): 
  tkmutrig = ""
  eltrig = "HLT_Ele35_WPTight_Gsf || HLT_Photon200"
if(year == "2016" or year == "2016APV"):
  tkmutrig = " || HLT_TkMu50"
  eltrig = "HLT_Ele115_CaloIdVT_GsfTrkIdT || HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 || HLT_Photon175"
if(year == "2016APV" and (era == "A" or era == "B")):
    tkmutrig = ""
    eltrig = "HLT_Ele115_CaloIdVT_GsfTrkIdT || HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 || HLT_Photon175"
ROOT.gInterpreter.Declare('string tkmutrig = "'+tkmutrig+'"; string eltrig = "'+eltrig+'"; ')

#lVars.Add("isMu", "Form(\"(nMuon>0) && (HLT_Mu50%s) && (nSignalIsoMu==1) && (nVetoIsoLep==0) && (nElectron == 0 || nSignalIsoEl == 0)\",tkmutrig.c_str())")
#lVars.Add("isEl", "Form(\"(nElectron>0) && (%s) && (nSignalIsoEl==1) && (nVetoIsoLep==0) && (nMuon == 0 || nSignalIsoMu == 0)\",eltrig.c_str())")
lVars.Add("isMu", "(nMuon>0) && (HLT_Mu50"+ tkmutrig +") && (nSignalIsoMu==1) && (nVetoIsoLep==0) && (nElectron == 0 || nSignalIsoEl == 0)")
lVars.Add("isEl", "(nElectron>0) && ("+ eltrig +") && (nSignalIsoEl==1) && (nVetoIsoLep==0) && (nMuon == 0 || nSignalIsoMu == 0)")
#lVars.Add("isMu","nMuon > 0 && nTightMu == 1 && (nElectron == 0 || nTightEl == 0) && nVetoLep == 0 && (HLT_Mu50 == 1 || HLT_Mu15_IsoVVVL_PFHT450 == 1)") 
#lVars.Add("isEl","nElectron > 0 && nTightEl == 1 && (nMuon == 0 || nTightMu == 0) && nVetoLep==0 && (HLT_Ele35_WPTight_Gsf == 1 || HLT_Ele15_IsoVVVL_PFHT450 == 1)") 
        
        # filter lepton
lCuts = CutGroup('LeptonCuts')
lCuts.Add("Event is either muon or electron", "isMu || isEl")
        # assign lepton
#lVars.Add("assignleps","assign_leps(isMu,isEl,TightMu,TightEl,Muon_pt,Muon_eta,Muon_phi,Muon_mass,Muon_miniPFRelIso_all,Electron_pt,Electron_eta,Electron_phi,Electron_mass,Electron_miniPFRelIso_all)")
lVars.Add("assignleps", "assign_leps(isMu,isEl,SignalIsoMu,SignalIsoEl,Muon_pt,Muon_eta,Muon_phi,Muon_mass,Muon_miniPFRelIso_all,Electron_pt,Electron_eta,Electron_phi,Electron_mass,Electron_miniPFRelIso_all)")
lVars.Add("lepton_pt","assignleps[0]")
lVars.Add("lepton_eta","assignleps[1]")	
lVars.Add("lepton_phi","assignleps[2]")
lVars.Add("lepton_mass","assignleps[3]")
lVars.Add("lepton_miniIso","assignleps[4]")


# ------------------ JET Cleaning and JERC ------------------
jVars = VarGroup('JetCleaningVars')

jVars.Add("Jet_P4", "fVectorConstructor(Jet_pt,Jet_eta,Jet_phi,Jet_mass)")
jVars.Add("FatJet_P4", "fVectorConstructor(FatJet_pt,FatJet_eta,FatJet_phi,FatJet_mass)")
jVars.Add("Jet_EmEF","Jet_neEmEF + Jet_chEmEF")
#a.Define("DummyZero","0.0")
jVars.Add("DummyZero","float(0.0)")
        # Clean Jets

###TODO sily things: signal X is = tight X ???
#jVars.Add("SMuon_P4","TMuon_P4")
#jVars.Add("SMuon_jetIdx","TMuon_jetIdx")
#jVars.Add("SElectron_P4","TElectron_P4")
#jVars.Add("SElectron_jetIdx","TElectron_jetIdx")

#gc.disable()

if isMC:
  jVars.Add("GenJet_P4","fVectorConstructor(GenJet_pt,GenJet_eta,GenJet_phi,GenJet_mass)")
  jVars.Add("cleanedJets", "cleanJets(debug,jesvar,isMC,ak4corr,ak4corrL1,ak4corrUnc,ak4ptres,ak4jer,ak8corr,ak8corrUnc,Jet_P4,Jet_rawFactor,Jet_muonSubtrFactor,Jet_area,Jet_EmEF,Jet_jetId,GenJet_P4,Jet_genJetIdx,SMuon_P4,SMuon_jetIdx,SElectron_P4,SElectron_jetIdx,fixedGridRhoFastjetAll,DummyZero,DummyZero)") # muon and EM factors unused in this call
  jVars.Add("cleanMets", "cleanJets(debug,jesvar,isMC,ak4corr,ak4corrL1,ak4corrUnc,ak4ptres,ak4jer,ak8corr,ak8corrUnc,Jet_P4,Jet_rawFactor,Jet_muonSubtrFactor,Jet_area,Jet_EmEF,Jet_jetId,GenJet_P4,Jet_genJetIdx,SMuon_P4,SMuon_jetIdx,SElectron_P4,SElectron_jetIdx,fixedGridRhoFastjetAll,RawMET_pt,RawMET_phi)") # lepton args are unused in this call
  jVars.Add("GenJetAK8_P4", "fVectorConstructor(GenJetAK8_pt,GenJetAK8_eta,GenJetAK8_phi,GenJetAK8_mass)")
  jVars.Add("cleanFatJets", "cleanJets(debug,jesvar,isMC,ak4corr,ak4corrL1,ak4corrUnc,ak4ptres,ak4jer,ak8corr,ak8corrUnc,FatJet_P4,FatJet_rawFactor,FatJet_rawFactor,FatJet_area,FatJet_area,FatJet_jetId,GenJetAK8_P4,FatJet_genJetAK8Idx,SMuon_P4,SMuon_jetIdx,SElectron_P4,SElectron_jetIdx,fixedGridRhoFastjetAll,DummyZero,DummyZero)") # args 12 and 14 are dummies
else:
    # Replace all the GenJet arguments with fakes here for data. 
  jVars.Add("cleanedJets", "cleanJets(debug,jesvar,isMC,ak4corr,ak4corrL1,ak4corrUnc,ak4ptres,ak4jer,ak8corr,ak8corrUnc,Jet_P4,Jet_rawFactor,Jet_muonSubtrFactor,Jet_area,Jet_EmEF,Jet_jetId,Jet_P4,Jet_jetId,SMuon_P4,SMuon_jetIdx,SElectron_P4,SElectron_jetIdx,fixedGridRhoFastjetAll,DummyZero,DummyZero)") # muon and EM factors unused in this call, args 16-17 are dummies
  jVars.Add("cleanMets", "cleanJets(debug,jesvar,isMC,ak4corr,ak4corrL1,ak4corrUnc,ak4ptres,ak4jer,ak8corr,ak8corrUnc,Jet_P4,Jet_rawFactor,Jet_muonSubtrFactor,Jet_area,Jet_EmEF,Jet_jetId,Jet_P4,Jet_jetId,Muon_P4,Muon_jetIdx,SElectron_P4,SElectron_jetIdx,fixedGridRhoFastjetAll,RawMET_pt,RawMET_phi)") # lepton args unused in this call, args 16-17 are dummies
  jVars.Add("cleanFatJets", "cleanJets(debug,jesvar,isMC,ak4corr,ak4corrL1,ak4corrUnc,ak4ptres,ak4jer,ak8corr,ak8corrUnc,FatJet_P4,FatJet_rawFactor,FatJet_rawFactor,FatJet_area,FatJet_area,FatJet_jetId,FatJet_P4,FatJet_jetId,SMuon_P4,SMuon_jetIdx,SElectron_P4,SElectron_jetIdx,fixedGridRhoFastjetAll,DummyZero,DummyZero)") # args 12, 14, 16, 17 are dummies
        # Jet Assign

jVars.Add("cleanJet_pt", "return Map(cleanedJets, [&](const TLorentzVector& vec) { return vec.Pt(); });")
jVars.Add("cleanJet_eta", "return Map(cleanedJets, [&](const TLorentzVector& vec) { return vec.Eta(); });")
jVars.Add("cleanJet_phi", "return Map(cleanedJets, [&](const TLorentzVector& vec) { return vec.Phi(); });")
jVars.Add("cleanJet_mass", "return Map(cleanedJets, [&](const TLorentzVector& vec) { return vec.M(); });")

jVars.Add("cleanFatJet_pt", "return Map(cleanFatJets, [&](const TLorentzVector& vec) { return vec.Pt(); });")
jVars.Add("cleanFatJet_eta", "return Map(cleanFatJets, [&](const TLorentzVector& vec) { return vec.Eta(); });")
jVars.Add("cleanFatJet_phi", "return Map(cleanFatJets, [&](const TLorentzVector& vec) { return vec.Phi(); });")
jVars.Add("cleanFatJet_mass", "return Map(cleanFatJets, [&](const TLorentzVector& vec) { return vec.M(); });")

#print(a.GetCollectionNames())
#collectNames = a.GetCollectionNames()
#for cn in collectNames:
    #print(cn)
#    print(a._collectionOrg.GetCollectionAttributes(cn))

#jVars.Add("cleanJet_pt", "cleanedJets[0]")
#jVars.Add("cleanJet_eta", "cleanedJets[1]")
#jVars.Add("cleanJet_phi", "cleanedJets[2]")
#jVars.Add("cleanJet_mass", "cleanedJets[3]")
#jVars.Add("cleanJet_rawFactor", "cleanedJets[4]")
#jVars.Add("cleanFatJet_pt", "cleanFatJets[0]")
#jVars.Add("cleanFatJet_eta", "cleanFatJets[1]")
#jVars.Add("cleanFatJet_phi", "cleanFatJets[2]")
#jVars.Add("cleanFatJet_mass", "cleanFatJets[3]")
#jVars.Add("cleanFatJet_rawFactor", "cleanFatJets[4]")


# ------------------ MET Selection ------------------
metVars = VarGroup('METVars')

metVars.Add("corrMETnoxy_pt","return Map(cleanMets, [&](const TLorentzVector& vec) { return vec.Pt(); });")
metVars.Add("corrMETnoxy_phi","return Map(cleanMets, [&](const TLorentzVector& vec) { return vec.Phi(); });")
#metVars.Add("corrMETnoxy_pt","cleanMets[5][0]")
#metVars.Add("corrMETnoxy_phi","cleanMets[5][1]")

metVars.Add("metxyoutput", "metfunc(metptcorr, metphicorr, corrMETnoxy_pt, corrMETnoxy_phi, PV_npvs, run)")

metVars.Add("corrMET_pt","metxyoutput[0]")
metVars.Add("corrMET_phi","metxyoutput[1]")

metCuts.Add("Pass corr MET > 60", "corrMET_pt > 60")
metCuts.Add("Electron Triangle Cut", "isMu || corrMET_pt>((130/1.5)*DeltaPhi(lepton_phi, corrMET_phi)-130)")

# ------------------ HT Calculation and N Jets cuts ------------------
#TODO UNCOMMENT jVars.Add("DR_lepJets","DeltaR_VecAndFloat(cleanedJet_eta,cleanedJet_phi,lepton_eta,lepton_phi)")
#               jVars.Add("ptrel_lepJets","ptRel(cleanedJet_pt,cleanedJet_eta,cleanedJet_phi,cleanedJet_mass,lepton_pt,lepton_eta,lepton_phi,lepton_mass)") 
#               jVars.Add("goodcleanJets", "cleanedJet_pt > 30 && abs(cleanedJet_eta) < 2.4 && Jet_jetId > 1 && (DR_lepJets > 0.4 || ptrel_lepJets > 20)")

jVars.Add("DR_lepFatJets","DeltaR_VecAndFloat(FatJet_eta,FatJet_phi,lepton_eta,lepton_phi)")
jVars.Add("ptrel_lepFatJets","ptRel(FatJet_pt,FatJet_eta,FatJet_phi,FatJet_mass,lepton_pt,lepton_eta,lepton_phi,lepton_mass)") 
jVars.Add("goodcleanFatJets", "FatJet_pt > 200 && abs(FatJet_eta) < 2.4 && FatJet_jetId > 1 && (DR_lepFatJets > 0.8 || ptrel_lepFatJets > 20)")
jVars.Add("NFatJets_central", "(int) Sum(goodcleanFatJets)")

#jVars.Add('AK4HT', 'Sum(gcJet_pt)')    

jCuts = CutGroup('JetCuts')
#jCuts.Add('AK4 HT Pass', 'AK4HT > 510')    
jCuts.Add('3 AK8s Pass', 'NFatJets_central > 2')    # need to ensure three jets exist


# ------------------ Jet pt ordering, counting, lepton association ------------------
# requires clean jet things for this:
'''jVars.Add("gcJet_pt_unsort", "cleanJet_pt[goodcleanJets == true]")
jVars.Add("gcJet_ptargsort","ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(gcJet_pt_unsort))")
jVars.Add("gcJet_pt","reorder(gcJet_pt_unsort,gcJet_ptargsort)")
jVars.Add("gcJet_eta", "reorder(cleanJet_eta[goodcleanJets == true],gcJet_ptargsort)")
jVars.Add("gcJet_phi", "reorder(cleanJet_phi[goodcleanJets == true],gcJet_ptargsort)")
jVars.Add("gcJet_mass", "reorder(cleanJet_mass[goodcleanJets == true],gcJet_ptargsort)")
jVars.Add("gcJet_vetomap", "jetvetofunc(jetvetocorr, gcJet_eta, gcJet_phi)")'''
        #fatjet vars
jVars.Add("gcFatJet_pt_unsort", "FatJet_pt[goodcleanFatJets == true]")
jVars.Add("gcFatJet_ptargsort","ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(gcFatJet_pt_unsort))")
jVars.Add("gcFatJet_pt","reorder(gcFatJet_pt_unsort,gcFatJet_ptargsort)")
jVars.Add("gcFatJet_eta", "reorder(FatJet_eta[goodcleanFatJets == true],gcFatJet_ptargsort)")
jVars.Add("gcFatJet_phi", "reorder(FatJet_phi[goodcleanFatJets == true],gcFatJet_ptargsort)")
jVars.Add("gcFatJet_mass", "reorder(FatJet_mass[goodcleanFatJets == true],gcFatJet_ptargsort)")
jVars.Add("gcFatJet_sdmass", "reorder(FatJet_msoftdrop[goodcleanFatJets == true],gcFatJet_ptargsort)")
jVars.Add("gcFatJet_vetomap", "jetvetofunc(jetvetocorr, gcFatJet_eta, gcFatJet_phi)")
	#condition to only return isolated jets (not inside a fatjet), needs gcJet_* and gcFatJet_*
jVars.Add("Isolated_AK4","standalone_Jet(gcJet_eta, gcJet_phi, gcFatJet_eta, gcFatJet_phi)")

# ------------------ Add scale factors and MC jet-based calcs ------------------
#TODO could be a fatJetVar group
if isMC:
  jVars.Add("leptonRecoSF", "recofunc(electroncorr, muoncorr, yrstr, lepton_pt, lepton_eta, isEl)")
  jVars.Add("leptonIDSF", "idfunc(muonidcorr,elid_pts,elid_etas,elecidsfs,elecidsfuncs,yrstr, lepton_pt, lepton_eta, isEl)") #at(0) 
  jVars.Add("leptonIsoSF", "isofunc(muiso_pts,muiso_etas,muonisosfs,muonisosfunc,elid_pts,elid_etas,elecisosfs,elecisosfunc, lepton_pt, lepton_eta, isEl)")
  jVars.Add("leptonHLTSF", "hltfunc(muonhltcorr,elhlt_pts,elhlt_etas,elechltsfs,elechltuncs,yrstr, lepton_pt, lepton_eta, isEl)")


# ------------------ Results ------------------
rframeVars = VarGroup('restFrameVars')
rframeVars.Add('VLQ_mass', 'rfc.compute_mass(rdfslot_, lepton_pt, lepton_eta, lepton_phi, lepton_mass, gcFatJet_pt, gcFatJet_eta, gcFatJet_phi, gcFatJet_mass, MET_pt, MET_phi)')
rframeVars.Add('VLQ_mass_T', 'VLQ_mass[0]')
rframeVars.Add('VLQ_mass_Tbar', 'VLQ_mass[1]')
rframeVars.Add('VLQ_mass_T_r', 'VLQ_mass[2]')
rframeVars.Add('VLQ_mass_Tbar_r', 'VLQ_mass[3]')
rframeVars.Add('VLQ_mass_ratio', 'VLQ_mass_T/VLQ_mass_Tbar')
rframeVars.Add('VLQ_mass_avg', '(VLQ_mass_T+VLQ_mass_Tbar)*0.5')


# -------------------------------------


nodeToPlot = a.Apply([gjsonVars, gjsonCuts, lVars, lCuts]) #, jVars, jCuts, metVars, metCuts, rframeVars])
#nodeToPlot = a.Apply(metCuts)
#a.Apply(lVars)
#a.Apply(lCuts)
#a.ActiveNode.Apply([jVars, jCuts, rframeVars]) ## where problem is

# Solution to cleanJets() problem:
#       The analyzer .Apply() calls the analyzer .Define().  This .Define() calls self._collectionOrg.CollectionDefCheck(var, newNode).
#  This method executes this line: if re.search(r"\b" + re.escape(c+'s') + r"\b", action_str) and (c+'s' not in self._builtCollections):
#                                       print ('MAKING %ss for %s'%(c,action_str))
#  Apparently somethings get discarded from this _collectionOrg?

newNode = a.ActiveNode.Apply(jVars)
a.SetActiveNode(newNode)

a.Apply([jCuts, rframeVars])

allColumns = a.GetColumnNames()
columns = [] #allColumns

#i = 0
for col in allColumns:
    #i = i + 1
    #if i > 49: continue
    if col == "run": break # lets just skip all the original branches?

    if ("P4" in col) or ("cleanedJets" in col) or ("cleanFatJets" in col) or ("cleanMets" in col) or ("Dummy" in col): continue 
    if ("LHE" in col) and ("Weight" not in col) and (col != "LHE_HT") and (col != "LHE_Vpt") and (col != "gcHTCorr_WjetLHE"): continue
    if col.startswith("Muon") and ("_tightId" not in col) and ("_isPF" not in col) and ("tunep" not in col) and ("genPartFlav" not in col): continue
    if col.startswith("Electron") and ("genPartFlav" not in col): continue
    if col.startswith("Jet") and ("rawFactor" not in col): continue
    if col.startswith("FatJet") and ("rawFactor" not in col): continue
    if col.startswith("PPS") or col.startswith("Proton") or col.startswith("L1_"): continue
    if col.startswith("Gen") or col.startswith("Soft") or col.startswith("fixed"): continue
    if col.startswith("Sub") or col.startswith("RawPuppi") or col.startswith("Calo") or col.startswith("Chs"): continue
    if col.startswith("Corr") or col.startswith("Fsr") or col.startswith("Iso") or col.startswith("Tau"): continue
    if col.startswith("SV") or col.startswith("Puppi") or col.startswith("Photon") or col.startswith("Low"): continue
    if col.startswith("HLT") or col.startswith("HT") or col.startswith("boosted") or col.startswith("Deep"): continue
    if col.startswith("Flag") or col == "Bprime_gen_info" or col == "t_gen_info" or col == "W_gen_info" or col == "metxyoutput": continue
    if col == "assignleps" or col == "pnetoutput" or col == "t_output" or col == "Bprime_output" or col.startswith("Other"): continue
    if col.startswith("PS") or col.startswith("PV") or col.startswith("Tk") or col.startswith("Trig"): continue
    if col.startswith("nCorr") or col.startswith("nFsr"): continue
    if col.startswith("nGen") or col.startswith("nIso") or col.startswith("nLow"): continue
    if col.startswith("nOther") or col.startswith("nPS") or col.startswith("nPhoton"): continue
    if col.startswith("nSV") or col.startswith("nSub") or col.startswith("nTau") or col.startswith("nTrig"): continue
    if col.startswith("nboosted"): continue
    #TODO need to figure out how to exclude the things related to nSub and Sub.
    columns.append(col)

#TODO think do we really want to recreate this everytime?  or just create?
a.Snapshot(columns, "out_Tprime.root", "Events", lazy=False, openOption='RECREATE', saveRunChain=False)

myHist1 = a.GetActiveNode().DataFrame.Histo1D(('m_T_lab', 'Mass of T lab', 25, 500, 2000), 'VLQ_mass_T')
myHist2 = a.GetActiveNode().DataFrame.Histo1D(('m_Tbar_lab', 'Mass of Tbar lab', 25, 500, 2000), 'VLQ_mass_Tbar')
myHist1a = a.GetActiveNode().DataFrame.Histo1D(('m_T', 'Mass of T', 25, 500, 2000), 'VLQ_mass_T_r')
myHist2a = a.GetActiveNode().DataFrame.Histo1D(('m_Tbar', 'Mass of Tbar', 25, 500, 2000), 'VLQ_mass_Tbar_r')
myHist3 = a.GetActiveNode().DataFrame.Histo1D(('m_T/m_Tbar', 'Mass ratio of the two particles', 25, 0, 2), 'VLQ_mass_ratio')
myHist4 = a.GetActiveNode().DataFrame.Histo1D(('m_avg', 'Mass average of the two', 25, 500, 2000), 'VLQ_mass_avg')

out = ROOT.TFile.Open('test_Tprime_out.root','RECREATE') #'UPDATE')
myHist1.Write()
myHist2.Write()
myHist1a.Write()
myHist2a.Write()
myHist3.Write()
myHist4.Write()

print("--------- Analysis End ---------")

out.Close()


a.Close()

