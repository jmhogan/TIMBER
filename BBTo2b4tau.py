from TIMBER.Analyzer import *
from TIMBER.Tools.Common import *

# python3 BBTo2b4tau.py testfile_SIGNAL_2022.txt 1 2 2022
# python3 BBTo2b4tau.py full_SIGNAL_2022.txt 1 2 2022


# X python3 BBTo2b4tau.py full_SIGNAL_2022.txt 0 64 2022
# X python3 BBTo2b4tau.py full_SIGNAL_2022EE.txt 0 72 2022EE
# X python3 BBTo2b4tau.py full_SIGNAL_2023.txt 0 40 2023
# X python3 BBTo2b4tau.py full_SIGNAL_2023.txt 41 60 2023
# X python3 BBTo2b4tau.py full_SIGNAL_2023.txt 61 70 2023
# python3 BBTo2b4tau.py full_SIGNAL_2023.txt 70 75 2023 Gave me errors
# X python3 BBTo2b4tau.py full_SIGNAL_2023.txt 75 80 2023
# python3 BBTo2b4tau.py full_SIGNAL_2023BPix.txt 0 26 2023BPix

# M400 prompts
# python3 BBTo2b4tau.py fullsigtxts/full_SIGNAL_2022M400.txt 0 64 2022
# python3 BBTo2b4tau.py fullsigtxts/full_SIGNAL_2022EEM400.txt 0 72 2022EE
# python3 BBTo2b4tau.py fullsigtxts/full_SIGNAL_2023M400.txt 0 80 2023
# python3 BBTo2b4tau.py fullsigtxts/full_SIGNAL_2023BPixM400.txt 0 26 2023BPix
        #this one has one .root file

# M1000 prompts
# python3 BBTo2b4tau.py fullsigtxts/full_SIGNAL_2022M1000.txt 0 64 2022
# python3 BBTo2b4tau.py fullsigtxts/full_SIGNAL_2022EEM1000.txt 0 72 2022EE
# python3 BBTo2b4tau.py fullsigtxts/full_SIGNAL_2023M1000.txt 0 50 2023 
# python3 BBTo2b4tau.py fullsigtxts/full_SIGNAL_2023BPixM1000.txt 0 2 2023BPix

# M1600 prompts
# python3 BBTo2b4tau.py fullsigtxts/full_SIGNAL_2022M1600.txt 0 64 2022
# python3 BBTo2b4tau.py fullsigtxts/full_SIGNAL_2022EEM1600.txt 0 72 2022EE
# python3 BBTo2b4tau.py fullsigtxts/full_SIGNAL_2023M1600.txt 0 50 2023 
# python3 BBTo2b4tau.py fullsigtxts/full_SIGNAL_2023BPixM1600.txt 0 2 2023BPix

# QCD prompt
# python3 BBTo2b4tau.py fullsigtxts/full_SIGNAL_2023BPixQCD.txt 0 50 2023BPix

'''
source /cvmfs/cms.cern.ch/cmsset_default.sh 
voms-proxy-init --voms cms --valid 168:00 # will be valid for a week, only needed when it's needed
cd nobackup/BBto2b4tau/CMSSW_13_2_10
cmsenv
cd ../
source timber-env/bin/activate
cd TIMBER/
python3 BBTo2b4tau.py testfile_SIGNAL_2022.txt 0 0 2022
'''


import ROOT
from ROOT import TFile
import sys, os
import gc

gc.disable()

#from TIMBER.Tools.RestFramesHandler import load_restframes
import correctionlib
correctionlib.register_pyroot_binding()

sys.path.append('../../')
sys.path.append('../../../')

# ------------------ Command Line Arguments and Parsing -------------------
inputFiles = sys.argv[1] #fileList
# Are all the files running if 0 and 8?
testNum1 = sys.argv[2]   #first file in the list to use 
testNum2 = sys.argv[3]   #last file in the list to use
year = sys.argv[4]       #2022, 2022EE, 2023, 2023BPix

# Make the New .txt file from line testNum1 to testNum2 because TIMBER can handle .txt of .root's files   
print(f"Input File Path: {inputFiles}")
with open(inputFiles) as fp:
  lines = fp.readlines()

start = int(testNum1)
end = int(testNum2)

print(f"TestNum 1: {start} and TestNum 2: {end}")

print("Adding files to trimmed_input.txt")
#file_name = 'trimmed_input_'+inputFiles.replace('.txt','')+'_'+str(testNum1)+'.txt'
#listFiles = open(file_name, 'w')
  #if (listFiles.is_open())
filelist = []
for i, line in enumerate(lines): #, start=1):
  #if i in [start, end]:
  if i in range(start,end+1):
    #listFiles.write(line)    
    #print(line)
    filelist.append(line.strip())
#listFiles.close()
print(filelist)

print("Number of Entries:",len(filelist))
print("list contents:",filelist)
#sampleName = lines[start]
sampleName = filelist[0]
print(f"Sample Name: {sampleName}")

# Parse the incoming file names to assign labels  
isSig = ("Bprime" in sampleName)
isMadgraphBkg = (("QCD" in sampleName) or ("madgraphMLM" in sampleName))
isTOP = (("Mtt" in sampleName) or ("ST" in sampleName) or ("ttZ" in sampleName) or ("ttW" in sampleName) or ("ttH" in sampleName) or ("TTTo" in sampleName))
isTT = (("TT_Tune" in sampleName) or ("Mtt" in sampleName) or ("TTTo" in sampleName))
isVV = (("WW_" in sampleName) or ("WZ_" in sampleName) or ("ZZ_" in sampleName))
isSM = ("Muon" in sampleName)
isSE = (("SingleElectron" in sampleName) or ("EGamma" in sampleName))
isMC = not (("Single" in sampleName) or ("Muon" in sampleName) or ("EGamma" in sampleName) or ("MuonEG" in sampleName) or ("Tau/" in sampleName))
dataset = int(-1)
if "MuonEG" in sampleName: dataset = 2
elif ("Muon" in sampleName) or ("Muon0" in sampleName) or ("Muon1" in sampleName): dataset = 1
elif ("EGamma" in sampleName) or ("EGamma0" in sampleName) or ("EGamma1" in sampleName): dataset = 3
elif "Tau/" in sampleName: dataset = 4


#'root://cms-xrd-global.cern.ch//store/data/Run2018A/SingleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v2/2550000/28FF17A8-95EB-FD41-A55B-2EFAF2D6AF91.root' 
tokens = sampleName.split("/")
sample = tokens[7] # was 5
era = ''
ver = ''
if not isMC:
  runera = tokens[6] # was 4
  process = tokens[9] # was 7
  era = runera[-1] # last char
  print(process)
  ver = process[process.find('_')+1:process.find('_')+3]
del tokens

jecera = ''
if not isMC:
  jecera = era
  if year == '2022':
    jecera = 'CD'
  elif year == '2023':
    if(ver != 'v4'):
      jecera = 'Cv123'
    else:
      jecera = 'Cv4'
    
if isMC:
  if (("_ext1" in sampleName)): era = "ext1"
  if (("_ext2" in sampleName)): era = "ext2"
  if (("_ext3" in sampleName)): era = "ext3"

region = "Signal"
if isTT:
  region = "TTbar" # TPrimeTPrime or BPrimeBPrime
elif not isSig:
  region = "DataBkg"

print('============== RUN SETTINGS ===============')
print('isMC = ',isMC)
print('isSig = ',isSig)
print('sample = ',sample)
print('era = ',era)
print('jecera = ',jecera)
print('ver = ',ver)
print('region = ',region)
print('dataset = ', dataset)

# ------------------ TIMBER Analyzer inputs ------------------

num_threads = 1

# Import the C++
CompileCpp('TIMBER/Framework/include/common.h') # Compile (via gInterpreter) commonly used c++ code
CompileCpp('TIMBER/Framework/Tprime1lep/cleanjet.cc') # Compile Our vlq c++ code
CompileCpp('TIMBER/Framework/Tprime1lep/utilities.cc') # Compile Our vlq c++ code
CompileCpp('TIMBER/Framework/Tprime1lep/lumiMask.cc')
CompileCpp('TIMBER/Framework/Tprime1lep/selfDerived_corrs.cc')
CompileCpp('TIMBER/Framework/Tprime1lep/corr_funcs.cc') 
CompileCpp('TIMBER/Framework/Tprime1lep/topographInput.cc') 
CompileCpp('TIMBER/Framework/Tprime1lep/manualreco.cc') 
ROOT.gInterpreter.ProcessLine('#include "TString.h"')

# Enable using 4 threads
ROOT.ROOT.EnableImplicitMT(num_threads)

# load rest frames handler
#handler_name = 'Bprime_handler.cc'
#class_name = 'gc_Container'
#load_restframes(num_threads, handler_name, class_name, 'rfc')

# ------------------ Important Variables ------------------
debug = False

ROOT.gInterpreter.Declare("""
  string year = \"""" + year + """\"; 
  string sample = \"""" + sample + """\";
  string jecera = \"""" + jecera + """\";
  string region = \"""" + region + """\";
  string ver = \"""" + ver + """\";

  bool isMC = """+str(isMC).lower()+"""; 
  bool debug = """+str(debug).lower()+""";
  bool isSig = """+str(isSig).lower()+""";
  int dataset = """+str(dataset)+""";
  """)

def analyze(jesvar):
  ROOT.gInterpreter.ProcessLine('string jesvar = "' + jesvar + '"; ')

  # Create analyzer instance
  # is filelist still what you want it be here
  a = analyzer(filelist)
  
  print('==========================INITIALIZED ANALYZER========================')
  
  # ------------------ Golden JSON Data ------------------
  # change the jsonfile path to somewhere they have it in TIMBER
  jsonfile = "./TIMBER/data/LumiJSON/"
  if '2022' in year:
    jsonfile = jsonfile + "Cert_Collisions2022_355100_362760_Golden.json"
  elif '2023' in year:
    jsonfile = jsonfile + "Cert_Collisions2023_366442_370790_Golden.json"
  else:
    print(f'ERROR: Can\'t parse the year to assign a golden json file. Expected 2022(EE) or 2023(BPix). Got: {year}\n')
    
  ROOT.gInterpreter.Declare("""
    const auto myLumiMask = lumiMask::fromJSON(\"""" + jsonfile + """\");
  """)

  print('========= loaded lumimask ============')
  
  # ------------------ Self-derived corrections ------------------

  #TODO when we determine what Run 3 scale factors need to be computed by us, load them up here

  ROOT.gInterpreter.ProcessLine('initialize(year);')


  # ------------------ correctionsLib corrections ------------------

  # mutrig = "OldMu100_or_TkMu100"
  #deepjetL = {'2022':0.0583,'2022EE':0.0614,'2023':0.0479,'2023BPix':0.0683} #ParT
  PNetL = {'2022':0.047,'2022EE':0.0499,'2023':0.0358,'2023BPix':0.0359} #PN
  #deepjetL = {'2022':0.0583,'2022EE':0.0614,'2023':0.0479,'2023BPix':0.048} #DeepFlav
  deepjetM = {'2022':0.3086,'2022EE':0.3196,'2023':0.2431,'2023BPix':0.2435}
  yrstr = {'2022':"2022_Summer22",'2022EE':"2022_Summer22EE",'2023':"2023_Summer23",'2023BPix':"2023_Summer23BPix"}
  jecyr = {'2022':"Summer22_22Sep2023",'2022EE':"Summer22EE_22Sep2023",'2023':"Summer23Prompt23",'2023BPix':"Summer23BPixPrompt23"}
  if not isMC:
      jecyr = {'2022':"Summer22_22Sep2023_RunCD",'2022EE':"Summer22EE_22Sep2023_Run"+jecera,'2023':"Summer23Prompt23",'2023BPix':"Summer23BPixPrompt23"}
  jeryr = {'2022':"Summer22_22Sep2023",'2022EE':"Summer22EE_22Sep2023",'2023':"Summer23Prompt23_RunCv1234",'2023BPix':"Summer23BPixPrompt23_RunD"}
  jecver = {'2022':"V2",'2022EE':"V2",'2023':"V2",'2023BPix':"V3"}
  puname = {'2022':"Collisions2022_355100_357900_eraBCD_GoldenJson",'2022EE':"Collisions2022_359022_362760_eraEFG_GoldenJson",'2023':"Collisions2023_366403_369802_eraBC_GoldenJson",'2023BPix':"Collisions2023_369803_370790_eraD_GoldenJson"}
  jetvetoname = {'2022':"Summer22_23Sep2023_RunCD_V1",'2022EE':"Summer22EE_23Sep2023_RunEFG_V1",'2023':"Summer23Prompt23_RunC_V1",'2023BPix':"Summer23BPixPrompt23_RunD_V1"}
  elecyr = {'2022':"2022Re-recoBCD",'2022EE':"2022Re-recoE+PromptFG",'2023':"2023PromptC",'2023BPix':"2023PromptD"}
  tauyr = {'2022':"2022_preEE",'2022EE':"2022_postEE",'2023':"2023_preBPix",'2023BPix':"2023_postBPix"} 
  METyr = {'2022':"2022",'2022EE':"2022EE",'2023':"2023",'2023BPix':"2023BPix"} 
  METsimpleyr = {'2022':"2022",'2022EE':"2022",'2023':"2023",'2023BPix':"2023"} 
  btagname = {'2022':"particleNet_comb",'2022EE':"particleNet_comb",'2023':"deepJet_comb",'2023BPix':"deepJet_comb"}

  print(jecyr[year]+"_"+jecver[year]+"_DATA_L1L2L3Res_AK4PFPuppi")
 
  ROOT.gInterpreter.Declare("""
  float PNetL = """+str(PNetL[year])+""";
  float deepjetM = """+str(deepjetM[year])+""";
  string yrstr = \""""+yrstr[year]+"""\";
  string jecyr = \""""+jecyr[year]+"""\";
  string jeryr = \""""+jeryr[year]+"""\";
  string jecver = \""""+jecver[year]+"""\";
  string puname = \""""+puname[year]+"""\";
  string jetvetoname = \""""+jetvetoname[year]+"""\";
  string elecyr = \""""+elecyr[year]+"""\";
  string tauyr = \""""+tauyr[year]+"""\";
  string METyr = \""""+METyr[year]+"""\";
  string METsimpleyr = \""""+METsimpleyr[year]+"""\";
  string btagname = \""""+btagname[year]+"""\";
  """)

  # *************** muonisocorr does not match the muon definition !!!!! (change to mediumPFIso hopefully) *****************
  ROOT.gInterpreter.Declare("""
  auto pileupcorrset = correction::CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/LUM/"+yrstr+"/puWeights.json.gz");
  auto btagcorrset = correction::CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/"+yrstr+"/btagging.json.gz");
  auto jetvetocorrset = correction::CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/"+yrstr+"/jetvetomaps.json.gz");
  auto electroncorrset = correction::CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/"+yrstr+"/electron.json.gz");
  auto muoncorrset = correction::CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/"+yrstr+"/muon_Z.json.gz");
  auto taucorrset = correction::CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/TAU/"+yrstr+"/tau_DeepTau2018v2p5_"+tauyr+".json.gz");
  auto METcorrset = correction::CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/"+yrstr+"/met_xyCorrections_"+METsimpleyr+"_"+METyr+".json.gz");

  auto pileupcorr = pileupcorrset->at(puname);
  auto btagwpbccorr = btagcorrset->at(btagname);
  auto btagwplcorr = btagcorrset->at("particleNet_light");
  auto jetvetocorr = jetvetocorrset->at(jetvetoname);
  auto electroncorr = electroncorrset->at("Electron-ID-SF");
  auto muonidcorr = muoncorrset->at("NUM_MediumID_DEN_TrackerMuons");
  auto muonisocorr = muoncorrset->at("NUM_LoosePFIso_DEN_MediumID"); 
  auto tauidVSecorr = taucorrset->at("DeepTau2018v2p5VSe"); 
  auto tauidVSmucorr = taucorrset->at("DeepTau2018v2p5VSmu");
  auto tauidVSjetcorr = taucorrset->at("DeepTau2018v2p5VSjet");
  auto METcorr = METcorrset->at("met_xy_corrections");
  """)

  if not isMC:
    ROOT.gInterpreter.Declare("""
    auto ak4corrset = correction::CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/"+yrstr+"/jet_jerc.json.gz"); 
    auto ak8corrset = correction::CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/"+yrstr+"/fatJet_jerc.json.gz"); 

    auto ak4corr = ak4corrset->compound().at(jecyr+"_"+jecver+"_DATA_L1L2L3Res_AK4PFPuppi");
    auto ak4corrL1 = ak4corrset->at(jecyr+"_"+jecver+"_DATA_L1FastJet_AK4PFPuppi");
    auto ak8corr = ak8corrset->compound().at(jecyr+"_"+jecver+"_DATA_L1L2L3Res_AK8PFPuppi");
    """)
  else:
    ROOT.gInterpreter.Declare("""
    auto ak4corrset = correction::CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/"+yrstr+"/jet_jerc.json.gz"); 
    auto ak8corrset = correction::CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/"+yrstr+"/fatJet_jerc.json.gz"); 

    auto ak4corr = ak4corrset->compound().at(jecyr+"_"+jecver+"_MC_L1L2L3Res_AK4PFPuppi");
    auto ak4corrL1 = ak4corrset->at(jecyr+"_"+jecver+"_MC_L1FastJet_AK4PFPuppi");
    auto ak4corrUnc = ak4corrset->at(jecyr+"_"+jecver+"_MC_Total_AK4PFPuppi");
    auto ak4ptres = ak4corrset->at(jeryr+"_JRV1_MC_PtResolution_AK4PFPuppi");
    auto ak4jer = ak4corrset->at(jeryr+"_JRV1_MC_ScaleFactor_AK4PFPuppi");
    auto ak8corr = ak8corrset->compound().at(jecyr+"_"+jecver+"_MC_L1L2L3Res_AK8PFPuppi");
    auto ak8corrUnc = ak8corrset->at(jecyr+"_"+jecver+"_MC_Total_AK8PFPuppi");
    """)
    
 
 # ------------------ Flag Cuts ------------------
  flagCuts = CutGroup('FlagCuts')
  flagCuts.Add('Bad Event Filters', 'Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && Flag_goodVertices == 1 && Flag_eeBadScFilter == 1 && Flag_globalSuperTightHalo2016Filter == 1 && Flag_BadPFMuonFilter == 1 && Flag_BadPFMuonDzFilter == 1')
  flagCuts.Add('Event has jets', 'nJet > 0') # need jets   && nFatJet > 0
  
  # ------------------ Golden JSON (Data) || GEN Info (MC) ------------------
  gjsonVars = VarGroup('GoldenJsonVars')
  gjsonCuts = CutGroup('GoldenJsonCuts')
  if not isMC: # apply golden json to data
    gjsonVars.Add("passesJSON", "goldenjson(myLumiMask, run, luminosityBlock)") # function name, parameters are branches in the file or defined objects
    gjsonCuts.Add("Data passes Golden JSON", "passesJSON == 1") 
  else:
    gjsonVars.Add("PileupWeights", "pufunc(pileupcorr, Pileup_nTrueInt)")


  #### Let's keep this around in case we decide to bring back the idea of leptons in the future
  #lVars.Add("Electron_cutBasedIdNoIso_tight", "Electron_cutBasedIdNoIso_tight(nElectron, Electron_vidNestedWPBitmap)")
  # lVars.Add("TPassMu", "abs(Muon_eta)<2.4 && Muon_mediumId==1 && Muon_miniIsoId>=3 && abs(Muon_dz) < 0.5 && Muon_dxy < 0.2")
  # lVars.Add("TPassEl", "(abs(Electron_eta)<1.442 || (abs(Electron_eta)>1.566 && abs(Electron_eta)<2.5)) && Electron_cutBasedIdNoIso_tight==1 && Electron_miniPFRelIso_all<0.1")
  # lVars.Add("VetoMu", "TPassMu && (Muon_pt>25)")
  # lVars.Add("VetoEl", "TPassEl && (Electron_pt>25)")
  # lVars.Add("SignalIsoMu", "TPassMu && (Muon_pt>=55)")
  # lVars.Add("SignalIsoEl", "TPassEl && (Electron_pt>=55)")
  # lVars.Add("nVetoLep", "(int) (Sum(VetoMu)+Sum(VetoEl))")
  # lVars.Add("Muon_P4", "fVectorConstructor(Muon_pt,Muon_eta,Muon_phi,Muon_mass)")
  # lVars.Add("SMuon_P4", "fVectorConstructor(SMuon_pt,SMuon_eta,SMuon_phi,SMuon_mass)")
  # lVars.Add("SElectron_P4", "fVectorConstructor(SElectron_pt,SElectron_eta,SElectron_phi,SElectron_mass)")
  # lVars.Add("SMuon_jetIdx", "Muon_jetIdx[SignalIsoMu == true]")
  # lVars.Add("SElectron_jetIdx", "Electron_jetIdx[SignalIsoEl]")
  # lVars.Add("nSignalIsoMu", "(int) Sum(SignalIsoMu)")
  # lVars.Add("nSignalIsoEl", "(int) Sum(SignalIsoEl)")
  # lVars.Add("VetoIsoMu", "(VetoMu == true && Muon_pt < 55)")
  # lVars.Add("VetoIsoEl", "(VetoEl == true && Electron_pt < 55)")
  # lVars.Add("nVetoIsoLep", "(int) (Sum(VetoIsoMu)+Sum(VetoIsoEl))")

  # lVars.Add("isMu", "(nMuon>0) && (HLT_Mu50 || HLT_HighPtTkMu100) && (nSignalIsoMu==1) && (nVetoIsoLep==0) && (nElectron == 0 || nSignalIsoEl == 0)")
  # lVars.Add("isEl", "(nElectron>0) && (HLT_Ele115_CaloIdVT_GsfTrkIdT || HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 || HLT_Photon200) && (nSignalIsoEl==1) && (nVetoIsoLep==0) && (nMuon == 0 || nSignalIsoMu == 0)")

  # lCuts = CutGroup('LeptonCuts')
  # lCuts.Add("Event is either muon or electron", "isMu || isEl")
  # lVars.Add("assignleps", "assign_leps(isMu,isEl,SignalIsoMu,SignalIsoEl,Muon_pt,Muon_eta,Muon_phi,Muon_mass,Muon_miniPFRelIso_all,Electron_pt,Electron_eta,Electron_phi,Electron_mass,Electron_miniPFRelIso_all)")
  # lVars.Add("lepton_pt","assignleps[0]")
  # lVars.Add("lepton_eta","assignleps[1]") 
  # lVars.Add("lepton_phi","assignleps[2]")
  # lVars.Add("lepton_mass","assignleps[3]")
  # lVars.Add("lepton_miniIso","assignleps[4]")
  

  # ------------------ TAU Definitions ------------------
  # Count N taus per event of di:fferent types
  # (Command line? figure out typical tau momentum range)
  # What % of signal events pass criteria like "I have N taus of type A with momentum >= B". change N, change A, change B (within reason)
  #
  # Goal: maximize the signal efficiency without choose like N >= 0, A >= loosest thing, B >= 0
  # Goal: provide a ROOT file (or stage in this TIMBER analyzer) that has the NN inputs we would want for a test

  tVars = VarGroup('TauVars')
  
  #tVars.Add('matchedTau','Tau_pt > 20 && abs(Tau_eta) < 2.5 && abs(Tau_dz) < 0.2 && Tau_genPartFlav == 5')  
  #tVars.Add('NmatchedTau', 'Sum(matchedTau)')
 
  # these comments were the first choices, before the scale factor considerations
  #vs mu >= 1 
  tVars.Add('goodTaumu', 'Tau_pt > 20 && abs(Tau_eta) < 2.5 && abs(Tau_dz) < 0.2 && Tau_idDeepTau2018v2p5VSmu >= 1')
  tVars.Add('NgoodTaumu', 'Sum(goodTaumu)')
  #vs e >= 3
  tVars.Add('goodTaue', 'goodTaumu == 1 && Tau_idDeepTau2018v2p5VSe >= 2')
  tVars.Add('NgoodTaue', 'Sum(goodTaue)')
  #vs jet >= 2
  tVars.Add('goodTau', 'goodTaue == 1 && Tau_idDeepTau2018v2p5VSjet >= 2')
  tVars.Add('NgoodTau', 'Sum(goodTau)')

  tVars.Add('GoodTau_eta', 'Tau_eta[goodTau == true]')
  tVars.Add('GoodTau_phi', 'Tau_phi[goodTau == true]')
  tVars.Add('GoodTau_pt', 'Tau_pt[goodTau == true]')
  tVars.Add('GoodTau_mass', 'Tau_mass[goodTau == true]')
  tVars.Add('GoodTau_ID', 'ROOT::VecOps::RVec<int>(NgoodTau, 15)')
  tVars.Add('GoodTau_DM', 'Tau_decayMode[goodTau == true]')
  if (isMC):
     tVars.Add('GoodTau_genmch', 'Tau_genPartFlav[goodTau == true]')
  tVars.Add('GoodTau_charge', 'Tau_charge[goodTau == true]')

  #Tau Cuts
  tCuts = CutGroup('TauCuts')
  #tCuts.Add('genHadrTauCount is 1', 'genHadrTauCount[0] == 1')
  #tCuts.Add('matched Tau', 'NmatchedTau > 0')
  #tCuts.Add('good Taumu', 'NgoodTaumu > 0')
  #tCuts.Add('good Taue', 'NgoodTaue > 0')
  #tCuts.Add('good Tau (jet)', 'NgoodTau > 0')
  
  # -------------------- EL and MUON Definitions --------------------
  eandmuVars = VarGroup('ElandMuVars')

  #Good Electrons
  eandmuVars.Add('goodElectrons', 'Electron_pt > 10 && abs(Electron_eta) < 2.5 && (Electron_mvaIso_WP80 == 1)')
  eandmuVars.Add('NgoodElecs', 'Sum(goodElectrons)')
  
  #eandmuVars.Add('matchedElec','Electron_pt > 10 && abs(Electron_eta) < 2.5') # && abs(Electron_genPartFlav) == 15')  
  #eandmuVars.Add('NmatchedElec', 'Sum(matchedElec)')
  #eandmuVars.Add('goodElectrons', 'matchedElec == 1 && (Electron_cutBased >= 3 || Electron_cutBased >= 4) && ((abs(Electron_eta) <= 1.479 && abs(Electron_dz) < 0.10 && abs(Electron_dxy) < 0.05 ) || (abs(Electron_eta) > 1.479 && abs(Electron_dz) < 0.20 && abs(Electron_dxy) < 0.10))') 
   
  eandmuVars.Add('GoodEl_pt','Electron_pt[goodElectrons == true]')
  eandmuVars.Add('GoodEl_eta','Electron_eta[goodElectrons == true]')
  eandmuVars.Add('GoodEl_phi','Electron_phi[goodElectrons == true]')
  eandmuVars.Add('GoodEl_mass','Electron_mass[goodElectrons == true]')
  eandmuVars.Add('GoodEl_ID', 'ROOT::VecOps::RVec<int>(NgoodElecs, 11)')
  eandmuVars.Add('GoodEl_mtforTau', 'ROOT::VecOps::RVec<int>(NgoodElecs, -1)')
  eandmuVars.Add('GoodEl_charge', 'Electron_charge[goodElectrons == true]')

  #Good Muons
  eandmuVars.Add('goodMuons', 'Muon_pt >= 15 && abs(Muon_eta) < 2.4 && Muon_mediumId == 1 && abs(Muon_dxy) < 0.2 && abs(Muon_dz) < 0.5 && Muon_pfIsoId >= 3')
  eandmuVars.Add('NgoodMuons', 'Sum(goodMuons)')
  
  #eandmuVars.Add('matchedMuon','Muon_pt > 10 && abs(Muon_eta) < 2.4') # && abs(Muon_genPartFlav) == 15')  
  #eandmuVars.Add('NmatchedMuon', 'Sum(matchedMuon)')
  #eandmuVars.Add('goodMuons', 'matchedMuon == 1 && Muon_mvaTTH > 0.85') #note: will become Muon_promptMVA in Nanov15):
 
  eandmuVars.Add('GoodMu_pt','Muon_pt[goodMuons == true]')
  eandmuVars.Add('GoodMu_eta','Muon_eta[goodMuons == true]')
  eandmuVars.Add('GoodMu_phi','Muon_phi[goodMuons == true]')
  eandmuVars.Add('GoodMu_mass','Muon_mass[goodMuons == true]')
  eandmuVars.Add('GoodMu_ID', 'ROOT::VecOps::RVec<int>(NgoodMuons, 13)')
  eandmuVars.Add('GoodMu_mtforTau', 'ROOT::VecOps::RVec<int>(NgoodMuons, -1)')
  eandmuVars.Add('GoodMu_charge', 'Muon_charge[goodMuons == true]')

  #test filters
  testCuts = CutGroup('Test Cuts')
  #testCuts.Add('NmatchedElec filter','NmatchedElec > 0')
  #testCuts.Add('NgoodElecs filter', 'NgoodElecs > 0')
  #testCuts.Add('NmatchedMuon filter', 'NmatchedMuon > 0')
  #testCuts.Add('NgoodMuons filter', 'NgoodMuons > 0')
  
  #Gen Cuts
  GenCuts = CutGroup('Gen Cuts')
  #GenCuts.Add('genElecCount is 1 and pT>10', 'genElecCount == 1')
  #GenCuts.Add('genMuonCount is 1 and pT>10', 'genMuonCount == 1')
 
  # ------------------------ LEPTON Definitions and Selection -------------------------
  lVars = VarGroup('LeptonVars')

  lVars.Add('igLepton', 'ROOT::VecOps::Concatenate(goodElectrons, goodMuons)')
  lVars.Add('gLepton', 'ROOT::VecOps::Concatenate(igLepton, goodTau)')
  lVars.Add('NgoodLeptons', 'Sum(gLepton)')  

  lVars.Add('iLepton_pt', 'ROOT::VecOps::Concatenate(GoodEl_pt, GoodMu_pt)')
  lVars.Add('Lepton_pt', 'ROOT::VecOps::Concatenate(iLepton_pt, GoodTau_pt)')
 
  lVars.Add('iLepton_eta', 'ROOT::VecOps::Concatenate(GoodEl_eta, GoodMu_eta)')
  lVars.Add('Lepton_eta', 'ROOT::VecOps::Concatenate(iLepton_eta, GoodTau_eta)')

  lVars.Add('iLepton_phi', 'ROOT::VecOps::Concatenate(GoodEl_phi, GoodMu_phi)')
  lVars.Add('Lepton_phi', 'ROOT::VecOps::Concatenate(iLepton_phi, GoodTau_phi)')

  lVars.Add('iLepton_mass', 'ROOT::VecOps::Concatenate(GoodEl_mass, GoodMu_mass)')
  lVars.Add('Lepton_mass', 'ROOT::VecOps::Concatenate(iLepton_mass, GoodTau_mass)')

  lVars.Add('iLepton_ID', 'ROOT::VecOps::Concatenate(GoodEl_ID, GoodMu_ID)')
  lVars.Add('Lepton_ID', 'ROOT::VecOps::Concatenate(iLepton_ID, GoodTau_ID)')
  
  lVars.Add('iLepton_charge', 'ROOT::VecOps::Concatenate(GoodEl_charge, GoodMu_charge)')
  lVars.Add('Lepton_charge', 'ROOT::VecOps::Concatenate(iLepton_charge, GoodTau_charge)')
  
  lVars.Add('iLepton_mtElMu', 'ROOT::VecOps::Concatenate(GoodEl_mtforTau, GoodMu_mtforTau)')
  lVars.Add('Lepton_TauDM', 'ROOT::VecOps::Concatenate(iLepton_mtElMu, GoodTau_DM)')
  if (isMC):
     lVars.Add('Lepton_TauMch', 'ROOT::VecOps::Concatenate(iLepton_mtElMu, GoodTau_genmch)')
  
  lVars.Add("Lepton_ptargsort","ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(Lepton_pt))")
  lVars.Add("GoodLepton_pt","reorder(Lepton_pt,Lepton_ptargsort)")
  lVars.Add("GoodLepton_eta","reorder(Lepton_eta,Lepton_ptargsort)")
  lVars.Add("GoodLepton_phi","reorder(Lepton_phi,Lepton_ptargsort)")
  lVars.Add("GoodLepton_mass","reorder(Lepton_mass,Lepton_ptargsort)")
  lVars.Add("GoodLepton_ID","reorder(Lepton_ID,Lepton_ptargsort)")
  lVars.Add("GoodLepton_TauDM","reorder(Lepton_TauDM,Lepton_ptargsort)")
  if (isMC):
     lVars.Add("GoodLepton_TauMch","reorder(Lepton_TauMch,Lepton_ptargsort)")
  lVars.Add("GoodLepton_charge","reorder(Lepton_charge,Lepton_ptargsort)")
  
  lVars.Add("Good4Lepton_pt", "NgoodLeptons >= 4 ? ROOT::VecOps::Take(GoodLepton_pt, 4) : ROOT::VecOps::Take(GoodLepton_pt, 3)")
  lVars.Add("Good4Lepton_eta", "NgoodLeptons >= 4 ? ROOT::VecOps::Take(GoodLepton_eta, 4) : ROOT::VecOps::Take(GoodLepton_eta, 3)")
  lVars.Add("Good4Lepton_phi", "NgoodLeptons >= 4 ? ROOT::VecOps::Take(GoodLepton_phi, 4) : ROOT::VecOps::Take(GoodLepton_phi, 3)")
  lVars.Add("Good4Lepton_mass", "NgoodLeptons >= 4 ? ROOT::VecOps::Take(GoodLepton_mass, 4) : ROOT::VecOps::Take(GoodLepton_mass, 3)")
  lVars.Add("Good4Lepton_ID", "NgoodLeptons >= 4 ? ROOT::VecOps::Take(GoodLepton_ID, 4) : ROOT::VecOps::Take(GoodLepton_ID, 3)")
  lVars.Add("Good4Lepton_TauDM", "NgoodLeptons >= 4 ? ROOT::VecOps::Take(GoodLepton_TauDM, 4) : ROOT::VecOps::Take(GoodLepton_TauDM, 3)")
  if (isMC):
     lVars.Add("Good4Lepton_TauMch", "NgoodLeptons >= 4 ? ROOT::VecOps::Take(GoodLepton_TauMch, 4) : ROOT::VecOps::Take(GoodLepton_TauMch, 3)")
  lVars.Add("Good4Lepton_charge", "NgoodLeptons >= 4 ? ROOT::VecOps::Take(GoodLepton_charge, 4) : ROOT::VecOps::Take(GoodLepton_charge, 3)")

  lVars.Add("hasBosonishMass", "hasBosonishMfunc(Good4Lepton_pt, Good4Lepton_eta, Good4Lepton_phi, Good4Lepton_mass, Good4Lepton_ID, Good4Lepton_charge)")

  lCuts = CutGroup('Lepton Cuts')
  lCuts.Add('NgoodLeptons >= 3', 'NgoodLeptons >= 3')
  lCuts.Add('hasBosonishMass == 0', 'hasBosonishMass == 0')
  #lCuts.Add('NgoodLeptons >= 4', 'NgoodLeptons >= 4')
 
  # ----------------- corrections/SFs for el/mu/tau ----------------------
  lepSFs = VarGroup('Lepton Scale Factors')
  if (isMC):
    lepSFs.Add('elrecoSF', 'elrecofunc(electroncorr, elecyr, Good4Lepton_pt, Good4Lepton_eta, Good4Lepton_phi, Good4Lepton_ID)')
    lepSFs.Add('elidSF', 'elidfunc(electroncorr, elecyr, Good4Lepton_pt, Good4Lepton_eta, Good4Lepton_phi, Good4Lepton_ID)')
    lepSFs.Add('muonidSF', 'muidfunc(muonidcorr, Good4Lepton_pt, Good4Lepton_eta, Good4Lepton_ID)')
    lepSFs.Add('muonisoSF', 'muisofunc(muonisocorr, Good4Lepton_pt, Good4Lepton_eta, Good4Lepton_ID)')
    lepSFs.Add('tauidVSeSF', 'tauefunc(tauidVSecorr, Good4Lepton_eta, Good4Lepton_TauDM, Good4Lepton_TauMch, Good4Lepton_ID)')
    lepSFs.Add('tauidVSmuSF', 'taumufunc(tauidVSmucorr, Good4Lepton_eta, Good4Lepton_TauMch, Good4Lepton_ID)')
    lepSFs.Add('tauidVSjetSF', 'taujetfunc(tauidVSjetcorr, Good4Lepton_pt, Good4Lepton_TauDM, Good4Lepton_TauMch, Good4Lepton_ID)') 

  # -------------------------- TRIGGERS -------------------------------
  TrigVars = VarGroup('Triggers Vars')

  #should double check
  TrigVars.Add('NMu20', 'Sum(GoodMu_pt > 20)')
  TrigVars.Add('NMu10LT20', 'Sum(GoodMu_pt > 10 && GoodMu_pt <= 20)')
  TrigVars.Add('NMu23', 'Sum(GoodMu_pt > 25)')
  TrigVars.Add('NMu12LT23', 'Sum(GoodMu_pt > 15 && GoodMu_pt <= 25)')
  TrigVars.Add('NEle23', 'Sum(GoodEl_pt > 25)')
  TrigVars.Add('NEle12LT23', 'Sum(GoodEl_pt > 15 && GoodEl_pt <= 25)')
  TrigVars.Add('NEle30', 'Sum(GoodEl_pt > 33 && Electron_cutBased[goodElectrons == true] >= 4)')
  TrigVars.Add('NMu20forTau', 'Sum(GoodMu_pt > 23 && abs(GoodMu_eta) < 2.1)')
  TrigVars.Add('NEle24forTau', 'Sum(GoodEl_pt > 27 && Electron_cutBased[goodElectrons == true] >= 4 && abs(GoodEl_eta) < 2.1)') # cutbased >= 4 is tight
  TrigVars.Add('NTauMuTrig', 'Sum(GoodTau_pt > 30 && Tau_idDeepTau2018v2p5VSjet[goodTau == true] >= 4 && abs(GoodTau_eta) < 2.1)')
  TrigVars.Add('NTauEleTrig', 'Sum(GoodTau_pt > 33 && Tau_idDeepTau2018v2p5VSjet[goodTau == true] >= 4 && abs(GoodTau_eta) < 2.1)')
  TrigVars.Add('NTauMedTrig', 'Sum(GoodTau_pt > 38 && Tau_idDeepTau2018v2p5VSjet[goodTau == true] >= 5 && abs(GoodTau_eta) < 2.1)')
  TrigVars.Add('NTauLooseTrig', 'Sum(GoodTau_pt > 185 && Tau_idDeepTau2018v2p5VSjet[goodTau == true] >= 4 && abs(GoodTau_eta) < 2.1)')

  TrigVars.Add('isMuMuTrig', '(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 == 1) && (NMu20 >= 2 || (NMu20 >= 1 && NMu10LT20 >= 1))') # in muon dataset
  TrigVars.Add('isMuElTrig', '''(((HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ == 1) && (NMu23 >= 1 && (NEle23 >= 1 || NEle12LT23 >= 1))) 
                                || ((HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ == 1) && (NEle23 >= 1 && (NMu23 >= 1 || NMu12LT23 >= 1))))''') # in muonEG dataset
  TrigVars.Add('isMuTauTrig', '(HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1 == 1) && (NMu20forTau >= 1 && NTauMuTrig >= 1)') # in muon dataset
  TrigVars.Add('isElElTrig', '(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL == 1) && (NEle23 >= 2 || (NEle23 >= 1 && NEle12LT23 >= 1))') # in Egamma dataset
  TrigVars.Add('isElTauTrig', '(HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1 == 1) && (NEle24forTau >= 1 && NTauEleTrig >= 1)') # in Egamma dataset
  TrigVars.Add('isTauTauTrig', '(HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1 == 1) && (NTauMedTrig >= 2)') # in tau dataset
  TrigVars.Add('isSingleMuTrig', '(HLT_IsoMu24 == 1) && (NMu23 >= 1)') # in muon dataset
  TrigVars.Add('isSingleEleTrig', '(HLT_Ele30_WPTight_Gsf == 1) && (NEle30 >= 1)') # in Egamma dataset
  TrigVars.Add('isSingleTauTrig', '(HLT_LooseDeepTauPFTauHPS180_L2NN_eta2p1 == 1) && (NTauLooseTrig >= 1)') # in tau dataset

  TrigVars.Add('NlessRestrict', 'isMuMuTrig == 1 || isMuElTrig == 1 || isElElTrig == 1')
  TrigVars.Add('NwithElMu_Tau', 'isMuMuTrig == 1 || isMuElTrig == 1 || isElElTrig == 1 || isMuTauTrig == 1 || isElTauTrig == 1')
  TrigVars.Add('NwithDiTau', 'isMuMuTrig == 1 || isMuElTrig == 1 || isElElTrig == 1 || isTauTauTrig == 1')
  TrigVars.Add('NwithAll', 'isMuMuTrig == 1 || isMuElTrig == 1 || isElElTrig == 1 || isMuTauTrig == 1 || isElTauTrig == 1 || isTauTauTrig == 1')
  TrigVars.Add('Nnone', 'isMuMuTrig == 0 && isMuElTrig == 0 && isElElTrig == 0 && isMuTauTrig == 0 && isElTauTrig == 0 && isTauTauTrig == 0')
  
  TrigVars.Add('passesMuPD', '(isMuMuTrig == 1 || isSingleMuTrig == 1 || isMuTauTrig == 1) && (isMC == 1 || dataset == 1)')
  TrigVars.Add('passesMuEGPD', '(isMuElTrig == 1) && (isMC == 1 || dataset == 2)')
  TrigVars.Add('passesEGPD', '(isElElTrig == 1 || isSingleEleTrig == 1 || isElTauTrig == 1) && (isMC == 1 || dataset == 3)')
  TrigVars.Add('passesTauPD', '(isTauTauTrig == 1 || isSingleTauTrig == 1) && (isMC == 1 || dataset == 4)')

  TrigCuts = CutGroup('Trigger Cuts')
  TrigCuts.Add('Dataset Filter', 'passesMuPD == 1 || passesMuEGPD == 1 || passesEGPD == 1 || passesTauPD == 1') #getting errors?
  #TrigCuts.Add('NlessRestrict', 'isMuMuTrig == 1 || isMuElTrig == 1 || isElElTrig == 1')
  #TrigCuts.Add('NwithElMu_Tau', 'isMuMuTrig == 1 || isMuElTrig == 1 || isElElTrig == 1 || isMuTauTrig == 1 || isElTauTrig == 1')
  #TrigCuts.Add('NwithDiTau', 'isMuMuTrig == 1 || isMuElTrig == 1 || isElElTrig == 1 || isTauTauTrig == 1')
  #TrigCuts.Add('NwithAll', 'isMuMuTrig == 1 || isMuElTrig == 1 || isElElTrig == 1 || isMuTauTrig == 1 || isElTauTrig == 1 || isTauTauTrig == 1')
  #TrigCuts.Add('Nnone', 'isMuMuTrig == 0 && isMuElTrig == 0 && isElElTrig == 0 && isMuTauTrig == 0 && isElTauTrig == 0 && isTauTauTrig == 0')


  # ------------------ JET Cleaning and JERC ------------------
  jVars = VarGroup('JetCleaningVars')
  
  jVars.Add("Jet_P4", "fVectorConstructor(Jet_pt,Jet_eta,Jet_phi,Jet_mass)")
  #jVars.Add("FatJet_P4", "fVectorConstructor(FatJet_pt,FatJet_eta,FatJet_phi,FatJet_mass)")
  jVars.Add("Jet_EmEF","Jet_neEmEF + Jet_chEmEF")
  jVars.Add("DummyZero","float(0.0)")
  
  if isMC:          #TODO Can we use overloading even more to remove data/MC distinction?
    jVars.Add("GenJet_P4","fVectorConstructor(GenJet_pt,GenJet_eta,GenJet_phi,GenJet_mass)")
    jVars.Add("cleanedJets", "cleanJetsMC(debug,year,jesvar,ak4corr,ak4corrL1,ak4corrUnc,ak4ptres,ak4jer,ak8corr,ak8corrUnc,Jet_P4,Jet_rawFactor,Jet_muonSubtrFactor,Jet_area,Jet_EmEF,Jet_jetId,GenJet_P4,Jet_genJetIdx,Rho_fixedGridRhoFastjetAll,DummyZero,DummyZero)") # muon and EM factors unused in this call
    jVars.Add("cleanMets", "cleanJetsMC(debug,year,jesvar,ak4corr,ak4corrL1,ak4corrUnc,ak4ptres,ak4jer,ak8corr,ak8corrUnc,Jet_P4,Jet_rawFactor,Jet_muonSubtrFactor,Jet_area,Jet_EmEF,Jet_jetId,GenJet_P4,Jet_genJetIdx,Rho_fixedGridRhoFastjetAll,RawPuppiMET_pt,RawPuppiMET_phi)") # lepton args are unused in this call
    #jVars.Add("GenJetAK8_P4", "fVectorConstructor(GenJetAK8_pt,GenJetAK8_eta,GenJetAK8_phi,GenJetAK8_mass)")
    #jVars.Add("cleanFatJets", "cleanJetsMC(debug,year,jesvar,ak4corr,ak4corrL1,ak4corrUnc,ak4ptres,ak4jer,ak8corr,ak8corrUnc,FatJet_P4,FatJet_rawFactor,FatJet_rawFactor,FatJet_area,FatJet_area,FatJet_jetId,GenJetAK8_P4,FatJet_genJetAK8Idx,Rho_fixedGridRhoFastjetAll,DummyZero,DummyZero)") # args 12 and 14 are dummies
  else:
    jVars.Add("cleanedJets", "cleanJetsData(run,debug,year,ak4corr,ak4corrL1,ak8corr,Jet_P4,Jet_rawFactor,Jet_muonSubtrFactor,Jet_area,Jet_EmEF,Jet_jetId,Jet_P4,Jet_jetId,Rho_fixedGridRhoFastjetAll,DummyZero,DummyZero)") # muon and EM factors unused in this call, args 16-17 are dummies
    jVars.Add("cleanMets", "cleanJetsData(run,debug,year,ak4corr,ak4corrL1,ak8corr,Jet_P4,Jet_rawFactor,Jet_muonSubtrFactor,Jet_area,Jet_EmEF,Jet_jetId,Jet_P4,Jet_jetId,Rho_fixedGridRhoFastjetAll,RawPuppiMET_pt,RawPuppiMET_phi)") # lepton args unused in this call, args 16-17 are dummies
    #jVars.Add("cleanFatJets", "cleanJetsData(debug,year,ak4corr,ak4corrL1,ak8corr,FatJet_P4,FatJet_rawFactor,FatJet_rawFactor,FatJet_area,FatJet_area,FatJet_jetId,FatJet_P4,FatJet_jetId,Rho_fixedGridRhoFastjetAll,DummyZero,DummyZero)") # args 12, 14, 16, 17 are dummies
  #jVars.Add("NcleanJets", "Sum(cleanedJets)")
  jVars.Add("cleanJet_pt", "cleanedJets[0]")
  jVars.Add("cleanJet_eta", "cleanedJets[1]")
  jVars.Add("cleanJet_phi", "cleanedJets[2]")
  jVars.Add("cleanJet_mass", "cleanedJets[3]")
  #jVars.Add("cleanFatJet_pt", "cleanFatJets[0]")
  #jVars.Add("cleanFatJet_eta", "cleanFatJets[1]")
  #jVars.Add("cleanFatJet_phi", "cleanFatJets[2]")
  #jVars.Add("cleanFatJet_mass", "cleanFatJets[3]")

  # ------------------ MET Selection ------------------
  metVars = VarGroup('METVars')
  metVars.Add("corrMET_pt","cleanMets[4][0]")
  metVars.Add("corrMET_phi","cleanMets[4][1]")
  #metVars.Add("corrMET_dPhiLep","DeltaPhi(lepton_phi, corrMET_phi)")
  metVars.Add("MET_ptcorr", "METptfunc(METcorr, METyr, isMC, corrMET_pt, corrMET_phi, PV_npvs)")
  metVars.Add("MET_phicorr", "METphifunc(METcorr, METyr, isMC, corrMET_pt, corrMET_phi, PV_npvs)")

  metCuts = CutGroup('METCuts')
  metCuts.Add("Pass PuppiMET pt > 50", "corrMET_pt > 50")

  # ------------------ HT Calculation and N Jets cuts ------------------ 
  jVars.Add("minDR_lepsJets","minDR_jets_gtau(nJet, cleanJet_eta, cleanJet_phi, NgoodLeptons, Good4Lepton_eta, Good4Lepton_phi)")

  #jVars.Add("minDR_tauJets","minDR_jets_gtau(nJet, cleanJet_eta, cleanJet_phi, NgoodTau, GoodTau_eta, GoodTau_phi)") # Some function to get minimum DeltaR(cleanJet_eta,cleanJet_phi,goodtau_eta,goodtau_phi)") #tau_eta/phi are lists of *say* 4 entries
  #jVars.Add("ptrel_lepJets","ptRel(cleanJet_pt,cleanJet_eta,cleanJet_phi,cleanJet_mass,lepton_pt,lepton_eta,lepton_phi,lepton_mass)")

  ### ADD HERE a calculation of "minDR" to one of our hadronic taus for each jet. Require below that DR(jet, any tau) > 0.4
  jVars.Add("goodcleanJets", "cleanJet_pt > 30 && abs(cleanJet_eta) < 2.5 && Jet_jetId > 1 && minDR_lepsJets > 0.4 ")
  jVars.Add("NgoodcleanJets", "Sum(goodcleanJets)")
  jVars.Add("gcJet_HT","Sum(cleanJet_pt[goodcleanJets == true])")
  #jVars.Add("DR_lepFatJets","DeltaR_VecAndFloat(FatJet_eta,FatJet_phi,lepton_eta,lepton_phi)")
  #jVars.Add("ptrel_lepFatJets","ptRel(FatJet_pt,FatJet_eta,FatJet_phi,FatJet_mass,lepton_pt,lepton_eta,lepton_phi,lepton_mass)")  
  #jVars.Add("goodcleanFatJets", "FatJet_pt > 200 && abs(FatJet_eta) < 2.4 && FatJet_jetId > 1")
  #jVars.Add("NFatJets", "(int) Sum(goodcleanFatJets)")
  
  # ------------------ Jet pt ordering, counting, lepton association ------------------
  jVars.Add("gcJet_pt_unsort", "cleanJet_pt[goodcleanJets == true]")
  jVars.Add("gcJet_ptargsort","ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(gcJet_pt_unsort))")
  
  jVars.Add("gcJet_pt","reorder(gcJet_pt_unsort,gcJet_ptargsort)")
  jVars.Add("gcJet_eta", "reorder(cleanJet_eta[goodcleanJets == true],gcJet_ptargsort)")
  jVars.Add("gcJet_phi", "reorder(cleanJet_phi[goodcleanJets == true],gcJet_ptargsort)")
  jVars.Add("gcJet_mass", "reorder(cleanJet_mass[goodcleanJets == true],gcJet_ptargsort)")
  
  jVars.Add("gcJet_vetomap", "jetvetofunc(jetvetocorr, gcJet_eta, gcJet_phi)")
  if (isMC):
      jVars.Add("gcJet_hflav", "reorder(Jet_hadronFlavour[goodcleanJets == true],gcJet_ptargsort)")
 
  jVars.Add("gcJet_PNet", "reorder(Jet_btagPNetB[goodcleanJets == true],gcJet_ptargsort)")
  jVars.Add("gcJet_PNetL", "gcJet_PNet > PNetL") 
  #jVars.Add("gcJet_PNetM", "gcJet_PNet > deepjetM")
  jVars.Add("NJets_PNetL", "Sum(gcJet_PNetL)")
  #jVars.Add("NJets_PNetM", "Sum(gcJet_PNetM)")
  
  jVars.Add("gcBJet_eta", "gcJet_eta[gcJet_PNetL]")
  jVars.Add("gcBJet_phi", "gcJet_phi[gcJet_PNetL]")
  jVars.Add("gcBJet_pt", "gcJet_pt[gcJet_PNetL]")
  jVars.Add("gcBJet_mass", "gcJet_mass[gcJet_PNetL]")

  jVars.Add("gcJet_ht", "Sum(gcJet_pt)")

  jCuts = CutGroup('JetCuts')  
  #jCuts.Add('Pass HT > 510', 'gcJet_HT > 510') ## will this be helpful? Not sure...
  
  # Go with medium
  # jCuts.Add('2 GCJets Pass', 'gcJet_HT >= 2')
  jCuts.Add('NgoodcleanJets >= 2', 'NgoodcleanJets >= 2')
  jCuts.Add('2 B Jets Pass (Loose)', 'NJets_PNetL >= 2')
  #jCuts.Add('2 B Jets Pass (Medium)', 'NJets_PNetM >= 2')

  
  GenVars = VarGroup('GenVars')
  GenVars.Add('decayType', 'decayType(isSig, nGenPart, GenPart_pdgId, GenPart_mass, GenPart_pt, GenPart_phi, GenPart_eta, GenPart_genPartIdxMother, GenPart_status)')
  GenVars.Add('decayTypeSum', 'Sum(decayType)')

  GenVars.Add('GenInfo', 'GenInfo(isSig, nGenPart, GenPart_pdgId, GenPart_mass, GenPart_pt, GenPart_phi, GenPart_eta, GenPart_genPartIdxMother)')
  GenVars.Add('genElecPt', 'GenInfo[0]')
  GenVars.Add('genElecDRb', 'GenInfo[1]')
  GenVars.Add('genElecCount', 'Sum(genElecPt>10)')
  GenVars.Add('genMuonPt', 'GenInfo[2]')
  GenVars.Add('genMuonDRb', 'GenInfo[3]')
  GenVars.Add('genMuonCount', 'Sum(genMuonPt>10)')
  GenVars.Add('genTauPt', 'GenInfo[4]')
  GenVars.Add('genTauEta', 'GenInfo[8]')
  GenVars.Add('genTauDRb', 'GenInfo[5]')
  GenVars.Add('visBPt', 'GenInfo[6][0]')
  GenVars.Add('visBEta', 'GenInfo[6][1]')
  GenVars.Add('visBPhi', 'GenInfo[6][2]')
  GenVars.Add('visBM', 'GenInfo[6][3]')
  GenVars.Add('visBbarPt', 'GenInfo[6][4]')
  GenVars.Add('visBbarEta', 'GenInfo[6][5]')
  GenVars.Add('visBbarPhi', 'GenInfo[6][6]')
  GenVars.Add('visBbarM', 'GenInfo[6][7]')
  GenVars.Add('genNeutPt', 'GenInfo[7][0]')
  GenVars.Add('genNeutPhi', 'GenInfo[7][1]')
 
  GenVars.Add('tauBUG', 'tauBUG(isSig, nGenPart, GenPart_pdgId, GenPart_mass, GenPart_pt, GenPart_phi, GenPart_eta, GenPart_genPartIdxMother, GenVisTau_pt, GenVisTau_eta, GenVisTau_phi, GenVisTau_mass, GenVisTau_genPartIdxMother, GenVisTau_status, nGenVisTau)')
  GenVars.Add('piplusTauR', 'tauBUG[0]')
  GenVars.Add('pinegTauR', 'tauBUG[1]')

  # GoodTau = reco taus
  recoGenVars = VarGroup('recoGenVars')
  recoGenVars.Add('Matching', 'recoGenMatch(isSig, NgoodLeptons, Good4Lepton_eta, Good4Lepton_phi, NJets_PNetL, gcJet_PNet, gcBJet_eta, gcBJet_phi, gcBJet_pt, nGenPart, GenPart_pdgId, GenPart_mass, GenPart_pt, GenPart_phi, GenPart_eta, GenPart_genPartIdxMother, GenPart_status)')
  #recoGenVars.Add('MatchingSum', 'Sum(Matching)')
  recoGenVars.Add('ObjectList_indices', 'Matching[0]') # Topograph thing
  recoGenVars.Add('matchabilityArr', 'Matching[1]') # Topograph thing
  recoGenVars.Add('nObjects', 'Matching[2]') # Topograph thing
  recoGenVars.Add('nbjets', 'Matching[3]') # Topograph thing
  recoGenVars.Add('tauArr', 'Matching[4]')
  recoGenVars.Add('bjetArr', 'Matching[5]')
  recoGenVars.Add('vecbArr', 'Matching[6]')
  recoGenVars.Add('matchability', 'convertMatchToInt(matchabilityArr)')
  recoGenVars.Add('genHadrTauCount', 'Matching[7]')
  
  recoGenVars.Add('ObjectList', 'ObjectList(isSig, NgoodLeptons, Good4Lepton_pt, Good4Lepton_eta, Good4Lepton_phi, Good4Lepton_mass, NJets_PNetL, gcBJet_pt, gcBJet_eta, gcBJet_phi, gcBJet_mass)') # Topograph thing
  recoGenVars.Add('PtListObject', 'ObjectList[0]') # PtListObject, EtaListObject, PhiListObject, EnergyListObject, TaggedListObject
  recoGenVars.Add('EtaListObject', 'ObjectList[1]')
  recoGenVars.Add('PhiListObject', 'ObjectList[2]')
  recoGenVars.Add('EnergyListObject', 'ObjectList[3]')
  recoGenVars.Add('TaggedListObject', 'ObjectList[4]')
  
  recoGenVars.Add('GenList', 'GenList(isSig, tauArr, bjetArr, vecbArr, nGenPart, GenPart_pdgId, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)') # Topograph thing
  recoGenVars.Add('pdgIdListGen', 'GenList[0]') # pdgIdListGen, ptListGen, etaListGen, phiListGen, massListGen
  recoGenVars.Add('ptListGen', 'GenList[1]')
  recoGenVars.Add('etaListGen', 'GenList[2]')
  recoGenVars.Add('phiListGen', 'GenList[3]')
  recoGenVars.Add('massListGen', 'GenList[4]')
  
  # Shoudl isSig be included?
  manualVars = VarGroup('manualVars')
  manualVars.Add('manual', 'funcmassdiff(Good4Lepton_pt, Good4Lepton_eta, Good4Lepton_phi, Good4Lepton_mass, gcBJet_pt, gcBJet_eta, gcBJet_phi, gcBJet_mass)')
  manualVars.Add('B1finalPx', 'manual[0]')
  manualVars.Add('B1finalPy', 'manual[1]')
  manualVars.Add('B1finalPz', 'manual[2]')
  manualVars.Add('B1finalM', 'manual[3]')
  manualVars.Add('B2finalPx', 'manual[4]')
  manualVars.Add('B2finalPy', 'manual[5]')
  manualVars.Add('B2finalPz', 'manual[6]')
  manualVars.Add('B2finalM', 'manual[7]')

  # ---------------------------------------------------------------------------

  #jVars.Add("gcFatJet_pt_unsort", "FatJet_pt[goodcleanFatJets == true]")
  #jVars.Add("gcFatJet_ptargsort","ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(gcFatJet_pt_unsort))")
  #jVars.Add("gcFatJet_pt","reorder(gcFatJet_pt_unsort,gcFatJet_ptargsort)")  
  #jVars.Add("gcFatJet_eta", "reorder(FatJet_eta[goodcleanFatJets == true],gcFatJet_ptargsort)")
  #jVars.Add("gcFatJet_phi", "reorder(FatJet_phi[goodcleanFatJets == true],gcFatJet_ptargsort)")
  #jVars.Add("gcFatJet_mass", "reorder(FatJet_mass[goodcleanFatJets == true],gcFatJet_ptargsort)")
  #jVars.Add("gcFatJet_sdmass", "reorder(FatJet_msoftdrop[goodcleanFatJets == true],gcFatJet_ptargsort)")
  #jVars.Add("gcFatJet_vetomap", "jetvetofunc(jetvetocorr, gcFatJet_eta, gcFatJet_phi)")
  #WORK ON THIS MORE -- need to just be isolated from the 3 highest-pt fat jets, not any of them...
  #jVars.Add("Isolated_AK4","standalone_Jet(gcJet_eta, gcJet_phi, gcFatJet_eta, gcFatJet_phi)")

  # ------------------ Add scale factors and MC jet-based calcs ------------------
  if isMC:
    #jVars.Add("leptonRecoSF", "recofunc(electroncorr, muonidcorr, yrstr, lepton_pt, lepton_eta, isEl)") ## this is not right, but we'll figure out what corrections we need later
    #jVars.Add("leptonIDSF", "idfunc(muonidcorr,elid_pts,elid_etas,elecidsfs,elecidsfuncs,yrstr, lepton_pt, lepton_eta, isEl)") #at(0) 
    #jVars.Add("leptonIsoSF", "isofunc(muiso_pts,muiso_etas,muonisosfs,muonisosfunc,elid_pts,elid_etas,elecisosfs,elecisosfunc, lepton_pt, lepton_eta, isEl)")
    #jVars.Add("leptonHLTSF", "hltfunc(muonhltcorr,elhlt_pts,elhlt_etas,elechltsfs,elechltuncs,yrstr, lepton_pt, lepton_eta, isEl)")
    jVars.Add("btagWeights","btagshapefunc(yrstr, jesvar, btagwpbccorr, btagwplcorr, btagpts, btageffs, PNetL, gcJet_pt, gcJet_eta, gcJet_PNet, gcJet_hflav)")
    #### WORK ON THESE! Ethan can compute the efficiencies we need.


  # # ------------------ Results ------------------
  # # rframeVars = VarGroup('restFrameVars')
  # # rframeVars.Add('VLQ_mass', 'rfc.compute_mass(rdfslot_, lepton_pt, lepton_eta, lepton_phi, lepton_mass, gcFatJet_pt, gcFatJet_eta, gcFatJet_phi, gcFatJet_mass, MET_pt, MET_phi)')
  # # rframeVars.Add('VLQ_mass_T', 'VLQ_mass[0]')
  # # rframeVars.Add('VLQ_mass_Tbar', 'VLQ_mass[1]')
  # # rframeVars.Add('VLQ_mass_T_r', 'VLQ_mass[2]')
  # # rframeVars.Add('VLQ_mass_Tbar_r', 'VLQ_mass[3]')
  # # rframeVars.Add('VLQ_mass_ratio', 'VLQ_mass_T/VLQ_mass_Tbar')
  # # rframeVars.Add('VLQ_mass_avg', '(VLQ_mass_T+VLQ_mass_Tbar)*0.5')
  
  
  # # -------------------------------------

  nodeToPlot = a.Apply([flagCuts, gjsonVars, gjsonCuts, tVars, eandmuVars, lVars, lCuts, TrigVars])
  
  if isSig:
    a.Apply([GenVars, GenCuts])

  # # We want the BW decays that go to l + nu
  # ## These will be meaningless for non-signal, but shouldn't crash...
  # #a.Define("decayMODE", "decayModeSelection(region, nGenPart,GenPart_pdgId,GenPart_mass,GenPart_pt,GenPart_phi,GenPart_eta,GenPart_genPartIdxMother,GenPart_status)")	
  # #a.Define("isLeptWdecay", "(decayMODE == 101 || decayMODE == 105 || decayMODE == 106)") 
  # #a.Define("isLeptTdecay", "(1002 <= decayMODE && decayMODE <= 1007) || decayMODE == 1011 || decayMODE == 1012") 

  # # Solution to cleanJets() problem:
  # #       The analyzer .Apply() calls the analyzer .Define().  This .Define() calls self._collectionOrg.CollectionDefCheck(var, newNode).
  # #  This method executes this line: if re.search(r"\b" + re.escape(c+'s') + r"\b", action_str) and (c+'s' not in self._builtCollections):
  # #                                       print ('MAKING %ss for %s'%(c,action_str))
  # #  Apparently somethings get discarded from this _collectionOrg?
  # #  Instead force the .Apply() from the ActiveNode because the node .Apply() is better.
  
  newNode = a.ActiveNode.Apply(jVars)
  a.SetActiveNode(newNode)
  if isSig:
      a.Apply([recoGenVars])
  a.Apply([jCuts, metVars, manualVars, tCuts, metCuts, lepSFs, TrigCuts])  #,  testCuts, rframeVars
  
  allColumns = a.GetColumnNames()
  
  #old column list before append operation worked
  #columns = ["cleanJet_pt", "gcJet_HT"] #,'NJets_PNetM','corrMET_pt', 'decayTypeSum', 'ObjectList_indices', 'matchability', 'nObjects', 'nbjets', 'PtListObject', 'EtaListObject', 'PhiListObject', 'EnergyListObject', 'TaggedListObject', 'pdgIdListGen', 'ptListGen', 'etaListGen', 'phiListGen', 'massListGen', 'PtListObject', 'EtaListObject', 'PhiListObject', 'EnergyListObject', 'TaggedListObject', 'B1finalPx', 'B1finalPy', 'B1finalPz', 'B1finalM', 'B2finalPx', 'B2finalPy', 'B2finalPz','B2finalM', 'GenPart_pdgId', 'GenPart_mass', 'GenPart_pt', 'GenPart_phi', 'GenPart_eta', 'GenPart_genPartIdxMother', 'genElecPt', 'genElecDRb', 'genMuonPt', 'genMuonDRb', 'genTauPt', 'genTauEta', 'genTauDRb', 'visBPt', 'visBEta', 'visBPhi', 'visBM', 'visBbarPt', 'visBbarEta', 'visBbarPhi', 'visBbarM', 'genNeutPt', 'genNeutPhi', 'NgoodElecs', 'NgoodMuons', 'NgoodTau', 'NgoodLeptons', 'genElecCount', 'genHadrTauCount', 'Good4Lepton_pt', 'Good4Lepton_eta', 'Good4Lepton_phi', 'Good4Lepton_mass', 'Good4Lepton_ID', 'piplusTauR', 'pinegTauR', 'elrecoSF', 'elidSF', 'muonidSF', 'muonisoSF', 'tauidVSeSF', 'tauidVSmuSF', 'tauidVSjetSF', 'btagWeights', 'PileupWeights', 'NlessRestrict', 'NwithElMu_Tau', 'NwithDiTau' ,'NwithAll' ,'Nnone', 'MET_ptcorr', 'MET_phicorr', 'passesMuPD', 'passesMuEGPD', 'passesEGPD', 'passesTauPD'] 
  
  #entire list made from appendages; for data
  #columns = ['B1finalM', 'B1finalPx', 'B1finalPy', 'B1finalPz', 'B2finalM', 'B2finalPx', 'B2finalPy', 'B2finalPz', 'FatJet_rawFactor', 'Good4Lepton_ID', 'Good4Lepton_TauDM', 'Good4Lepton_charge', 'Good4Lepton_eta', 'Good4Lepton_mass', 'Good4Lepton_phi', 'Good4Lepton_pt', 'GoodEl_ID', 'GoodEl_charge', 'GoodEl_eta', 'GoodEl_mass', 'GoodEl_mtforTau', 'GoodEl_phi', 'GoodEl_pt', 'GoodLepton_ID', 'GoodLepton_TauDM', 'GoodLepton_charge', 'GoodLepton_eta', 'GoodLepton_mass', 'GoodLepton_phi', 'GoodLepton_pt', 'GoodMu_ID', 'GoodMu_charge', 'GoodMu_eta', 'GoodMu_mass', 'GoodMu_mtforTau', 'GoodMu_phi', 'GoodMu_pt', 'GoodTau_DM', 'GoodTau_ID', 'GoodTau_charge', 'GoodTau_eta', 'GoodTau_mass', 'GoodTau_phi', 'GoodTau_pt', 'Jet_rawFactor', 'L1Reco_step', 'Muon_isPFcand', 'Muon_tightId', 'Muon_tunepRelPt', 'NEle12LT23', 'NEle23', 'NEle24forTau', 'NEle30', 'NJets_PNetL', 'NMu10LT20', 'NMu12LT23', 'NMu20', 'NMu20forTau', 'NMu23', 'NTauEleTrig', 'NTauLooseTrig', 'NTauMedTrig', 'NTauMuTrig', 'NgoodElecs', 'NgoodLeptons', 'NgoodMuons', 'NgoodTau', 'NgoodTaue', 'NgoodTaumu', 'NgoodcleanJets', 'NlessRestrict', 'Nnone', 'NwithAll', 'NwithDiTau', 'NwithElMu_Tau', 'PV_chi2', 'PV_ndof', 'PV_npvs', 'PV_npvsGood', 'PV_score', 'PV_x', 'PV_y', 'PV_z', 'PuppiMET_phi', 'PuppiMET_phiJERDown', 'PuppiMET_phiJERUp', 'PuppiMET_phiJESDown', 'PuppiMET_phiJESUp', 'PuppiMET_phiUnclusteredDown', 'PuppiMET_phiUnclusteredUp', 'PuppiMET_pt', 'PuppiMET_ptJERDown', 'PuppiMET_ptJERUp', 'PuppiMET_ptJESDown', 'PuppiMET_ptJESUp', 'PuppiMET_ptUnclusteredDown', 'PuppiMET_ptUnclusteredUp', 'PuppiMET_sumEt', 'RawPuppiMET_phi', 'RawPuppiMET_pt', 'RawPuppiMET_sumEt', 'Rho_fixedGridRhoAll', 'Rho_fixedGridRhoFastjetAll', 'Rho_fixedGridRhoFastjetCentral', 'Rho_fixedGridRhoFastjetCentralCalo', 'Rho_fixedGridRhoFastjetCentralChargedPileUp', 'Rho_fixedGridRhoFastjetCentralNeutral', 'bunchCrossing', 'cleanJet_eta', 'cleanJet_mass', 'cleanJet_phi', 'cleanJet_pt', 'corrMET_phi', 'corrMET_pt', 'event', 'gLepton', 'gcBJet_eta', 'gcBJet_mass', 'gcBJet_phi', 'gcBJet_pt', 'gcJet_HT', 'gcJet_PNet', 'gcJet_PNetL', 'gcJet_eta', 'gcJet_ht', 'gcJet_mass', 'gcJet_phi', 'gcJet_pt', 'gcJet_pt_unsort', 'gcJet_ptargsort', 'gcJet_vetomap', 'goodElectrons', 'goodMuons', 'goodTau', 'goodTaue', 'goodTaumu', 'goodc:qleanJets', 'hasBosonishMass', 'igLepton', 'isElElTrig', 'isElTauTrig', 'isMuElTrig', 'isMuMuTrig', 'isMuTauTrig', 'isSingleEleTrig', 'isSingleMuTrig', 'isSingleTauTrig', 'isTauTauTrig', 'luminosityBlock', 'minDR_lepsJets', 'nElectron', 'nFatJet', 'nJet', 'nMuon', 'nPPSLocalTrack', 'nProton_multiRP', 'nProton_singleRP', 'nSoftActivityJet', 'passesEGPD', 'passesJSON', 'passesMuEGPD', 'passesMuPD', 'passesTauPD', 'run']

  #debugging list: taking out MET corrs
 
  #removed: corrMET_pt, corrMET_phi

   
  columns = []
  for col in allColumns:
     if ("P4" in col) or ("cleanedJets" in col) or ("cleanFatJets" in col) or ("cleanMets" in col) or ("Dummy" in col): continue 
     if ("LHE" in col) and ("Weight" not in col) and (col != "LHE_HT") and (col != "LHE_Vpt") and (col != "gcHTCorr_WjetLHE"): continue
     if col.startswith("Muon") and ("_tightId" not in col) and ("_isPF" not in col) and ("tunep" not in col) and ("genPartFlav" not in col): continue
     if col.startswith("Electron") and ("genPartFlav" not in col): continue
     if col.startswith("Jet") and ("rawFactor" not in col): continue
     if col.startswith("FatJet") and ("rawFactor" not in col): continue
     if col.startswith("PPS") or col.startswith("Proton") or col.startswith("L1_"): continue
     if col.startswith("Gen") or col.startswith("Soft") or col.startswith("fixed"): continue
     if col.startswith("Sub")  or col.startswith("Calo") or col.startswith("Chs"): continue
     if col.startswith("Corr") or col.startswith("Fsr") or col.startswith("Iso") or col.startswith("Tau"): continue
     if col.startswith("SV") or col.startswith("Photon") or col.startswith("Low"): continue
     if col.startswith("HLT") or col.startswith("HT") or col.startswith("boosted") or col.startswith("Deep"): continue
     if col.startswith("Flag") or col == "Bprime_gen_info" or col == "t_gen_info" or col == "W_gen_info" or col == "metxyoutput": continue
     if col == "assignleps" or col == "pnetoutput" or col == "t_output" or col == "Bprime_output" or col.startswith("Other"): continue
     if col.startswith("PS") or col.startswith("Tk") or col.startswith("Trig"): continue
     if col.startswith("nCorr") or col.startswith("nFsr"): continue
     if col.startswith("nGen") or col.startswith("nIso") or col.startswith("nLow"): continue
     if col.startswith("nOther") or col.startswith("nPS") or col.startswith("nPhoton"): continue
     if col.startswith("nSV") or col.startswith("nSub") or col.startswith("nTau") or col.startswith("nTrig"): continue
     if col.startswith("nboosted"): continue
     if col == "tauBUG": continue
     if col == "Matching": continue
     if col == "ObjectList": continue
     if col == "manual": continue
     if col.startswith("BeamSpot"): continue
     if col.startswith("Lepton"): continue
     if col.startswith("iLepton"): continue
     if col.startswith("MET"): continue
     if col.startswith("RawMET"): continue
     
     columns.append(col)
     

  finalFile = sample + era + "_" + year + "_" + str(testNum1) + ".root"
  if not isMC:
    finalFile = sample + era + ver + "_" + year + "_" + str(testNum1) + ".root";
 
  mode = 'RECREATE'
  if jesvar != "Nominal":
    mode = 'UPDATE'
   
  a.Snapshot(columns, finalFile, "Events_"+jesvar, lazy=False, openOption=mode, saveRunChain=True)
    
  if jesvar == "Nominal":
    print("Cut statistics:")
    rep = a.DataFrame.Report()
    rep.Print()

  print("--------- Analysis End ---------")
    
  a.Close()

if not isMC:
  analyze("Nominal")
else:
  analyze("Nominal")
  #TODO fix why this not work?  shifts = ["Nominal","JECup","JECdn","JERup","JERdn"]
  #for shift in shifts:
  #  analyze(shift)
