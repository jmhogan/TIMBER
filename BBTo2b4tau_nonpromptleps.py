from TIMBER.Analyzer import *
from TIMBER.Tools.Common import *
from rates import *
import ROOT
from ROOT import TFile
import sys, os
import gc

gc.disable()

from TIMBER.Tools.RestFramesHandler import load_restframes
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
filelist = []
for i, line in enumerate(lines):
  if i in range(start,end+1):
    filelist.append(line.strip())
print(filelist)

print("Number of Entries:",len(filelist))
print("list contents:",filelist)
sampleName = filelist[0]
print(f"Sample Name: {sampleName}")

# Parse the incoming file names to assign labels  
isMC = not (("Single" in sampleName) or ("Muon" in sampleName) or ("EGamma" in sampleName) or ("MuonEG" in sampleName) or ("Tau/" in sampleName))
dataset = int(-1)
if "MuonEG" in sampleName: dataset = 2
elif ("Muon" in sampleName) or ("Muon0" in sampleName) or ("Muon1" in sampleName): dataset = 1
elif ("EGamma" in sampleName) or ("EGamma0" in sampleName) or ("EGamma1" in sampleName): dataset = 3
elif "Tau/" in sampleName: dataset = 4


#'root://cms-xrd-global.cern.ch//store/data/Run2018A/SingleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v2/2550000/28FF17A8-95EB-FD41-A55B-2EFAF2D6AF91.root' 
tokens = sampleName.split("/")
sample = tokens[7] # was 5
runera = tokens[6] # was 4
process = tokens[9] # was 7
era = runera[-1] # last char
print(process)
ver = process[process.find('_')+1:process.find('_')+3]
del tokens

jecera = era
if year == '2022':
  jecera = 'CD'
elif year == '2023':
  if(ver != 'v4'):
    jecera = 'Cv123'
  else:
    jecera = 'Cv4'
    
region = "DataBkg"

print('============== RUN SETTINGS ===============')
print('isMC = ',isMC)
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
CompileCpp('TIMBER/Framework/Tprime1lep/nonpromptweight.cc')
ROOT.gInterpreter.ProcessLine('#include "TString.h"')

#ROOT.gInterpreter.ProcessLine('#include "TIMBER/Framework/Tprime1lep/neononpromptweight.cc"')


# Enable using 4 threads
ROOT.ROOT.EnableImplicitMT(num_threads)

# load rest frames handler
handler_name = 'Bprime_handler_new.cc'
class_name = 'Bprime_RestFrames_Container_new'
load_restframes(num_threads, handler_name, class_name, 'B_rfc')

#CompileCpp('bin/restframes/helper.cc')
# Essentially just called return doubles
# Now we implement directly

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
  
  ROOT.gInterpreter.ProcessLine('initialize(year);')


  # ------------------ correctionsLib corrections ------------------

  BTagL = {'2022':0.047,'2022EE':0.0499,'2023':0.0358,'2023BPix':0.0359, '2024':0.0246, '2025':0.0246} #PNet for 22-23BPix, UParT for 2024-2025
  yrstr = {'2022':"Run3-22CDSep23-Summer22-NanoAODv12",'2022EE':"Run3-22EFGSep23-Summer22EE-NanoAODv12",'2023':"Run3-23CSep23-Summer23-NanoAODv12",'2023BPix':"Run3-23DSep23-Summer23BPix-NanoAODv12",'2024':"Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15",'2025':"Run3-25Prompt-Summer24-NanoAODv15"}
  jmetags = {'2022':'2026-06-05','2022EE':'2026-06-05','2023':'2026-06-05','2023BPix':'2026-06-05','2024':'2026-06-05', '2025':'2026-06-05'} 
  jecver = {'2022':"V4",'2022EE':"V4",'2023':"V4",'2023BPix':"V4",'2024':"V3",'2025':"V3"}   
  jmeyrstr = {'2022':yrstr['2022'],'2022EE':yrstr['2022EE'],'2023':yrstr['2023'],'2023BPix':yrstr['2023BPix'],'2024':yrstr['2024'],'2025':'Run3-25Prompt-Winter25-NanoAODv15'} # yes, really use 24 for 25
  jecyr = {'2022':"Summer22_22Sep2023",'2022EE':"Summer22EE_22Sep2023",'2023':"Summer23Prompt23",'2023BPix':"Summer23BPixPrompt23",'2024':"Summer24Prompt24", '2025':"Winter25Prompt25"}
  jeryr = {'2022':"Summer22_22Sep2023_JRV2",'2022EE':"Summer22EE_22Sep2023_JRV2",'2023':"Summer23Prompt23_RunCv1234_JRV2",'2023BPix':"Summer23BPixPrompt23_RunD_JRV2",'2024':"Summer24Prompt24_JRV1",'2025':"Summer24Prompt25_JRV1"}
  jetvetoname = {'2022':"Summer22_23Sep2023_RunCD_V1",'2022EE':"Summer22EE_23Sep2023_RunEFG_V1",'2023':"Summer23Prompt23_RunC_V1",'2023BPix':"Summer23BPixPrompt23_RunD_V1",'2024':"Summer24Prompt24_RunBCDEFGHI_V1",'2025':"Winter25Prompt25_RunCDEFG_V1"    }      
 
  ROOT.gInterpreter.Declare("""
  float BTagL = """+str(BTagL[year])+""";
  string yrstr = \""""+yrstr[year]+"""\";
  string jmetag = \""""+jmetags[year]+"""\";
  string jecyr = \""""+jecyr[year]+"""\";
  string jecver = \""""+jecver[year]+"""\";
  string jetvetoname = \""""+jetvetoname[year]+"""\";
  """)

  # *************** muonisocorr does not match the muon definition !!!!! (change to mediumPFIso hopefully) *****************
  ROOT.gInterpreter.Declare("""
  auto jetvetocorrset = correction::CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/"+yrstr+"/"+jmetag+"/jetvetomaps.json.gz");

  auto jetvetocorr = jetvetocorrset->at(jetvetoname);
  """)

  print('Finished declaring corrections')
  
  ROOT.gInterpreter.Declare("""
  auto ak4corrset = correction::CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/"+yrstr+"/"+jmetag+"/jet_jerc.json.gz"); 
  auto ak8corrset = correction::CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/"+yrstr+"/"+jmetag+"/fatJet_jerc.json.gz"); 

  auto ak4corr = ak4corrset->compound().at(jecyr+"_"+jecver+"_DATA_L1L2L3Res_AK4PFPuppi");
  auto ak4corrL1 = ak4corrset->at(jecyr+"_"+jecver+"_DATA_L1FastJet_AK4PFPuppi");
  auto ak8corr = ak8corrset->compound().at(jecyr+"_"+jecver+"_DATA_L1L2L3Res_AK8PFPuppi");
  """)

  ROOT.gInterpreter.Declare("""
  std::vector<float> etabinsel = {-2.5,-2.0,-1.566,-1.444,-0.8,0.0,0.8,1.444,1.566,2.0,2.5};
  std::vector<float> etabinsmu = {0.0,0.9,1.2,2.1,2.4};
  std::vector<float> ptbins = {10.0,20.0,30.0,40.0,50.0,70.0,100.0,200.0,9999.0};
  std::vector<float> ptbins_datafake = {10.0,15.0,20.0,25.0,30.0,40.0,50.0,70.0,9999.0};

  std::vector<std::vector<float>> ePromptRatesMC = """ + to_cpp_vec2d(prompt_ttbar_el[year]) + """;
  std::vector<std::vector<float>> mPromptRatesMC = """ + to_cpp_vec2d(prompt_ttbar_mu[year]) + """;
  std::vector<std::vector<float>> tPromptRatesMC = """ + to_cpp_vec2d(prompt_ttbar_tau[year]) + """;
  std::vector<std::vector<float>> eFakeRatesMC = """ + to_cpp_vec2d(fake_QCD_el[year]) + """;
  std::vector<std::vector<float>> mFakeRatesMC = """ + to_cpp_vec2d(fake_QCD_mu[year]) + """;
  std::vector<std::vector<float>> tFakeRatesMC = """ + to_cpp_vec2d(fake_QCD_tau[year]) + """;
  std::vector<std::vector<float>> ePromptRatesData = """ + to_cpp_vec2d(prompt_Data_el[year]) + """;
  std::vector<std::vector<float>> mPromptRatesData = """ + to_cpp_vec2d(prompt_Data_mu[year]) + """;
  std::vector<std::vector<float>> eFakeRatesData = """ + to_cpp_vec2d(fake_Data_el[year]) + """;
  std::vector<std::vector<float>> mFakeRatesData = """ + to_cpp_vec2d(fake_Data_mu[year]) + """;
  """)
  
  # ------------------ Flag Cuts ------------------
  flagCuts = CutGroup('FlagCuts')
  flagCuts.Add('Bad Event Filters', 'Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && Flag_goodVertices == 1 && Flag_eeBadScFilter == 1 && Flag_globalSuperTightHalo2016Filter == 1 && Flag_BadPFMuonFilter == 1 && Flag_BadPFMuonDzFilter == 1')
  flagCuts.Add('Event has jets', 'nJet > 0') # need jets   && nFatJet > 0
  
  # ------------------ Golden JSON (Data) || GEN Info (MC) ------------------
  gjsonVars = VarGroup('GoldenJsonVars')
  gjsonCuts = CutGroup('GoldenJsonCuts')
  gjsonVars.Add("passesJSON", "goldenjson(myLumiMask, run, luminosityBlock)") # function name, parameters are branches in the file or defined objects
  gjsonCuts.Add("Data passes Golden JSON", "passesJSON == 1") 

  # ------------------ Tau selection criteria ----------------
  tVars = VarGroup('TauVars')
  tVars.Add('goodTaumu', 'Tau_pt > 20 && abs(Tau_eta) < 2.5 && abs(Tau_dz) < 0.2 && Tau_idDeepTau2018v2p5VSmu >= 1')
  tVars.Add('NgoodTaumu', 'Sum(goodTaumu)')

  tVars.Add('goodTaue', 'goodTaumu == 1 && Tau_idDeepTau2018v2p5VSe >= 2')
  tVars.Add('NgoodTaue', 'Sum(goodTaue)')

  #tVars.Add('goodTau', 'goodTaue == 1 && Tau_idDeepTau2018v2p5VSjet >= 4')
  #tVars.Add('NgoodTau', 'Sum(goodTau)')
  tVars.Add('looseTau', 'goodTaue == 1 && Tau_idDeepTau2018v2p5VSjet >= 2')
  tVars.Add('NlooseTau', 'Sum(looseTau)')

  tVars.Add('LooseTau_eta', 'Tau_eta[looseTau == true]')
  tVars.Add('LooseTau_pt', 'Tau_pt[looseTau == true]')
  tVars.Add('LooseTau_phi', 'Tau_phi[looseTau == true]')
  tVars.Add('LooseTau_mass', 'Tau_mass[looseTau == true]')
  tVars.Add('LooseTau_charge', 'Tau_charge[looseTau == true]')
  tVars.Add('LooseTau_ID', 'ROOT::VecOps::RVec<int>(NlooseTau, 15)')
  tVars.Add('LooseTau_isTight','Tau_idDeepTau2018v2p5VSjet[looseTau == true] >= 4')

  # tVars.Add('GoodTau_eta', 'Tau_eta[goodTau == true]')
  # tVars.Add('GoodTau_phi', 'Tau_phi[goodTau == true]')
  # tVars.Add('GoodTau_pt', 'Tau_pt[goodTau == true]')
  # tVars.Add('GoodTau_mass', 'Tau_mass[goodTau == true]')
  # tVars.Add('GoodTau_ID', 'ROOT::VecOps::RVec<int>(NgoodTau, 15)')
  # tVars.Add('GoodTau_DM', 'Tau_decayMode[goodTau == true]')
  # if (isMC):
  #    tVars.Add('GoodTau_genmch', 'Tau_genPartFlav[goodTau == true]')
  # tVars.Add('GoodTau_charge', 'Tau_charge[goodTau == true]')
  
  # -------------------- EL and MUON Definitions --------------------
  eandmuVars = VarGroup('ElandMuVars')

  #Good Electrons
  eandmuVars.Add('looseElectrons', 'Electron_pt > 10 && abs(Electron_eta) < 2.5 && (Electron_mvaIso_WP90 == 1)')
  eandmuVars.Add('NlooseElecs', 'Sum(looseElectrons)')
  #eandmuVars.Add('goodElectrons', 'Electron_pt > 10 && abs(Electron_eta) < 2.5 && (Electron_mvaIso_WP80 == 1)')
  #eandmuVars.Add('NgoodElecs', 'Sum(goodElectrons)')
  
  eandmuVars.Add('LooseEl_pt','Electron_pt[looseElectrons == true]')
  eandmuVars.Add('LooseEl_eta','Electron_eta[looseElectrons == true]')
  eandmuVars.Add('LooseEl_phi', 'Electron_phi[looseElectrons == true]')
  eandmuVars.Add('LooseEl_mass', 'Electron_mass[looseElectrons == true]')
  eandmuVars.Add('LooseEl_charge', 'Electron_charge[looseElectrons == true]')
  eandmuVars.Add('LooseEl_ID', 'ROOT::VecOps::RVec<int>(NlooseElecs, 11)')
  eandmuVars.Add('LooseEl_isTight','Electron_mvaIso_WP80[looseElectrons == true] == 1')

  # eandmuVars.Add('GoodEl_pt','Electron_pt[goodElectrons == true]')
  # eandmuVars.Add('GoodEl_eta','Electron_eta[goodElectrons == true]')
  # eandmuVars.Add('GoodEl_phi','Electron_phi[goodElectrons == true]')
  # eandmuVars.Add('GoodEl_mass','Electron_mass[goodElectrons == true]')
  # eandmuVars.Add('GoodEl_ID', 'ROOT::VecOps::RVec<int>(NgoodElecs, 11)')
  # eandmuVars.Add('GoodEl_mtforTau', 'ROOT::VecOps::RVec<int>(NgoodElecs, -1)')
  # eandmuVars.Add('GoodEl_charge', 'Electron_charge[goodElectrons == true]')

  #Good Muons
  eandmuVars.Add('looseMuons', 'Muon_pt >= 15 && abs(Muon_eta) < 2.4 && Muon_looseId == 1 && abs(Muon_dxy) < 0.2 && abs(Muon_dz) < 0.5 && Muon_pfIsoId >= 2')
  eandmuVars.Add('NlooseMuons', 'Sum(looseMuons)')
  #eandmuVars.Add('goodMuons', 'Muon_pt >= 15 && abs(Muon_eta) < 2.4 && Muon_mediumId == 1 && abs(Muon_dxy) < 0.2 && abs(Muon_dz) < 0.5 && Muon_pfIsoId >= 3')
  #eandmuVars.Add('NgoodMuons', 'Sum(goodMuons)')
  
  eandmuVars.Add('LooseMu_pt','Muon_pt[looseMuons == true]')
  eandmuVars.Add('LooseMu_eta','Muon_eta[looseMuons == true]')
  eandmuVars.Add('LooseMu_phi', 'Muon_phi[looseMuons == true]')
  eandmuVars.Add('LooseMu_mass', 'Muon_mass[looseMuons == true]')
  eandmuVars.Add('LooseMu_charge', 'Muon_charge[looseMuons == true]')
  eandmuVars.Add('LooseMu_ID', 'ROOT::VecOps::RVec<int>(NlooseMuons, 13)')
  eandmuVars.Add('LooseMu_isTight','Muon_mediumId[looseMuons == true] == 1 && Muon_pfIsoId[looseMuons == true] >= 3')
  
  # eandmuVars.Add('GoodMu_pt','Muon_pt[goodMuons == true]')
  # eandmuVars.Add('GoodMu_eta','Muon_eta[goodMuons == true]')
  # eandmuVars.Add('GoodMu_phi','Muon_phi[goodMuons == true]')
  # eandmuVars.Add('GoodMu_mass','Muon_mass[goodMuons == true]')
  # eandmuVars.Add('GoodMu_ID', 'ROOT::VecOps::RVec<int>(NgoodMuons, 13)')
  # eandmuVars.Add('GoodMu_mtforTau', 'ROOT::VecOps::RVec<int>(NgoodMuons, -1)')
  # eandmuVars.Add('GoodMu_charge', 'Muon_charge[goodMuons == true]')

 
  # ------------------------ LEPTON Definitions and Selection -------------------------
  lVars = VarGroup('LeptonVars')

  lVars.Add('ilLepton', 'ROOT::VecOps::Concatenate(looseElectrons, looseMuons)')
  lVars.Add('lLepton', 'ROOT::VecOps::Concatenate(ilLepton, looseTau)')
  lVars.Add('NlooseLeptons', 'Sum(lLepton)')  
  lVars.Add('ilLepton_pt', 'ROOT::VecOps::Concatenate(LooseEl_pt, LooseMu_pt)')
  lVars.Add('lLepton_pt', 'ROOT::VecOps::Concatenate(ilLepton_pt, LooseTau_pt)') 
  lVars.Add('ilLepton_eta', 'ROOT::VecOps::Concatenate(LooseEl_eta, LooseMu_eta)')
  lVars.Add('lLepton_eta', 'ROOT::VecOps::Concatenate(ilLepton_eta, LooseTau_eta)')
  lVars.Add('ilLepton_phi', 'ROOT::VecOps::Concatenate(LooseEl_phi, LooseMu_phi)')
  lVars.Add('lLepton_phi', 'ROOT::VecOps::Concatenate(ilLepton_phi, LooseTau_phi)')
  lVars.Add('ilLepton_mass', 'ROOT::VecOps::Concatenate(LooseEl_mass, LooseMu_mass)')
  lVars.Add('lLepton_mass', 'ROOT::VecOps::Concatenate(ilLepton_mass, LooseTau_mass)')
  lVars.Add('ilLepton_charge', 'ROOT::VecOps::Concatenate(LooseEl_charge, LooseMu_charge)')
  lVars.Add('lLepton_charge', 'ROOT::VecOps::Concatenate(ilLepton_charge, LooseTau_charge)')
  lVars.Add('ilLepton_ID', 'ROOT::VecOps::Concatenate(LooseEl_ID, LooseMu_ID)')
  lVars.Add('lLepton_ID', 'ROOT::VecOps::Concatenate(ilLepton_ID, LooseTau_ID)')
  lVars.Add('ilLepton_isTight', 'ROOT::VecOps::Concatenate(LooseEl_isTight, LooseMu_isTight)')
  lVars.Add('lLepton_isTight', 'ROOT::VecOps::Concatenate(ilLepton_isTight, LooseTau_isTight)')

  lVars.Add("lLepton_ptargsort","ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(lLepton_pt))")
  lVars.Add("LooseLepton_pt","reorder(lLepton_pt,lLepton_ptargsort)")
  lVars.Add("LooseLepton_eta","reorder(lLepton_eta,lLepton_ptargsort)")
  lVars.Add("LooseLepton_phi","reorder(lLepton_phi,lLepton_ptargsort)")
  lVars.Add("LooseLepton_mass","reorder(lLepton_mass,lLepton_ptargsort)")
  lVars.Add("LooseLepton_charge","reorder(lLepton_charge,lLepton_ptargsort)")
  lVars.Add("LooseLepton_ID","reorder(lLepton_ID,lLepton_ptargsort)")
  lVars.Add("LooseLepton_isTight","reorder(lLepton_isTight,lLepton_ptargsort)")

  lVars.Add("Loose4Lepton_pt", "NlooseLeptons >= 4 ? ROOT::VecOps::Take(LooseLepton_pt, 4) : ROOT::VecOps::Take(LooseLepton_pt, 3)")
  lVars.Add("Loose4Lepton_eta", "NlooseLeptons >= 4 ? ROOT::VecOps::Take(LooseLepton_eta, 4) : ROOT::VecOps::Take(LooseLepton_eta, 3)")
  lVars.Add("Loose4Lepton_phi", "NlooseLeptons >= 4 ? ROOT::VecOps::Take(LooseLepton_phi, 4) : ROOT::VecOps::Take(LooseLepton_phi, 3)")
  lVars.Add("Loose4Lepton_mass", "NlooseLeptons >= 4 ? ROOT::VecOps::Take(LooseLepton_mass, 4) : ROOT::VecOps::Take(LooseLepton_mass, 3)")
  lVars.Add("Loose4Lepton_charge", "NlooseLeptons >= 4 ? ROOT::VecOps::Take(LooseLepton_charge, 4) : ROOT::VecOps::Take(LooseLepton_charge, 3)")
  lVars.Add("Loose4Lepton_ID", "NlooseLeptons >= 4 ? ROOT::VecOps::Take(LooseLepton_ID, 4) : ROOT::VecOps::Take(LooseLepton_ID, 3)")
  lVars.Add("Loose4Lepton_isTight", "NlooseLeptons >= 4 ? ROOT::VecOps::Take(LooseLepton_isTight, 4) : ROOT::VecOps::Take(LooseLepton_isTight, 3)")

  lVars.Add("hasQuarkonia", "hasQuarkoniafunc(Loose4Lepton_pt, Loose4Lepton_eta, Loose4Lepton_phi, Loose4Lepton_mass, Loose4Lepton_ID, Loose4Lepton_charge)")
  lVars.Add("Loose4Lepton_fromZ", "hasZfunc(Loose4Lepton_pt, Loose4Lepton_eta, Loose4Lepton_phi, Loose4Lepton_mass, Loose4Lepton_ID, Loose4Lepton_charge)")
  lVars.Add("hasZ", "Sum(Loose4Lepton_fromZ) > 0")

  lCuts = CutGroup('Lepton Cuts')
  lCuts.Add('NlooseLeptons >= 3', 'NlooseLeptons >= 3')
  lCuts.Add('hasQuarkonia == 0', 'hasQuarkonia == 0') 

  # ------------------ JET Cleaning and JERC ------------------
  jVars = VarGroup('JetCleaningVars')
  
  jVars.Add("Jet_P4", "fVectorConstructor(Jet_pt,Jet_eta,Jet_phi,Jet_mass)")
  #jVars.Add("FatJet_P4", "fVectorConstructor(FatJet_pt,FatJet_eta,FatJet_phi,FatJet_mass)")
  jVars.Add("Jet_EmEF","Jet_neEmEF + Jet_chEmEF")
  jVars.Add("DummyZero","float(0.0)")
  
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
  # This function needs work to understand it, and the output is a vector.
  # I'm not sure it's for PuppiMET, it might actually be for PFMET...
  #metVars.Add("corrMET_pt", "METptfunc(METcorr, METyr, isMC, cleanMET_pt, cleanMET_phi, PV_npvs)")
  #metVars.Add("corrMET_phi", "METphifunc(METcorr, METyr, isMC, cleanMET_pt, cleanMET_phi, PV_npvs)")

  metCuts = CutGroup('METCuts')
  metCuts.Add("Pass PuppiMET pt > 50", "corrMET_pt > 50")

  # -------------------------- TRIGGERS -------------------------------
  TrigVars = VarGroup('Triggers Vars')

  #should double check
  TrigVars.Add('NMu20', 'Sum(LooseMu_pt > 20)')
  TrigVars.Add('NMu10LT20', 'Sum(LooseMu_pt > 10 && LooseMu_pt <= 20)')
  TrigVars.Add('NMu25', 'Sum(LooseMu_pt > 25)')
  TrigVars.Add('NMu12LT23', 'Sum(LooseMu_pt > 15 && LooseMu_pt <= 25)')
  TrigVars.Add('NEle23', 'Sum(LooseEl_pt > 25)')
  TrigVars.Add('NEle12LT23', 'Sum(LooseEl_pt > 15 && LooseEl_pt <= 25)')
  TrigVars.Add('NEle30', 'Sum(LooseEl_pt > 33 && Electron_cutBased[looseElectrons == true] >= 4)')
  TrigVars.Add('NMu20forTau', 'Sum(LooseMu_pt > 23 && abs(LooseMu_eta) < 2.1)')
  TrigVars.Add('NEle24forTau', 'Sum(LooseEl_pt > 27 && Electron_cutBased[looseElectrons == true] >= 4 && abs(LooseEl_eta) < 2.1)') # cutbased >= 4 is tight
  TrigVars.Add('NTauMuTrig', 'Sum(LooseTau_pt > 30 && Tau_idDeepTau2018v2p5VSjet[looseTau == true] >= 4 && abs(LooseTau_eta) < 2.1)')
  TrigVars.Add('NTauEleTrig', 'Sum(LooseTau_pt > 33 && Tau_idDeepTau2018v2p5VSjet[looseTau == true] >= 4 && abs(LooseTau_eta) < 2.1)')
  TrigVars.Add('NTauMedTrig', 'Sum(LooseTau_pt > 38 && Tau_idDeepTau2018v2p5VSjet[looseTau == true] >= 5 && abs(LooseTau_eta) < 2.1)')
  TrigVars.Add('NTauLooseTrig', 'Sum(LooseTau_pt > 185 && Tau_idDeepTau2018v2p5VSjet[looseTau == true] >= 4 && abs(LooseTau_eta) < 2.1)')

  TrigVars.Add('isMuMuTrig', '(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 == 1) && (NMu20 >= 2 || (NMu20 >= 1 && NMu10LT20 >= 1))') # in muon dataset
  TrigVars.Add('isMuElTrig', '''(((HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ == 1) && (NMu25 >= 1 && (NEle23 >= 1 || NEle12LT23 >= 1))) 
                                || ((HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ == 1) && (NEle23 >= 1 && (NMu25 >= 1 || NMu12LT23 >= 1))))''') # in muonEG dataset
  TrigVars.Add('isMuTauTrig', '(HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1 == 1) && (NMu20forTau >= 1 && NTauMuTrig >= 1)') # in muon dataset
  TrigVars.Add('isElElTrig', '(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL == 1) && (NEle23 >= 2 || (NEle23 >= 1 && NEle12LT23 >= 1))') # in Egamma dataset
  TrigVars.Add('isElTauTrig', '(HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1 == 1) && (NEle24forTau >= 1 && NTauEleTrig >= 1)') # in Egamma dataset
  TrigVars.Add('isTauTauTrig', '(HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1 == 1) && (NTauMedTrig >= 2)') # in tau dataset
  TrigVars.Add('isSingleMuTrig', '(HLT_IsoMu24 == 1) && (NMu25 >= 1)') # in muon dataset
  TrigVars.Add('isSingleEleTrig', '(HLT_Ele30_WPTight_Gsf == 1) && (NEle30 >= 1)') # in Egamma dataset
  TrigVars.Add('isSingleTauTrig', '(HLT_LooseDeepTauPFTauHPS180_L2NN_eta2p1 == 1) && (NTauLooseTrig >= 1)') # in tau dataset

  TrigVars.Add('SortMuPD','isMuMuTrig == 1 || isSingleMuTrig == 1 || isMuTauTrig == 1')
  TrigVars.Add('SortMuEGPD','SortMuPD == 0 && isMuElTrig == 1')
  TrigVars.Add('SortEGPD','SortMuPD == 0 && SortMuEGPD == 0 && (isElElTrig == 1 || isSingleEleTrig == 1 || isElTauTrig == 1)')
  TrigVars.Add('SortTauPD','SortMuPD == 0 && SortMuEGPD == 0 && SortEGPD == 0 && (isTauTauTrig == 1 || isSingleTauTrig == 1)')
  
  TrigVars.Add('passesMuPD', 'SortMuPD && (isMC == 1 || dataset == 1)')
  TrigVars.Add('passesMuEGPD', 'SortMuEGPD && (isMC == 1 || dataset == 2)')
  TrigVars.Add('passesEGPD', 'SortEGPD && (isMC == 1 || dataset == 3)')
  TrigVars.Add('passesTauPD', 'SortTauPD && (isMC == 1 || dataset == 4)')

  TrigCuts = CutGroup('Trigger Cuts')
  TrigCuts.Add('Dataset Filter', 'passesMuPD == 1 || passesMuEGPD == 1 || passesEGPD == 1 || passesTauPD == 1')

  # ------------------ Jet Calculations ------------------ 
  jVars.Add("minDR_lepsJets","minDR_jets_gtau(nJet, cleanJet_eta, cleanJet_phi, NlooseLeptons, Loose4Lepton_eta, Loose4Lepton_phi)")
  jVars.Add("goodcleanJets", "cleanJet_pt > 30 && abs(cleanJet_eta) < 2.5 && Jet_jetId > 1 && minDR_lepsJets > 0.4 ")
  jVars.Add("NgoodcleanJets", "Sum(goodcleanJets)")

  jVars.Add("gcJet_HT","Sum(cleanJet_pt[goodcleanJets == true])")
  jVars.Add("gcJet_pt_unsort", "cleanJet_pt[goodcleanJets == true]")
  jVars.Add("gcJet_ptargsort","ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(gcJet_pt_unsort))")
  
  jVars.Add("gcJet_pt","reorder(gcJet_pt_unsort,gcJet_ptargsort)")
  jVars.Add("gcJet_eta", "reorder(cleanJet_eta[goodcleanJets == true],gcJet_ptargsort)")
  jVars.Add("gcJet_phi", "reorder(cleanJet_phi[goodcleanJets == true],gcJet_ptargsort)")
  jVars.Add("gcJet_mass", "reorder(cleanJet_mass[goodcleanJets == true],gcJet_ptargsort)")
  
  jVars.Add("gcJet_vetomap", "jetvetofunc(jetvetocorr, gcJet_eta, gcJet_phi)")
  if (isMC):
      jVars.Add("gcJet_hflav", "reorder(Jet_hadronFlavour[goodcleanJets == true],gcJet_ptargsort)")
 
  if year != '2024' and year != '2025':
    jVars.Add("gcJet_BTag", "reorder(Jet_btagPNetB[goodcleanJets == true],gcJet_ptargsort)")
  else:
    jVars.Add("gcJet_BTag", "reorder(Jet_btagUParTAK4B[goodcleanJets == true],gcJet_ptargsort)")
  jVars.Add("gcJet_BTagL", "gcJet_BTag > BTagL") 
  jVars.Add("NJets_BTagL", "Sum(gcJet_BTagL)")
  
  jVars.Add("gcBJet_eta", "gcJet_eta[gcJet_BTagL]")
  jVars.Add("gcBJet_phi", "gcJet_phi[gcJet_BTagL]")
  jVars.Add("gcBJet_pt", "gcJet_pt[gcJet_BTagL]")
  jVars.Add("gcBJet_mass", "gcJet_mass[gcJet_BTagL]")

  jVars.Add("gcJet_ht", "Sum(gcJet_pt)")

  jCuts = CutGroup('JetCuts')
  jCuts.Add('Event has no vetoed jets', 'Sum(gcJet_vetomap) == 0')
  jCuts.Add('NgoodcleanJets >= 2', 'NgoodcleanJets >= 2')
  jCuts.Add('2 B Jets Pass (Loose)', 'NJets_BTagL >= 2')

  # ---------------- Nonprompt Weight ----------------
  npVars = VarGroup('NonpromptVars')
  npVars.Add("nonpromptWeightMCMC", "neononpromptweight(ptbins, ptbins, etabinsel, etabinsmu, ePromptRatesMC, mPromptRatesMC, tPromptRatesMC, eFakeRatesMC, mFakeRatesMC, tFakeRatesMC, Loose4Lepton_pt, Loose4Lepton_eta, Loose4Lepton_ID, Loose4Lepton_isTight)")
  npVars.Add("nonpromptWeightDD", "neononpromptweight(ptbins, ptbins_datafake, etabinsel, etabinsmu, ePromptRatesData, mPromptRatesData, tPromptRatesMC, eFakeRatesData, mFakeRatesData, tFakeRatesMC, Loose4Lepton_pt, Loose4Lepton_eta, Loose4Lepton_ID, Loose4Lepton_isTight)")
  npVars.Add("nonpromptWeightDMC", "neononpromptweight(ptbins, ptbins, etabinsel, etabinsmu, ePromptRatesData, mPromptRatesData, tPromptRatesMC, eFakeRatesMC, mFakeRatesMC, tFakeRatesMC, Loose4Lepton_pt, Loose4Lepton_eta, Loose4Lepton_ID, Loose4Lepton_isTight)")

  # ------------------ Results ------------------
  manualVars = VarGroup('manualVars')
  manualVars.Add('manual', 'funcmassdiff(Loose4Lepton_pt, Loose4Lepton_eta, Loose4Lepton_phi, Loose4Lepton_mass, gcBJet_pt, gcBJet_eta, gcBJet_phi, gcBJet_mass)')
  manualVars.Add('B1finalPx', 'manual[0]')
  manualVars.Add('B1finalPy', 'manual[1]')
  manualVars.Add('B1finalPz', 'manual[2]')
  manualVars.Add('B1finalM', 'manual[3]')
  manualVars.Add('B2finalPx', 'manual[4]')
  manualVars.Add('B2finalPy', 'manual[5]')
  manualVars.Add('B2finalPz', 'manual[6]')
  manualVars.Add('B2finalM', 'manual[7]')

  rframeVars = VarGroup('restFrameVars')
  
  rframeVars.Add('VLQ', 'B_rfc.return_doubles(rdfslot_, Loose4Lepton_pt, Loose4Lepton_eta, Loose4Lepton_phi, Loose4Lepton_mass, Loose4Lepton_charge, gcBJet_pt, gcBJet_eta, gcBJet_phi, gcBJet_mass, MET_pt, MET_phi)')
  rframeVars.Add('VLQ_BBbar_mass', 'VLQ[0]')
  rframeVars.Add('VLQ_BBbar_cosDecayAngle', 'VLQ[1]')
  rframeVars.Add('VLQ_BBbar_deltaPhiDecayAngle', 'VLQ[2]')
  rframeVars.Add('VLQ_B_mass', 'VLQ[3]')
  rframeVars.Add('VLQ_B_cosDecayAngle', 'VLQ[4]')
  rframeVars.Add('VLQ_B_deltaPhiDecayAngle', 'VLQ[5]')
  rframeVars.Add('VLQ_Bbar_mass', 'VLQ[6]')
  rframeVars.Add('VLQ_Bbar_cosDecayAngle', 'VLQ[7]')
  rframeVars.Add('VLQ_Bbar_deltaPhiDecayAngle', 'VLQ[8]')
  rframeVars.Add('VLQ_tau11_mass', 'VLQ[9]')
  rframeVars.Add('VLQ_tau12_mass', 'VLQ[10]')
  rframeVars.Add('VLQ_tau21_mass', 'VLQ[11]')
  rframeVars.Add('VLQ_tau22_mass', 'VLQ[12]')
  rframeVars.Add('VLQ_BBbar_DeltaPhiVisible', 'VLQ[13]')
  rframeVars.Add('VLQ_BBbar_DeltaPhiDecayVisible', 'VLQ[14]')
  rframeVars.Add('VLQ_BBbar_DeltaPhiBoostVisible', 'VLQ[15]')
  rframeVars.Add('VLQ_BBbar_VisibleShape', 'VLQ[16]')

  # # -------------------------------------

  nodeToPlot = a.Apply([flagCuts, gjsonVars, gjsonCuts, tVars, eandmuVars, lVars, lCuts, TrigVars, TrigCuts])

  # # Solution to cleanJets() problem:
  # #       The analyzer .Apply() calls the analyzer .Define().  This .Define() calls self._collectionOrg.CollectionDefCheck(var, newNode).
  # #  This method executes this line: if re.search(r"\b" + re.escape(c+'s') + r"\b", action_str) and (c+'s' not in self._builtCollections):
  # #                                       print ('MAKING %ss for %s'%(c,action_str))
  # #  Apparently somethings get discarded from this _collectionOrg?
  # #  Instead force the .Apply() from the ActiveNode because the node .Apply() is better.
  
  newNode = a.ActiveNode.Apply(jVars)
  a.SetActiveNode(newNode)
  a.Apply([jCuts, metVars, metCuts, npVars, manualVars])#, rframeVars])
  
  allColumns = a.GetColumnNames()
     
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
     if col == "assignleps" or col == "pnetoutput" or col == "t_output" or col == "Bprime_output" or col == "hasZoutput": continue
     if col.startswith("PS") or col.startswith("Tk") or col.startswith("Trig"): continue
     if col.startswith("nCorr") or col.startswith("nFsr") or col.startswith("Other"): continue
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
     

  finalFile = 'Nonprompt'+ sample + era + ver + "_" + year + "_" + str(testNum1) + ".root";
 
  mode = 'RECREATE'
  if jesvar != "Nominal":
    mode = 'UPDATE'
  #print('\n(1)\n') 
  #sys.setprofile(trace_calls)
  a.Snapshot(columns, finalFile, "Events_"+jesvar, lazy=False, openOption=mode, saveRunChain=True)
  #print('\n(2)\n')
  if jesvar == "Nominal":
    print("Cut statistics:")
    rep = a.DataFrame.Report()
    rep.Print()

  print("--------- Analysis End ---------")
    
  a.Close()

analyze("Nominal")
