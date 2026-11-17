from TIMBER.Analyzer import *
from TIMBER.Tools.Common import *
import ROOT
from ROOT import TFile
import sys, os
import gc
import array


gc.disable()

#from TIMBER.Tools.RestFramesHandler import load_restframes
import correctionlib
correctionlib.register_pyroot_binding()

sys.path.append('../../')
sys.path.append('../../../')

ROOT.gInterpreter.ProcessLine('#pragma GCC diagnostic ignored "-Wdeprecated-declarations"') #Command to ignore certain warning messages

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

ROOT.gInterpreter.Declare("""
ROOT::VecOps::RVec<float> DeltaR(
    const ROOT::VecOps::RVec<float>& jet_eta,
    const ROOT::VecOps::RVec<float>& jet_phi,
    float lep_eta,
    float lep_phi
) {
    return ROOT::VecOps::Map(
        jet_eta,
        jet_phi,
        [&](float eta, float phi) {
            float dphi = ROOT::VecOps::DeltaPhi(phi, lep_phi);
            float deta = eta - lep_eta;
            return std::sqrt(dphi*dphi + deta*deta);
        }
    );
}
""")

ROOT.gInterpreter.Declare("""
ROOT::VecOps::RVec<int> zero_and_return(
    const ROOT::VecOps::RVec<int>& vector,
    int index
) {
    ROOT::VecOps::RVec<int> new_vec = vector;
    new_vec[index] = 0;
    return new_vec;
}
""")

ROOT.gInterpreter.Declare("""
ROOT::VecOps::RVec<float> zero_and_return(
    const ROOT::VecOps::RVec<float>& vector,
    float index
) {
    ROOT::VecOps::RVec<float> new_vec = vector;
    new_vec[index] = 0.0;
    return new_vec;
}
""")

# Enable using 4 threads
ROOT.ROOT.EnableImplicitMT(num_threads)

# load rest frames handler
#handler_name = 'Bprime_handler_new.cc'
#class_name = 'Bprime_RestFrames_Container_new'
#load_restframes(num_threads, handler_name, class_name, 'B_rfc')

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
  bool isSM = """+str(isSM).lower()+"""; 
  bool isSE = """+str(isSE).lower()+"""; 
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
  elif '2024' in year:
      jsonfile = jsonfile + "Cert_Collisions2024_378981_386951_Golden.json"
  elif '2025' in year:
      jsonfile = jsonfile + "Cert_Collisions2025_391658_398903_Golden.json"
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
  jecyr = {'2022':"Summer22_22Sep2023",'2022EE':"Summer22EE_22Sep2023",'2023':"Summer23Prompt23",'2023BPix':"Summer23BPixPrompt23",'2024':"Summer24Prompt24", '2025':"Summer24Prompt24"}
  jecver = {'2022':"V4",'2022EE':"V4",'2023':"V4",'2023BPix':"V4",'2024':"V3",'2025':"V3"}   
  jetvetoname = {'2022':"Summer22_23Sep2023_RunCD_V1",'2022EE':"Summer22EE_23Sep2023_RunEFG_V1",'2023':"Summer23Prompt23_RunC_V1",'2023BPix':"Summer23BPixPrompt23_RunD_V1",'2024':"Summer24Prompt24_RunBCDEFGHI_V1",'2025':"Summer24Prompt24_RunBCDEFGHI_V1"}      
  jmeyrstr = {'2022':yrstr['2022'],'2022EE':yrstr['2022EE'],'2023':yrstr['2023'],'2023BPix':yrstr['2023BPix'],'2024':yrstr['2024'],'2025':yrstr['2024']} # yes, really use 24 for 25

  if not isMC: #is DATA
    jmeyrstr['2025'] = 'Run3-25Prompt-Winter25-NanoAODv15'
    jecyr['2025'] = "Winter25Prompt25"
    jetvetoname['2025'] = "Winter25Prompt25_RunCDEFG_V1"    

  ROOT.gInterpreter.Declare("""
  float BTagL = """+str(BTagL[year])+""";
  string yrstr = \""""+yrstr[year]+"""\";
  string jmetag = \""""+jmetags[year]+"""\";
  string jecyr = \""""+jecyr[year]+"""\";
  string jecver = \""""+jecver[year]+"""\";
  string jmeyrstr = \""""+jmeyrstr[year]+"""\";
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

  auto path = (jmeyrstr == "Run3-25Prompt-Winter25-NanoAODv15")
            ? "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/2026-06-05/jetid.json.gz"
            : "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/" + jmeyrstr + "/" + jmetag + "/jetid.json.gz";

  auto jetidcorrset = correction::CorrectionSet::from_file(path);

  auto jetidAK4Tcorr = jetidcorrset->at("AK4PUPPI_Tight");
  auto jetidAK4TLcorr = jetidcorrset->at("AK4PUPPI_TightLeptonVeto");
  auto jetidAK8Tcorr = jetidcorrset->at("AK8PUPPI_Tight");
  auto jetidAK8TLcorr = jetidcorrset->at("AK8PUPPI_TightLeptonVeto");
  auto ak4corr = ak4corrset->compound().at(jecyr+"_"+jecver+"_DATA_L1L2L3Res_AK4PFPuppi");
  auto ak4corrL1 = ak4corrset->at(jecyr+"_"+jecver+"_DATA_L1FastJet_AK4PFPuppi");
  auto ak8corr = ak8corrset->compound().at(jecyr+"_"+jecver+"_DATA_L1L2L3Res_AK8PFPuppi");
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
  eandmuVars.Add('preselElectrons', 'Electron_pt > 10 && abs(Electron_eta) < 2.5')
  eandmuVars.Add('looseElectrons', 'preselElectrons == 1 && (Electron_cutBased >= 2)')
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
  eandmuVars.Add('preselMuons', 'Muon_pt >= 15 && abs(Muon_eta) < 2.4')
  eandmuVars.Add('looseMuons', 'preselMuons == 1 && Muon_looseId == 1 && Muon_pfIsoId >= 2')
  eandmuVars.Add('NlooseMuons', 'Sum(looseMuons)')
  #eandmuVars.Add('goodMuons', 'Muon_pt >= 15 && abs(Muon_eta) < 2.4 && Muon_mediumId == 1 && abs(Muon_dxy) < 0.2 && abs(Muon_dz) < 0.5 && Muon_pfIsoId >= 3')
  #eandmuVars.Add('NgoodMuons', 'Sum(goodMuons)')

  eandmuVars.Add('LooseMu_pt','Muon_pt[looseMuons == true]')
  eandmuVars.Add('LooseMu_eta','Muon_eta[looseMuons == true]')
  eandmuVars.Add('LooseMu_phi', 'Muon_phi[looseMuons == true]')
  eandmuVars.Add('LooseMu_mass', 'Muon_mass[looseMuons == true]')
  eandmuVars.Add('LooseMu_charge', 'Muon_charge[looseMuons == true]')
  eandmuVars.Add('LooseMu_ID', 'ROOT::VecOps::RVec<int>(NlooseMuons, 13)')
  eandmuVars.Add('LooseMu_isTight','Muon_mediumId[looseMuons == true] == 1 && Muon_pfIsoId[looseMuons == true] >= 3 && abs(Muon_dxy) < 0.2 && abs(Muon_dz) < 0.5')

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
  lVars.Add("LooseLepton_eta","reorder(lLepton_eta,lLepton_ptargsort)")
  lVars.Add("LooseLepton_phi","reorder(lLepton_phi,lLepton_ptargsort)")
  lVars.Add("LooseLepton_isTight","reorder(lLepton_isTight,lLepton_ptargsort)")


  lCuts = CutGroup('Lepton Cuts')
  lCuts.Add('NlooseLeptons >= 1', 'NlooseLeptons >= 1')


  # ------------------ JET Cleaning and JERC ------------------
  jVars = VarGroup('JetCleaningVars')

  jVars.Add("Jet_P4", "fVectorConstructor(Jet_pt,Jet_eta,Jet_phi,Jet_mass)")
  jVars.Add("Jet_EmEF","Jet_neEmEF + Jet_chEmEF")
  jVars.Add("DummyZero","float(0.0)")

  if year == '2024' or year == '2025' or year =='2023':
    jVars.Add("Jet_jetId","jetidfunc(jetidAK4Tcorr,jetidAK4TLcorr,Jet_eta,Jet_chHEF,Jet_neHEF,Jet_chEmEF,Jet_neEmEF,Jet_muEF,Jet_chMultiplicity,Jet_neMultiplicity)")
    jVars.Add("FatJet_jetId","fatjetidfunc(jetidAK8Tcorr,jetidAK8TLcorr,FatJet_eta,FatJet_chHEF,FatJet_neHEF,FatJet_chEmEF,FatJet_neEmEF,FatJet_muEF,FatJet_chMultiplicity,FatJet_neMultiplicity)")

  jVars.Add("cleanedJets", "cleanJetsData(run,debug,year,ak4corr,ak4corrL1,ak8corr,Jet_P4,Jet_rawFactor,Jet_muonSubtrFactor,Jet_area,Jet_EmEF,Jet_jetId,Jet_P4,Jet_jetId,Rho_fixedGridRhoFastjetAll,DummyZero,DummyZero)") # muon and EM factors unused in this call, args 16-17 are dummies
  jVars.Add("cleanMets", "cleanJetsData(run,debug,year,ak4corr,ak4corrL1,ak8corr,Jet_P4,Jet_rawFactor,Jet_muonSubtrFactor,Jet_area,Jet_EmEF,Jet_jetId,Jet_P4,Jet_jetId,Rho_fixedGridRhoFastjetAll,RawPuppiMET_pt,RawPuppiMET_phi)") # lepton args unused in this call, args 16-17 are dummies
  jVars.Add("cleanJet_pt", "cleanedJets[0]")
  jVars.Add("cleanJet_eta", "cleanedJets[1]")
  jVars.Add("cleanJet_phi", "cleanedJets[2]")
  jVars.Add("cleanJet_mass", "cleanedJets[3]")

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

  # ------------------ Jet Calculations ------------------
  jVars.Add("minDR_lepsJets","minDR_jets_gtau(nJet, cleanJet_eta, cleanJet_phi, NlooseLeptons, LooseLepton_eta, LooseLepton_phi)")
  jVars.Add("goodcleanJets", "cleanJet_pt > 30 && abs(cleanJet_eta) < 2.5 && Jet_jetId > 1 && minDR_lepsJets > 0.4 ")
  jVars.Add("NgoodcleanJets", "Sum(goodcleanJets)")

  jVars.Add("gcJet_pt_unsort", "cleanJet_pt[goodcleanJets == true]")
  jVars.Add("gcJet_ptargsort","ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(gcJet_pt_unsort))")

  jVars.Add("gcJet_pt","reorder(gcJet_pt_unsort,gcJet_ptargsort)")
  jVars.Add("gcJet_eta", "reorder(cleanJet_eta[goodcleanJets == true],gcJet_ptargsort)")
  jVars.Add("gcJet_phi", "reorder(cleanJet_phi[goodcleanJets == true],gcJet_ptargsort)")
  jVars.Add("gcJet_mass", "reorder(cleanJet_mass[goodcleanJets == true],gcJet_ptargsort)")

  jVars.Add("gcJet_vetomap", "jetvetofunc(jetvetocorr, gcJet_eta, gcJet_phi)")
  jVars.Add("gcJet_PNet", "reorder(Jet_btagPNetB[goodcleanJets == true],gcJet_ptargsort)")
  jVars.Add("gcJet_PNetL", "gcJet_PNet > BTagL")
  jVars.Add("NJets_PNetL", "Sum(gcJet_PNetL)")

  jVars.Add("gcBJet_eta", "gcJet_eta[gcJet_PNetL]")
  jVars.Add("gcBJet_phi", "gcJet_phi[gcJet_PNetL]")
  jVars.Add("gcBJet_pt", "gcJet_pt[gcJet_PNetL]")
  jVars.Add("gcBJet_mass", "gcJet_mass[gcJet_PNetL]")

  jCuts = CutGroup('JetCuts')
  jCuts.Add('Events has no vetoed jets','Sum(gcJet_vetomap) == 0')
  jCuts.Add('NgoodcleanJets >= 1', 'NgoodcleanJets >= 1')
  #jCuts.Add('1 B Jet Pass (Loose)', 'NJets_PNetL >= 1')



  # # --------------- Prompt Rate Calculations -------------
  preTagCuts = CutGroup("preTagCuts")
  preTagCuts.Add("Flavor_XOR_Cut", "(NlooseElecs > 0 && NlooseMuons > 0) || (NlooseElecs > 0 && NlooseTau > 0) || (NlooseMuons > 0 && NlooseTau > 0)")

  tagDefVars = VarGroup("tagTestVariables")
  tagDefVars.Add("TagMu","Muon_pt > 30 && abs(Muon_eta)<2.4 && Muon_tightId==1 && Muon_pfIsoId>=3 && abs(Muon_dz) < 0.5 && Muon_dxy < 0.2")
  tagDefVars.Add("ProbeMuNotTag","looseMuons == 1 && (Muon_tightId == 0 || Muon_pfIsoId < 3 || Muon_pt < 30)")
  tagDefVars.Add("nTagMu","(int) Sum(TagMu)")
  tagDefVars.Add("nProbeMuNotTag","(int) Sum(ProbeMuNotTag)")
  tagDefVars.Add("passTagMuTrig","(HLT_IsoMu27 == 1)")
  tagDefVars.Add("tagMu_idx", "ROOT::VecOps::Nonzero(TagMu)")
  tagDefVars.Add("probeMu_idx", "ROOT::VecOps::Nonzero(looseMuons)")
  tagDefVars.Add("probeTau_idx", "ROOT::VecOps::Nonzero(looseTau)")
  
  tagDefVars1 = VarGroup("split 1")
  tagDefVars1.Add("TagEl","Electron_pt > 40 && (abs(Electron_eta)<1.442 || (abs(Electron_eta)>1.566 && abs(Electron_eta)<2.5)) && Electron_cutBased>=4 ")
  tagDefVars1.Add("ProbeElNotTag","looseElectrons == 1 && (Electron_pt < 40 || Electron_cutBased < 4)")
  tagDefVars1.Add("nTagEl","(int) Sum(TagEl)")
  tagDefVars1.Add("nProbeElNotTag","(int) Sum(ProbeElNotTag)")
  tagDefVars1.Add("passTagElTrig","(HLT_Ele35_WPTight_Gsf == 1 || HLT_Ele38_WPTight_Gsf == 1 || HLT_Photon200)")
  tagDefVars1.Add("tagEl_idx", "ROOT::VecOps::Nonzero(TagEl)")
  tagDefVars1.Add("probeEl_idx", "ROOT::VecOps::Nonzero(looseElectrons)")

  tagDefVars2 = VarGroup("split 2")
  tagDefVars2.Add("is_tau_mu", "isSM && passTagMuTrig && NlooseTau == 1 && nTagMu == 1 && NlooseElecs == 0 && nProbeMuNotTag == 0")
  tagDefVars2.Add("is_tau_el", "isSE && passTagElTrig && NlooseTau == 1 && nTagEl == 1 && NlooseMuons == 0 && nProbeElNotTag == 0")
  tagDefVars2.Add("is_mu_el", "isSE && passTagElTrig && NlooseMuons == 1 && nTagEl == 1 && NlooseTau == 0 && nProbeElNotTag == 0")
  tagDefVars2.Add("is_el_mu", "isSM && passTagMuTrig && NlooseElecs == 1 && nTagMu == 1 && NlooseTau == 0 && nProbeMuNotTag == 0")
  tagDefVars2.Add('cnLooseTau_isTight','goodTaue == 1 && Tau_idDeepTau2018v2p5VSjet >= 4')
  tagDefVars2.Add('cnLooseEl_isTight','preselElectrons == 1 && Electron_mvaIso_WP80 == 1')
  tagDefVars2.Add('cnLooseMu_isTight','preselMuons == 1 && Muon_mediumId == 1 && Muon_pfIsoId >= 3')
  
  cndtagDefVars = VarGroup("Conditional tag variables")
  cndtagDefVars.Add("probe_idx", "(is_tau_mu || is_tau_el) ? probeTau_idx[0] : (is_mu_el ? probeMu_idx[0] : (is_el_mu ? probeEl_idx[0] : -1))")
  cndtagDefVars.Add("probe_pt", "(is_tau_mu || is_tau_el) ? Tau_pt[probe_idx] : (is_mu_el ? Muon_pt[probe_idx] : (is_el_mu ? Electron_pt[probe_idx] : -9))")
  cndtagDefVars.Add("probe_eta","(is_tau_mu || is_tau_el) ? Tau_eta[probe_idx] : (is_mu_el ? Muon_eta[probe_idx]: (is_el_mu ? Electron_eta[probe_idx] : -9))")
  cndtagDefVars.Add("probe_abseta", "abs(probe_eta)")
  cndtagDefVars.Add("probe_phi","(is_tau_mu || is_tau_el) ? Tau_phi[probe_idx] : (is_mu_el ? Muon_phi[probe_idx] : (is_el_mu ? Electron_phi[probe_idx] : -9))")
  cndtagDefVars.Add("probe_flavor","(is_tau_mu || is_tau_el) ? 15 : (is_mu_el ? 13: (is_el_mu ? 11 : -1))")

  cndtagDefVars.Add("tag_idx", "(is_tau_mu || is_el_mu) ? tagMu_idx[0] : ((is_tau_el || is_mu_el) ? tagEl_idx[0] : -1)")
  cndtagDefVars.Add("tag_eta","(is_tau_mu || is_el_mu) ? Muon_eta[tag_idx] : ((is_tau_el || is_mu_el) ? Electron_eta[tag_idx] : -9)")
  cndtagDefVars.Add("tag_phi","(is_tau_mu || is_el_mu) ? Muon_phi[tag_idx] : ((is_tau_el || is_mu_el) ? Electron_phi[tag_idx] : -9)")
  cndtagDefVars.Add("tag_flavor","(is_tau_mu || is_el_mu) ? 13 : ((is_tau_el || is_mu_el) ? 11 : -1)")

  cndtagDefVars.Add("dR_jet_tag", "DeltaR(gcJet_eta, gcJet_phi, tag_eta, tag_phi)")
  cndtagDefVars.Add("dR_jet_probe", "DeltaR(gcJet_eta, gcJet_phi, probe_eta, probe_phi)")
  cndtagDefVars.Add("avoidant_jet", "((dR_jet_tag > 0.4) && (dR_jet_probe > 0.4))")
  cndtagDefVars.Add("n_avoidant_jet", "(int) Sum(avoidant_jet)")
  cndtagDefVars.Add("n_valid_channels", "(int)(is_tau_mu + is_tau_el + is_mu_el + is_el_mu)")
  
  idxtagDefVars = VarGroup("Index-related tag variables")
  #idxtagDefVars.Add("n_tags_in_channel", "is_tau_mu  ? (int)tagMu_idx.size()  : (is_tau_el ? (int)tagEl_idx.size()  : (is_mu_el  ? (int)tagEl_idx.size()  : (is_el_mu  ? (int)tagMu_idx.size()  : 0)))")
  #idxtagDefVars.Add("n_probes_in_channel","is_tau_mu  ? (int)probeTau_idx.size() : (is_tau_el ? (int)probeTau_idx.size() : (is_mu_el  ? (int)probeMu_idx.size()  : (is_el_mu  ? (int)probeEl_idx.size()  : 0)))")
  idxtagDefVars.Add("probe_isTight","(is_tau_mu || is_tau_el) ? cnLooseTau_isTight[probe_idx] : (is_mu_el ? cnLooseMu_isTight[probe_idx]: (is_el_mu ? cnLooseEl_isTight[probe_idx] : 0))")

  postTagCuts = CutGroup("postTagCuts")
  postTagCuts.Add("has_tag", "(tag_flavor != -1)")
  postTagCuts.Add("has_avoidant_jet","n_avoidant_jet > 0")
  tagProbeUniquenessCuts = CutGroup("tagProbeUniquenessCuts")
  tagProbeUniquenessCuts.Add("unique_channel", "n_valid_channels == 1")
  #tagProbeUniquenessCuts.Add("unique_tag",     "n_tags_in_channel == 1")
  #tagProbeUniquenessCuts.Add("unique_probe",   "n_probes_in_channel == 1")

  #---------Z Boson and B Jet

  zBJTagDefVars = VarGroup("zBJTagVariables")

  # Channel flags — reuse existing idx vectors
  zBJTagDefVars.Add("tag_mu_idx","nTagMu > 0 ? tagMu_idx[0] : -1")
  zBJTagDefVars.Add("tag_el_idx","nTagEl > 0 ? tagEl_idx[0] : -1")
  # In the loose leptons pass/fail lists, remove the tag of the same flavor. 
  zBJTagDefVars.Add("nontag_mu_probes","zero_and_return(looseMuons,tag_mu_idx)")
  zBJTagDefVars.Add("nontag_el_probes","zero_and_return(looseElectrons,tag_el_idx)")
  zBJTagDefVars.Add("probeMuNonTag_idx", "ROOT::VecOps::Nonzero(nontag_mu_probes)")
  zBJTagDefVars.Add("probeElNonTag_idx", "ROOT::VecOps::Nonzero(nontag_el_probes)")

  # choose flavor
  zBJTagDefVars.Add("is_mu_mu", "isSM && passTagMuTrig && nTagMu == 1 && Sum(nontag_mu_probes) == 1 && NlooseElecs == 0 && nTagEl == 0")
  zBJTagDefVars.Add("is_el_el", "isSE && passTagElTrig && nTagEl == 1 && Sum(nontag_el_probes) == 1 && NlooseMuons == 0 && nTagMu == 0")
  zBJTagDefVars.Add("n_zBJ_valid_channels", "(int)(is_mu_mu + is_el_el)")

  # Index and flavor selection (mu-mu takes priority, mirroring opposite-flavor logic)
  zBJTagDefVars.Add("zBJ_tag_idx",    "is_mu_mu ? tag_mu_idx : (is_el_el ? tag_el_idx : -1)")
  zBJTagDefVars.Add("zBJ_probe_idx",  "is_mu_mu ? probeMuNonTag_idx[0] : (is_el_el ? probeElNonTag_idx[0] : -1)")
  zBJTagDefVars.Add("zBJ_tag_flavor", "is_mu_mu ? 13 : (is_el_el ? 11 : -1)")

  # Kinematics — branch on flavor once, reuse for mass
  zBJTagDefVars.Add("zBJ_tag_pt",    "(zBJ_tag_flavor == 13) ? Muon_pt[zBJ_tag_idx]    : (is_el_el ? Electron_pt[zBJ_tag_idx] : -9)")
  zBJTagDefVars.Add("zBJ_tag_eta",   "(zBJ_tag_flavor == 13) ? Muon_eta[zBJ_tag_idx]   : (is_el_el ? Electron_eta[zBJ_tag_idx] : -9)")
  zBJTagDefVars.Add("zBJ_tag_phi",   "(zBJ_tag_flavor == 13) ? Muon_phi[zBJ_tag_idx]   : (is_el_el ? Electron_phi[zBJ_tag_idx] : -9)")
  zBJTagDefVars.Add("zBJ_probe_pt",  "(zBJ_tag_flavor == 13) ? Muon_pt[zBJ_probe_idx]  : (is_el_el ? Electron_pt[zBJ_probe_idx] : -9)")
  zBJTagDefVars.Add("zBJ_probe_eta", "(zBJ_tag_flavor == 13) ? Muon_eta[zBJ_probe_idx] : (is_el_el ? Electron_eta[zBJ_probe_idx] : -9)")
  zBJTagDefVars.Add("zBJ_probe_abseta", "abs(zBJ_probe_eta)")
  zBJTagDefVars.Add("zBJ_probe_phi", "(zBJ_tag_flavor == 13) ? Muon_phi[zBJ_probe_idx] : (is_el_el ? Electron_phi[zBJ_probe_idx] : -9)")
  zBJTagDefVars.Add("zBJ_probe_isTight","(is_mu_mu) ? cnLooseMu_isTight[zBJ_probe_idx] : (is_el_el ? cnLooseEl_isTight[zBJ_probe_idx] : -1)")

  # Invariant mass via ROOT::Math::PtEtaPhiMVector with correct particle masses
  zBJTagDefVars.Add("zBJ_pair_mass", "(ROOT::Math::PtEtaPhiMVector(zBJ_tag_pt, zBJ_tag_eta, zBJ_tag_phi, (zBJ_tag_flavor==13 ? 0.10566 : 0.000511)) + ROOT::Math::PtEtaPhiMVector(zBJ_probe_pt, zBJ_probe_eta, zBJ_probe_phi, (zBJ_tag_flavor==13 ? 0.10566 : 0.000511))).M()")

  # B-jet avoidance — require at least one b-jet well separated from both tag and probe
  zBJTagDefVars.Add("zBJ_dR_jet_tag", "DeltaR(gcJet_eta, gcJet_phi, zBJ_tag_eta, zBJ_tag_phi)")
  zBJTagDefVars.Add("zBJ_dR_jet_probe", "DeltaR(gcJet_eta, gcJet_phi, zBJ_probe_eta, zBJ_probe_phi)")
  zBJTagDefVars.Add("zBJ_avoidant_jet", "(zBJ_dR_jet_tag > 0.4) && (zBJ_dR_jet_probe > 0.4)")
  zBJTagDefVars.Add("zBJ_n_avoidant_jet", "(int) Sum(zBJ_avoidant_jet)")

  # Uniqueness counts — same pattern as opposite-flavor
  #zBJTagDefVars.Add("n_zBJ_tags_in_channel", "has_mu_mu ? (int)tagMu_idx.size()   : (has_el_el ? (int)tagEl_idx.size()   : 0)")
  #zBJTagDefVars.Add("n_zBJ_probes_in_channel", "has_mu_mu ? (int)probeMu_idx.size() : (has_el_el ? (int)probeEl_idx.size() : 0)")


  zBJPostTagCuts = CutGroup("zBJPostTagCuts")
  zBJPostTagCuts.Add("zBJ_has_tag",           "zBJ_tag_flavor != -1")
  zBJPostTagCuts.Add("zBJ_unique_channel",    "n_zBJ_valid_channels == 1")
  #zBJPostTagCuts.Add("zBJ_unique_tag",        "n_zBJ_tags_in_channel == 1")
  #zBJPostTagCuts.Add("zBJ_unique_probe",      "n_zBJ_probes_in_channel == 1")
  #zBJPostTagCuts.Add("zBJ_distinct_objects",  "zBJ_tag_idx != zBJ_probe_idx")
  zBJPostTagCuts.Add("zBJ_has_avoidant_jet", "zBJ_n_avoidant_jet > 0")
  zBJPostTagCuts.Add("zBJ_Z_mass_window",     "zBJ_pair_mass > 81 && zBJ_pair_mass < 101")

  #---------------------------------------------------------

  ptbins  = [10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0, 200.0, 9999.0]
  etabinsmu = [0.0, 0.9, 1.2, 2.1, 2.4]
  etabinsel = [-2.5, -2.0, -1.566, -1.444, -0.8, 0.0, 0.8, 1.444, 1.566, 2.0, 2.5]

  nodeToPlot = a.Apply([flagCuts, gjsonVars, gjsonCuts, tVars, eandmuVars, lVars, lCuts])

  # # Solution to cleanJets() problem:
  # #       The analyzer .Apply() calls the analyzer .Define().  This .Define() calls self._collectionOrg.CollectionDefCheck(var, newNode).
  # #  This method executes this line: if re.search(r"\b" + re.escape(c+'s') + r"\b", action_str) and (c+'s' not in self._builtCollections):
  # #                                       print ('MAKING %ss for %s'%(c,action_str))
  # #  Apparently somethings get discarded from this _collectionOrg?
  # #  Instead force the .Apply() from the ActiveNode because the node .Apply() is better.

  newNode = a.ActiveNode.Apply(jVars)
  a.SetActiveNode(newNode)
  a.Apply([jCuts])
  #a.Apply(preTagCuts)
  a.Apply(tagDefVars)
  a.Apply(tagDefVars1)
  a.Apply(tagDefVars2)
  
  divpt = a.GetActiveNode()

  #zbj = a.Duplicate()
  #ttbar = a.Duplicate()
  a.Apply(cndtagDefVars)
  a.Apply(idxtagDefVars)
  a.Apply(postTagCuts)
  a.Apply(metVars)
  a.Apply(metCuts)
  a.Apply(tagProbeUniquenessCuts)
  ttbar_cp = a.GetActiveNode()
  

  a.SetActiveNode(divpt)
  a.Apply(zBJTagDefVars)
  a.Apply(zBJPostTagCuts)
  zbj_cp = a.GetActiveNode()
  
  #a_loose = ttbar
  #a_tight = a_loose.Cut("probe_tight","probe_isTight")
  #allColumns = a.GetColumnNames()

  a.SetActiveNode(ttbar_cp)
  
  ptbins_arr = array.array('d', ptbins)
  etabinsel_arr = array.array('d', etabinsel)
  etabinsmu_arr = array.array('d', etabinsmu)

  finalFile = "DataPromptRate_"+sample + era + ver + "_" + year + "_" + str(testNum1) + ".root";

  outFile = TFile.Open(finalFile, "UPDATE")

  if isSM:
    a.SetActiveNode(ttbar_cp)
    a.Cut("ElectronCut", "is_el_mu")
    hDen_el = a.DataFrame.Histo2D(("hDen_el", "(tt) Electron prompt rate denominator;p_{T} [GeV];#eta",len(ptbins)-1,ptbins_arr,len(etabinsel)-1, etabinsel_arr),"probe_pt","probe_eta")
    a.Cut("ElectronNumCut", "probe_isTight == 1")
    hNum_el = a.DataFrame.Histo2D(("hNum_el", "(tt) Electron prompt rate numerator;p_{T} [GeV];#eta",len(ptbins)-1, ptbins_arr,len(etabinsel)-1, etabinsel_arr),"probe_pt","probe_eta")

    print("Cut statistics ttbar electrons:")
    rep = a.DataFrame.Report()
    rep.Print()

    a.SetActiveNode(ttbar_cp)
    a.Cut("TauCut", "is_tau_mu")
    hDen_tau = a.DataFrame.Histo2D(("hDen_tau", "(tt) Tau prompt rate denominator;p_{T} [GeV];#eta",len(ptbins)-1, ptbins_arr,len(etabinsmu)-1, etabinsmu_arr),"probe_pt","probe_abseta")
    a.Cut("TauNumCut", "probe_isTight == 1")
    hNum_tau = a.DataFrame.Histo2D(("hNum_tau", "(tt) Tau prompt rate numerator;p_{T} [GeV];#eta",len(ptbins)-1, ptbins_arr,len(etabinsmu)-1, etabinsmu_arr),"probe_pt","probe_abseta")

    print("Cut statistics ttbar taus:")
    rep = a.DataFrame.Report()
    rep.Print()

    a.SetActiveNode(zbj_cp)
    a.Cut("MuonCut", "is_mu_mu")
    zhDen_mu = a.DataFrame.Histo2D(("zhDen_mu", "(Z) Muon prompt rate denominator;p_{T} [GeV];#eta",
                                    len(ptbins)-1, ptbins_arr,
                                    len(etabinsmu)-1, etabinsmu_arr),
                                   "zBJ_probe_pt","zBJ_probe_abseta")
    
    a.Cut("MuonNumCut", "zBJ_probe_isTight")
    zhNum_mu = a.DataFrame.Histo2D(("zhNum_mu", "(Z) Muon prompt rate numerator;p_{T} [GeV];#eta",
                                    len(ptbins)-1, ptbins_arr,
                                    len(etabinsmu)-1, etabinsmu_arr),
                                   "zBJ_probe_pt","zBJ_probe_abseta")
    
    print("Cut statistics Z muons:")
    rep = a.DataFrame.Report()
    rep.Print()

    hDen_el.GetValue().SetName("ttbar_h_den_el")
    hNum_el.GetValue().SetName("ttbar_h_num_el")
    hNum_el.GetValue().Write()
    hDen_el.GetValue().Write()
    hDen_tau.GetValue().SetName("ttbar_h_den_tau")
    hNum_tau.GetValue().SetName("ttbar_h_num_tau")
    hNum_tau.GetValue().Write()
    hDen_tau.GetValue().Write()
    zhDen_mu.GetValue().SetName("zbj_h_den_mu")
    zhNum_mu.GetValue().SetName("zbj_h_num_mu")
    zhDen_mu.GetValue().Write()
    zhNum_mu.GetValue().Write()

  elif isSE:

    a.SetActiveNode(ttbar_cp)
    a.Cut("MuonCut", "is_mu_el")
    hDen_mu = a.DataFrame.Histo2D(("hDen_mu", "(tt) Muon prompt rate denominator;p_{T} [GeV];#eta",len(ptbins)-1, ptbins_arr,len(etabinsmu)-1, etabinsmu_arr),"probe_pt","probe_abseta")
    a.Cut("MuonNumCut", "probe_isTight == 1")
    hNum_mu = a.DataFrame.Histo2D(("hNum_mu", "(tt) Muon prompt rate numerators;p_{T} [GeV];#eta",len(ptbins)-1, ptbins_arr,len(etabinsmu)-1, etabinsmu_arr),"probe_pt","probe_abseta")
    
    print("Cut statistics ttbar muons:")
    rep = a.DataFrame.Report()
    rep.Print()

    a.SetActiveNode(ttbar_cp)
    a.Cut("TauCut", "is_tau_el")
    hDen_tau = a.DataFrame.Histo2D(("hDen_tau", "(tt) Tau prompt rate denominator;p_{T} [GeV];#eta",len(ptbins)-1, ptbins_arr,len(etabinsmu)-1, etabinsmu_arr),"probe_pt","probe_abseta")
    a.Cut("TauNumCut", "probe_isTight == 1")
    hNum_tau = a.DataFrame.Histo2D(("hNum_tau", "(tt) Tau prompt rate numerator;p_{T} [GeV];#eta",len(ptbins)-1, ptbins_arr,len(etabinsmu)-1, etabinsmu_arr),"probe_pt","probe_abseta")
    
    print("Cut statistics ttbar taus:")
    rep = a.DataFrame.Report()
    rep.Print()

    a.SetActiveNode(zbj_cp)
    a.Cut("ElectronCut", "is_el_el")
    zhDen_el = a.DataFrame.Histo2D(("zhDen_el", "(Z) Electron prompt rate denominator;p_{T} [GeV];#eta",
                                    len(ptbins)-1, ptbins_arr,
                                    len(etabinsel)-1, etabinsel_arr),
                                   "zBJ_probe_pt","zBJ_probe_eta")

    a.Cut("ElectronNumCut", "zBJ_probe_isTight")
    zhNum_el = a.DataFrame.Histo2D(("zhNum_el", "(Z) Electron prompt rate numerator;p_{T} [GeV];#eta",
                                    len(ptbins)-1, ptbins_arr,
                                    len(etabinsel)-1, etabinsel_arr),
                                   "zBJ_probe_pt","zBJ_probe_eta")
    
    print("Cut statistics Z electrons:")
    rep = a.DataFrame.Report()
    rep.Print()

    hDen_mu.GetValue().SetName("ttbar_h_den_mu")
    hNum_mu.GetValue().SetName("ttbar_h_num_mu")
    hNum_mu.GetValue().Write()
    hDen_mu.GetValue().Write()
    hDen_tau.GetValue().SetName("ttbar_h_den_tau")
    hNum_tau.GetValue().SetName("ttbar_h_num_tau")
    hNum_tau.GetValue().Write()
    hDen_tau.GetValue().Write()
    zhDen_el.GetValue().SetName("zbj_h_den_el")
    zhNum_el.GetValue().SetName("zbj_h_num_el")
    zhDen_el.GetValue().Write()
    zhNum_el.GetValue().Write()
    
  outFile.Close()

  print("--------- Analysis End ---------")

  a.Close()

analyze("Nominal")
