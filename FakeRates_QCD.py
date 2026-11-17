from TIMBER.Analyzer import analyzer, VarGroup, CutGroup
from TIMBER.Tools.Common import CompileCpp
import ROOT
import array


CompileCpp("""
#include <cmath>

ROOT::VecOps::RVec<float> DeltaRVecFloat(
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

ROOT::VecOps::RVec<float> InvariantMassVecFloat(
    const ROOT::VecOps::RVec<float>& jet_pt,
    const ROOT::VecOps::RVec<float>& jet_eta,
    const ROOT::VecOps::RVec<float>& jet_phi,
    const ROOT::VecOps::RVec<float>& jet_mass,
    float lep_pt,
    float lep_eta,
    float lep_phi,
    float lep_mass
) {
    return ROOT::VecOps::Map(
        jet_pt,
        jet_eta,
        jet_phi,
        jet_mass,
        [&](float pt, float eta, float phi, float mass) {
            return (ROOT::Math::PtEtaPhiMVector(pt, eta, phi, mass)+ROOT::Math::PtEtaPhiMVector(lep_pt, lep_eta, lep_phi, lep_mass)).M();
        }
    );
}

float calcMT(float lep_pt, float lep_phi, float met_pt, float met_phi) {
    float dphi = TVector2::Phi_mpi_pi(lep_phi - met_phi);
    return std::sqrt(2.0 * lep_pt * met_pt * (1.0 - std::cos(dphi)));
}
""")
from TIMBER.Analyzer import *
from TIMBER.Tools.Common import *
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


CompileCpp('TIMBER/Framework/include/common.h') # Compile (via gInterpreter) commonly used c++ code
CompileCpp('TIMBER/Framework/Tprime1lep/cleanjet.cc') # Compile Our vlq c++ code
CompileCpp('TIMBER/Framework/Tprime1lep/utilities.cc') # Compile Our vlq c++ code
CompileCpp('TIMBER/Framework/Tprime1lep/lumiMask.cc')
CompileCpp('TIMBER/Framework/Tprime1lep/selfDerived_corrs.cc')
CompileCpp('TIMBER/Framework/Tprime1lep/corr_funcs.cc')
CompileCpp('TIMBER/Framework/Tprime1lep/topographInput.cc')
CompileCpp('TIMBER/Framework/Tprime1lep/manualreco.cc')

debug = False
BTagL = {'2022':0.047,'2022EE':0.0499,'2023':0.0358,'2023BPix':0.0359, '2024':0.0246, '2025':0.0246} #PNet for 22-23BPix, UParT for 2024-2025
yrstr = {'2022':"Run3-22CDSep23-Summer22-NanoAODv12",'2022EE':"Run3-22EFGSep23-Summer22EE-NanoAODv12",'2023':"Run3-23CSep23-Summer23-NanoAODv12",'2023BPix':"Run3-23DSep23-Summer23BPix-NanoAODv12",'2024':"Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15",'2025':"Run3-25Prompt-Summer24-NanoAODv15"}
jmetags = {'2022':'2026-06-05','2022EE':'2026-06-05','2023':'2026-06-05','2023BPix':'2026-06-05','2024':'2026-06-05', '2025':'2026-06-05'} 
jecyr = {'2022':"Summer22_22Sep2023",'2022EE':"Summer22EE_22Sep2023",'2023':"Summer23Prompt23",'2023BPix':"Summer23BPixPrompt23",'2024':"Summer24Prompt24", '2025':"Summer24Prompt24"}
jecver = {'2022':"V4",'2022EE':"V4",'2023':"V4",'2023BPix':"V4",'2024':"V3",'2025':"V3"}   
jetvetoname = {'2022':"Summer22_23Sep2023_RunCD_V1",'2022EE':"Summer22EE_23Sep2023_RunEFG_V1",'2023':"Summer23Prompt23_RunC_V1",'2023BPix':"Summer23BPixPrompt23_RunD_V1",'2024':"Summer24Prompt24_RunBCDEFGHI_V1",'2025':"Summer24Prompt24_RunBCDEFGHI_V1"}      
jmeyrstr = {'2022':yrstr['2022'],'2022EE':yrstr['2022EE'],'2023':yrstr['2023'],'2023BPix':yrstr['2023BPix'],'2024':yrstr['2024'],'2025':yrstr['2024']} # yes, really use 24 for 25

if not isMC: #is DATA
  jecyr['2025'] = "Winter25Prompt25"
  jetvetoname['2025'] = "Winter25Prompt25_RunCDEFG_V1"    
  jmeyrstr['2025'] = 'Run3-25Prompt-Winter25-NanoAODv15'

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


print('============== RUN SETTINGS ===============')
print('isMC = ',isMC)
print('isSM = ',isSM)
print('isSE = ',isSE)
print('isSig = ',isSig)
print('sample = ',sample)
print('era = ',era)
print('jecera = ',jecera)
print('ver = ',ver)
print('region = ',region)
print('dataset = ', dataset)

def analyze(jesvar):
  ROOT.gInterpreter.ProcessLine('string jesvar = "' + jesvar + '"; ')

  # Create analyzer instance
  # is filelist still what you want it be here
  a = analyzer(filelist)
  # Trigger
  trigCuts = CutGroup("trigCuts")
  trigCuts.Add("HLT", "HLT_Mu17 == 1 || HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30 == 1")
  a.Apply(trigCuts)

  #Golden JSON

  gjsonVars = VarGroup('GoldenJsonVars')
  gjsonCuts = CutGroup('GoldenJsonCuts')
  gjsonVars.Add("passesJSON", "goldenjson(myLumiMask, run, luminosityBlock)") # function name, parameters are branches in the file or defined objects
  gjsonCuts.Add("Data passes Golden JSON", "passesJSON == 1")
  a.Apply([gjsonVars,gjsonCuts])

  eandmuVars = VarGroup("eandmuVars")
  eandmuVars.Add('looseElectrons', 'Electron_pt > 10 && abs(Electron_eta) < 2.5 && (Electron_cutBased >= 2)')
  eandmuVars.Add('looseMuons', 'Muon_pt >= 15 && abs(Muon_eta) < 2.4 && Muon_looseId == 1 && Muon_pfIsoId >= 2')
  eandmuVars.Add('LooseEl_isTight','Electron_mvaIso_WP80[looseElectrons == 1] == 1')
  eandmuVars.Add('LooseMu_isTight','Muon_mediumId[looseMuons == 1] == 1 && Muon_pfIsoId[looseMuons == 1] >= 3')
  eandmuVars.Add("nLooseMuons",     "Sum(looseMuons)")
  eandmuVars.Add("nLooseElectrons", "Sum(looseElectrons)")
  eandmuVars.Add("nLooseLeptons",   "nLooseMuons + nLooseElectrons")
  eandmuVars.Add("isMuon", "nLooseMuons == 1 && HLT_Mu17 == 1 && isSM")
  eandmuVars.Add("isEl", "nLooseElectrons == 1 && HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30 == 1 && isSE")         
  a.Apply(eandmuVars)

  lepVeto = CutGroup("lepVeto")
  lepVeto.Add("single_loose_lepton", "nLooseLeptons == 1")
  lepVeto.Add("isVIPlepton", "isMuon || isEl")
  a.Apply(lepVeto)

  # Lepton kinematics
  # (safe to index [0], exactly one lepton survives)
  lepVars = VarGroup("lepVars")
  lepVars.Add("lep_pt",   "isMuon ? Muon_pt  [looseMuons][0] : Electron_pt  [looseElectrons][0]")
  lepVars.Add("lep_eta",  "isMuon ? Muon_eta [looseMuons][0] : Electron_eta [looseElectrons][0]")
  lepVars.Add("lep_abseta", "abs(lep_eta)")
  lepVars.Add("lep_phi",  "isMuon ? Muon_phi [looseMuons][0] : Electron_phi [looseElectrons][0]")
  lepVars.Add("lep_mass", "isMuon ? Muon_mass[looseMuons][0] : Electron_mass[looseElectrons][0]")
  a.Apply(lepVars)

  jVars = VarGroup('JetCleaningVars')

  jVars.Add("Jet_P4", "fVectorConstructor(Jet_pt,Jet_eta,Jet_phi,Jet_mass)")
  jVars.Add("Jet_EmEF","Jet_neEmEF + Jet_chEmEF")
  jVars.Add("DummyZero","float(0.0)")

  if year == '2024' or year == '2025':
    jVars.Add("Jet_jetId","jetidfunc(jetidAK4Tcorr,jetidAK4TLcorr,Jet_eta,Jet_chHEF,Jet_neHEF,Jet_chEmEF,Jet_neEmEF,Jet_muEF,Jet_chMultiplicity,Jet_neMultiplicity)"
            )
  jVars.Add("cleanedJets", "cleanJetsData(run,debug,year,ak4corr,ak4corrL1,ak8corr,Jet_P4,Jet_rawFactor,Jet_muonSubtrFactor,Jet_area,Jet_EmEF,Jet_jetId,Jet_P4,Jet_jetId,Rho_fixedGridRhoFastjetAll,DummyZero,DummyZero)") # muon and EM factors unused in this call, args 16-17 are dummies
  jVars.Add("cleanMets", "cleanJetsData(run,debug,year,ak4corr,ak4corrL1,ak8corr,Jet_P4,Jet_rawFactor,Jet_muonSubtrFactor,Jet_area,Jet_EmEF,Jet_jetId,Jet_P4,Jet_jetId,Rho_fixedGridRhoFastjetAll,RawPuppiMET_pt,RawPuppiMET_phi)") # lepton args unused in this call, args 16-17 are dummies
  jVars.Add("cleanJet_pt", "cleanedJets[0]")
  jVars.Add("cleanJet_eta", "cleanedJets[1]")
  jVars.Add("cleanJet_phi", "cleanedJets[2]")
  jVars.Add("cleanJet_mass", "cleanedJets[3]")

  # ------------------ Jet Calculations ------------------
  jVars.Add("DR_lepJets","DeltaRVecFloat(cleanJet_eta, cleanJet_phi, lep_eta, lep_phi)")
  jVars.Add("goodcleanJets", "cleanJet_pt > 30 && abs(cleanJet_eta) < 2.5 && Jet_jetId > 1 && DR_lepJets > 0.4 ")
  jVars.Add("NgoodcleanJets", "Sum(goodcleanJets)")

  jVars.Add("gcJet_pt_unsort", "cleanJet_pt[goodcleanJets == true]")
  jVars.Add("gcJet_ptargsort","ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(gcJet_pt_unsort))")

  jVars.Add("gcJet_pt","reorder(gcJet_pt_unsort,gcJet_ptargsort)")
  jVars.Add("gcJet_eta", "reorder(cleanJet_eta[goodcleanJets == true],gcJet_ptargsort)")
  jVars.Add("gcJet_phi", "reorder(cleanJet_phi[goodcleanJets == true],gcJet_ptargsort)")
  jVars.Add("gcJet_mass", "reorder(cleanJet_mass[goodcleanJets == true],gcJet_ptargsort)")

  jVars.Add("gcJet_vetomap", "jetvetofunc(jetvetocorr, gcJet_eta, gcJet_phi)")

  a.Apply(jVars)

  # MET and mT
  metVars = VarGroup("metVars")
  metVars.Add("corrMET_pt","cleanMets[4][0]")
  metVars.Add("corrMET_phi","cleanMets[4][1]")
  metVars.Add("mT", "calcMT(lep_pt, lep_phi, corrMET_pt, corrMET_phi)")
  a.Apply(metVars)

  metCuts = CutGroup("metCuts")
  metCuts.Add("Event has no vetoed jets","Sum(gcJet_vetomap) == 0")
  metCuts.Add("MET_cut", "corrMET_pt < 25.0")
  metCuts.Add("mT_cut",  "mT < 25.0")
  a.Apply(metCuts)

  # Jet mask + count
  # Lead-jet indexing is deferred until AFTER the
  # nGoodJets cut to avoid empty-vector segfaults
  jetMaskVars = VarGroup("jetMaskVars")
  jetMaskVars.Add("oppositeJet_mask", "goodcleanJets == 1 && DR_lepJets > 1.0")
  jetMaskVars.Add("nOppositeJets",    "Sum(oppositeJet_mask)")
  jetMaskVars.Add("lepJet_invM", "InvariantMassVecFloat(gcJet_pt, gcJet_eta, gcJet_phi, gcJet_mass, lep_pt, lep_eta, lep_phi, lep_mass)")
  jetMaskVars.Add("lepJet_onZ", "lepJet_invM > 81.1 && lepJet_invM < 101.1")
  jetMaskVars.Add("nBadJets","Sum(lepJet_onZ)")
  a.Apply(jetMaskVars)

  jetCuts = CutGroup("jetCuts")
  jetCuts.Add("1+ opposite-side jets", "nOppositeJets >= 1")
  jetCuts.Add("No Z-mass jet+lep", "nBadJets == 0")
  a.Apply(jetCuts)

  # Tight ID definitions
  # LooseMu_isTight / LooseEl_isTight are already
  # sub-vectors filtered to loose candidates only,
  # so [0] gives the single surviving lepton's result
  tightVars = VarGroup("tightVars")
  tightVars.Add("lepPassesTight","isMuon ? LooseMu_isTight[0] : LooseEl_isTight[0]")
  a.Apply(tightVars)

  pt_bins  = [10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0, 200.0, 9999.0]
  etabinsmu = [0.0, 0.9, 1.2, 2.1, 2.4]
  etabinsel = [-2.5, -2.0, -1.566, -1.444, -0.8, 0.0, 0.8, 1.444, 1.566, 2.0, 2.5]

  pt_binsarr = array.array('d', pt_bins)
  etabinsmuarr = array.array('d', etabinsmu)
  etabinselarr = array.array('d', etabinsel)

  checkpoint = a.GetActiveNode()

  finalFile = "DataFakeRate_"+sample + era + ver + "_" + year + "_" + str(testNum1) + ".root";
  outFile = TFile.Open(finalFile, "UPDATE")
    
  if isSM:
    a.SetActiveNode(checkpoint)
    a.Cut("isMuon_den", "(isMuon == true) && (HLT_Mu17 == 1)")
    hDen_mu = a.DataFrame.Histo2D(
      ("hDen_mu", "(QCD) Muon fake rate denominator;p_{T} [GeV];#eta",
       len(pt_bins)-1,  pt_binsarr,
       len(etabinsmu)-1, etabinsmuarr),
      "lep_pt", "lep_abseta"
    )

    a.Cut("isMuon_num", "lepPassesTight == true")
    hNum_mu = a.DataFrame.Histo2D(
      ("hNum_mu", "(QCD) Muon fake rate numerator;p_{T} [GeV];#eta",
       len(pt_bins)-1,  pt_binsarr,
       len(etabinsmu)-1, etabinsmuarr),
      "lep_pt", "lep_abseta"
    )

    print("Cut statistics QCD muons:")
    rep = a.DataFrame.Report()
    rep.Print()

    hDen_mu.GetValue()
    hNum_mu.GetValue()
    hDen_mu.Write()
    hNum_mu.Write()

  elif isSE:
    a.SetActiveNode(checkpoint)
    a.Cut("isEle_den", "isMuon == false && HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30 == 1")
    hDen_el = a.DataFrame.Histo2D(
      ("hDen_el", "(QCD) Electron fake rate denominator;p_{T} [GeV];#eta",
       len(pt_bins)-1,  pt_binsarr,
       len(etabinsel)-1, etabinselarr),
      "lep_pt", "lep_eta"
    )

    a.Cut("isEle_num", "lepPassesTight == true")
    hNum_el = a.DataFrame.Histo2D(
      ("hNum_el", "(QCD) Electron fake rate numerator;p_{T} [GeV];#eta",
       len(pt_bins)-1,  pt_binsarr,
       len(etabinsel)-1, etabinselarr),
      "lep_pt", "lep_eta"
    )

    print("Cut statistics QCD electrons:")
    rep = a.DataFrame.Report()
    rep.Print()

    hDen_el.GetValue()
    hNum_el.GetValue()
    hDen_el.Write()
    hNum_el.Write()

  outFile.Close()

  print("--------- Analysis End ---------")

  a.Close()
  
analyze("Nominal")
