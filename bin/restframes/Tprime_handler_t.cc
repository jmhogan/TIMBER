#include "include/RestFramesHandler.hh"

#include "TLorentzVector.h"
#include "TVector3.h"
#include "include/RestFrames.hh"
#include <ROOT/RVec.hxx>
#include <algorithm>
#include <iostream>
//#include <array> //TODO think about replacing std::array<float, 4> with RVec<float>?
#include <memory>
#include <mutex>
#include <random>

using namespace RestFrames;
using namespace ROOT::VecOps;

class Tprime_RestFrames_Handler_t : public RestFramesHandler {
    private:
	//(H/Z)t tree
        // Reconstruction frames
        std::unique_ptr<DecayRecoFrame> TTbar;
        std::unique_ptr<DecayRecoFrame> T;
        std::unique_ptr<VisibleRecoFrame> Tbar;

	std::unique_ptr<DecayRecoFrame> t;
	std::unique_ptr<VisibleRecoFrame> J0; // this is the Z/H boson
		
	std::unique_ptr<DecayRecoFrame> W;
        std::unique_ptr<VisibleRecoFrame> b;

        std::unique_ptr<VisibleRecoFrame> l;
        std::unique_ptr<InvisibleRecoFrame> nu;

        // Groups
        std::unique_ptr<CombinatoricGroup> JETS;

	std::unique_ptr<CombinatoricGroup> SMTOP;
       
       	std::unique_ptr<InvisibleGroup> INV;
        
        // Jigsaws
	std::unique_ptr<SetMassInvJigsaw> NuM;
	std::unique_ptr<SetRapidityInvJigsaw> NuR;
	//std::unique_ptr<MinMassDiffInvJigsaw> MinDeltaMt;
        std::unique_ptr<ContraBoostInvJigsaw> MinContraMt;
	
	//std::unique_ptr<MinMassesCombJigsaw> MinMJets;
	std::unique_ptr<MinMassDiffCombJigsaw> MinDiffJets;

	std::unique_ptr<MinMassChi2CombJigsaw> MinChi2;

	void define_tree() override;
	void define_groups_jigsaws() override;

    public:
        Tprime_RestFrames_Handler_t();
        RVec<double> calculate_doubles(TLorentzVector &lepton, TVector3 &met3, TLorentzVector &jet1, TLorentzVector &jet2, TLorentzVector &jet3, RVec<TLorentzVector> &AK4s);
	//std::tuple<float,float>
	
	RVec<TLorentzVector> calculate_vecs();
};

Tprime_RestFrames_Handler_t::Tprime_RestFrames_Handler_t() {
    initialize();
};

void Tprime_RestFrames_Handler_t::define_tree() {
    LAB.reset(new LabRecoFrame("LAB","LAB"));
    TTbar.reset(new DecayRecoFrame("TTbar", "T#bar{T}"));
    LAB->AddChildFrame(*TTbar);

    // Vector Like T quark particle production
    T.reset(new DecayRecoFrame("T", "T"));
    Tbar.reset(new VisibleRecoFrame("Tbar", "#bar{T}"));
    TTbar->AddChildFrame(*T);
    TTbar->AddChildFrame(*Tbar);
    
    // T -> t J0 (Z or H) 
    t.reset(new DecayRecoFrame("t","t"));
    J0.reset(new VisibleRecoFrame("J0", "J0_{AK8}"));
    T->AddChildFrame(*t);
    T->AddChildFrame(*J0);

/*    J1.reset(new VisibleRecoFrame("J1", "J1_{AK8}"));
    J2.reset(new VisibleRecoFrame("J2","J2_{AK8}"));
    Tbar->AddChildFrame(*J1);
    Tbar->AddChildFrame(*J2);	*/
    
    // t -> W b 
    W.reset(new DecayRecoFrame("W","W"));
    b.reset(new VisibleRecoFrame("b", "b"));
    t->AddChildFrame(*W);
    t->AddChildFrame(*b);

    // W -> l nu
    l.reset(new VisibleRecoFrame("l", "#it{l}"));
    nu.reset(new InvisibleRecoFrame("nu", "#nu"));
    W->AddChildFrame(*l);
    W->AddChildFrame(*nu);
}

void Tprime_RestFrames_Handler_t::define_groups_jigsaws() {
    // Combinatoric Group for jets
    JETS.reset(new CombinatoricGroup("JETS", "Jet Jigsaws"));
    JETS->AddFrame(*J0);
    JETS->AddFrame(*Tbar);
    //JETS->AddFrame(*J1);
    //JETS->AddFrame(*J2);

    // jet frames must have at least one element
    JETS->SetNElementsForFrame(*J0, 1);
    JETS->SetNElementsForFrame(*Tbar, 2);
    //JETS->SetNElementsForFrame(*J1, 1);
    //JETS->SetNElementsForFrame(*J2, 1);
    
    // Combinatoric Group for SM top quark reconstruction
    SMTOP.reset(new CombinatoricGroup("SMTOP", "Standard Model top Jigsaw"));
    SMTOP->AddFrame(*b);
    SMTOP->SetNElementsForFrame(*b, 1);  //TODO was 1!
   /*
    * So if I put it to 1 it errors cause some events don't have the good ak4 isolated b-tagged jets?
    * if I put it to 0 it never puts anything to the b
    *
    * Is there a way to switch trees for when there are no ak4 isolated b-tagged jets in this event?  revert to the 1st tree?
    * ah have 2 trees, but then switch which tree you put the stuff in (MET, jets, lepton, etc), and which tree you analyze() 
    *
    *
    */ 
    // Invisible Group for Neutrino
    INV.reset(new InvisibleGroup("INV", "MET Jigsaws"));
    INV->AddFrame(*nu);

    // -------------------- Define Jigsaws for reconstruction trees --------------
    std::string jigsaw_name;

    // 1 Minimize difference Mt jigsaws                not: Minimize equal (vector) top masses neutrino jigsaws
    jigsaw_name = "M_{#nu} = f(m_{b#it{l}J_{0}J_{1}} , m_{b#it{l}} , m_{J_{0}J_{1}})";
    NuM.reset(new SetMassInvJigsaw("NuM", jigsaw_name));
    INV->AddJigsaw(*NuM); 
    
    // 2 
    jigsaw_name = "#eta_{#nu} = #eta_{J_{0} b #it{l} Tbar}";
    NuR.reset(new SetRapidityInvJigsaw("NuR", jigsaw_name));
    INV->AddJigsaw(*NuR);
    NuR->AddVisibleFrame(*l); 
    NuR->AddVisibleFrame(*b); 
    NuR->AddVisibleFrame(*Tbar); 
    NuR->AddVisibleFrame(*J0);
    //+*b+*Tbar);	//TODO not sure about this line

    // 3

    jigsaw_name = "min M_{T}, M_{T} = M_{Tbar}";
    MinContraMt.reset(new ContraBoostInvJigsaw("MinContraMt", jigsaw_name));
    INV->AddJigsaw(*MinContraMt);
    MinContraMt->AddVisibleFrames(*J0+*l+*b, 0);
    MinContraMt->AddVisibleFrame(*Tbar, 1);
    MinContraMt->AddInvisibleFrame(*nu, 0);

    // MinMassDiffInv was ok, not best
    /*jigsaw_name = "min ( M_{T}- M_{Tbar} )^{2}";
    MinDeltaMt.reset(new MinMassDiffInvJigsaw("MinDeltaMt", jigsaw_name, 2));
    INV->AddJigsaw(*MinDeltaMt);
    MinDeltaMt->AddInvisibleFrame(*nu, 0);
    //MinDeltaMt.AddInvisibleFrame(Nb_R4, 1);
    MinDeltaMt->AddVisibleFrames(*l+*b, 0);
    MinDeltaMt->AddVisibleFrame(*Tbar, 1); //OR *J0+*J1, 1) ???
    MinDeltaMt->AddMassFrame(*T, 0);
    MinDeltaMt->AddMassFrame(*Tbar, 1);
    //MinDeltaMt.AddMassFrame(Lb_R4, 1); //???
*/
    // 4 Combinatoric Jigsaws 
    // MinMassesSqCombJigsaw worked but same problem as MinMassesCombJigsaw
    // MinMassDiffCombJigsaw Now Works!!!
    // MinMassChi2ComJigsaw works very well

    // MinMassesCombJigsaw, combinatoric jigsaws for everything else...
/*    jigsaw_name = "Minimize M(b #it{l} ) , M(Tbar)"; //M(J0 J1 )
    
    MinMJets.reset(new MinMassesCombJigsaw("MinCombJets", jigsaw_name));
    JETS->AddJigsaw(*MinMJets);
    MinMJets->AddFrames(*l+*b,0);
    MinMJets->AddFrame(*Tbar,1);
*/
    // MinMassDiffCombJigsaw
    jigsaw_name = "min ( M_{T}- M_{Tbar} )^{2}";
  
    MinDiffJets.reset(new MinMassDiffCombJigsaw("MinDiffJets", jigsaw_name, 2, 1)); // last param is the # of object frames that need to be calculated.  2 doesn't work
    JETS->AddJigsaw(*MinDiffJets);
    MinDiffJets->AddObjectFrames(*J0+*l+*b, 0);
    MinDiffJets->AddCombFrame(*J0, 0);
    MinDiffJets->AddObjectFrame(*Tbar, 1); 
    MinDiffJets->AddCombFrame(*Tbar, 1);
    
    // JET group 2 for SM t quark
    jigsaw_name = "Minimize Chi^2";
    MinChi2.reset(new MinMassChi2CombJigsaw("MinChi2", jigsaw_name, 1, 1)); // may need to change second number to 2
    SMTOP->AddJigsaw(*MinChi2);
    MinChi2->AddObjectFrame(*l, 0);
    MinChi2->AddObjectFrame(*b, 0);
    MinChi2->AddCombFrame(*b, 0);
    //MinChi2->AddObjectFrame(*t, 0);
    //MinChi2->AddCombFrame(*t, 0);
    MinChi2->SetMass(171.77, 0);
    MinChi2->SetSigma(0.38, 0);
};

RVec<double> Tprime_RestFrames_Handler_t::calculate_doubles(TLorentzVector &lepton, TVector3 &met3, TLorentzVector &jet1, TLorentzVector &jet2, TLorentzVector &jet3, RVec<TLorentzVector> &AK4s) {
    before_analysis();
    
    INV->SetLabFrameThreeVector(met3);	
    l->SetLabFrameFourVector(lepton);
    //TAUS->AddLabFrameFourVector(tau1);
    //TAUS->AddLabFrameFourVector(tau2);

    std::vector<RFKey> JETS_ID; // ID for tracking jets in tree
    JETS_ID.push_back(JETS->AddLabFrameFourVector(jet3));
    JETS_ID.push_back(JETS->AddLabFrameFourVector(jet1));
    JETS_ID.push_back(JETS->AddLabFrameFourVector(jet2));
    
    for (int i = 0; i < AK4s.size(); i++) {
      JETS_ID.push_back(SMTOP->AddLabFrameFourVector(AK4s[i]));
    } // if this loop is skipped that means there were no T -> t -> b W

    LAB->AnalyzeEvent(); // analyze the event


    RVec<double> observables; // = {T_mass, Tbar_mass};

    observables.push_back(TTbar->GetMass());
    observables.push_back(TTbar->GetCosDecayAngle());
    observables.push_back(TTbar->GetDeltaPhiDecayAngle());
    
    observables.push_back(T->GetMass());
    observables.push_back(T->GetCosDecayAngle());
    observables.push_back(T->GetDeltaPhiDecayAngle());
    
    observables.push_back(Tbar->GetMass());
    observables.push_back(Tbar->GetCosDecayAngle());
  //observables.push_back(Tbar->GetDeltaPhiDecayAngle());
    
    observables.push_back(W->GetMass());
    observables.push_back(W->GetCosDecayAngle());
    observables.push_back(W->GetDeltaPhiDecayAngle());
    
    observables.push_back(b->GetMass());
    observables.push_back(b->GetCosDecayAngle());
  //observables.push_back(b->GetDeltaPhiDecayAngle());		DON'T GET this, doesn't work, ERRORS
    /*
    observables.push_back(J0->GetMass());
    observables.push_back(J0->GetCosDecayAngle());
    observables.push_back(J0->GetDeltaPhiDecayAngle());

    observables.push_back(J1->GetMass());
    observables.push_back(J1->GetCosDecayAngle());
    observables.push_back(J1->GetDeltaPhiDecayAngle());
    */
//wasn't able to save any l nor nu things

    /*std::default_random_engine generator;
    std::bernoulli_distribution dist(0.5);
    bool which = dist(generator);
    
    if (which) return std::make_tuple(calc_mass_Tbar, calc_mass_T); */
    return observables; //std::make_tuple(calc_mass_T, calc_mass_Tbar);
};

// calculate_vecs() returns all the important four vectors of the frames in the tree
RVec<TLorentzVector> Tprime_RestFrames_Handler_t::calculate_vecs() {
    RVec<TLorentzVector> observables;

    observables.push_back(TTbar->GetFourVector());
    
    observables.push_back(T->GetFourVector());
    observables.push_back(T->GetFourVector(*TTbar));
    
    observables.push_back(Tbar->GetFourVector());
    observables.push_back(Tbar->GetFourVector(*TTbar));
    
    observables.push_back(W->GetFourVector());
    observables.push_back(W->GetFourVector(*TTbar));
    observables.push_back(W->GetFourVector(*T));
    
    observables.push_back(b->GetFourVector());
    observables.push_back(b->GetFourVector(*TTbar));
    observables.push_back(b->GetFourVector(*T));
    
    /*observables.push_back(J0->GetFourVector());
    observables.push_back(J0->GetFourVector(*TTbar));
    observables.push_back(J0->GetFourVector(*Tbar));
    observables.push_back(J1->GetFourVector());
    observables.push_back(J1->GetFourVector(*TTbar));
    observables.push_back(J1->GetFourVector(*Tbar)); */
   
    // no 4vecs for l and nu because it didn't work 

    after_analysis();
    
    return observables;
};



class Tprime_RestFrames_Container_t : public RestFramesContainer {
    public:
        Tprime_RestFrames_Container_t(int num_threads);
        RestFramesHandler *create_handler() override;

        RVec<double> return_doubles(int thread_index, float lepton_pt, float lepton_eta, float lepton_phi, float lepton_mass, RVec<float> fatjet_pt, RVec<float> fatjet_eta, RVec<float> fatjet_phi, RVec<float> fatjet_mass, RVec<float> fatjet_DeepFlav, float met_pt, float met_phi, RVec<TLorentzVector> jets);
	
	RVec<TLorentzVector> return_vecs(int thread_index);
};

Tprime_RestFrames_Container_t::Tprime_RestFrames_Container_t (int num_threads) : RestFramesContainer(num_threads){
    initialize();
};

RestFramesHandler * Tprime_RestFrames_Container_t::create_handler() {
    return new Tprime_RestFrames_Handler_t;
}

// return_doubles() returns all the masses, cos angles, and deltaPhi angles of the frames in the tree
RVec<double> Tprime_RestFrames_Container_t::return_doubles(int thread_index, float lepton_pt, float lepton_eta, float lepton_phi, float lepton_mass, RVec<float> fatjet_pt, RVec<float> fatjet_eta, RVec<float> fatjet_phi, RVec<float> fatjet_mass, RVec<float> fatjet_DeepFlav, float met_pt, float met_phi, RVec<TLorentzVector> jets) {

    // This pointer should explicitly not be deleted!
    Tprime_RestFrames_Handler_t *rfh = static_cast<Tprime_RestFrames_Handler_t *>(get_handler(thread_index));

    TLorentzVector fatjet_1;
    TLorentzVector fatjet_2;
    TLorentzVector fatjet_3;

    TLorentzVector lepton;

    TVector3 met3;
    
    lepton.SetPtEtaPhiM(lepton_pt, lepton_eta, lepton_phi, lepton_mass);  //TODO in the future we want to make sure the muon/lepton is good enough
    //lepton.SetPtEtaPhiM(lepton_pt[0], lepton_eta[0], lepton_phi[0], lepton_mass[0]);  //TODO in the future we want to make sure the muon/lepton is good enough

    fatjet_1.SetPtEtaPhiM(fatjet_pt[0], fatjet_eta[0], fatjet_phi[0], fatjet_mass[0]);
    fatjet_2.SetPtEtaPhiM(fatjet_pt[1], fatjet_eta[1], fatjet_phi[1], fatjet_mass[1]);
    fatjet_3.SetPtEtaPhiM(fatjet_pt[2], fatjet_eta[2], fatjet_phi[2], fatjet_mass[2]);

    double MET_px  = met_pt*std::cos(met_phi);
    double MET_py  = met_pt*std::sin(met_phi);
    met3  = TVector3(MET_px, MET_py, 0.0);
    //std::tuple<float, float> masses = rfh->calculate_doubles(lepton, met3, jet_1, jet_2, jet_3); //, jet_4);
   
    RVec<double> observables = rfh->calculate_doubles(lepton, met3, fatjet_1, fatjet_2, fatjet_3, jets); //jet_4);

    return observables;
}


// return_vecs() returns all the four vectors/TLorentzVectors of the frames in the tree
RVec<TLorentzVector> Tprime_RestFrames_Container_t::return_vecs(int thread_index) {
    // This pointer should explicitly not be deleted!
    Tprime_RestFrames_Handler_t *rfh = static_cast<Tprime_RestFrames_Handler_t *>(get_handler(thread_index));

    return rfh->calculate_vecs();
}
