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

class Tprime_RestFrames_Handler_W : public RestFramesHandler {
    private:

        // Reconstruction frames
        std::unique_ptr<DecayRecoFrame> TTbar;
        std::unique_ptr<DecayRecoFrame> T;
        std::unique_ptr<VisibleRecoFrame> Tbar;
        
	std::unique_ptr<DecayRecoFrame> W;
        std::unique_ptr<VisibleRecoFrame> b;
        //std::unique_ptr<VisibleRecoFrame> J0;
        //std::unique_ptr<VisibleRecoFrame> J1;

        std::unique_ptr<VisibleRecoFrame> l;
        std::unique_ptr<InvisibleRecoFrame> nu;

        // Groups
        std::unique_ptr<CombinatoricGroup> JETS;

        std::unique_ptr<InvisibleGroup> INV;
        
        // Jigsaws
	std::unique_ptr<SetMassInvJigsaw> NuM;
	std::unique_ptr<SetRapidityInvJigsaw> NuR;
	//std::unique_ptr<MinMassDiffInvJigsaw> MinDeltaMt;
        std::unique_ptr<ContraBoostInvJigsaw> MinContraMt;
	
	//std::unique_ptr<MinMassChi2CombJigsaw> MinChi2;
	//std::unique_ptr<MinMassesCombJigsaw> MinMJets;
	std::unique_ptr<MinMassDiffCombJigsaw> MinDiffJets;

	void define_tree() override;
	void define_groups_jigsaws() override;

    public:
        Tprime_RestFrames_Handler_W();
        RVec<double> calculate_doubles(TLorentzVector &lepton, TVector3 &met3, TLorentzVector &jet1, TLorentzVector &jet2, TLorentzVector &jet3); //, TLorentzVector &jet4);
	//std::tuple<float,float>
	
	RVec<TLorentzVector> calculate_vecs();
};

Tprime_RestFrames_Handler_W::Tprime_RestFrames_Handler_W() {
    initialize();
};

void Tprime_RestFrames_Handler_W::define_tree() {
    LAB.reset(new LabRecoFrame("LAB","LAB"));
    TTbar.reset(new DecayRecoFrame("TTbar", "T#bar{T}"));
    LAB->AddChildFrame(*TTbar);

    // Vector Like T quark particle production
    T.reset(new DecayRecoFrame("T", "T"));
    Tbar.reset(new VisibleRecoFrame("Tbar", "#bar{T}"));
    TTbar->AddChildFrame(*T);
    TTbar->AddChildFrame(*Tbar);
    // T -> W b
    W.reset(new DecayRecoFrame("W","W"));
    b.reset(new VisibleRecoFrame("b", "b"));
    T->AddChildFrame(*W);
    T->AddChildFrame(*b);

/*    J1.reset(new VisibleRecoFrame("J1", "J1_{AK8}"));
    J0.reset(new VisibleRecoFrame("J0","J0_{AK8}"));
    Tbar->AddChildFrame(*J1);
    Tbar->AddChildFrame(*J0);	*/
    
    // W -> l nu
    l.reset(new VisibleRecoFrame("l", "#it{l}"));
    nu.reset(new InvisibleRecoFrame("nu", "#nu"));
    W->AddChildFrame(*l);
    W->AddChildFrame(*nu);
}

void Tprime_RestFrames_Handler_W::define_groups_jigsaws() {
    // Combinatoric Group for jets
    JETS.reset(new CombinatoricGroup("JETS", "Jet Jigsaws"));
    JETS->AddFrame(*b);
    JETS->AddFrame(*Tbar);
    //JETS->AddFrame(*J1);
    //JETS->AddFrame(*J0);

    // jet frames must have at least one element
    JETS->SetNElementsForFrame(*b, 1);
    JETS->SetNElementsForFrame(*Tbar, 2);
    //JETS->SetNElementsForFrame(*J1, 1);
    //JETS->SetNElementsForFrame(*J0, 1);
    
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
    jigsaw_name = "#eta_{#nu} = #eta_{b #it{l} Tbar}";
    NuR.reset(new SetRapidityInvJigsaw("NuR", jigsaw_name));
    INV->AddJigsaw(*NuR);
    NuR->AddVisibleFrame(*l); 
    NuR->AddVisibleFrame(*b); 
    NuR->AddVisibleFrame(*Tbar); 
    //+*b+*Tbar);	//TODO not sure about this line

    // 3

    jigsaw_name = "min M_{T}, M_{T} = M_{Tbar}";
    MinContraMt.reset(new ContraBoostInvJigsaw("MinContraMt", jigsaw_name));
    INV->AddJigsaw(*MinContraMt);
    MinContraMt->AddVisibleFrames(*l+*b, 0);
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
    // MinMassDiffCombJigsaw Initialized Analysis but fell into some infinite loop
    // MinMassChi2ComJigsaw works very well
    /*jigsaw_name = "Minimize Chi^2";
    MinChi2.reset(new MinMassChi2CombJigsaw("MinChi2", jigsaw_name, 2, 2));
    JETS->AddJigsaw(*MinChi2);
    MinChi2->AddObjectFrame(*l, 0);
    MinChi2->AddObjectFrame(*b, 0);
    MinChi2->AddCombFrame(*b, 0);
    MinChi2->AddObjectFrame(*Tbar, 1);
    MinChi2->AddCombFrame(*Tbar, 1);
    MinChi2->SetMass(1435, 0);
    MinChi2->SetSigma(205.2, 0);
    MinChi2->SetMass(1456, 1);
    MinChi2->SetSigma(173.2, 1); */ 

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
    MinDiffJets->AddObjectFrames(*l+*b, 0);
    MinDiffJets->AddCombFrame(*b, 0);
    MinDiffJets->AddObjectFrame(*Tbar, 1); 
    MinDiffJets->AddCombFrame(*Tbar, 1);
    /*
    MinDiffJets->AddVisibleFrames(*l+*b, 0);
    MinDiffJets->AddVisibleFrame(*Tbar, 1); //OR *J0+*J1, 1) ???
    MinDiffJets->AddMassFrame(*T, 0);
    MinDiffJets->AddMassFrame(*Tbar, 1);
    */
};

RVec<double> Tprime_RestFrames_Handler_W::calculate_doubles(TLorentzVector &lepton, TVector3 &met3, TLorentzVector &jet1, TLorentzVector &jet2, TLorentzVector &jet3) { //, TLorentzVector &jet4) {
    before_analysis();
    
    INV->SetLabFrameThreeVector(met3);	
    l->SetLabFrameFourVector(lepton);

    std::vector<RFKey> JETS_ID; // ID for tracking jets in tree
    JETS_ID.push_back(JETS->AddLabFrameFourVector(jet3));
    JETS_ID.push_back(JETS->AddLabFrameFourVector(jet1));
    JETS_ID.push_back(JETS->AddLabFrameFourVector(jet2));
    //JETS_ID.push_back(JETS->AddLabFrameFourVector(jet4));
    //b->SetLabFrameFourVector(jet4));

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

    observables.push_back(TTbar->GetDeltaPhiVisible());
    observables.push_back(TTbar->GetDeltaPhiDecayVisible());
    observables.push_back(TTbar->GetDeltaPhiBoostVisible());
    observables.push_back(TTbar->GetVisibleShape());
    
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

    //TODO moved? after_analysis();

    /*std::default_random_engine generator;
    std::bernoulli_distribution dist(0.5);
    bool which = dist(generator);
    
    if (which) return std::make_tuple(calc_mass_Tbar, calc_mass_T); */
    return observables; //std::make_tuple(calc_mass_T, calc_mass_Tbar);
};

// calculate_vecs() returns all the important four vectors of the frames in the tree
RVec<TLorentzVector> Tprime_RestFrames_Handler_W::calculate_vecs() {
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



class Tprime_RestFrames_Container_W : public RestFramesContainer {
    public:
        Tprime_RestFrames_Container_W(int num_threads);
        RestFramesHandler *create_handler() override;

        RVec<double> return_doubles(int thread_index, float lepton_pt, float lepton_eta, float lepton_phi, float lepton_mass, RVec<float> fatjet_pt, RVec<float> fatjet_eta, RVec<float> fatjet_phi, RVec<float> fatjet_mass, RVec<float> fatjet_DeepFlav, float met_pt, float met_phi);
	
	RVec<TLorentzVector> return_vecs(int thread_index);
};

Tprime_RestFrames_Container_W::Tprime_RestFrames_Container_W (int num_threads) : RestFramesContainer(num_threads){
    initialize();
};

RestFramesHandler * Tprime_RestFrames_Container_W::create_handler() {
    return new Tprime_RestFrames_Handler_W;
}

// return_doubles() returns all the masses, cos angles, and deltaPhi angles of the frames in the tree
RVec<double> Tprime_RestFrames_Container_W::return_doubles(int thread_index, float lepton_pt, float lepton_eta, float lepton_phi, float lepton_mass, RVec<float> fatjet_pt, RVec<float> fatjet_eta, RVec<float> fatjet_phi, RVec<float> fatjet_mass, RVec<float> fatjet_DeepFlav, float met_pt, float met_phi) {

    // This pointer should explicitly not be deleted!
    Tprime_RestFrames_Handler_W *rfh = static_cast<Tprime_RestFrames_Handler_W *>(get_handler(thread_index));

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
    //std::tuple<float, float> masses = rfh->calculate_doubles(lepton, met3, jet_1, jet_2, jet_3);
    RVec<double> observables = rfh->calculate_doubles(lepton, met3, fatjet_1, fatjet_2, fatjet_3); 

    //std::cout << "in return_doubles, about to return" << std::endl;
    return observables;
}


// return_vecs() returns all the four vectors/TLorentzVectors of the frames in the tree
RVec<TLorentzVector> Tprime_RestFrames_Container_W::return_vecs(int thread_index) {
  // This pointer should explicitly not be deleted!
  //std::cout << "in return_vecs, about to get the handler" << std::endl;
  Tprime_RestFrames_Handler_W *rfh = static_cast<Tprime_RestFrames_Handler_W *>(get_handler(thread_index));
  //std::cout << "in return_vecs, about to call calculates_vecs" << std::endl;
  return rfh->calculate_vecs();
  //std::cout << "done with calculate vecs" << std::endl;

}
