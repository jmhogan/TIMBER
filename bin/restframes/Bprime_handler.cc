#inelude "include/RestFramesHandler.hh"

#include "TLorentzVector.h"
#include "TVector3.h"
#include "include/RestFrames.hh"
#include <ROOT/RVec.hxx>
#include <algorithm>
#include <iostream>
#include <memory>
#include <mutex>
#include <random>

using namespace RestFrames;
using namespace ROOT::VecOps;

class gc_Handler : public RestFramesHandler {
    private:

        // Reconstruction frames
        std::unique_ptr<DecayRecoFrame> BB;

        std::unique_ptr<DecayRecoFrame> B1;
	std::unique_ptr<VisibleRecoFrame> tau11;
        std::unique_ptr<VisibleRecoFrame> tau12;
	std::unique_ptr<VisibleRecoFrame> b1;
        
        std::unique_ptr<DecayRecoFrame> B2;
	std::unique_ptr<VisibleRecoFrame> tau21;
	std::unique_ptr<VisibleRecoFrame> tau22;
	std::unique_ptr<VisibleRecoFrame> b2;

        // Groups
        std::unique_ptr<CombinatoricGroup> taus;
	std::unique_ptr<CombinatoricGroup> bs;
        
        // Jigsaws
	
	//std::unique_ptr<SetMassInvJigsaw> NuM;
	//std::unique_ptr<SetRapidityInvJigsaw> NuR;
	//std::unique_ptr<MinMassDiffInvJigsaw> MinDeltaMt;
        //std::unique_ptr<ContraBoostInvJigsaw> MinContraMt;
	//std::unique_ptr<MinMassChi2CombJigsaw> MinChi2;
	//std::unique_ptr<MinMassesCombJigsaw> MinMJets;
	
	std::unique_ptr<MinMassDiffCombJigsaw> MiMaDifTau;
	std::unique_ptr<MinMassDiffCombJigsaw> MiMaDifb;

	void define_tree() override;
	void define_groups_jigsaws() override;

    public:
        gc_Handler();
        RVec<double> calculate_doubles(TLorentzVector &recotau11, TLorentzVector &recotau12, TLorentzVector &recob1, TLorentzVector &recotau21, TLorentzVector &recotau22, TLorentzVector &recob2);
	RVec<TLorentzVector> calculate_vecs();
};

gc_Handler::gc_Handler() {
    initialize();
}

void gc_Handler::define_tree() {
    LAB.reset(new LabRecoFrame("LAB","LAB"));
    BB.reset(new DecayRecoFrame("BB","BBar"));
    LAB->AddChildFrame(*BB);
    
    B1.reset(new DecayRecoFrame("B1", "B1")); //full names?? (# is like \ latex)
    BB->AddChildFrame(*B1);
    tau11.reset(new VisibleRecoFrame("tau11", "tau11")); 
    tau12.reset(new VisibleRecoFrame("tau12", "tau12"));
    b1.reset(new VisibleRecoFrame("b1", "b1"));
    B1->AddChildFrame(*tau11);
    B1->AddChildFrame(*tau12);
    B1->AddChildFrame(*b1);
    
    B2.reset(new DecayRecoFrame("B2", "B2"));
    BB->AddChildFrame(*B2);
    tau21.reset(new VisibleRecoFrame("tau21", "tau21"));
    tau22.reset(new VisibleRecoFrame("tau22", "tau22")); 
    b2.reset(new VisibleRecoFrame("b2", "b2"));
    B2->AddChildFrame(*tau21);
    B2->AddChildFrame(*tau22);
    B2->AddChildFrame(*b2);
}  

void gc_Handler::define_groups_jigsaws() {
    taus.reset(new CombinatoricGroup("taus", "tau jigsaw"));
    taus->AddFrame(*tau11); 
    taus->AddFrame(*tau12);
    taus->AddFrame(*tau21);
    taus->AddFrame(*tau22);

    taus->SetNElementsForFrame(*tau11, 1);
    taus->SetNElementsForFrame(*tau12, 1);
    taus->SetNElementsForFrame(*tau21, 1);
    taus->SetNElementsForFrame(*tau22, 1);

   // bs.reset(new CombinatoricGroup("bs", "b jet jigsaw"));
   // bs->AddFrame(*b1);
    //bs->AddFrame(*b2);

    //bs->SetNElementsForFrame(*b1, 1);
    //bs->SetNElementsForFrame(*b2, 1);
    
    // -------------------- Define Jigsaws for reconstruction trees --------------
    // MinMassDiffCombJigsaw
    MiMaDifTau.reset(new MinMassDiffCombJigsaw("MiMaDifTau", "MinMassDiffTaus Jigsaw", 2, 1));
    
    taus->AddJigsaw(*MiMaDifTau);
    MiMaDifTau->AddObjectFrames(*tau11+*tau12+*b1, 0);
    MiMaDifTau->AddObjectFrames(*tau21+*tau22+*b2, 1);

    MiMaDifTau->AddCombFrame(*tau11, 0);
    MiMaDifTau->AddCombFrame(*tau12, 0);
    MiMaDifTau->AddCombFrame(*tau21, 1);
    MiMaDifTau->AddCombFrame(*tau22, 1);

    /*MiMaDifb.reset(new MinMassDiffCombJigsaw("MiMaDifb", "MinMassDiffbs Jigsaw", 2, 2));
    bs->AddJigsaw(*MiMaDifb);
    MiMaDifb->AddObjectFrame(*b1, 0);
    MiMaDifb->AddObjectFrame(*b2, 1);

    MiMaDifb->AddCombFrame(*b1, 0);
    MiMaDifb->AddCombFrame(*b2, 1);
    */
    std::cout<<"print1"<<std::endl; 
}

RVec<double> gc_Handler::calculate_doubles(TLorentzVector &recotau11, TLorentzVector &recotau12, TLorentzVector &recob1, TLorentzVector &recotau21, TLorentzVector &recotau22, TLorentzVector &recob2) {
    
    std::cout<<"print2"<<std::endl; 
    before_analysis();
    std::cout<<"print3"<<std::endl; 
    tau11->SetLabFrameFourVector(recotau11);
    std::vector<RFKey> taus_ID; // ID for tracking jets in tree
    taus_ID.push_back(taus->AddLabFrameFourVector(recotau11));
    taus_ID.push_back(taus->AddLabFrameFourVector(recotau12));
    taus_ID.push_back(taus->AddLabFrameFourVector(recotau21));
    taus_ID.push_back(taus->AddLabFrameFourVector(recotau22));


    std::vector<RFKey> bs_ID; // ID for tracking jets in tree
    bs_ID.push_back(bs->AddLabFrameFourVector(recob1));
    bs_ID.push_back(bs->AddLabFrameFourVector(recob2));

    LAB->AnalyzeEvent(); // analyze the event


    RVec<double> observables;

    observables.push_back(B1->GetMass());
    observables.push_back(B1->GetCosDecayAngle());
    observables.push_back(B1->GetDeltaPhiDecayAngle());
    
    observables.push_back(B2->GetMass());
    observables.push_back(B2->GetCosDecayAngle());
    observables.push_back(B2->GetDeltaPhiDecayAngle());
    
    observables.push_back(tau11->GetMass());
    observables.push_back(tau11->GetCosDecayAngle());
    observables.push_back(tau11->GetDeltaPhiDecayAngle());

    observables.push_back(tau12->GetMass());
    observables.push_back(tau12->GetCosDecayAngle());
    observables.push_back(tau12->GetDeltaPhiDecayAngle());
    
    observables.push_back(tau21->GetMass());
    observables.push_back(tau21->GetCosDecayAngle());
    observables.push_back(tau21->GetDeltaPhiDecayAngle());
    
    observables.push_back(tau22->GetMass());
    observables.push_back(tau22->GetCosDecayAngle());
    observables.push_back(tau22->GetDeltaPhiDecayAngle());
    
    observables.push_back(b1->GetMass());
    observables.push_back(b1->GetCosDecayAngle());
    observables.push_back(b1->GetDeltaPhiDecayAngle());

    observables.push_back(b2->GetMass());
    observables.push_back(b2->GetCosDecayAngle());
    observables.push_back(b2->GetDeltaPhiDecayAngle());

    return observables;
}

//returns all the important four vectors of the frames in the tree
RVec<TLorentzVector> gc_Handler::calculate_vecs() {
    
    RVec<TLorentzVector> observables;


    observables.push_back(B1->GetFourVector());
    observables.push_back(B2->GetFourVector());
    
    observables.push_back(tau11->GetFourVector());
    observables.push_back(tau11->GetFourVector(*B1));
    
    observables.push_back(tau12->GetFourVector());
    observables.push_back(tau12->GetFourVector(*B1));

    observables.push_back(b1->GetFourVector());
    observables.push_back(b1->GetFourVector(*B1));
    
    observables.push_back(tau21->GetFourVector());
    observables.push_back(tau21->GetFourVector(*B2));
    
    observables.push_back(tau22->GetFourVector());
    observables.push_back(tau22->GetFourVector(*B2));

    observables.push_back(b2->GetFourVector());
    observables.push_back(b2->GetFourVector(*B2));
   
   
    after_analysis();
    
    return observables;
};



class gc_Container : public RestFramesContainer {
    public:
        gc_Container(int num_threads);
        RestFramesHandler *create_handler() override;

        RVec<double> return_doubles(int thread_index, RVec<float> GoodTau_pt, RVec<float> GoodTau_eta, RVec<float> GoodTau_phi, RVec<float> GoodTau_mass, RVec<float> gcBJet_pt, RVec<float> gcBJet_eta, RVec<float> gcBJet_phi, RVec<float> gcBJet_mass);
	
	RVec<TLorentzVector> return_vecs(int thread_index);
};

gc_Container::gc_Container (int num_threads) : RestFramesContainer(num_threads){
    initialize();
};

RestFramesHandler * gc_Container::create_handler() {
    return new gc_Handler;
}

// return_doubles() returns all the masses, cos angles, and deltaPhi angles of the frames in the tree
RVec<double> gc_Container::return_doubles(int thread_index, RVec<float> GoodTau_pt, RVec<float> GoodTau_eta, RVec<float> GoodTau_phi, RVec<float> GoodTau_mass, RVec<float> gcBJet_pt, RVec<float> gcBJet_eta, RVec<float> gcBJet_phi, RVec<float> gcBJet_mass) {

    // This pointer should explicitly not be deleted!
    gc_Handler *rfh = static_cast<gc_Handler *>(get_handler(thread_index));
    
    TLorentzVector tau_1;
    TLorentzVector tau_2;
    TLorentzVector tau_3;
    TLorentzVector tau_4;
    TLorentzVector jet_1;
    TLorentzVector jet_2;

    tau_1.SetPtEtaPhiM(GoodTau_pt[0], GoodTau_eta[0], GoodTau_phi[0], GoodTau_mass[0]);
    tau_2.SetPtEtaPhiM(GoodTau_pt[1], GoodTau_eta[1], GoodTau_phi[1], GoodTau_mass[1]);
    tau_3.SetPtEtaPhiM(GoodTau_pt[2], GoodTau_eta[2], GoodTau_phi[2], GoodTau_mass[2]);
    tau_4.SetPtEtaPhiM(GoodTau_pt[3], GoodTau_eta[3], GoodTau_phi[3], GoodTau_mass[3]);
    jet_1.SetPtEtaPhiM(gcBJet_pt[0], gcBJet_eta[0], gcBJet_phi[0], gcBJet_mass[0]);
    jet_2.SetPtEtaPhiM(gcBJet_pt[1], gcBJet_eta[1], gcBJet_phi[1], gcBJet_mass[1]);
    
    RVec<double> observables = rfh->calculate_doubles(tau_1, tau_2, tau_3, tau_4, jet_1, jet_2); 

    //std::cout << "in return_doubles, about to return" << std::endl;
    return observables;
}


// return_vecs() returns all the four vectors/TLorentzVectors of the frames in the tree
RVec<TLorentzVector> gc_Container::return_vecs(int thread_index) {
  // This pointer should explicitly not be deleted!
  //std::cout << "in return_vecs, about to get the handler" << std::endl;
  gc_Handler *rfh = static_cast<gc_Handler *>(get_handler(thread_index));
  //std::cout << "in return_vecs, about to call calculates_vecs" << std::endl;
  return rfh->calculate_vecs();
  //std::cout << "done with calculate vecs" << std::endl;

}
