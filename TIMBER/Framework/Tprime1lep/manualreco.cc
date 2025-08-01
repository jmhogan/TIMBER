//Manual Reco function for BBto2b4tau

#include "TLorentzVector.h"
#include <ROOT/RVec.hxx>
#include <iostream>

using namespace std;

RVec<double> funcmassdiff(RVec<double> goodLepton_pt, RVec<double> goodLepton_eta, RVec<double> goodLepton_phi, RVec<double> goodLepton_mass, RVec<double> goodbjet_pt, RVec<double> goodbjet_eta, RVec<double> goodbjet_phi, RVec<double> goodbjet_mass) {

	TLorentzVector B1, B2;
	TLorentzVector B1final, B2final;
	TLorentzVector b1, b2;

	b1.SetPtEtaPhiM(goodbjet_pt[0], goodbjet_eta[0], goodbjet_phi[0], goodbjet_mass[0]);
	b2.SetPtEtaPhiM(goodbjet_pt[1], goodbjet_eta[1], goodbjet_phi[1], goodbjet_mass[1]);

	//cout << "goodbjet stuff: " << goodbjet_pt[0] << "," << goodbjet_eta[0] << "," << goodbjet_phi[0] << "," << goodbjet_mass[0] << endl;
	//cout << "b1 stuff: " << b1.Pt() << "," << b1.Eta() << "," << b1.Phi() << "," << b1.M() << endl;

	double massdiff = 9999999999999.0;
	double massdifftemp;

	TLorentzVector LeptonA, LeptonB, LeptonC, LeptonD;

	LeptonA.SetPtEtaPhiM(goodLepton_pt[0], goodLepton_eta[0], goodLepton_phi[0], goodLepton_mass[0]);
	LeptonB.SetPtEtaPhiM(goodLepton_pt[1], goodLepton_eta[1], goodLepton_phi[1], goodLepton_mass[1]);
	LeptonC.SetPtEtaPhiM(goodLepton_pt[2], goodLepton_eta[2], goodLepton_phi[2], goodLepton_mass[2]);
	
	if (goodLepton_pt.size() == 4) {
		LeptonD.SetPtEtaPhiM(goodLepton_pt[3], goodLepton_eta[3], goodLepton_phi[3], goodLepton_mass[3]);
	} else if (goodLepton_pt.size() == 3) {
		LeptonD.SetPtEtaPhiM(0, 0, 0, 0);
	}	

	//'looping' through the options
	//First
	B1 = b1+LeptonA+LeptonB;
	B2 = b2+LeptonC+LeptonD;
	massdifftemp = abs(B1.M() - B2.M());

	if(massdifftemp < massdiff) {
		massdiff = massdifftemp;
		B1final = B1;
		B2final = B2;
	}

	//Second
	B1 = b1+LeptonA+LeptonC;
	B2 = b2+LeptonB+LeptonD;
	massdifftemp = abs(B1.M() - B2.M());

	if(massdifftemp < massdiff) {
		massdiff = massdifftemp;
		B1final = B1;
		B2final = B2;
	}

	//Third
	B1 = b1+LeptonA+LeptonD;
	B2 = b2+LeptonB+LeptonC;
	massdifftemp = abs(B1.M() - B2.M());

	if(massdifftemp < massdiff) {
		massdiff = massdifftemp;
		B1final = B1;
		B2final = B2;
	}

	//Fourth
	B1 = b1+LeptonB+LeptonC;
	B2 = b2+LeptonA+LeptonD;
	massdifftemp = abs(B1.M() - B2.M());

	if(massdifftemp < massdiff) {
		massdiff = massdifftemp;
		B1final = B1;
		B2final = B2;
	}

	//Fifth
	B1 = b1+LeptonB+LeptonD;
	B2 = b2+LeptonA+LeptonC;
	massdifftemp = abs(B1.M() - B2.M());

	if(massdifftemp < massdiff) {
		massdiff = massdifftemp;
		B1final = B1;
		B2final = B2;
	}

	//Sixth
	B1 = b1+LeptonC+LeptonD;
	B2 = b2+LeptonA+LeptonB;
	massdifftemp = abs(B1.M() - B2.M());

	if(massdifftemp < massdiff) {
		massdiff = massdifftemp;
		B1final = B1;
		B2final = B2;
	}

	//return the remaining B1/B2final

	RVec<double> returnVec = {B1final.Px(), B1final.Py(), B1final.Pz(), B1final.M(), B2final.Px(), B2final.Py(), B2final.Pz(), B2final.M()};


	return returnVec;

}
