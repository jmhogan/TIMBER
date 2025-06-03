//Manual Reco function for BBto2b4tau

#include "TLorentzVector.h"
#include <ROOT/RVec.hxx>
#include <iostream>

using namespace std;

RVec<double> funcmassdiff(RVec<double> goodtau_pt, RVec<double> goodtau_eta, RVec<double> goodtau_phi, RVec<double> goodtau_mass, RVec<double> goodbjet_pt, RVec<double> goodbjet_eta, RVec<double> goodbjet_phi, RVec<double> goodbjet_mass) {

TLorentzVector B1, B2;
TLorentzVector B1final, B2final;
TLorentzVector b1, b2;

b1.SetPtEtaPhiM(goodbjet_pt[0], goodbjet_eta[0], goodbjet_phi[0], goodbjet_mass[0]);
b2.SetPtEtaPhiM(goodbjet_pt[1], goodbjet_eta[1], goodbjet_phi[1], goodbjet_mass[1]);

//cout << "goodbjet stuff: " << goodbjet_pt[0] << "," << goodbjet_eta[0] << "," << goodbjet_phi[0] << "," << goodbjet_mass[0] << endl;
//cout << "b1 stuff: " << b1.Pt() << "," << b1.Eta() << "," << b1.Phi() << "," << b1.M() << endl;

double massdiff = 9999999999999.0;
double massdifftemp;

TLorentzVector tauA, tauB, tauC, tauD;

tauA.SetPtEtaPhiM(goodtau_pt[0], goodtau_eta[0], goodtau_phi[0], goodtau_mass[0]);
tauB.SetPtEtaPhiM(goodtau_pt[1], goodtau_eta[1], goodtau_phi[1], goodtau_mass[1]);
tauC.SetPtEtaPhiM(goodtau_pt[2], goodtau_eta[2], goodtau_phi[2], goodtau_mass[2]);
tauD.SetPtEtaPhiM(goodtau_pt[3], goodtau_eta[3], goodtau_phi[3], goodtau_mass[3]);

//'looping' through the options
//First
B1 = b1+tauA+tauB;
B2 = b2+tauC+tauD;
massdifftemp = abs(B1.M() - B2.M());

if(massdifftemp < massdiff) {
	massdiff = massdifftemp;
	B1final = B1;
	B2final = B2;
}

//Second
B1 = b1+tauA+tauC;
B2 = b2+tauB+tauD;
massdifftemp = abs(B1.M() - B2.M());

if(massdifftemp < massdiff) {
	massdiff = massdifftemp;
	B1final = B1;
	B2final = B2;
}

//Third
B1 = b1+tauA+tauD;
B2 = b2+tauB+tauC;
massdifftemp = abs(B1.M() - B2.M());

if(massdifftemp < massdiff) {
	massdiff = massdifftemp;
	B1final = B1;
	B2final = B2;
}

//Fourth
B1 = b1+tauB+tauC;
B2 = b2+tauA+tauD;
massdifftemp = abs(B1.M() - B2.M());

if(massdifftemp < massdiff) {
	massdiff = massdifftemp;
	B1final = B1;
	B2final = B2;
}

//Fifth
B1 = b1+tauB+tauD;
B2 = b2+tauA+tauC;
massdifftemp = abs(B1.M() - B2.M());

if(massdifftemp < massdiff) {
	massdiff = massdifftemp;
	B1final = B1;
	B2final = B2;
}

//Sixth
B1 = b1+tauC+tauD;
B2 = b2+tauA+tauB;
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
