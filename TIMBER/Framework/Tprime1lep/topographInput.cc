// Methods in this file:
// decayType, recoGenMatch

using namespace std;
using namespace ROOT::VecOps;

#include <cmath>
#include <iostream>
#include <string>
#include "TLorentzVector.h"
#include <ROOT/RVec.hxx>

// -----------------------------------------------------------------------
// 		fxn for filtering events with bosonish mass
// -----------------------------------------------------------------------
bool hasBosonishMfunc(RVec<float> G4L_pt, RVec<float> G4L_eta, RVec<float> G4L_phi, RVec<float> G4L_mass, RVec<int> G4L_ID, RVec<int> G4L_charge) {

	for(int i = 0; i < G4L_ID.size(); i++) {
		for(int j = i+1; j < G4L_ID.size(); j++) {
			if (G4L_ID[i] != G4L_ID[j]) continue;
			if (G4L_charge[i] == G4L_charge[j]) continue;

			TLorentzVector Lep1, Lep2;
			Lep1.SetPtEtaPhiM(G4L_pt[i], G4L_eta[i], G4L_phi[i], G4L_mass[i]);
			Lep2.SetPtEtaPhiM(G4L_pt[j], G4L_eta[j], G4L_phi[j], G4L_mass[j]);
		
			if ((Lep1 + Lep2).M() < 12) return 1;
			else if ((Lep1 + Lep2).M() > 85 && (Lep1 + Lep2).M() < 97) return 1;
			//else if (G4L_ID[i] == 15 && (Lep1 + Lep2).M() > 115 && (Lep1 + Lep2).M() < 135) return 1; //taking out

		}
	}
	return 0;	
}

// ------------------------------------------------------------------------
//               convert matchibility from RVec to bitstring
// ------------------------------------------------------------------------

int convertMatchToInt(RVec<int> match) {
  int decimalValue = 0;
  int arraySize = 8;

  for (int i = 2; i < arraySize; ++i) {
      if (match[i] == 1) {
          decimalValue += static_cast<int>(std::pow(2, i-2));
      }
  }
  return decimalValue;
}


// ------------------------------------------------------------------------
//               convert matchibility from RVec to bitstring
// ------------------------------------------------------------------------

RVec<float> minDR_jets_gtau(unsigned int NcleanJets, RVec<float> &cleanJet_eta, RVec<float> &cleanJet_phi, unsigned int NisGood, RVec<float> &GoodTau_eta, RVec<float> &GoodTau_phi) {
  RVec<float> minDR(NcleanJets);
  for (int i = 0; i < NcleanJets; i++) {
    minDR[i] = DeltaR(GoodTau_eta[0], cleanJet_eta[i], GoodTau_phi[0], cleanJet_phi[i]);
    for (int j = 1; j < NisGood; j++) {
      float tempMinDR = DeltaR(GoodTau_eta[j], cleanJet_eta[i], GoodTau_phi[j], cleanJet_phi[i]);
      if (tempMinDR < minDR[i]) {
        minDR[i] = tempMinDR;
      }
    }
  }
  return minDR;
}



// -------------------------------------------------------
//               figure out this event’s taus
// -------------------------------------------------------

RVec<int> decayType(bool isSig, unsigned int nGenPart, RVec<int> &GenPart_pdgId, RVec<float> &GenPart_mass, RVec<float> &GenPart_pt, RVec<float> &GenPart_phi, RVec<float> &GenPart_eta, RVec<short> &GenPart_genPartIdxMother, RVec<int> &GenPart_status)
{
  if (isSig) {
    // Find the taus and their indices
    int tauArr[4];
    int tau = 0;
    for(unsigned int i = 0; i < nGenPart; i++){
      if (abs(GenPart_pdgId[i]) == 15 && abs(GenPart_pdgId[GenPart_genPartIdxMother[i]]) == 9000005) { // && GenPart_pdgId[GenPart_genPartIdxMother[i]] >= 6000000 && GenPart_pdgId[GenPart_genPartIdxMother[i]] <= 6000009
        // It is a tau

        int igen = i;
        for(unsigned int j = i; j < nGenPart; j++){
          if (GenPart_pdgId[j] != GenPart_pdgId[i]) {
            continue;
          }
          if (GenPart_genPartIdxMother[j] != igen) {
            continue;
          }
          igen = j;
        }
        tauArr[tau] = igen;
        tau++;
      }
    }

    // See if they decay to E or Mu
    int countTauE = 0;
    int countTauMu = 0;
    for(unsigned int i = 0; i < nGenPart; i++){
      if (abs(GenPart_pdgId[i]) == 11) {
        // It is a electron
        for(int j = 0; j < 4; j++) {
          if (GenPart_genPartIdxMother[i] == tauArr[j]) { //see if mother is tau
            countTauE++;
          } 
        }
      }

      if (abs(GenPart_pdgId[i]) == 13) {
        // It is a muon
        for(int j = 0; j < 4; j++) {
          if (GenPart_genPartIdxMother[i] == tauArr[j]) {
            countTauMu++;
          } 
        }
      }
    }
    
    RVec<int> counts{countTauE, countTauMu};
    return counts;
  }
  RVec<int> counts{0};
  return counts;
};

// ---------------------------------------------------------------
//                  decays including e,mu,pions
// ---------------------------------------------------------------

RVec<RVec<double>> GenInfo(bool isSig, unsigned int nGenPart, RVec<int> &GenPart_pdgId, RVec<float> &GenPart_mass, RVec<float> &GenPart_pt, RVec<float> &GenPart_phi, RVec<float> &GenPart_eta, RVec<short> &GenPart_genPartIdxMother)
{	
	RVec<double> ElecPt;
	RVec<double> ElecDRb;
	RVec<double> MuonPt;
	RVec<double> MuonDRb;
	RVec<double> TauPt;
	RVec<double> TauEta;
	RVec<double> TauDRb;
	RVec<double> Bvisible = {0,0,0,0,0,0,0,0};
	RVec<double> Neutrino;

	if(isSig)
	{
		//cout << "-----------------------------------------------------" << endl;
		int tauArrB[2] = {-1, -1};
		int tauArrBbar[2] = {-1, -1};
		int tauB = 0;
		int tauBbar = 0;

		TLorentzVector B;
		TLorentzVector Bbar;

		TLorentzVector b1;
		TLorentzVector b2;
		bool haveb1 = false;
		bool haveb2 = false;
		
		TLorentzVector neut;

		for(unsigned int i = 0; i < nGenPart; i++)
		{
			//first finding the b jets (1)		
			if (abs(GenPart_pdgId[i]) == 5 && GenPart_pdgId[GenPart_genPartIdxMother[i]] == 9000005 && haveb1 == false) 
			{
				b1.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], 4.18);
				B += b1;

				//cout << "b1: " << b1.Pt() << " " << b1.M() << endl;
				//cout << "B is " << B.Pt() << " " << B.M() << endl;

				haveb1 = true;
			}
		
			//here is the second b jet
			if (abs(GenPart_pdgId[i]) == 5 && GenPart_pdgId[GenPart_genPartIdxMother[i]] == -9000005 && haveb2 == false) 
			{
				b2.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], 4.18);
				Bbar += b2;
				
				//cout << "b2: " << b2.Pt() << " " << b2.M() << endl;
				//cout << "Bbar is " << Bbar.Pt() << " " << Bbar.M() << endl;
			
				haveb2 = true;
			}	

			//now finding the taus from B
			if (abs(GenPart_pdgId[i]) == 15 && GenPart_pdgId[GenPart_genPartIdxMother[i]] == 9000005) 
			{
				double mass = 1.77686;
				
				TLorentzVector taus1, taus1A;
				taus1.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], mass);

				TauPt.push_back(taus1.Pt()); 
				TauEta.push_back(taus1.Eta());
				TauDRb.push_back(taus1.DeltaR(b1));

				//cout << "genID " << GenPart_pdgId[i] << " genMother " << GenPart_genPartIdxMother[i] << " genIdofMother " << GenPart_pdgId[GenPart_genPartIdxMother[i]] << endl;

				int igen = i;
				//cout << "found tau " << igen << " " << i << endl;
				for(unsigned int j = i; j < nGenPart; j++)
				{
					if (GenPart_pdgId[j] != GenPart_pdgId[i]) {continue;} //is different particle?
					if (GenPart_genPartIdxMother[j] != igen) {continue;}  
					
					igen = j;
				}

				tauArrB[tauB] = igen;
				//cout << "tauEB " << tauEsB[tauB] << " mass " << taus1.M() << " pt " << taus1.Pt() << " eta " << taus1.Eta() << " phi " << taus1.Phi() << endl;
				
				taus1A.SetPtEtaPhiM(GenPart_pt[igen], GenPart_eta[igen], GenPart_phi[igen], mass);
				//cout << "tauEB igen mass " << taus1A.M() << " pt " << taus1A.Pt() << " eta " << taus1A.Eta() << " phi " << taus1A.Phi() << " E " << taus1A.E() << " otrE " << tauEsB[tauB] << endl;
				
				//cout << "stored tau " << igen << " " << tauB << " " << tauArrB[tauB] << endl;
        			tauB++;
			}
		
			//finding the taus from Bbar
			if (abs(GenPart_pdgId[i]) == 15 && GenPart_pdgId[GenPart_genPartIdxMother[i]] == -9000005) 
			{
				double mass = 1.77686;

				TLorentzVector taus2, taus2A;
				taus2.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], mass);

				TauPt.push_back(taus2.Pt()); 
				TauEta.push_back(taus2.Eta());
				TauDRb.push_back(taus2.DeltaR(b2));

				int igen = i;
				//cout << "found tau " << igen << " " << i << endl;
				for(unsigned int j = i; j < nGenPart; j++)
				{
					if (GenPart_pdgId[j] != GenPart_pdgId[i]) {continue;}
					if (GenPart_genPartIdxMother[j] != igen) {continue;}
					
					igen = j;
				}
			
				tauArrBbar[tauBbar] = igen;
				//cout << "tauEBbar " << tauEsB[tauBbar] << " mass " << taus2.M() << " pt " << taus2.Pt() << " eta " << taus2.Eta() << " phi " << taus2.Phi() << endl;
				
				taus2A.SetPtEtaPhiM(GenPart_pt[igen], GenPart_eta[igen], GenPart_phi[igen], mass);
				//cout << "tauEBbar igen mass " << taus2A.M() << " pt " << taus2A.Pt() << " eta " << taus2A.Eta() << " phi " << taus2A.Phi() << " E " << taus2A.E() << " otrE " << tauEsBbar[tauBbar]<< endl;
				
				//cout << "stored tau " << igen << " " << tauBbar << " " << tauArrBbar[tauBbar] << endl;
        			tauBbar++;
			}

		}
		
		//cout << "tauArrB: " << tauArrB[0] << " " << tauArrB[1] << "tauArrBbar:" << tauArrBbar[0] << " " << tauArrBbar[1] << endl;
				
		//matches of e, mu, pis, in tauArrB or tauArrBbar; and now adding neutrinos
		for(unsigned int i = 0; i < nGenPart; i++)
		{
      			if (abs(GenPart_pdgId[i]) == 11 || abs(GenPart_pdgId[i]) == 13 || abs(GenPart_pdgId[i]) == 111 || abs(GenPart_pdgId[i]) == 211)  
			{
				double mass = 0.000511;
				if (abs(GenPart_pdgId[i]) == 13) { mass = 0.1057; }
				else if (abs(GenPart_pdgId[i]) == 111) { mass = 0.135; }
				else if (abs(GenPart_pdgId[i]) == 211) { mass = 0.139570; }

        			for(int j = 0; j < 2; j++) 
				{
          				if (GenPart_genPartIdxMother[i] == tauArrB[j]) 
					{
            					TLorentzVector theVec;
						theVec.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], mass);
						
						B += theVec;
							
						if (abs(GenPart_pdgId[i]) == 11) {	
							ElecPt.push_back(theVec.Pt()); 
							ElecDRb.push_back(theVec.DeltaR(b1));
						} else if (abs(GenPart_pdgId[i]) == 13) {
							MuonPt.push_back(theVec.Pt()); 
							MuonDRb.push_back(theVec.DeltaR(b1));
						}
					
						//cout << "theVec: " << theVec.Pt() << " " << theVec.M() << " " << GenPart_pdgId[i] << " " << j << endl;
						//cout << "B is " << B.Pt() << " " << B.M() << endl;
          				} 

          				if (GenPart_genPartIdxMother[i] == tauArrBbar[j]) 
					{
            					TLorentzVector theVec;
						theVec.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], mass);
						
						Bbar += theVec;
						
						if (abs(GenPart_pdgId[i]) == 11) {	
							ElecPt.push_back(theVec.Pt()); 
							ElecDRb.push_back(theVec.DeltaR(b2));
						} else if (abs(GenPart_pdgId[i]) == 13) {
							MuonPt.push_back(theVec.Pt()); 
							MuonDRb.push_back(theVec.DeltaR(b2));
						}
	
						//cout << "theVec: " << theVec.Pt() << " " << theVec.M() << " " << GenPart_pdgId[i] << " " << j << endl;
						//cout << "Bbar is " << Bbar.Pt() << " " << Bbar.M() << endl;

          				} 
        			}
      			}
    		
      			if (abs(GenPart_pdgId[i]) == 12 || abs(GenPart_pdgId[i]) == 14 || abs(GenPart_pdgId[i]) == 16)
			{
				double mass = 0.0;
        			
				for(int j = 0; j < 2; j++) 
				{
          				if (GenPart_genPartIdxMother[i] == tauArrB[j]) 
					{
						TLorentzVector theNuVec;
						theNuVec.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], mass);
						
						neut += theNuVec;
						//cout << "theNuVec: " << theNuVec.Pt() << " " << theNuVec.Phi() << " " << theNuVec.M() << endl;
						//cout << "neut: " << neut.Pt() << " " << neut.Phi() << " " << neut.M() << endl;					
					} 
          				
					if (GenPart_genPartIdxMother[i] == tauArrBbar[j]) 
					{
						TLorentzVector theNuVec;
						theNuVec.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], mass);
						
						neut += theNuVec;
						//cout << "theNuVec: " << theNuVec.Pt() << " " << theNuVec.Phi() << " " << theNuVec.M() << endl;
						//cout << "neut: " << neut.Pt() << " " << neut.Phi() << " " << neut.M() << endl;					
					
					}
        			}
			}	
      			
		}

		Bvisible = {B.Pt(), B.Eta(), B.Phi(), B.M(), Bbar.Pt(), Bbar.Eta(), Bbar.Phi(), Bbar.M()};

		Neutrino = {neut.Pt(), neut.Phi()};	
	} 

	RVec<RVec<double>> theOutput = {ElecPt, ElecDRb, MuonPt, MuonDRb, TauPt, TauDRb, Bvisible, Neutrino, TauEta};

		return theOutput;
}

// -------------------------------------------------------------
//                understanding the tau bug
// -------------------------------------------------------------

RVec<RVec<double>> tauBUG(bool isSig, unsigned int nGenPart, RVec<int> &GenPart_pdgId, RVec<float> &GenPart_mass, RVec<float> &GenPart_pt, RVec<float> &GenPart_phi, RVec<float> &GenPart_eta, RVec<short> &GenPart_genPartIdxMother, RVec<float> &GenVisTau_pt, RVec<float> &GenVisTau_eta, RVec<float> &GenVisTau_phi, RVec<float> &GenVisTau_mass, RVec<short> &GenVisTau_genPartIdxMother, RVec<unsigned char> &GenVisTau_status, unsigned int nGenVisTau) {

	RVec<double> piplusTauR;
	RVec<double> pinegTauR;
	float tauMass = 1.77686;
	float pionMass = 0.13957;
	

	if (isSig) {
 	 		
		for (int i = 0; i < nGenVisTau; i++) {
			
			if (GenVisTau_status[i] != 0) { continue; }
			TLorentzVector genVisTau, genPartTau;
			genVisTau.SetPtEtaPhiM(GenVisTau_pt[i], GenVisTau_eta[i], GenVisTau_phi[i], tauMass);
			
			float minDR = 999999;
			int minDRi = -1;
			
			for (int j = 0; j < nGenPart; j++) {
				
				if (abs(GenPart_pdgId[j]) == 15 && abs(GenPart_pdgId[GenPart_genPartIdxMother[j]]) == 9000005) {
						
					int igen = j;

					for(unsigned int ii = j; ii < nGenPart; ii++) {
						if (GenPart_pdgId[ii] != GenPart_pdgId[j]) {continue;}
						if (GenPart_genPartIdxMother[ii] != igen) {continue;}
						igen = ii;
					}
							
					genPartTau.SetPtEtaPhiM(GenPart_pt[igen], GenPart_eta[igen], GenPart_phi[igen], tauMass);
					
					float tempDR = genVisTau.DeltaR(genPartTau);
					if (tempDR < minDR) {
						minDR = tempDR;
						minDRi = igen;
					}	
				}

			}	

			TLorentzVector PionNeg, PionPos;


			for (int k = 0; k < nGenPart; k++) {
				
				//cout << "----------------------------------------------------" << endl;
				if (GenPart_pdgId[k] == 211 && GenPart_genPartIdxMother[k] == minDRi) {

					PionPos.SetPtEtaPhiM(GenPart_pt[k], GenPart_eta[k], GenPart_phi[k], pionMass);
					piplusTauR.push_back( (PionPos.E() / genPartTau.E()) );

					if ( (PionPos.E() / genPartTau.E()) > 5 ) {
						//cout << "Pion Pos stuff " << PionPos.Pt() << " " << PionPos.E() << " " << PionPos.Eta() << endl;
						//cout << "Pion Pos from Gen " << GenPart_pt[k] << "          " << GenPart_eta[k] << endl;		
					}

				}
				
				if (GenPart_pdgId[k] == -211 && GenPart_genPartIdxMother[k] == minDRi) {

					PionNeg.SetPtEtaPhiM(GenPart_pt[k], GenPart_eta[k], GenPart_phi[k], pionMass);
					pinegTauR.push_back( (PionNeg.E() / genPartTau.E()) );
					
					if ( (PionNeg.E() / genPartTau.E()) > 5 ) {
						//cout << "Pion Neg stuff " << PionNeg.Pt() << " " << PionNeg.E() << " " << PionNeg.Eta() << endl;		
						//cout << "Pion Neg from Gen " << GenPart_pt[k] << "          " << GenPart_eta[k] << endl;		
					}

				}

			}



			//cout << "minDR: " << minDR << endl;			
			//cout << "Vis tau: " << genVisTau.Pt() << " " << genVisTau.Eta() << " " << genVisTau.Phi() << " " << genVisTau.M() << " " << genVisTau.E() <<endl;
			//cout << "Part tau: " << genPartTau.Pt() << " " << genPartTau.Eta() << " " << genPartTau.Phi() << " " << genPartTau.M() << " " << genVisTau.E() << endl;
		}	

	}
	
	RVec<RVec<double>> returnVec = {piplusTauR, pinegTauR};
	return returnVec;
}

// -------------------------------------------------------
//               figure out this event’s taus and bjets
// -------------------------------------------------------

RVec<RVec<int>> recoGenMatch(bool isSig, unsigned int NisGood, RVec<float> &GoodTau_eta, RVec<float> &GoodTau_phi, unsigned int NJets_DeepFlavM, RVec<float> &gcJet_DeepFlav, RVec<float> &gcBJet_eta, RVec<float> &gcBJet_phi, RVec<float> &gcBJet_pt, unsigned int nGenPart, RVec<int> &GenPart_pdgId, RVec<float> &GenPart_mass, RVec<float> &GenPart_pt, RVec<float> &GenPart_phi, RVec<float> &GenPart_eta, RVec<short> &GenPart_genPartIdxMother, RVec<int> &GenPart_status) {

  // Variables
  int matched = 0;
  int matchedB = 0;
  int tau = 0;
  int tauH = 0;
  int bjet = 0;
  int bjetcount = 0;
  int tauIndex = 0;
  int vecb = 0;
  int bjetIndex = NisGood;
  int HadrTau = 0;

  // ---------------------- Variables ----------------------
  RVec<int> tauArr(4);
  RVec<int> tauArrOriginal(4);
  RVec<int> bjetArr(2);
  RVec<int> bjetArrOriginal(2);
  RVec<int> vecbArr(2);
  // For tau numbers
  int b = 1;
  int antiB = 4;

  // Things for Topograph
  int numOfObjects = NisGood + NJets_DeepFlavM;
  RVec<int> ObjectList_indices(numOfObjects);
  // Bits from 0 to 5: b1 t1 t1 b2 t2 t2 (1 if matched, 0 if not) (1= from bottom, 2=from antibottom)
  RVec<int> matchability = {0, -1, 0, 0, 0, 0, 0, 0};
  

  if (isSig) {
    // ---------------------- Find Gens ----------------------

    // Find the gen taus and their indices
    for(unsigned int i = 0; i < nGenPart; i++){
      
      if (abs(GenPart_pdgId[i]) == 15 && abs(GenPart_pdgId[GenPart_genPartIdxMother[i]]) == 9000005) { 
        // It is a tau
        if (vecb == 0) {
          vecbArr[vecb] = GenPart_genPartIdxMother[i];
          vecb++;
        } else if (vecb == 1) {
          if (GenPart_genPartIdxMother[i] != vecbArr[0]) {
            vecbArr[vecb] = GenPart_genPartIdxMother[i];
            vecb++;
          }
        }
        
        int igen = i;
        for(unsigned int j = i; j < nGenPart; j++){
          if (GenPart_pdgId[j] != GenPart_pdgId[i]) {
            continue;
          }
          if (GenPart_genPartIdxMother[j] != igen) {
            continue;
          }
          igen = j;
        }

        tauArr[tau] = igen;
        tauArrOriginal[tau] = i;
        tau++;
        tauH++;
      }
    }

    // See if they come from E or Mu. If so, remove them from the list
    for(unsigned int i = 0; i < nGenPart; i++){
      if (abs(GenPart_pdgId[i]) == 11 || abs(GenPart_pdgId[i]) == 13) {
        // It is a electron
        for(int j = 0; j < 4; j++) {
          if (GenPart_genPartIdxMother[i] == tauArr[j]) {
            tauArr[j] = -1;
            tauArrOriginal[j] = -1;
            tauH--;
          }
        }
      }
    }

    // Find the bjets and their indices
    for(unsigned int i = 0; i < nGenPart; i++){
      if (abs(GenPart_pdgId[i]) == 5 && abs(GenPart_pdgId[GenPart_genPartIdxMother[i]]) == 9000005) { // && GenPart_pdgId[GenPart_genPartIdxMother[i]] >= 6000000 && GenPart_pdgId[GenPart_genPartIdxMother[i]] <= 6000009
        // It is a bjet
        int igen = i;
        for(unsigned int j = i; j < nGenPart; j++){
          if (GenPart_pdgId[j] != GenPart_pdgId[i]) {
            continue;
          }
          if (GenPart_genPartIdxMother[j] != igen) {
            continue;
          }
          igen = j;
        }
        bjetArr[bjet] = igen;
        bjetArrOriginal[bjet] = i;
        bjet++;
      }
    }

    HadrTau = tauH;	    

    // ---------------------- Matching ----------------------


    // Match Taus: Loop through Gen Taus
    for(unsigned int i = 0; i < tau; i++){
      if (tauArr[i] > -1) {
        tauIndex = 0;
        // Loop through reco taus if gen tau is hadronic
        for(unsigned int j = 0; j < NisGood; j++) {
          // Delta R (eta, phi)
          double delta = DeltaR(GoodTau_eta[j], GenPart_eta[tauArr[i]], GoodTau_phi[j], GenPart_phi[tauArr[i]]);
          if (delta < 0.2) { 
            matched++;
            if (GenPart_pdgId[GenPart_genPartIdxMother[tauArrOriginal[i]]] > 0) {
              ObjectList_indices[tauIndex] = b;
              matchability[b+2] = 1; //index 3,4
              b++;
              tauH++;
            } else {
              ObjectList_indices[tauIndex] = antiB;
              matchability[antiB+2] = 1; //index 6,7
              antiB++;
              tauH++;
            }
          } else {
            if (i == 0) {
              ObjectList_indices[tauIndex] = -1;
            } // Sets all to -1 initially except those matched with first tau
          }
          tauIndex++;
        }
      }
    }
    // Match B Jets
    for(unsigned int i = 0; i < bjet; i++){
      double delta = DeltaR(gcBJet_eta[0], GenPart_eta[bjetArr[i]], gcBJet_phi[0], GenPart_phi[bjetArr[i]]);
      bjetIndex = NisGood;
      for(unsigned int j = 0; j < NJets_DeepFlavM; j++) {
        // Delta R (eta, phi)
        double delta = DeltaR(gcBJet_eta[j], GenPart_eta[bjetArr[i]], gcBJet_phi[j], GenPart_phi[bjetArr[i]]);
        
        if (delta < 0.4) { 
          matchedB++;

          if (GenPart_pdgId[GenPart_genPartIdxMother[bjetArrOriginal[i]]] > 0) {
            ObjectList_indices[bjetIndex] = 0;
            matchability[2] = 1;
            bjetcount++;
          } else {
            ObjectList_indices[bjetIndex] = 3;
            matchability[5] = 1;
            bjetcount++;
          }
        } else {
          if (i == 0) {
            ObjectList_indices[bjetIndex] = -1;
          }
        }
        bjetIndex++;
      }
    }


  }

  //cout << "My Htau: " << HadrTau << " otr tau " << tauH << endl;

  // int nObjects = tauH + bjet + vecb;
  int nObjects = tauH + bjetcount;
  return {ObjectList_indices, matchability, {nObjects}, {bjet}, tauArr, bjetArr, vecbArr, {HadrTau}};
}


// ------------------------------------------------------------------------
//               Make Object list for Topograph
// ------------------------------------------------------------------------


RVec<RVec<float>> ObjectList(bool isSig, unsigned int NisGood, RVec<float> &GoodTau_pt, RVec<float> &GoodTau_eta, RVec<float> &GoodTau_phi, RVec<float> &GoodTau_energy, unsigned int NJets_DeepFlavM, RVec<float> &gcBJet_pt, RVec<float> &gcBJet_eta, RVec<float> &gcBJet_phi, RVec<float> &gcBJet_energy) {
  int numOfObjects = NisGood + NJets_DeepFlavM;
  RVec<float> PtList(numOfObjects);
  RVec<float> EtaList(numOfObjects);
  RVec<float> PhiList(numOfObjects);
  RVec<float> EnergyList(numOfObjects);
  RVec<float> TaggedList(numOfObjects);
  if (isSig) {
    // rearrange so pt is RVec 1 and rearrange to correct format later
    // transpose
    for(int i = 0; i < NisGood; i++) {
      PtList[i] = GoodTau_pt[i];
      EtaList[i] = GoodTau_eta[i];
      PhiList[i] = GoodTau_phi[i];
      EnergyList[i] = GoodTau_energy[i];
      TaggedList[i] = 0;
    } 

    for(int i = NisGood; i < numOfObjects; i++) {
      PtList[i] = gcBJet_pt[i];
      EtaList[i] = gcBJet_eta[i];
      PhiList[i] = gcBJet_phi[i];
      EnergyList[i] = gcBJet_energy[i];
      TaggedList[i] = 1;
    }
  }
  RVec<RVec<float>> ObjectList = {PtList, EtaList, PhiList, EnergyList, TaggedList};
  return ObjectList;
}

// ------------------------------------------------------------------------
//               Make Parton list for Topograph
// ------------------------------------------------------------------------

// Transpose here as well
RVec<RVec<float>> GenList(bool isSig, RVec<int> tauArr, RVec<int> bjetArr, RVec<int> vecbArr, unsigned int nGenPart, RVec<int> &GenPart_pdgId, RVec<float> &GenPart_pt, RVec<float> &GenPart_eta, RVec<float> &GenPart_phi, RVec<float> &GenPart_mass) {
  int num = 8;
  for (int i = 0; i < 4; i++) {
    if (tauArr[i] == -1) {
      //cout << "-1" << endl;
      num--;
    }
  }
  RVec<float> pdgIdList(num);
  RVec<float> ptList(num);
  RVec<float> etaList(num);
  RVec<float> phiList(num);
  RVec<float> massList(num);
  RVec<int> objectArrTemp = Concatenate(vecbArr, bjetArr); 
  RVec<int> objectArr = Concatenate(objectArrTemp, tauArr); 

  // tau = 1.777
  // b = 4.18

  if (isSig) {
    float mass = 1.777;
    for (int i = 0; i < num; i++) {
      if (objectArr[i] > -1) {
        mass = GenPart_mass[objectArr[i]];
        if (i >= 2) {
          mass = 4.18;
        }
        if (i >= 4) {
          mass = 1.777;
        }
        // cout << "Index: " << objectArr[i] << " ID: " << GenPart_pdgId[objectArr[i]] << " Mass: " << mass << endl;
        pdgIdList[i] = static_cast<float>(GenPart_pdgId[objectArr[i]]);
        ptList[i] = GenPart_pt[objectArr[i]];
        etaList[i] = GenPart_eta[objectArr[i]];
        phiList[i] = GenPart_phi[objectArr[i]];
        massList[i] = mass;
      }
    }
  }
  RVec<RVec<float>> GenList = {pdgIdList, ptList, etaList, phiList, massList};
  return GenList;
}

