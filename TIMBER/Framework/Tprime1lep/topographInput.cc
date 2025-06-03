// Methods in this file:
// decayType, recoGenMatch

using namespace std;
using namespace ROOT::VecOps;

#include <cmath>
#include <iostream>
#include <string>
#include "TLorentzVector.h"
#include <ROOT/RVec.hxx>

// ------------------------------------------------------------------------
//               convert matchibility from RVec to bitstring
// ------------------------------------------------------------------------

int convertMatchToInt(RVec<int> match) {
  int decimalValue = 0;
  int arraySize = 8;

  for (int i = 2; i < arraySize; ++i) {
      if (match[i] == 1) {
          decimalValue += static_cast<int>(std::pow(2, arraySize - 1 - i));
      }
  }
  return decimalValue;
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

    // See if they come from E or Mu
    int countTauE = 0;
    int countTauMu = 0;
    for(unsigned int i = 0; i < nGenPart; i++){
      if (abs(GenPart_pdgId[i]) == 11) {
        // It is a electron
        for(int j = 0; j < 4; j++) {
          if (GenPart_genPartIdxMother[i] == tauArr[j]) {
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
  int tauIndex = 0;
  int vecb = 0;
  int bjetIndex = NisGood;

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
            } else {
              ObjectList_indices[tauIndex] = antiB;
              matchability[antiB+2] = 1; //index 6,7
              antiB++;
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
          } else {
            ObjectList_indices[bjetIndex] = 3;
            matchability[5] = 1;
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

  int nObjects = tauH + bjet + vecb;
  return {ObjectList_indices, matchability, {nObjects}, {bjet}, tauArr, bjetArr, vecbArr};
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
RVec<RVec<float>> GenList(bool isSig, RVec<int> nObjects, RVec<int> nbjets, RVec<int> tauArr, RVec<int> bjetArr, RVec<int> vecbArr, unsigned int nGenPart, RVec<int> &GenPart_pdgId, RVec<float> &GenPart_pt, RVec<float> &GenPart_eta, RVec<float> &GenPart_phi, RVec<float> &GenPart_mass) {
  RVec<float> pdgIdList(nObjects[0]);
  RVec<float> ptList(nObjects[0]);
  RVec<float> etaList(nObjects[0]);
  RVec<float> phiList(nObjects[0]);
  RVec<float> massList(nObjects[0]);
  RVec<int> objectArrTemp = Concatenate(vecbArr, tauArr); 
  RVec<int> objectArr = Concatenate(objectArrTemp, bjetArr); 

  int Ntaus = nObjects[0] - nbjets[0];  
  //tau = 1.777
  //b = 4.18

  if (isSig) {
    float mass = 1.777;
    for (int i = 0; i < (nObjects[0]); i++) {
      if (i == Ntaus) {
        mass = 4.18;
      }
      pdgIdList[i] = static_cast<float>(GenPart_pdgId[objectArr[i]]);
      ptList[i] = GenPart_pt[objectArr[i]];
      etaList[i] = GenPart_eta[objectArr[i]];
      phiList[i] = GenPart_phi[objectArr[i]];
      massList[i] = mass;
    }
  }
  RVec<RVec<float>> GenList = {pdgIdList, ptList, etaList, phiList, massList};
  return GenList;
}





// ------------------manual reco fxn----------------------
//mass difference function--updated version is in Tprime1lep/topographinput.cc
/*
RVec<float> funcmassdiff(RVec<float> goodtau_pt, RVec<float> goodtau_eta, RVec<float> goodtau_phi, RVec<float> goodtau_mass, RVec<float> goodbjet_pt, RVec<float> goodbjet_eta, RVec<float> goodbjet_phi, RVec<float> goodbjet_mass) {
  TLorentzVector B1, B2;
  TLorentzVector B1final, B2final;
  TLorentzVector b1, b2;
  TLorentzVector tauA, tauB, tauC, tauD;

  float massdiff = 999999999.0;
  float massdifftemp;

  b1.SetPtEtaPhiM(goodbjet_pt[0], goodbjet_eta[0], goodbjet_phi[0], goodbjet_mass[0]);
  b2.SetPtEtaPhiM(goodbjet_pt[1], goodbjet_eta[1], goodbjet_phi[1], goodbjet_mass[1]);

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
  RVec<float> returnVec = {static_cast<float>(B1final.Px()), static_cast<float>(B1final.Py()), static_cast<float>(B1final.Pz()), static_cast<float>(B2final.Px()), static_cast<float>(B2final.Py()), static_cast<float>(B2final.Pz())}; 
  return returnVec;
}*/
