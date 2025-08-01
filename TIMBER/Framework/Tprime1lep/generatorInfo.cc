// Methods in this file:
// genttbarMassCalc() decayModeSelection() 	//TODO BtoTW generatorInfo.cc has a lot more funcs
#include <iostream>
using namespace ROOT::VecOps;
using namespace std;

// ttbar background mass CALCULATOR:
int genttbarMassCalc(string sample, unsigned int nGenPart, RVec<int> &GenPart_pdgId, RVec<float> &GenPart_mass, RVec<float> &GenPart_pt, RVec<float> &GenPart_phi, RVec<float> &GenPart_eta, RVec<int> &GenPart_genPartIdxMother, RVec<int> &GenPart_status)
{
  int returnVar = 0;
  if (sample.find("TTTo") != std::string::npos || sample.find("Mtt") != std::string::npos)
    {
      int genTTbarMass = -999;
      double topPtWeight = 1.0;
      TLorentzVector top, antitop;
      bool gottop = false;
      bool gotantitop = false;
      bool gottoppt = false;
      bool gotantitoppt = false;
      float toppt, antitoppt;
      for (unsigned int p = 0; p < nGenPart; p++)
      {
          int id = GenPart_pdgId[p];
          if (abs(id) != 6)
          {
              continue;
          }
          if (GenPart_mass[p] < 10)
          {
              continue;
          }
          int motherid = GenPart_pdgId[GenPart_genPartIdxMother[p]];
          if (abs(motherid) != 6)
          {
              if (!gottop && id == 6)
              {
                  top.SetPtEtaPhiM(GenPart_pt[p], GenPart_eta[p], GenPart_phi[p], GenPart_mass[p]);
                  gottop = true;
              }
              if (!gotantitop && id == -6)
              {
                  antitop.SetPtEtaPhiM(GenPart_pt[p], GenPart_eta[p], GenPart_phi[p], GenPart_mass[p]);
                  gotantitop = true;
              }
          }
          if (GenPart_status[p] == 62)
          {
              if (!gottoppt && id == 6)
              {
                  toppt = GenPart_pt[p];
                  gottoppt = true;
              }
              if (!gotantitoppt && id == -6)
              {
                  antitoppt = GenPart_pt[p];
                  gotantitoppt = true;
              }
          }
      }
      if (gottop && gotantitop)
      {
          genTTbarMass = (top + antitop).M();
      }
      if (gottoppt && gotantitoppt)
      {
          float SFtop = TMath::Exp(0.0615 - 0.0005 * toppt);
          float SFantitop = TMath::Exp(0.0615 - 0.0005 * antitoppt);
          topPtWeight = TMath::Sqrt(SFtop * SFantitop);
      }
      returnVar = genTTbarMass;
  }
  return returnVar;
};

// ----------------------------------------------------
//   		DECAY CALCULATOR:
// ----------------------------------------------------
int decayModeSelection(string region, unsigned int nGenPart, ROOT::VecOps::RVec<int>& GenPart_pdgId, ROOT::VecOps::RVec<float>& GenPart_mass, ROOT::VecOps::RVec<float>& GenPart_pt, ROOT::VecOps::RVec<float>& GenPart_phi, ROOT::VecOps::RVec<float>& GenPart_eta, ROOT::VecOps::RVec<int>& GenPart_genPartIdxMother, ROOT::VecOps::RVec<int>& GenPart_status)
{
  //  std::cout << "Hello! Made it to decayModeSelection!" << std::endl;
  std::vector<int> tPrimeID;
  std::vector<int> bPrimeID;
  std::vector<int> listofQuarkIDs;
  std::vector<int> listofBosonIDs;
  std::vector<unsigned int> quarks;
  std::vector<unsigned int> bosons;
  std::vector<int> wIDs;
  std::vector<int> tIDs;

  bool isBWBW = false;
  bool isTZTZ = false;
  bool isTHTH = false;
  bool isTZTH = false;
  bool isTZBW = false;
  bool isTHBW = false;
  
  bool isTWTW = false;
  bool isBZBZ = false;
  bool isBHBH = false;
  bool isBZBH = false;
  bool isBZTW = false;
  bool isBHTW = false;
  
  int decayMode = 0;
  
  tPrimeID.clear();
  bPrimeID.clear();
  listofQuarkIDs.clear();
  listofBosonIDs.clear();
  quarks.clear();
  bosons.clear();
  wIDs.clear();
  tIDs.clear();

  // format of decayMode: [0-2 leptonic t's] [0-2 leptonic W's] [0-11 Bprime and Tprime decay modes]

  for(unsigned int p = 0; p < nGenPart; p++)
    {
      int id=GenPart_pdgId[p];
      // find T' and B' particles
      if(abs(id) != 6000006 && abs(id) != 6000005){continue;}
      
      bool hasTdaughter = false;
      vector<unsigned int> daughters;
      daughters.clear();
      for(unsigned int  dau = 0; dau < nGenPart; dau++)
  	{
  	if(GenPart_genPartIdxMother[dau]!=p){continue;}
  	daughters.push_back(dau);
  	if(abs(id) == 6000006 && abs(GenPart_pdgId[dau]) == 6000006){hasTdaughter = true;}
  	if(abs(id) == 6000005 && abs(GenPart_pdgId[dau]) == 6000005){hasTdaughter = true;}
  	}
      if(hasTdaughter){continue;}
      int mother = GenPart_genPartIdxMother[p];
      int mother_id = GenPart_pdgId[mother];
      if(abs(id) == 6000006)
  	{
  	if(abs(mother_id) == 6000006){tPrimeID.push_back(GenPart_pdgId[mother]);}
  	else{tPrimeID.push_back(GenPart_pdgId[p]);}
  	}
      if(abs(id) == 6000005)
  	{
  	if(abs(mother_id) == 6000005){bPrimeID.push_back(GenPart_pdgId[mother]);}
  	else{bPrimeID.push_back(GenPart_pdgId[p]);}
  	}
      for(unsigned int j = 0; j < daughters.size(); j++)
  	{
  	unsigned int d = daughters.at(j);
  	int dauId = GenPart_pdgId[d];
  	if(abs(dauId) == 5 || abs(dauId) == 6)
  	  {
  	    quarks.push_back(d);
  	    listofQuarkIDs.push_back(dauId);
  	  }
  	else if(abs(dauId) > 22 && abs(dauId) < 26)
  	  {
  	    bosons.push_back(d);
  	    listofBosonIDs.push_back(dauId);
  	  }
  	else{continue;}
  	}
    }
  
  if(tPrimeID.size() > 0 && bPrimeID.size() > 0) {std::cout << "Found both T' and B' " << std::endl;}
  if(listofQuarkIDs.size() != 0 && listofQuarkIDs.size() != 2)
    {
      std::cout << "More/less than 2 quarks stored: " << listofQuarkIDs.size() << std::endl;
      for(unsigned int i = 0; i < listofQuarkIDs.size(); i++){std::cout << "quark " << i << " = " << listofQuarkIDs.at(i) << std::endl;}
      int test = listofQuarkIDs.at(0)*listofQuarkIDs.at(1);
      int sign = -1;
      if(test > 0){sign = 1;}
      if(sign > 0)
  	{
  	if(listofQuarkIDs.size() == 4)
  	  {
  	    std::swap(listofQuarkIDs.at(2),listofQuarkIDs.at(3));
  	    std::swap(quarks.at(2),quarks.at(3));
  	  }
  	std::swap(listofQuarkIDs.at(1),listofQuarkIDs.at(2));
  	std::swap(quarks.at(1),quarks.at(2));
  	test = listofQuarkIDs.at(0)*listofQuarkIDs.at(1);
  	sign = -1;
  	if(test > 0){sign = 1;}
  	if(sign < 0){std::cout << "Signs are fixed!" << std::endl;}
  	}
      if(listofQuarkIDs.size() > 3 && abs(listofQuarkIDs.at(3)) == 6)
  	{
  	std::swap(listofQuarkIDs.at(2),listofQuarkIDs.at(3));
  	std::swap(quarks.at(2),quarks.at(3));
  	}
      if(listofQuarkIDs.size() > 2 && abs(listofQuarkIDs.at(2)) == 6)
  	{
  	std::swap(listofQuarkIDs.at(1),listofQuarkIDs.at(2));
  	std::swap(quarks.at(1),quarks.at(2));
  	}
    }
  if(listofBosonIDs.size() != 0 && listofBosonIDs.size() != 2)
    {
      std::cout << "More/less than 2 bosons stored: " << listofBosonIDs.size() << std::endl;
    }
  // tag the decay chains according to ID'd quarks and bosons.
  
  // TPrime Decay Mode Selector
  if(tPrimeID.size() > 1 && bPrimeID.size() == 0)
    {
      if(abs(listofQuarkIDs.at(0)) == 5 && abs(listofQuarkIDs.at(1)) == 5)
  	{
  	if(abs(listofBosonIDs.at(0)) == 24 && abs(listofBosonIDs.at(1)) == 24)
  	  {
  	    isBWBW = true;
  	    decayMode = 1; // BWBW ID!
	    wIDs.push_back(bosons.at(0));
	    wIDs.push_back(bosons.at(1));
  	  }
  	}
      // 2 t quarks, check for Z's and H's
      else if(abs(listofQuarkIDs.at(0)) == 6 && abs(listofQuarkIDs.at(1)) == 6)
  	{
	tIDs.push_back(quarks.at(0));
	tIDs.push_back(quarks.at(1));
  	if(listofBosonIDs.at(0) == 23 && listofBosonIDs.at(1) == 23)
  	  {
  	    isTZTZ = true;
  	    decayMode = 2; // TZTZ ID!
  	  }
  	else if(listofBosonIDs.at(0) == 25 && listofBosonIDs.at(1) == 25)
  	  {
  	    isTHTH = true;
  	    decayMode = 3; // THTH ID!
  	  }
  	else if(listofBosonIDs.at(0) == 25 && listofBosonIDs.at(1) == 23)
  	  {
  	    isTZTH = true;
  	    decayMode = 4; //TZTH ID!
  	  }
  	else if(listofBosonIDs.at(0) == 23 && listofBosonIDs.at(1) == 25)
  	  {
  	    isTZTH = true;
  	    decayMode = 4; // TZTH ID!
  	  }
  	else
  	  {
  	    std::cout << "2 t daughters didn't match tZtZ, tHtH, or tZtH" << listofBosonIDs.at(0) << ", " << listofBosonIDs.at(1) << std::endl;
  	  }
  	}
      // t-b pairs, check for correlating bosons in the right spots
      else if(abs(listofQuarkIDs.at(0)) == 6 && abs(listofQuarkIDs.at(1)) == 5)
  	{
	tIDs.push_back(quarks.at(0));
  	if(listofBosonIDs.at(0) == 23 && abs(listofBosonIDs.at(1)) == 24)
  	  {
  	    isTZBW = true;
  	    decayMode = 5; // TZBW ID!
	    wIDs.push_back(bosons.at(1));
  	  }
  	else if(listofBosonIDs.at(0) == 25 && abs(listofBosonIDs.at(1)) == 24)
  	  {
  	    isTHBW = true;
  	    decayMode = 6; // THBW ID!
	    wIDs.push_back(bosons.at(1));
  	  }
  	else{std::cout<< "t - b pair didn't match Z/H - W pair" << listofBosonIDs.at(0)<<", "<<listofBosonIDs.at(1) << std::endl;}
  	}
      // b-t pairs, check for correlating bosons in the right spots
      else if(abs(listofQuarkIDs.at(1)) == 6 && abs(listofQuarkIDs.at(0)) == 5)
  	{
	tIDs.push_back(quarks.at(1));
  	if(listofBosonIDs.at(1) == 23 && abs(listofBosonIDs.at(0)) == 24)
  	  {
  	    isTZBW = true;
  	    decayMode = 5; // TZBW ID!
	    wIDs.push_back(bosons.at(0));
  	  }
  	else if(listofBosonIDs.at(1) == 25 && abs(listofBosonIDs.at(0)) == 24)
  	  {
  	    isTHBW = true;
  	    decayMode = 6; //THBW ID!
	    wIDs.push_back(bosons.at(0));
  	  }
  	else{std::cout<< "b - t pair didn't match W - Z/H pair" << listofBosonIDs.at(0)<<", "<<listofBosonIDs.at(1) << std::endl;}
  	}
      // error messages if we found something else entirely
      else
  	{
  	std::cout << "T' daughters didn't match a recognized pattern" << std::endl;
  	for(size_t i = 0; i < listofQuarkIDs.size(); i++)
  	  {
  	    std::cout << "quark " << i << " = " << listofQuarkIDs.at(i) << std::endl;
  	  }
  	for(size_t i = 0; i < listofBosonIDs.size(); i++)
  	  {
  	    std::cout << "boson " << i << " = " << listofBosonIDs.at(i) << std::endl;
  	  }
  	decayMode = -1;
  	}
    }
  // BPrime Decay Mode Selector
  if(bPrimeID.size() > 1 && tPrimeID.size() == 0)
    {
      // 2 t quarks, check for matching W's
      if(abs(listofQuarkIDs.at(0)) == 6 && abs(listofQuarkIDs.at(1)) == 6)
  	{
	tIDs.push_back(quarks.at(0));
	tIDs.push_back(quarks.at(1));
  	if(abs(listofBosonIDs.at(0)) == 24 && abs(listofBosonIDs.at(1)) == 24)
  	  {
  	    isTWTW = true;
  	    decayMode = 7; // TWTW ID!
	    wIDs.push_back(bosons.at(0));
	    wIDs.push_back(bosons.at(1));
  	  }
  	else{std::cout<< "2 t daughters didn't match tWtW: " <<listofBosonIDs.at(0)<<", "<<listofBosonIDs.at(1) << std::endl;}
  	}
      // 2 b quarks, check for Z's and H's
      else if(abs(listofQuarkIDs.at(0)) == 5 && abs(listofQuarkIDs.at(1)) == 5)
  	{
  	if(listofBosonIDs.at(0) == 23 && listofBosonIDs.at(1) == 23)
  	  {
  	    isBZBZ = true;
  	    decayMode = 8; // BZBZ ID!
  	  }
  	else if(listofBosonIDs.at(0) == 25 && listofBosonIDs.at(1) == 25)
  	  {
  	    isBHBH = true;
  	    decayMode = 9; // BHBH ID!
  	  }
  	else if(listofBosonIDs.at(0) == 25 && listofBosonIDs.at(1) == 23)
  	  {
  	    isBZBH = true;
  	    decayMode = 10; // BZBH ID!
  	  }
  	else if(listofBosonIDs.at(0) == 23 && listofBosonIDs.at(1) == 25)
  	  {
  	    isBZBH = true;
  	    decayMode = 10; //BZBH ID!
  	  }
  	else
  	  {
  	    std::cout << "2 b daughters didn't match bZbZ, bHbH, or bZbH" << listofBosonIDs.at(0) << ", " << listofBosonIDs.at(1) << std::endl;
  	  }
  	}
      // b-t pairs, check for correlating bosons in the right spots
      else if(abs(listofQuarkIDs.at(0)) == 5 && abs(listofQuarkIDs.at(1)) == 6)
  	{
	tIDs.push_back(quarks.at(1));
  	if(listofBosonIDs.at(0) == 23 && abs(listofBosonIDs.at(1)) == 24)
  	  {
  	    isBZTW = true;
  	    decayMode = 11; // BZTW ID!
	    wIDs.push_back(bosons.at(1));
  	  }
  	else if(listofBosonIDs.at(0) == 25 && abs(listofBosonIDs.at(1)) == 24)
  	  {
  	    isBHTW = true;
  	    decayMode = 12; // BHTW ID!
	    wIDs.push_back(bosons.at(1));
  	  }
  	else{std::cout<< "b - t pair didn't match Z/H - W pair" << listofBosonIDs.at(0)<<", "<<listofBosonIDs.at(1) << std::endl;}
  	}
      // t-b pairs, check for correlating bosons in the right spots
      else if(abs(listofQuarkIDs.at(1)) == 5 && abs(listofQuarkIDs.at(0)) == 6)
  	{
	tIDs.push_back(quarks.at(0));
  	if(listofBosonIDs.at(1) == 23 && abs(listofBosonIDs.at(0)) == 24)
  	  {
  	    isBZTW = true;
  	    decayMode = 11; // BZTW ID!
	    wIDs.push_back(bosons.at(0));
  	  }
  	else if(listofBosonIDs.at(1) == 25 && abs(listofBosonIDs.at(0)) == 24)
  	  {
  	    isBHTW = true;
  	    decayMode = 12; // BHTW ID!
	    wIDs.push_back(bosons.at(0));
  	  }
  	else{std::cout<< "t - b pair didn't match W - Z/H pair" << listofBosonIDs.at(0)<<", "<<listofBosonIDs.at(1) << std::endl;}
      }
      // error messages if we found something else entirely
      else
      {
  	std::cout << "B' daughters didn't match a recognized pattern" << std::endl;
  	for(size_t i = 0; i < listofQuarkIDs.size(); i++)
  	  {
  	    std::cout << "quark " << i << " = " << listofQuarkIDs.at(i) << std::endl;
  	  }
  	for(size_t i = 0; i < listofBosonIDs.size(); i++)
  	  {
  	    std::cout << "boson " << i << " = " << listofBosonIDs.at(i) << std::endl;
  	  }
  	decayMode = -1;
     }
  }

  if (decayMode < 0) { return decayMode; }
  
  // Look for leptons from W's and t's
  int index, prev, child;
  // W's
  for (int i = 0; i < wIDs.size(); i++) {
    index = prev = wIDs.at(i);
    
    for (++index;index < nGenPart; index++) {
      if (GenPart_genPartIdxMother[index] == prev) { // Found a child of W!
        child = abs(GenPart_pdgId[index]);

	if (child == 24) { // it's a W
	  prev = index;
	} else if (11 <= child && child <= 16) { // it's a lepton!
          decayMode += 100;
	  break;
	} else { break; } // non-leptonic decay
      }
    }
  }

  // t's
  for (int i = 0; i < tIDs.size(); i++) {
    index = prev = tIDs.at(i);
    
    for (++index; index < nGenPart; index++) {
      if (GenPart_genPartIdxMother[index] == prev) { // Found a child of t!
        child = abs(GenPart_pdgId[index]);

	if (child == 6) { // it's a t
	  prev = index;
	} else if (child == 24) { // it's a W
	  prev = index;
	} else if (11 <= child && child <= 16) { // it's a lepton!
          decayMode += 1000;
	  break;
	} else { break; } // non-leptonic decay
      }
    }
  }

  //  std::cout << "Returning decayMode = " << decayMode << std::endl;
  return decayMode;
}

// Wtolnu() makes sure there is a W that decays to a lepton and neutrino
bool Wtolnu(unsigned int nGenPart, ROOT::VecOps::RVec<int>& GenPart_pdgId, ROOT::VecOps::RVec<int>& GenPart_genPartIdxMother, ROOT::VecOps::RVec<int>& GenPart_status)
{
  for (int i = 0; i < 10 && i < nGenPart; i++) {
    cout << "Gen part index: " << i << " has ID " << GenPart_pdgId[i] << " and mother at index " << GenPart_genPartIdxMother[i] << endl;
    if (11 <= GenPart_pdgId[i] && GenPart_pdgId[i] <= 16) { // this is a lepton
      if (abs( GenPart_pdgId[GenPart_genPartIdxMother[i]] ) == 24) { 
        //mother is a W
      }
    }
  }
  return false; //TODO actully write this function
}
