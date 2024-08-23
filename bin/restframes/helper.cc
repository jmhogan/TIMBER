 // processDecayTree(), returnVectors()
// These two functions help to differentiate between the bW and (H/Z)t trees
#include <iostream>

RVec<double> processDecayTree(Tprime_RestFrames_Container_W * W_rfc, Tprime_RestFrames_Container_t * t_rfc, int thread_index, float lepton_pt, float lepton_eta, float lepton_phi, float lepton_mass, RVec<float> fatjet_pt, RVec<float> fatjet_eta, RVec<float> fatjet_phi, RVec<float> fatjet_mass, RVec<float> fatjet_DeepFlav, float met_pt, float met_phi, RVec<float> jet_pt, RVec<float> jet_eta, RVec<float> jet_phi, RVec<float> jet_mass, RVec<float> jet_DeepFlav, RVec<float> isoAK4) {
  RVec<TLorentzVector> jets;
  int i = 0; 
  // Make an RVec of valid jets for the possible b.
  TLorentzVector jet;
  for (; i < isoAK4.size(); i++) {
    //  stand alone         b-tagged
    if (isoAK4[i] == 1 && jet_DeepFlav[i] > 0.9) {
      jet.SetPtEtaPhiM(jet_pt[i], jet_eta[i], jet_phi[i], jet_mass[i]);
      jets.push_back(jet);
    } 
  }   

  RVec<double> result;
  if (jets.size() == 0) { // analyze bW tree
    result = W_rfc->return_doubles(thread_index, lepton_pt, lepton_eta, lepton_phi, lepton_mass, fatjet_pt, fatjet_eta, fatjet_phi, fatjet_mass, fatjet_DeepFlav, met_pt, met_phi);
    //cout << "\tW tree!\n";

    result.push_back(0.0); // 0 is for W tree
  } else { // analyze (H/Z)t tree
    result = t_rfc->return_doubles(thread_index, lepton_pt, lepton_eta, lepton_phi, lepton_mass, fatjet_pt, fatjet_eta, fatjet_phi, fatjet_mass, fatjet_DeepFlav, met_pt, met_phi, jets);

    //cout << "\tt tree!\n";
    result.push_back(1.0); // 1 is for t tree
  }
  return result;
  // do they need: gcJet_DeepFlav ?
}

RVec<TLorentzVector> returnVectors(Tprime_RestFrames_Container_W * W_rfc, Tprime_RestFrames_Container_t * t_rfc, int thread_index, RVec<float> jet_DeepFlav, RVec<float> isoAK4) {
  int i = 0;
  for (; (isoAK4[i] != 1 || jet_DeepFlav[i] <= 0.9) && i < isoAK4.size(); i++); // find the first stand alone ak4 b-tagged jet
  
  RVec<TLorentzVector> result;
  if (i == isoAK4.size()) { // didn't find a good b jet so we're at the bW decay tree
    result = W_rfc->return_vecs(thread_index);
  } else { // we're at the (H/Z)t tree
    result = t_rfc->return_vecs(thread_index);
  }
  return result;
}

