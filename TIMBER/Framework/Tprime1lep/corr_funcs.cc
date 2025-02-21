// Methods in this file:  		There are several corrlib functions and some non corrlib functions.
// goldenjson() pufunc() recofunc() idfunc() isofunc() metfunc() hltfunc() jetvetofunc() 

using namespace ROOT::VecOps;
using namespace std;

// Make sure the run passes the json? //TODO better explanation
bool goldenjson(lumiMask myLumiMask, const unsigned int &run, const unsigned int &luminosityBlock)
{
  return myLumiMask.accept(run, luminosityBlock);
}; 

// Pile Up Function
RVec<double> pufunc(correction::Correction::Ref& pileupcorr, const float &numTrueInt) 
{
  RVec<double> pu = {pileupcorr->evaluate({numTrueInt, "nominal"}), pileupcorr->evaluate({numTrueInt, "up"}), pileupcorr->evaluate({numTrueInt, "down"})};
  return pu;
};

// Reconstruct the lepton?
RVec<double> recofunc(correction::Correction::Ref& electroncorr, correction::Correction::Ref& muoncorr, string yrstr, const float &pt, const float &eta, const bool &isEl)
{
  RVec<double> reco;
  if(isEl == 0) { 
    reco = {muoncorr->evaluate({yrstr+"_UL",abs(eta),pt,"sf"}), 
      muoncorr->evaluate({yrstr+"_UL",abs(eta),pt,"systup"}), 
      muoncorr->evaluate({yrstr+"_UL",abs(eta),pt,"systdown"})};
  }else{
    reco = {electroncorr->evaluate({yrstr,"sf","RecoAbove20",eta,pt}), 
      electroncorr->evaluate({yrstr,"sfup","RecoAbove20",eta,pt}), 
      electroncorr->evaluate({yrstr,"sfdown","RecoAbove20",eta,pt})};
  }
  return reco;
}; 

// Get the ID of leptons?
RVec<float> idfunc(correction::Correction::Ref& muonidcorr, vector<float> &elid_pts, vector<float> &elid_etas, vector<vector<float>> &elecidsfs, vector<vector<float>> &elecidsfuncs, string &yrstr, const float &pt, const float &eta, const bool &isEl)
{
  RVec<float> id;
  if(isEl > 0){
    int ptbin = (std::upper_bound(elid_pts.begin(), elid_pts.end(), pt) - elid_pts.begin())-1;
    int etabin = (std::upper_bound(elid_etas.begin(), elid_etas.end(), eta) - elid_etas.begin())-1;
    
    id = {elecidsfs[ptbin][etabin], elecidsfuncs[ptbin][etabin]}; //PTL
     
  }else{  
    id = {static_cast<float>(muonidcorr->evaluate({yrstr+"_UL",abs(eta),pt,"sf"})), 
      static_cast<float>(muonidcorr->evaluate({yrstr+"_UL",abs(eta),pt,"systup"})), 
      static_cast<float>(muonidcorr->evaluate({yrstr+"_UL",abs(eta),pt,"systdown"}))};
  }
  return id;
}; 

// iso function
RVec<double> isofunc(vector<float> muiso_pts, vector<float> muiso_etas, vector<vector<float>> muonisosfs, float muonisosfunc, vector<float> elid_pts, vector<float> elid_etas, vector<vector<float>> elecisosfs, float elecisosfunc, const float &pt, const float &eta, const bool &isEl)
{
  RVec<double> iso;
  if(isEl > 0){
    int ptbin = (std::upper_bound(elid_pts.begin(), elid_pts.end(), pt) - elid_pts.begin())-1;
    int etabin = (std::upper_bound(elid_etas.begin(), elid_etas.end(), eta) - elid_etas.begin())-1;
    iso = {elecisosfs[ptbin][etabin], elecisosfunc};
  }else{
    int ptbin = (std::upper_bound(muiso_pts.begin(), muiso_pts.end(), pt) - muiso_pts.begin())-1;
    int etabin = (std::upper_bound(muiso_etas.begin(), muiso_etas.end(), abs(eta)) - muiso_etas.begin())-1;
    iso = {muonisosfs[ptbin][etabin], muonisosfunc};
  }
  return iso;
};

// What is MET? = missing transverse energy = neutrinos. It returns the corrections for met's transverse momentum and phi angle. 
RVec<float> metfunc(correction::Correction::Ref& metptcorr, correction::Correction::Ref& metphicorr, const float &met, const float &phi, const int &npvs, const unsigned int &run)
//assigned types as per BtoTW git analyzer_RDF.cc
{
   float floatrun = run;
   float floatnpvs = npvs;
   float tmpmet = met;
      if(tmpmet > 6500) tmpmet = 6499;
         RVec<float> corrmet = {static_cast<float>(metptcorr->evaluate({tmpmet, phi, floatnpvs, floatrun})), static_cast<float>(metphicorr->evaluate({tmpmet, phi, floatnpvs, floatrun}))};
  return corrmet;
};

// IDK what this function does. What's HLT? seems to be calling a bunch of variables under HLT 
// (electrons? and a muon correction)
RVec<double> hltfunc(correction::Correction::Ref& muonhltcorr, vector<float> &elhlt_pts, vector<float> &elhlt_etas, vector<vector<float>> &elechltsfs, vector<vector<float>> &elechltuncs, string &yrstr, const float &pt, const float &eta, const bool &isEl)
// I assumed hltfunc is a double type vector because of the definition of its return value below.
{
    RVec<double> hlt;
    if(isEl > 0){
      int ptbin = (std::upper_bound(elhlt_pts.begin(), elhlt_pts.end(), pt) - elhlt_pts.begin())-1;
      int etabin = (std::upper_bound(elhlt_etas.begin(), elhlt_etas.end(), 
	            abs(eta)) - elhlt_etas.begin())-1;
      hlt = {elechltsfs[ptbin][etabin], elechltuncs[ptbin][etabin]};
    }
    else {
      //we have lookup tables we might put it into correction format.	      
      
      hlt = {1,1,1};  
      /*hlt = {muonhltcorr->evaluate({yrstr+"_UL",abs(eta),pt,"sf"}), 
	     muonhltcorr->evaluate({yrstr+"_UL",abs(eta),pt,"systup"}), 
	     muonhltcorr->evaluate({yrstr+"_UL",abs(eta),pt,"systdown"})};*/
    }
    return hlt;
 }; 

// Jet veto function
RVec<double> jetvetofunc(correction::Correction::Ref& jetvetocorr, const RVec<float> &eta, const RVec<float> &phi){
  RVec<double> map;
  for(unsigned int ijet = 0; ijet < eta.size(); ijet++){
    float phitemp = phi.at(ijet);
    if(phitemp < -3.14159) phitemp = -3.14159;
    else if(phitemp > 3.14159) phitemp = 3.14159;
    map.push_back(jetvetocorr->evaluate({"jetvetomap",eta.at(ijet),phitemp}));
  }
  return map;
};

// WORK ON THIS JULIE
// RVec<float> btagshapefunc(string year, string jesvar, correction::Correction::Ref& btagwpbccorr, correction::Correction::Ref& btagwplcorr, btagpts,btageffs,nominal,deepjetL]const RVec<float> &pt, const RVec<float> &eta, const RVec<float> &disc, const RVec<unsigned char> &flav){

//   std::string nominal = "central";
//   if(jesvar == "JECup") nominal = "up_jes";
//   else if (jesvar == "JECdn") nominal = "down_jes";

//   RVec<float> weights(9, 1.0); // collect product of SFs over jets
//   for(unsigned int ijet = 0; ijet < eta.size(); ijet++){
//       int ptbin = (std::upper_bound(btagpts.begin(), btagpts.end(), pt.at(ijet)) - btagpts.begin())-1;
//       correction::Correction::Ref wpcorr;
//       float eff;
//       string downcorrelated = "down_correlated"; string upcorrelated = "up_correlated";
//       std::vector<int> shift = {0,4};
//       if(flav.at(ijet) < 4){wpcorr = btagwplcorr; eff = btageffs[ptbin][2]; shift = {4,0}; downcorrelated = "down"; upcorrelated = "up";}
//       else if(flav.at(ijet) == 4){wpcorr = btagwpbccorr; eff = btageffs[ptbin][1];}
//       else{wpcorr = btagwpbccorr; eff = btageffs[ptbin][0];}

//       if(disc.at(ijet) > deepjetL){
// 	weights[0] *= wpcorr->evaluate({"central","L",flav.at(ijet),abs(eta.at(ijet)), pt.at(ijet)}); // seemed like P(Data)/P(MC) reduces to SF
// 	weights[1+shift[0]] *= wpcorr->evaluate({upcorrelated,"L",flav.at(ijet),abs(eta.at(ijet)), pt.at(ijet)}); 
// 	weights[2+shift[0]] *= wpcorr->evaluate({downcorrelated,"L",flav.at(ijet),abs(eta.at(ijet)), pt.at(ijet)});
// 	if(year == "2023" or year == "2023BPix"){
// 	  weights[3+shift[0]] *= wpcorr->evaluate({"up_uncorrelated","L",flav.at(ijet),abs(eta.at(ijet)), pt.at(ijet)});
// 	  weights[4+shift[0]] *= wpcorr->evaluate({"down_uncorrelated","L",flav.at(ijet),abs(eta.at(ijet)), pt.at(ijet)});
// 	}
// 	weights[1+shift[1]] *= wpcorr->evaluate({"central","L",flav.at(ijet),abs(eta.at(ijet)), pt.at(ijet)});
// 	weights[2+shift[1]] *= wpcorr->evaluate({"central","L",flav.at(ijet),abs(eta.at(ijet)), pt.at(ijet)});
// 	weights[3+shift[1]] *= wpcorr->evaluate({"central","L",flav.at(ijet),abs(eta.at(ijet)), pt.at(ijet)});
// 	weights[4+shift[1]] *= wpcorr->evaluate({"central","L",flav.at(ijet),abs(eta.at(ijet)), pt.at(ijet)});
//       }else{
// 	weights[0] *= (1.0-eff*wpcorr->evaluate({"central","L",flav.at(ijet),abs(eta.at(ijet)), pt.at(ijet)}))/(1.0-eff);
// 	weights[1+shift[0]] *= (1.0-eff*wpcorr->evaluate({upcorrelated,"L",flav.at(ijet),abs(eta.at(ijet)), pt.at(ijet)}))/(1.0-eff);
// 	weights[2+shift[0]] *= (1.0-eff*wpcorr->evaluate({downcorrelated,"L",flav.at(ijet),abs(eta.at(ijet)), pt.at(ijet)}))/(1.0-eff);
// 	if(year == "2023" or year == "2023BPix"){
// 	  weights[3+shift[0]] *= (1.0-eff*wpcorr->evaluate({"up_uncorrelated","L",flav.at(ijet),abs(eta.at(ijet)), pt.at(ijet)}))/(1.0-eff);
// 	  weights[4+shift[0]] *= (1.0-eff*wpcorr->evaluate({"down_uncorrelated","L",flav.at(ijet),abs(eta.at(ijet)), pt.at(ijet)}))/(1.0-eff);
// 	}
// 	weights[1+shift[1]] *= (1.0-eff*wpcorr->evaluate({"central","L",flav.at(ijet),abs(eta.at(ijet)), pt.at(ijet)}))/(1.0-eff);
// 	weights[2+shift[1]] *= (1.0-eff*wpcorr->evaluate({"central","L",flav.at(ijet),abs(eta.at(ijet)), pt.at(ijet)}))/(1.0-eff);
// 	weights[3+shift[1]] *= (1.0-eff*wpcorr->evaluate({"central","L",flav.at(ijet),abs(eta.at(ijet)), pt.at(ijet)}))/(1.0-eff);
// 	weights[4+shift[1]] *= (1.0-eff*wpcorr->evaluate({"central","L",flav.at(ijet),abs(eta.at(ijet)), pt.at(ijet)}))/(1.0-eff);
//       }
//     }
//     return weights;      
//   };
