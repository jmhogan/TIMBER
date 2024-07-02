// Methods in this file:
// genttbarMassCalc()		//TODO BtoTW generatorInfo.cc has a lot more funcs

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
