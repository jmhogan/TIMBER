import os, sys, math, re
from ROOT import TFile, TTree, TH1D, TH2D, TCanvas, gStyle, gPad, TLegend, kBlue, kRed, kGreen, TLatex, kRainBow
from array import array

callAlgoEff = False
callMassEff = True
dir_str = "root://cmseos.fnal.gov//store/user/lpchtop/TTBB_Jun2026_Run3/"
sample_files = ["RDF_TprimeTprime_Par-M-1200_TuneCP5_13p6TeV_amcatnlo-pythia8_2024_0.root","RDF_TprimeTprime_Par-M-1300_TuneCP5_13p6TeV_amcatnlo-pythia8_2024_0.root","RDF_TprimeTprime_Par-M-1400_TuneCP5_13p6TeV_amcatnlo-pythia8_2024_0.root","RDF_TprimeTprime_Par-M-1600_TuneCP5_13p6TeV_amcatnlo-pythia8_2024_0.root","RDF_TprimeTprime_Par-M-1700_TuneCP5_13p6TeV_amcatnlo-pythia8_2024_0.root","RDF_TprimeTprime_Par-M-1800_TuneCP5_13p6TeV_amcatnlo-pythia8_2024_0.root","RDF_TprimeTprime_Par-M-1900_TuneCP5_13p6TeV_amcatnlo-pythia8_2024_0.root","RDF_TprimeTprime_Par-M-2000_TuneCP5_13p6TeV_amcatnlo-pythia8_2024_0.root","RDF_TprimeTprime_Par-M-2100_TuneCP5_13p6TeV_amcatnlo-pythia8_2024_0.root","RDF_TprimeTprime_Par-M-2200_TuneCP5_13p6TeV_amcatnlo-pythia8_2024_0.root"]

test_files = ["RDF_TprimeTprime_Par-M-1200_TuneCP5_13p6TeV_amcatnlo-pythia8_2024_0.root","RDF_TprimeTprime_Par-M-1300_TuneCP5_13p6TeV_amcatnlo-pythia8_2024_0.root"]

files = [dir_str + x for x in sample_files]

ptbins = array('d',[200,380,470,560,650,740,1000,1280,1640,2000])
nbins = len(ptbins) - 1
print(nbins)

def algoEff(file_str,ptbins,nbins):
    print(f"opening {file_str}")
    inFile = TFile.Open(file_str)
    t = inFile.Get("Events_Nominal")

    pattern = r"RDF_([TB]).*?Par-M-(\d+)"
    mass = re.search(pattern, file_str)
    mass_name = f"{mass.group(1)}_{mass.group(2)[0]}p{mass.group(2)[1]}TeV"
    
    print(mass_name)
    
    for matched in {"matched","mismatched"}:
        for part in {"t_quark","W_boson","Z_boson","H_boson"}:
            if part == "t_quark":
                ID = 6
            elif part == "W_boson":
                ID = 24
            elif part == "Z_boson":
                ID = 23
            else:
                ID = 25

            axis_str = ";" + part + "jet pt GeV;N(tagged " + matched + " jets)/N(" + matched + " jets)"
            truth = TH1D(f"truth_{part}_{matched}",axis_str,nbins,ptbins)
            PNWM = TH1D(f"{part}_{matched}_PNWM",axis_str,nbins,ptbins)
            GPT = TH1D(f"{part}_{matched}_GPT",axis_str,nbins,ptbins)
            GPTWM = TH1D(f"{part}_{matched}_GPTWM",axis_str,nbins,ptbins)
            
            truth.Sumw2() #Sum weights squared -> account for uncertainties. in future we divide this
            PNWM.Sumw2()
            GPT.Sumw2()
            GPTWM.Sumw2()
    
            for ievent in range(t.GetEntries()):
                t.GetEntry(ievent)
                nJets = t.nFatJet
                
                for ijet in range(nJets):
                    if t.gcFatJet_truth[ijet] == ID and matched == "matched":
                        truth.Fill(t.gcFatJet_pt[ijet])
            
                        if t.gcFatJet_PNWMtags[ijet] == ID:
                            PNWM.Fill(t.gcFatJet_pt[ijet])

                        if t.gcFatJet_GPTtags[ijet] == ID:
                            GPT.Fill(t.gcFatJet_pt[ijet])

                        if t.gcFatJet_GPTWMtags[ijet] == ID:
                            GPTWM.Fill(t.gcFatJet_pt[ijet])

                    elif t.gcFatJet_truth[ijet] == ID and matched == "mismatched":
                        truth.Fill(t.gcFatJet_pt[ijet])
            
                        if t.gcFatJet_PNWMtags[ijet] != ID:
                            PNWM.Fill(t.gcFatJet_pt[ijet])

                        if t.gcFatJet_GPTtags[ijet] != ID:
                            GPT.Fill(t.gcFatJet_pt[ijet])

                        if t.gcFatJet_GPTWMtags[ijet] != ID:
                            GPTWM.Fill(t.gcFatJet_pt[ijet])
                        
            PNWM.Divide(PNWM, truth, 1, 1, "B")
            GPT.Divide(GPT, truth, 1, 1, "B")
            GPTWM.Divide(GPTWM, truth, 1, 1, "B")

            print("divided done")
            
            PNWM.SetMarkerStyle(20)
            GPT.SetMarkerStyle(20)
            GPTWM.SetMarkerStyle(20)
            
            PNWM.SetMarkerColor(kBlue)
            GPT.SetMarkerColor(kRed)
            GPTWM.SetMarkerColor(kGreen)
            
            PNWM.SetLineColor(kBlue)
            GPT.SetLineColor(kRed)
            GPTWM.SetLineColor(kGreen)
            
            #text = "M(T) = 1.7TeV"
            #tex = TLatex()
            #tex.DrawLatex(0.2, 0.8, text)

            canvas_name = f"canv_{part}_{matched}"            
            print(f"create {canvas_name} below")
            canv1 = TCanvas(canvas_name,part,800,600)
            
            gStyle.SetOptStat(0)
            canv1.SetLeftMargin(0.15)
            
            PNWM.Draw("pe")    
            GPT.Draw("pe same")
            GPTWM.Draw("pe same")
            
            png_name = f"taggerEff/{part}_{matched}_{mass_name}_taggerEff.png"
            canv1.BuildLegend()
            canv1.SaveAs(png_name)

            print(f"taggerEff/{mass_name}_taggerEff.root")
            
            histFile = TFile.Open(f"taggerEff/{mass_name}_taggerEff.root", "recreate")
            GPT.Write()
            histFile.Write()
            histFile.Close()

            
            truth.Reset()
            PNWM.Reset()
            GPT.Reset()
            GPTWM.Reset()
            canv1.Close()
            del canv1


def massEff(file_strs,ptbins,nbins):
    mass_branches = []
    for i,file_str in enumerate(file_strs): #RDF_TprimeTprime_Par-M-1700_TuneCP5_13p6TeV_amcatnlo-pythia8_2024_0.root
        print(f"opening {file_str}")
        pattern = r"RDF_([TB]).*?Par-M-(\d+)"
        mass = re.search(pattern, file_str)
        mass_name = f"{mass.group(1)}_{mass.group(2)[0]}p{mass.group(2)[1]}TeV"
        print(f"processing {mass_name}")
        inFile = TFile.Open(file_str)
        t = inFile.Get("Events_Nominal")
        
        for matched in {"matched","mismatched"}:
            for part in {"t_quark","W_boson","Z_boson","H_boson"}:
                if part == "t_quark":
                    ID = 6
                elif part == "W_boson":
                    ID = 24
                elif part == "Z_boson":
                    ID = 23
                else:
                    ID = 25

                axis_str = ";" + part + "jet pt GeV;N(tagged " + matched + " jets)/N(" + matched + " jets)"
                truth = TH1D(f"truth_{part}_{matched}_{mass_name}",axis_str,nbins,ptbins)
                PNWM = TH1D(f"{part}_{matched}_PNWM_{mass_name}", axis_str,nbins,ptbins)
                GPT = TH1D(f"{part}_{matched}_GPT_{mass_name}", axis_str,nbins,ptbins)
                GPTWM = TH1D(f"{part}_{matched}_GPTWM_{mass_name}", axis_str,nbins,ptbins)
                
                truth.Sumw2() #Sum weights squared -> account for uncertainties. in future we divide this
                PNWM.Sumw2()
                GPT.Sumw2()
                GPTWM.Sumw2()
    
                for ievent in range(t.GetEntries()):
                    t.GetEntry(ievent)
                    nJets = t.nFatJet
                
                    for ijet in range(nJets):
                        if t.gcFatJet_truth[ijet] == ID and matched == "matched":
                            truth.Fill(t.gcFatJet_pt[ijet])
            
                            if t.gcFatJet_PNWMtags[ijet] == ID:
                                PNWM.Fill(t.gcFatJet_pt[ijet])

                            if t.gcFatJet_GPTtags[ijet] == ID:
                                GPT.Fill(t.gcFatJet_pt[ijet])

                            if t.gcFatJet_GPTWMtags[ijet] == ID:
                                GPTWM.Fill(t.gcFatJet_pt[ijet])

                        elif t.gcFatJet_truth[ijet] == ID and matched == "mismatched":
                            truth.Fill(t.gcFatJet_pt[ijet])
            
                            if t.gcFatJet_PNWMtags[ijet] != ID:
                                PNWM.Fill(t.gcFatJet_pt[ijet])

                            if t.gcFatJet_GPTtags[ijet] != ID:
                                GPT.Fill(t.gcFatJet_pt[ijet])

                            if t.gcFatJet_GPTWMtags[ijet] != ID:
                                GPTWM.Fill(t.gcFatJet_pt[ijet])
                        
                PNWM.Divide(PNWM, truth, 1, 1, "B")
                GPT.Divide(GPT, truth, 1, 1, "B")
                GPTWM.Divide(GPTWM, truth, 1, 1, "B")

                if i == 0:
                    mode = "recreate"
                else:
                    mode = "update"
                    
                histFile = TFile.Open(f"massEff/{part}_{matched}_PNWM_massEff.root", mode)
                PNWM.Write()
                histFile.Write()
                histFile.Close()
                PNWM.Reset()
                
                histFile = TFile.Open(f"massEff/{part}_{matched}_GPT_massEff.root", mode)
                GPT.Write()
                histFile.Write()
                histFile.Close()
                GPT.Reset()
                
                histFile = TFile.Open(f"massEff/{part}_{matched}_GPTWM_massEff.root", mode)
                GPTWM.Write()
                histFile.Write()
                histFile.Close()
                GPTWM.Reset()

                truth.Reset()
                
        mass_branches.append(mass_name)
    
    for matched in {"matched","mismatched"}:
        for part in {"t_quark","W_boson","Z_boson","H_boson"}:
            for tagger in {"PNWM","GPT","GPTWM"}:
                canvas_name = f"canv_{part}_{matched}_{tagger}"
                canv1 = TCanvas(canvas_name,part,800,600)
                gStyle.SetOptStat(0)
                canv1.SetLeftMargin(0.15)
                gStyle.SetPalette(kRainBow)

                print("opening tagger file, branches below")
                tagger_file = TFile.Open(f"massEff/{part}_{matched}_{tagger}_massEff.root")
                tagger_file.ls()
                print(f"mass_branches: {mass_branches}")
                for i,mass_name in enumerate(mass_branches):
                    print(f"{part}_{matched}_{tagger}_{mass_name}")
                    mass_hist = tagger_file.Get(f"{part}_{matched}_{tagger}_{mass_name}")
                    color_index = int(i * (gStyle.GetNumberOfColors() / len(mass_branches)))
                    color = gStyle.GetColorPalette(color_index)
                    
                    mass_hist.SetMarkerStyle(20)
                    mass_hist.SetMarkerColor(color)
                    mass_hist.SetLineColor(color)
                    mass_hist.GetYaxis().SetRangeUser(0,1)
                    
                    if i == 0:
                        mass_hist.Draw("pe")    
                    else:
                        mass_hist.Draw("pe same")

                    #mass_hist.Reset()
            
                png_name = f"massEff/{part}_{matched}_{tagger}_massEff.png"
                canv1.BuildLegend()
                canv1.SaveAs(png_name)
                tagger_file.Close()
                canv1.Close()
                del canv1

            

if callAlgoEff:
    for f in files:
        algoEff(f,ptbins,nbins)
if callMassEff:
    massEff(files,ptbins,nbins)
