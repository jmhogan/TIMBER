import os, sys, math, re
from ROOT import TFile, TTree, TH1D, TH2D, TCanvas, gStyle, gPad, TLegend, kBlue, kRed, kGreen, TLatex, kRainBow
from array import array

callAlgoEff = False
callMassEff = True
dir_str = "root://cmseos.fnal.gov//store/user/lpchtop/TTBB_Jun2026_Run3/"
sample_files = ["RDF_TprimeTprime_Par-M-1200_TuneCP5_13p6TeV_amcatnlo-pythia8_2024_0.root","RDF_TprimeTprime_Par-M-1300_TuneCP5_13p6TeV_amcatnlo-pythia8_2024_0.root","RDF_TprimeTprime_Par-M-1400_TuneCP5_13p6TeV_amcatnlo-pythia8_2024_0.root","RDF_TprimeTprime_Par-M-1600_TuneCP5_13p6TeV_amcatnlo-pythia8_2024_0.root","RDF_TprimeTprime_Par-M-1700_TuneCP5_13p6TeV_amcatnlo-pythia8_2024_0.root","RDF_TprimeTprime_Par-M-1800_TuneCP5_13p6TeV_amcatnlo-pythia8_2024_0.root","RDF_TprimeTprime_Par-M-1900_TuneCP5_13p6TeV_amcatnlo-pythia8_2024_0.root","RDF_TprimeTprime_Par-M-2000_TuneCP5_13p6TeV_amcatnlo-pythia8_2024_0.root","RDF_TprimeTprime_Par-M-2100_TuneCP5_13p6TeV_amcatnlo-pythia8_2024_0.root","RDF_TprimeTprime_Par-M-2200_TuneCP5_13p6TeV_amcatnlo-pythia8_2024_0.root"]

test_files = ["RDF_TprimeTprime_Par-M-1200_TuneCP5_13p6TeV_amcatnlo-pythia8_2024_0.root","RDF_TprimeTprime_Par-M-1300_TuneCP5_13p6TeV_amcatnlo-pythia8_2024_0.root"]

files = [dir_str + x for x in test_files]

ptbins = array('d',[200,300,400,480,560,650,740,1000,1280,1640,2000])
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
        if matched == "matched":
            isMatched = True
        else:
            isMatched = False
            
        for part in {"t_quark","W_boson","Z_boson","H_boson"}:
            if part == "t_quark":
                ID = 6
            elif part == "W_boson":
                ID = 24
            elif part == "Z_boson":
                ID = 23
            else:
                ID = 25

            xaxis_str = f" jet pt (GeV)"
            yaxis_str = f" {part} "
            yaxis_str_d = f"tagged {part}"
            mistagged = "tot"
            if isMatched:
                xaxis_str = part + xaxis_str
                yaxis_str = " "
                yaxis_str_d = matched
                mistagged = ""
                
            truth = TH1D(f"truth_{part}_{matched}", f";{xaxis_str};",nbins,ptbins)
            PNWM = TH1D(f"{part}_{matched}_{tot}PNWM",f";{xaxis_str};N(PNWM{part}{matched} jets)/N({yaxis_str_d} jets)",nbins,ptbins)
            GPT = TH1D(f"{part}_{matched}_{tot}GPT",f";{xaxis_str};N(GPT{part}{matched} jets)/N({yaxis_str_d} jets)",nbins,ptbins)
            GPTWM = TH1D(f"{part}_{matched}_{tot}GPTWM",f";{xaxis_str};N(GPTWM{part}{matched} jets)/N({yaxis_str_d} jets)",nbins,ptbins)
            if not isMatched:
                misPNWM = TH1D(f"{part}_{matched}_mistaggedPNWM",f";{xaxis_str};N(PNWM{part}{matched} jets)/N({yaxis_str_d} jets)",nbins,ptbins)
                misGPT = TH1D(f"{part}_{matched}_mistaggedGPT",f";{xaxis_str};N(GPT{part}{matched} jets)/N({yaxis_str_d} jets)",nbins,ptbins)
                misGPTWM = TH1D(f"{part}_{matched}_mistaggedGPTWM",f";{xaxis_str};N(GPTWM{part}{matched} jets)/N({yaxis_str_d} jets)",nbins,ptbins)                    
                
            truth.Sumw2() #Sum weights squared -> account for uncertainties. in future we divide this
            PNWM.Sumw2()
            GPT.Sumw2()
            GPTWM.Sumw2()
            if not isMatched:
                misPNWM.Sumw2()
                misGPT.Sumw2()
                misGPTWM.Sumw2()
            
            for ievent in range(t.GetEntries()):
                t.GetEntry(ievent)
                nJets = t.nFatJet
                
                for ijet in range(nJets):
                    if t.gcFatJet_truth[ijet] == ID:
                        truth.Fill(t.gcFatJet_pt[ijet])
            
                        if t.gcFatJet_PNWMtags[ijet] == ID:
                            PNWM.Fill(t.gcFatJet_pt[ijet])

                        if t.gcFatJet_GPTtags[ijet] == ID:
                            GPT.Fill(t.gcFatJet_pt[ijet])

                        if t.gcFatJet_GPTWMtags[ijet] == ID:
                            GPTWM.Fill(t.gcFatJet_pt[ijet])

                    if not isMatched:
                        if t.gcFatJet_PNWMtags[ijet] == ID and t.gcFatJet_truth[ijet] != ID:
                            PNWM.Fill(t.gcFatJet_pt[ijet])
                            misPNWM.Fill(t.gcFatJet_pt[ijet])
                            
                        if t.gcFatJet_GPTtags[ijet] != ID:
                            GPT.Fill(t.gcFatJet_pt[ijet])
                            misGPT.Fill(t.gcFatJet_pt[ijet])

                        if t.gcFatJet_GPTWMtags[ijet] != ID:
                            GPTWM.Fill(t.gcFatJet_pt[ijet])
                            misGPTWM.Fill(t.gcFatJet_pt[ijet])
                        

            if isMatched:                
                PNWM.Divide(PNWM, truth, 1, 1, "B")
                GPT.Divide(GPT, truth, 1, 1, "B")
                GPTWM.Divide(GPTWM, truth, 1, 1, "B")

                PNWM.SetMarkerStyle(20)
                GPT.SetMarkerStyle(20)
                GPTWM.SetMarkerStyle(20)
            
                PNWM.SetMarkerColor(kBlue)
                GPT.SetMarkerColor(kRed)
                GPTWM.SetMarkerColor(kGreen)
            
                PNWM.SetLineColor(kBlue)
                GPT.SetLineColor(kRed)
                GPTWM.SetLineColor(kGreen)
            else:
                misPNWM.Divide(misPNWM,PNWM,1,1,"B")
                misGPT.Divide(misGPT,GPT,1,1,"B")
                misGPTWM.Divide(misGPTWM,GPTWM,1,1,"B")

                misPNWM.SetMarkerStyle(20)
                misGPT.SetMarkerStyle(20)
                misGPTWM.SetMarkerStyle(20)
            
                misPNWM.SetMarkerColor(kBlue)
                misGPT.SetMarkerColor(kRed)
                misGPTWM.SetMarkerColor(kGreen)
            
                misPNWM.SetLineColor(kBlue)
                misGPT.SetLineColor(kRed)
                misGPTWM.SetLineColor(kGreen)            
            
            canvas_name = f"canv_{part}_{matched}"            
            print(f"create {canvas_name} below")
            canv1 = TCanvas(canvas_name,part,800,600)
            
            gStyle.SetOptStat(0)
            canv1.SetLeftMargin(0.15)

            if isMatched:
                PNWM.Draw("pe")    
                GPT.Draw("pe same")
                GPTWM.Draw("pe same")
            else:
                misPNWM.Draw("pe")    
                misGPT.Draw("pe same")
                misGPTWM.Draw("pe same")

            png_name = f"taggerEff/{part}_{matched}_{mass_name}_taggerEff.png"
            canv1.BuildLegend()
            canv1.SaveAs(png_name)

            print(f"taggerEff/{mass_name}_taggerEff.root")

            if isMatched:
                histFile = TFile.Open(f"taggerEff/{mass_name}_taggerEff.root", "recreate")
                PNWM.Write()
                GPT.Write()
                GPTWM.Write()
                histFile.Write()
                histFile.Close()
            else:
                histFile = TFile.Open(f"taggerEff/{mass_name}_taggerEff.root", "update")
                misPNWM.Write()
                misGPT.Write()
                misGPTWM.Write()
                histFile.Write()
                histFile.Close()

                misPNWM.Reset()
                misGPT.Reset()
                misGPTWM.Reset()
                
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
        
        for matched in {"tagged","mistagged"}:
            if matched == "tagged":
                isMatched = True
            else:
                isMatched = False
                
            for part in {"t_quark","W_boson","Z_boson","H_boson"}:
                if part == "t_quark":
                    ID = 6
                elif part == "W_boson":
                    ID = 24
                elif part == "Z_boson":
                    ID = 23
                else:
                    ID = 25

                
                if not isMatched:
                    yaxis_str_d = f"{part} tagged"
                    tot = "tot"
                else:
                    yaxis_str_d = part
                    tot = ""
                    
                truth = TH1D(f"truth_{part}_{matched}", f";x;",nbins,ptbins)
                PNWM = TH1D(f"{part}_{matched}_{tot}PNWM_{mass_name}",f";Jet pt (GeV);N({part} {matched} PNWM jets)/N({yaxis_str_d} jets)",nbins,ptbins)
                GPT = TH1D(f"{part}_{matched}_{tot}GPT_{mass_name}",f";Jet pt (GeV);N({part} {matched} GPT jets)/N({yaxis_str_d} jets)",nbins,ptbins)
                GPTWM = TH1D(f"{part}_{matched}_{tot}GPTWM_{mass_name}",f";Jet pt (GeV);N({part} {matched} GPTWM jets)/N({yaxis_str_d} jets)",nbins,ptbins)

                if not isMatched:
                    misPNWM = TH1D(f"{part}_{matched}_PNWM_{mass_name}",f";Jet pt (GeV);N({part} {matched} PNWM jets)/N({yaxis_str_d} jets)",nbins,ptbins)
                    misGPT = TH1D(f"{part}_{matched}_GPT_{mass_name}",f";Jet pt (GeV);N({part} {matched} GPT jets)/N({yaxis_str_d} jets)",nbins,ptbins)
                    misGPTWM = TH1D(f"{part}_{matched}_GPTWM_{mass_name}",f";Jet pt (GeV);N({part} {matched} GPTWM jets)/N({yaxis_str_d} jets)",nbins,ptbins)                    
                    
                    misPNWM.Sumw2()
                    misGPT.Sumw2()
                    misGPTWM.Sumw2()
                    
                truth.Sumw2() #Sum weights squared -> account for uncertainties. in future we divide this
                PNWM.Sumw2()
                GPT.Sumw2()
                GPTWM.Sumw2()
    
                for ievent in range(t.GetEntries()):
                    t.GetEntry(ievent)
                    nJets = t.nFatJet
                
                    for ijet in range(nJets):
                        if t.gcFatJet_truth[ijet] == ID:
                            truth.Fill(t.gcFatJet_pt[ijet])
            
                            if t.gcFatJet_PNWMtags[ijet] == ID:
                                PNWM.Fill(t.gcFatJet_pt[ijet])

                            if t.gcFatJet_GPTtags[ijet] == ID:
                                GPT.Fill(t.gcFatJet_pt[ijet])

                            if t.gcFatJet_GPTWMtags[ijet] == ID:
                                GPTWM.Fill(t.gcFatJet_pt[ijet])

                    if not isMatched:
                        if t.gcFatJet_PNWMtags[ijet] == ID and t.gcFatJet_truth[ijet] != ID:
                            PNWM.Fill(t.gcFatJet_pt[ijet])
                            misPNWM.Fill(t.gcFatJet_pt[ijet])
                            
                        if t.gcFatJet_GPTtags[ijet] != ID:
                            GPT.Fill(t.gcFatJet_pt[ijet])
                            misGPT.Fill(t.gcFatJet_pt[ijet])

                        if t.gcFatJet_GPTWMtags[ijet] != ID:
                            GPTWM.Fill(t.gcFatJet_pt[ijet])
                            misGPTWM.Fill(t.gcFatJet_pt[ijet])

                if isMatched:
                    print(GPTWM.GetBinContent(1))

                    PNWM.Divide(PNWM, truth, 1, 1, "B")
                    GPT.Divide(GPT, truth, 1, 1, "B")
                    GPTWM.Divide(GPTWM, truth, 1, 1, "B")
                else:
                    misPNWM.Divide(misPNWM,PNWM,1,1,"B")
                    misGPT.Divide(misGPT,GPT,1,1,"B")
                    misGPTWM.Divide(misGPTWM,GPTWM,1,1,"B")
                    
                if i == 0:
                    mode = "recreate"
                else:
                    mode = "update"
                    
                histFile = TFile.Open(f"massEff/{part}_{matched}_PNWM_massEff.root", mode)
                if isMatched:
                    PNWM.Write()
                else:
                    misPNWM.Write()
                    misPNWM.Reset()
                PNWM.Reset()
                histFile.Write()
                histFile.Close()
                
                histFile = TFile.Open(f"massEff/{part}_{matched}_GPT_massEff.root", mode)
                if isMatched:
                    GPT.Write()
                else:
                    misGPT.Write()
                    misGPT.Reset()
                GPT.Reset()
                histFile.Write()
                histFile.Close()
                                
                histFile = TFile.Open(f"massEff/{part}_{matched}_GPTWM_massEff.root", mode)
                if isMatched:
                    GPTWM.Write()
                else:
                    misGPTWM.Write()
                    misGPTWM.Reset()
                GPTWM.Reset()
                histFile.Write()
                histFile.Close()
                
                truth.Reset()
                
        mass_branches.append(mass_name)
    
    for matched in {"tagged","mistagged"}:
        for part in {"t_quark","W_boson","Z_boson","H_boson"}:
            for tagger in {"PNWM","GPT","GPTWM"}:
                canvas_name = f"canv_{part}_{matched}_{tagger}"
                canv1 = TCanvas(canvas_name,part,800,600)
                gStyle.SetOptStat(0)
                canv1.SetLeftMargin(0.15)
                gStyle.SetPalette(kRainBow)
                if matched == "tagged":
                    legend = TLegend(0.25,0.15,0.38,0.33)
                    legend.SetTextFont(43)
                    legend.SetTextSize(24)
                    legend.SetHeader(f"{part} {matched} {tagger}","R")
                else: #mistagged
                    legend = TLegend(0.65,0.65,0.80,0.80)
                    legend.SetTextFont(43)
                    legend.SetTextSize(24)
                    legend.SetHeader(f"{part} {matched} {tagger}","C")
                
                print("opening tagger file, branches below")
                tagger_file = TFile.Open(f"massEff/{part}_{matched}_{tagger}_massEff.root")
                tagger_file.ls()
                
                print(f"mass_branches: {mass_branches}")
                for i,mass_name in enumerate(mass_branches):
                    print(f"getting {part}_{matched}_{tagger}_{mass_name}")
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

                    print(f"legend name from {mass_name} is ")
                    
                    legend.AddEntry(mass_hist, f"{mass_name[0]}{mass_name[0]} {mass_name[2]}.{mass_name[4]}TeV")
                    #mass_hist.Reset()
            
                png_name = f"massEff/{part}_{matched}_{tagger}_massEff.png"
                #canv1.BuildLegend(title=, option=)
                legend.SetBorderSize(0)
                legend.SetFillStyle(0)
                legend.Draw()
                canv1.SaveAs(png_name)
                tagger_file.Close()
                canv1.Close()
                legend.Clear()
                del canv1

if callAlgoEff:
    for f in files:
        algoEff(f,ptbins,nbins)
if callMassEff:
    massEff(files,ptbins,nbins)
