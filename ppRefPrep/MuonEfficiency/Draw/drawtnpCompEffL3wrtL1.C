#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>

#include <TH1.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TF1.h>
#include <TTree.h>
#include <TTreeCache.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TTimeStamp.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLine.h>

#define DEBUG (0)
#include "PlotTools.h"

using namespace std;

void printRunTime(TStopwatch timer_)
{
  Double_t cpuTime = timer_.CpuTime();
  Double_t realTime = timer_.RealTime();

  cout << endl;
  cout << "************************************************" << endl;
  cout << "Total real time: " << realTime << " (seconds)" << endl;
  cout << "Total CPU time:  " << cpuTime << " (seconds)" << endl;
  cout << "  CPU time / real time = " << cpuTime / realTime << endl;
  cout << "************************************************" << endl;
}

void drawtnpCompEffL3wrtL1(
  // TString efftag0 = "PPRefL1SingleMu0_Cosmics", TString efftag1 = "PPRefL1SingleMu7", TString efftag2 = "PPRefL1SingleMu12",
  TString efftag0 = "PPRefL2SingleMu7", TString efftag1 = "PPRefL2SingleMu12", TString efftag2 = "PPRefL2SingleMu15", TString efftag3 = "PPRefL2SingleMu20",
  // TString efftag0 = "PPRefL3SingleMu3", TString efftag1 = "PPRefL3SingleMu5", TString efftag2 = "PPRefL3SingleMu7",
  // TString efftag0 = "PPRefL3SingleMu12", TString efftag1 = "PPRefL3SingleMu15", TString efftag2 = "PPRefL3SingleMu20",
  // TString efftag0 = "PPRefL1DoubleMu0_Open",  TString efftag1 = "PPRefL1DoubleMu0",
  // TString efftag0 = "PPRefL2DoubleMu0_Open", TString efftag1 = "PPRefL2DoubleMu0",
  // TString efftag0 = "PPRefL3DoubleMu0_Open", TString efftag1 = "PPRefL3DoubleMu0",
  // TString ver = "vRun3_JPsi_PU4_CMSSW_14_0_13", TString SAMPLE = "Run2024 MC prep", TString tag = "Muon",
  TString ver = "vRun3_Zmm_PU4_CMSSW_14_0_13", TString SAMPLE = "Run2024 MC prep", TString tag = "Muon",
  TString L1tag = "L1DQ2", TString L1str = "L1 qual > 7, p_{T}^{L1} > 2 GeV",
  //TString L1tag = "L1DQ8", TString L1str = "L1 qual > 7, p_{T}^{L1} > 8 GeV",
  bool isLogy = false  // HERE
) {
  TStopwatch timer_total;
  timer_total.Start();

  gStyle->SetPalette(kRainBow);
  TH1::SetDefaultSumw2(kTRUE);

  TString Dir = "./plots_"+ver+"/"+tag+"/Eff_"+efftag1+"/"+L1tag+"/";
  if (gSystem->mkdir(Dir,kTRUE) != -1)
    gSystem->mkdir(Dir,kTRUE);

  vector<TString> v_var = {"pt_zoom", "pt", "eta", "phi", "nvtx"};//, "pu"};
  vector< vector<double> > range = {
    {1, 10, 100},  // pt zoom
    {1, 0, 100},  // pt
    {1, -2.4, 2.4},  // eta
    {1, -TMath::Pi(), TMath::Pi()},
    {1, 1, 10}  // nvtx
    // {1, 10, 75}  // PU

  };
  if (tag == "JPsi" || tag == "Bs") {
    range.at(0) = {1, 0, 40};
    range.at(1) = {1, 0, 40};
  }

  if (tag == "Zmm" || tag == "Bs") {
    range.at(0) = {1, 0, 40};
    range.at(1) = {1, 0, 40};
  }
  /// for offline pt cut 6
  // const int n_pt_bins = 24-1;
  // double pt_bins[n_pt_bins+1] = {
  //   0, 1, 2, 3, 4,
  //   5, 5.9, 6, 6.1, 6.3, 
  //   6.7, 7, 8, 9, 10, 
  //   11, 12, 13, 14, 15, 
  //   17, 21, 26, 30
  // };

  // /// for offline pt cut 8
  // const int n_pt_bins = 22-1;
  // double pt_bins[n_pt_bins+1] = {
  //   0, 1, 2, 3, 4,
  //   5, 6, 7, 8, 9, 
  //   10, 11, 12, 13, 14, 
  //   15, 17, 21, 26, 31, 
  //   36, 40
  // };

  /// for offline pt cut 8
  const int n_pt_bins = 24-1;
  double pt_bins[n_pt_bins+1] = {
    10, 11, 12, 13, 14, 
    15, 17, 20, 25, 30, 
    35, 40, 45, 50, 55,
    60, 65, 70, 75, 80,
    85, 90, 95, 100
  };

  // /// for offline pt cut 12
  // const int n_pt_bins = 11-1;
  // double pt_bins[n_pt_bins+1] = {
  //   10, 11, 12, 12.2, 12.8, 
  //   13.8, 16, 19, 22, 25, 
  //   30
  // };

  // /// for offline pt cut 15
  // const int n_pt_bins = 13-1;
  // double pt_bins[n_pt_bins+1] = {
  //    10, 11, 12, 13, 14, 
  //    15, 15.1, 15.6, 17, 
  //    20, 23, 26, 30
  // };

  // /// for offline pt cut 20
  // const int n_pt_bins = 14-1;
  // double pt_bins[n_pt_bins+1] = {
  //    10, 11, 12, 13, 14, 
  //    15, 16, 17, 20, 20.1, 
  //    20.6, 23, 26, 30
  // };

  // int n_eta_bins = 23-1;
  // double eta_bins[23] = {
  //   -2.4, -2.1, -1.9, -1.7, -1.6, -1.5, -1.3, -1.2, -0.9,
  //   -0.3, -0.2,  0.0,  0.2,  0.3,
  //    0.9,  1.2, 1.3, 1.5, 1.6, 1.7, 1.9, 2.1,  2.4
  // };

  // int n_eta_bins = 15-1;
  // double eta_bins[15] = {
  //   -2.4, -2.1, -1.6, -1.2, -0.9,
  //   -0.3, -0.2,  0.0,  0.2,  0.3,
  //    0.9,  1.2,  1.6,  2.1,  2.4
  // };
  int n_eta_bins = 22-1;
  double eta_bins[22] = {
    -2.4, -2.1, -1.9, -1.7, -1.6, 
    -1.3, -1.2, -0.9, -0.3, -0.2,  
     0.0,  0.2,  0.3,  0.9,  1.2,  
     1.3,  1.5,  1.6,  1.7,  1.9,
     2.1,  2.4
  };


  vector<TString> etas_str = {"I"};//, "BB", "BE", "EB", "EE"};
  vector<TString> etas_str_long = {"|#eta^{offline}| < 2.4"};//, "|#eta^{offline}| < 0.9", "0.9 < |#eta^{offline}| < 1.2", "1.2 < |#eta^{offline}| < 2.1", "2.1 < |#eta^{offline}| < 2.4"};

  vector<Color_t> v_color = {
    // kBlack,
    kRed,
    // kOrange,
    // kGreen+2,
    kCyan+2,
    // kPink+1,
    kBlue-2,
    kGray+2,
    // kMagenta,
  };
  vector<int> v_marker = {
    20,
    22,
    23,
    24,
    25,
    26,
    23,
    //22,
    //26,
    //23,
    //32,
  };
  // -- input file names
  vector<TString> files = {
    // "../Analyzer/hist-v00-TEST-Eff_ppRefMCJPsi_PU4_CMSSW14_0_13.root",
    // "../Analyzer/hist-v00-TEST-Eff_ppRefMCJPsi_PU4_CMSSW14_0_13.root",
    // "../Analyzer/hist-v00-TEST-Eff_ppRefMCJPsi_PU4_CMSSW14_0_13.root",
    // "../Analyzer/hist-v00-TEST-Eff_ppRefMCJPsi_PU4_CMSSW14_0_13.root",
    
    // "../Analyzer/hist-v00-TEST-Eff_ppRefMCJPsi_PU10_CMSSW14_0_13.root",
    // "../Analyzer/hist-v00-TEST-Eff_ppRefMCJPsi_PU10_CMSSW14_0_13.root",
    // "../Analyzer/hist-v00-TEST-Eff_ppRefMCJPsi_PU10_CMSSW14_0_13.root",
    // "../Analyzer/hist-v00-TEST-Eff_ppRefMCJPsi_PU10_CMSSW14_0_13.root",
    // "../Analyzer/hist-v00-TEST-Eff_ppRefMCJPsi_PU10_CMSSW14_0_13.root",
    // "../Analyzer/hist-v00-TEST-Eff_ppRefMCJPsi_PU10_CMSSW14_0_13.root",

    // "../Analyzer/hist-v00-TEST-Eff_ppRefMCZmm_PU4_CMSSS14_0_13.root"

    "../Analyzer/hist-v00-TEST-Eff_ppRefMCZmm_PU4_CMSSW14_0_13.root",
    "../Analyzer/hist-v00-TEST-Eff_ppRefMCZmm_PU4_CMSSW14_0_13.root",
    "../Analyzer/hist-v00-TEST-Eff_ppRefMCZmm_PU4_CMSSW14_0_13.root",
    "../Analyzer/hist-v00-TEST-Eff_ppRefMCZmm_PU4_CMSSW14_0_13.root",
  };

  vector<TString> types = {
    //TString("Eff/"+efftag+"/num_Eff_"+L1tag+"_"+efftag+"_RunAll").ReplaceAll("my", ""),
    "Eff/"+efftag0+"/num_Eff_"+L1tag+"_"+efftag0+"_RunAll",
    "Eff/"+efftag1+"/num_Eff_"+L1tag+"_"+efftag1+"_RunAll",
    "Eff/"+efftag2+"/num_Eff_"+L1tag+"_"+efftag2+"_RunAll",
    "Eff/"+efftag3+"/num_Eff_"+L1tag+"_"+efftag3+"_RunAll",
    // "Eff/"+efftag4+"/num_Eff_"+L1tag+"_"+efftag4+"_RunAll",
    // "Eff/"+efftag5+"/num_Eff_"+L1tag+"_"+efftag5+"_RunAll",
    // "Eff/"+efftag+"/num_Eff_"+L1tag+"_"+efftag+"_RunAll",
    // "Eff/"+efftag+"/num_Eff_"+L1tag+"_"+efftag+"_RunAll",
  };
  vector<TString> types_den = {
    //TString("Eff/"+efftag+"/den_Eff_"+L1tag+"_"+efftag+"_RunAll").ReplaceAll("my", ""),
    "Eff/"+efftag0+"/den_Eff_"+L1tag+"_"+efftag0+"_RunAll",
    "Eff/"+efftag1+"/den_Eff_"+L1tag+"_"+efftag1+"_RunAll",
    "Eff/"+efftag2+"/den_Eff_"+L1tag+"_"+efftag2+"_RunAll",
    "Eff/"+efftag3+"/den_Eff_"+L1tag+"_"+efftag3+"_RunAll",
    // "Eff/"+efftag4+"/den_Eff_"+L1tag+"_"+efftag4+"_RunAll",
    // "Eff/"+efftag5+"/den_Eff_"+L1tag+"_"+efftag5+"_RunAll",
    // "Eff/"+efftag+"/den_Eff_"+L1tag+"_"+efftag+"_RunAll",
    // "Eff/"+efftag+"/den_Eff_"+L1tag+"_"+efftag+"_RunAll",
  };
  vector<TString> types_str = {
    // "Run2022G Data",
    // //"Run3Winer23 DY",
    // "Run2023B Data",
    // "Run2023C Data",
    // "Run2023D Data",
    // "Before fix",
    // "After fix"
    // "2024 ppRef MC hltL1sDoubleMuOpen  "
    // "2024 ppRef MC L3SingleMu12  ",
    // "2024 ppRef MC hltL3FromL2Merged",
    // "2024 ppRef MC hltL3Merged"
    // "2024 ppRef MC L1SingleMu7"
    // "2024 ppRef MC L1SingleMu12"
    // "2024 ppRef MC L2DoubleMu0"
    efftag0,
    efftag1,
    efftag2,
    efftag3,
    // efftag4,
    // efftag5,
    // Form("2024 ppRef MC %s", (char)efftag)
  };

  vector<TString> v_pts = {
    "genpt0",
    "genpt2",
    "genpt4",
    "genpt6",
    "genpt8",
    "genpt10",
    "genpt12",
    "genpt15",
    "genpt20",
    "genpt30",
    "genpt57",
    "genpt100",
  };

  vector<TString> v_pts_str = {
    "p_{T}^{offline} > 0 GeV",
    "p_{T}^{offline} > 2 GeV",
    "p_{T}^{offline} > 4 GeV",
    "p_{T}^{offline} > 6 GeV",
    "p_{T}^{offline} > 8 GeV",
    "p_{T}^{offline} > 10 GeV",
    "p_{T}^{offline} > 12 GeV",
    "p_{T}^{offline} > 15 GeV",
    "p_{T}^{offline} > 20 GeV",
    "p_{T}^{offline} > 30 GeV",
    "p_{T}^{offline} > 57 GeV",
    "p_{T}^{offline} > 100 GeV",
  };

  for(unsigned i_eta=0; i_eta<etas_str.size(); i_eta++){
    for(int ipt=0; ipt<(int)v_pts.size(); ++ipt) {
      for(int ivar=0; ivar<(int)v_var.size(); ++ivar) {

        double xmin = range[ivar][1];
        double xmax = range[ivar][2];
        double ymin = 0.0;
        double ymax = 1.6;

	if(!v_var[ivar].Contains("pt") || v_var[ivar] == "pt_zoom") {
          ymin = 0.;//0.82;//0.5;//0.6;//0.85;
          ymax = 1.6;//1.1;//1.25;//1.2;//1.1;
        }

        TString canvasName = TString::Format("Eff_%s_%s_%s_%s_%s_%s_%s_%s",
                                             efftag0.Data(),
                                             efftag1.Data(),
                                             efftag2.Data(),
                                             efftag3.Data(),
                                             // efftag4.Data(),
                                             // efftag5.Data(),
                                             L1tag.Data(),
                                             etas_str.at(i_eta).Data(),
                                             v_pts[ipt].Data(),
                                             v_var[ivar].Data());
        canvasName.ReplaceAll(".","p").ReplaceAll("-","_").ReplaceAll("my", "");
        TCanvas *c;
        SetCanvas_Square( c, canvasName, kFALSE, kFALSE, 900, 900 );
        c->cd();
        if(isLogy) c->SetLogy();
        if(tag == "Zprime" && v_var[ivar].Contains("pt")) c->SetLogx();

        TLegend *legend;
        SetLegend( legend, 0.17, 0.69, 0.94, 0.82, -1);

        bool isFirst = true;
        for(int i = 0; i<(int)files.size(); ++i) {
          TString fileName = files.at(i);
          cout << "fileName: "<< fileName << endl;

          TString the_type_num = types[i];
          TString the_type_den = types_den[i];
          TString the_type_str = types_str[i].ReplaceAll("my","");
          //cout<<the_type_num<<", " <<the_type_den<<endl;

          TString hist_var = v_var[ivar];
          hist_var.ReplaceAll("_zoom", "");

          TString titleX = GetTitleX(hist_var+"_reco");
          TString titleY = "HLT efficiency"; //"L3/L1 efficiency";
          if(efftag1.Contains("L2Muon")) titleY.ReplaceAll("L3", "L2");
          if(efftag1.Contains("PixelTracks")) titleY.ReplaceAll("L3", "PixelTrack");

          TString den_name = TString::Format("%s_%s_%s_%s", the_type_den.Data(), etas_str.at(i_eta).Data(), v_pts[ipt].Data(), hist_var.Data());
          TString num_name = TString::Format("%s_%s_%s_%s", the_type_num.Data(), etas_str.at(i_eta).Data(), v_pts[ipt].Data(), hist_var.Data());

          cout << "den_name: " << den_name << endl;
          cout << "num_name: " << num_name << endl;

          TH1F* den = Get_Hist( fileName, den_name );
          TH1F* num = Get_Hist( fileName, num_name );

          int nbins0 = den->GetNbinsX();
          // cout << nbins0<< endl;

          if(v_var[ivar] == "eta") {
            for(int ip=1; ip<=nbins0; ++ip) {
              cout << "ip, den: " << ip << ", "<< den -> GetBinContent(ip) << endl;
              cout << "ip, num: " << ip << ", " << num -> GetBinContent(ip) << endl;
            }
          }
          if(v_var[ivar] == "phi") {
            for(int ip=1; ip<=nbins0; ++ip) {
              cout << "ip, den: " << ip << ", "<< den -> GetBinContent(ip) << endl;
              cout << "ip, num: " << ip << ", " << num -> GetBinContent(ip) << endl;
            }
          }

          if(v_var[ivar].Contains("pt")) {
            den = (TH1F*)den->Rebin(n_pt_bins, den_name+"_rb", pt_bins);
            num = (TH1F*)num->Rebin(n_pt_bins, num_name+"_rb", pt_bins);
            for(int ip=1; ip<=nbins0; ++ip) {
              cout << "ip, den: " << ip << ", "<< den -> GetBinContent(ip) << endl;
              cout << "ip, num: " << ip << ", " << num -> GetBinContent(ip) << endl;
            }
          }
          else if(v_var[ivar] == "eta") {
            den = (TH1F*)den->Rebin(n_eta_bins, den_name+"_rb", eta_bins);
            num = (TH1F*)num->Rebin(n_eta_bins, num_name+"_rb", eta_bins);
          }
          else{
            den = (TH1F*)den->Rebin(range[ivar][0]);
            num = (TH1F*)num->Rebin(range[ivar][0]);
          }

          int nbins = den->GetNbinsX();

          if(v_var[ivar] == "eta") {
            for(int ip=0; ip<nbins; ++ip) {
              cout << "ip, den: " << ip << ", "<< den -> GetBinContent(ip+1) << endl;
              cout << "ip, num: " << ip << ", " << num -> GetBinContent(ip+1) << endl;
            }
          }

          c->cd();

          TGraphAsymmErrors* g = new TGraphAsymmErrors(nbins);
          g->Divide(num, den, "n e0");
          //g->Divide(num, den, "pois");

          for(int ip=0; ip<nbins; ++ip) {
            if(g->GetPointY(ip) == 0.)  g->SetPointEYhigh(ip, 0.0);
          }

          g->SetTitle("");
          g->SetMarkerSize(1.5);
          g->SetMarkerStyle(v_marker[i]);
          g->SetMarkerColor(v_color[i]);
          g->SetLineColor(  v_color[i]);
          g->SetLineWidth(1);

          g->GetXaxis()->SetLimits( xmin, xmax );
          g->GetXaxis()->SetRangeUser( xmin, xmax );
          g->GetYaxis()->SetRangeUser( ymin, ymax );

          SetAxis_SinglePad( g->GetXaxis(), g->GetYaxis(), titleX, titleY );

          if(isFirst) {
            g->Draw("APE");
            isFirst = false;
          }
          else {
            g->Draw("PE same");
          }

          legend->AddEntry( g, TString::Format("%s", the_type_str.Data()), "lep" );
        }

        legend->Draw();

        TString L3str = "";
        if(efftag1 == "L2Muon") L3str = "L2 Muon";
        else if(efftag1 == "hltOI") L3str = "Outside-in L3 MuonTrack";
        else if(efftag1 == "hltPixelTracksInRegionL2") L3str = "PixelTrack near L2";
        else if(efftag1 == "hltPixelTracksInRegionL1") L3str = "PixelTrack near L1";
        else if(efftag1 == "hltPixelTracks") L3str = "PixelTrack";
        else if(efftag1 == "hltIter0FromL1") L3str = "Inside-out L3 MuonTrack from L1";
        else if(efftag1 == "hltL3FromL2Merged") L3str = "L3 MuonTrack from L2";
        else if(efftag1 == "hltL3Merged") L3str = "L3 MuonTrack";
        else if(efftag1.Contains("hltIterL3MuonNoID")) L3str = "L3 Muon before Trigger ID";
        else if(efftag1 == "hltIterL3Muon") L3str = "L3 Muon after Trigger ID";
        
        else if(efftag1.Contains("PPRefL1SingleMu")) L3str = "PPRef L1 Single Mu";

        else if(efftag1.Contains("PPRefL2SingleMu")) L3str = "PPRef L2 Single Mu";

        else if(efftag1.Contains("PPRefL3SingleMu")) L3str = "PPRef L3 Single Mu";

        else if(efftag1.Contains("PPRefL1DoubleMu")) L3str = "PPRef L1 Double Mu";

        else L3str = efftag1;
        // else if(efftag.Contains("IsoMu24")) L3str = "Isolated muon with p_{T} > 24 GeV";
        // else if(efftag.Contains("Mu24")) L3str = "Non-isolated muon with p_{T} > 24 GeV";
        // else if(efftag.Contains("Mu50OrOldMu100OrTkMu100")) L3str = "Non-isolated muon with p_{T} > 50 GeV";
        // else if(efftag.Contains("Mu50")) L3str = "Non-isolated muon with p_{T} > 50 GeV";

        // // lines
        // TLine* horline = new TLine(0, 1, 20, 1);

        // horline->SetLineStyle(kDashed);
        // horline->SetLineColor(kGray);
        // horline->SetLineWidth(2);
        // horline->Draw("SAME");

        TLatex latex;
	      // Latex_Preliminary_13p6TeV( latex );
        Latex_Preliminary_5p36TeV( latex );
        latex.DrawLatexNDC(0.45,0.96, "#scale[0.8]{#font[42]{"+SAMPLE+"}}");
        latex.DrawLatexNDC(0.16, 0.90, "#font[42]{#scale[0.6]{"+L3str+"}}");
        latex.DrawLatexNDC(0.16, 0.85, "#font[42]{#scale[0.6]{"+L1str+"}}");
        latex.DrawLatexNDC((i_eta==2?0.66:0.70), 0.89, "#font[42]{#scale[0.8]{"+etas_str_long.at(i_eta)+"}}");
        if(v_var[ivar] != "" ) latex.DrawLatexNDC(0.68, 0.84, "#font[42]{#scale[0.8]{"+v_pts_str.at(ipt)+"}}");

        TString logy_tag = isLogy ? "_log" : "";
        // CMS_lumi(c, 98, 11);
        c->Modified();  c->Update();  c->RedrawAxis();
        gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;");
        //c->SaveAs(Dir+canvasName+logy_tag+".pdf","pdf");
        c->SaveAs(Dir+canvasName+logy_tag+".png","png");
        gROOT->ProcessLine( "gErrorIgnoreLevel = kPrint;");

        c->Close();
      }
    }
  }
  printRunTime(timer_total);
}
