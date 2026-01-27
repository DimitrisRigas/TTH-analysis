// overlay_signal_keybkgs.C
// Mode-2 only overlays: Signal vs ONE background with S/B ratio pads
//
// Run:
//   .L overlay_signal_keybkgs.C
//   overlay_signal_keybkgs(2,"tth12gev",true,true);  // tt dilepton, normalized, keep windows
//   overlay_signal_keybkgs(3,"tth12gev",true,true);  // tt semilep
//   overlay_signal_keybkgs(1,"tth12gev",true,true);  // ttH->bb
//
// Arguments:
//   bkgCode: 1=ttH->bb, 2=tt dilep (output_ttbar.root), 3=tt semilep (output_TTtoLNu2Q.root)
//   signalTag: e.g. "tth12gev"
//   normalize: true => unit area shapes; false => raw entries
//   keepWindowsOpen: true keeps canvases alive after macro ends

#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TString.h>

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>

// Keep canvases alive
static std::vector<TCanvas*> gCanvases;

static TH1* CloneToROOT(TFile* f, const std::string& keyPrefix, const char* hname) {
  static std::unordered_map<std::string, TH1*> cache;
  static long long ctr = 0;

  const std::string key = keyPrefix + "::" + hname;
  auto it = cache.find(key);
  if (it != cache.end()) return it->second;

  TH1* h = (TH1*)f->Get(hname);
  if (!h) {
    std::cout << "WARNING: Missing histogram " << hname
              << " in file " << f->GetName() << std::endl;
    return nullptr;
  }

  const TString uniqName = Form("%s__keep__%s__%lld", hname, keyPrefix.c_str(), ctr++);
  TH1* hc = (TH1*)h->Clone(uniqName);
  hc->SetDirectory(gROOT);

  cache[key] = hc;
  return hc;
}

static void Prep(TH1* h, bool normalize) {
  if (!h) return;
  h->SetStats(0);
  if (normalize) {
    const double I = h->Integral();
    if (I > 0) h->Scale(1.0 / I);
  }
}

// Draw one histogram overlay inside a pad region of a canvas: top (main) + bottom (ratio)
static void DrawOverlayWithRatio(TVirtualPad* parentPad,
                                 TH1* hs_in, TH1* hb_in,
                                 const char* sigLabel, const char* bkgLabel,
                                 bool normalize) {
  if (!parentPad || !hs_in || !hb_in) return;

  parentPad->cd();

  // Create split pads inside this sub-pad
  // Top pad
  TPad* pTop = new TPad(Form("%s_top", parentPad->GetName()), "", 0.0, 0.30, 1.0, 1.0);
  pTop->SetBottomMargin(0.02);
  pTop->SetLeftMargin(0.12);
  pTop->SetRightMargin(0.05);
  pTop->SetTopMargin(0.06);
  pTop->Draw();

  // Bottom pad (ratio)
  TPad* pBot = new TPad(Form("%s_bot", parentPad->GetName()), "", 0.0, 0.00, 1.0, 0.30);
  pBot->SetTopMargin(0.02);
  pBot->SetBottomMargin(0.35);
  pBot->SetLeftMargin(0.12);
  pBot->SetRightMargin(0.05);
  pBot->Draw();

  // Clone for per-draw styling
  TH1* hs = (TH1*)hs_in->Clone(Form("%s__draw_sig", hs_in->GetName()));
  TH1* hb = (TH1*)hb_in->Clone(Form("%s__draw_bkg", hb_in->GetName()));
  hs->SetDirectory(gROOT);
  hb->SetDirectory(gROOT);

  Prep(hs, normalize);
  Prep(hb, normalize);

  // Style
  hs->SetLineWidth(3);
  hb->SetLineWidth(3);
  hs->SetLineColor(kRed+1);
  hb->SetLineColor(kBlue+1);

  // TOP DRAW
  pTop->cd();
  pTop->SetGrid();

  // axis cosmetics (top)
  hs->GetXaxis()->SetLabelSize(0.0); // hide x labels on top pad
  hs->GetXaxis()->SetTitleSize(0.0);

  const double ymax = std::max(hs->GetMaximum(), hb->GetMaximum());
  hs->SetMaximum(1.30 * ymax);

  hs->Draw("HIST");
  hb->Draw("HIST SAME");

  TLegend* leg = new TLegend(0.55, 0.75, 0.90, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hs, sigLabel, "l");
  leg->AddEntry(hb, bkgLabel, "l");
  leg->Draw();

  // RATIO DRAW: S/B
  pBot->cd();
  pBot->SetGrid();

  TH1* hr = (TH1*)hs->Clone(Form("%s__ratio", hs->GetName()));
  hr->SetDirectory(gROOT);
  hr->Divide(hb);

  hr->SetLineColor(kBlack);
  hr->SetLineWidth(2);
  hr->SetMarkerStyle(20);
  hr->SetMarkerSize(0.6);

  hr->SetTitle("");
  hr->GetYaxis()->SetTitle("S/B");
  hr->GetYaxis()->CenterTitle(true);

  hr->GetYaxis()->SetNdivisions(505);
  hr->GetYaxis()->SetTitleSize(0.12);
  hr->GetYaxis()->SetLabelSize(0.10);
  hr->GetYaxis()->SetTitleOffset(0.50);

  hr->GetXaxis()->SetTitleSize(0.14);
  hr->GetXaxis()->SetLabelSize(0.12);
  hr->GetXaxis()->SetTitleOffset(1.10);

  // A reasonable y-range for ratios; adjust as you like
  hr->SetMinimum(0.0);
  hr->SetMaximum(2.0);

  hr->Draw("EP");

  parentPad->Modified();
  parentPad->Update();
}

void overlay_signal_keybkgs(int bkgCode = 2,
                            const char* signalTag = "tth12gev",
                            bool normalize = true,
                            bool keepWindowsOpen = true) {
  gStyle->SetOptStat(0);

  // Signal file
  const std::string fSig = std::string("output_signal_") + signalTag + ".root";

  // Background selection
  std::string fBkg, bkgLabel, bkgTag;
  if (bkgCode == 1) {
    fBkg     = "output_TTH_Hbb.root";
    bkgLabel = "ttH #rightarrow bb";
    bkgTag   = "TTH_Hbb";
  } else if (bkgCode == 2) {
    fBkg     = "output_ttbar.root";
    bkgLabel = "t#bar{t} dilepton";
    bkgTag   = "TT_dilep";
  } else if (bkgCode == 3) {
    fBkg     = "output_TTtoLNu2Q.root";
    bkgLabel = "t#bar{t} semileptonic";
    bkgTag   = "TT_semilep";
  } else {
    std::cout << "ERROR: Unknown bkgCode = " << bkgCode << "\n"
              << "Use: 1=ttH->bb, 2=tt dilep, 3=tt semilep\n";
    return;
  }

  TFile* fs = TFile::Open(fSig.c_str(), "READ");
  if (!fs || fs->IsZombie()) {
    std::cout << "ERROR: cannot open signal file " << fSig << "\n";
    return;
  }
  TFile* fb = TFile::Open(fBkg.c_str(), "READ");
  if (!fb || fb->IsZombie()) {
    std::cout << "ERROR: cannot open background file " << fBkg << "\n";
    fs->Close(); delete fs;
    return;
  }

  const std::string outTag = std::string("_sig_") + signalTag + "_vs_" + bkgTag
                           + (normalize ? "_norm" : "_raw");

  const std::string sigLabel = std::string("Signal ") + signalTag;

  // =========================
  // CANVAS 1: Higgs (3 pads)
  // =========================
  {
    TCanvas* c = new TCanvas(Form("c_higgs%s", outTag.c_str()), "Higgs (Mode-2)", 1600, 900);
    c->Divide(3, 1);

    const char* hnames[3] = {"h_Hdbb_mass", "h_Hdbb_pt", "h_Hdbb_eta"};
    for (int i = 0; i < 3; ++i) {
      TH1* hs = CloneToROOT(fs, std::string("SIG_") + signalTag, hnames[i]);hs->Rebin(2);
      TH1* hb = CloneToROOT(fb, std::string("BKG_") + bkgTag,   hnames[i]);hb->Rebin(2);
      if (!hs || !hb) continue;
      DrawOverlayWithRatio(c->cd(i+1), hs, hb, sigLabel.c_str(), bkgLabel.c_str(), normalize);
    }

    c->SaveAs(Form("Higgs_reco%s.png", outTag.c_str()));
    if (keepWindowsOpen) gCanvases.push_back(c);
    else delete c;
  }


  // =========================================
  // CANVAS 2: Double-b jet kinematics (8 pads)
  // =========================================
  {
    TCanvas* c = new TCanvas(Form("c_dbjets%s", outTag.c_str()), "Double-b jets (Mode-2)", 2000, 1000);
    c->Divide(4, 2);

    const char* hnames[8] = {
      "h_dbj1_mass", "h_dbj1_pt", "h_dbj1_eta", "h_dbj1_phi_final",
      "h_dbj2_mass", "h_dbj2_pt", "h_dbj2_eta", "h_dbj2_phi_final"
    };

    for (int i = 0; i < 8; ++i) {
      TH1* hs = CloneToROOT(fs, std::string("SIG_") + signalTag, hnames[i]);
      TH1* hb = CloneToROOT(fb, std::string("BKG_") + bkgTag,   hnames[i]);
      if (!hs || !hb) continue;
      DrawOverlayWithRatio(c->cd(i+1), hs, hb, sigLabel.c_str(), bkgLabel.c_str(), normalize);
    }

    c->SaveAs(Form("doubleb_all_kinematics%s.png", outTag.c_str()));
    if (keepWindowsOpen) gCanvases.push_back(c);
    else delete c;
  }

  // =========================================
  // CANVAS 3: Pair vars + MET/HT (2x3 = 6)
  // =========================================
  {
    TCanvas* c = new TCanvas(Form("c_pair_met_ht%s", outTag.c_str()), "Pair/MET/HT (Mode-2)", 1800, 900);
    c->Divide(3, 2);

    const char* hnames[6] = {
      "h_dR_dbj12_final", "h_dM_bj12_final", "h_MET_pt_final",
      "h_MET_phi_final",  "h_HT",            "h_mll"
    };

    for (int i = 0; i < 6; ++i) {
      TH1* hs = CloneToROOT(fs, std::string("SIG_") + signalTag, hnames[i]);
      TH1* hb = CloneToROOT(fb, std::string("BKG_") + bkgTag,   hnames[i]);
      if (!hs || !hb) continue;
      DrawOverlayWithRatio(c->cd(i+1), hs, hb, sigLabel.c_str(), bkgLabel.c_str(), normalize);
    }

    c->SaveAs(Form("pair_MET_HT%s.png", outTag.c_str()));
    if (keepWindowsOpen) gCanvases.push_back(c);
    else delete c;
  }

  // =========================================
  // CANVAS 4: Final b-jets (2x2 = 4)
  // =========================================
  {
    TCanvas* c = new TCanvas(Form("c_bjets%s", outTag.c_str()), "Final b-jets (Mode-2)", 1400, 900);
    c->Divide(2, 2);

    const char* hnames[4] = {
      "h_bj1_pt_final", "h_bj1_eta_final",
      "h_bj2_pt_final", "h_bj2_eta_final"
    };

    for (int i = 0; i < 4; ++i) {
      TH1* hs = CloneToROOT(fs, std::string("SIG_") + signalTag, hnames[i]);
      TH1* hb = CloneToROOT(fb, std::string("BKG_") + bkgTag,   hnames[i]);
      if (!hs || !hb) continue;
      DrawOverlayWithRatio(c->cd(i+1), hs, hb, sigLabel.c_str(), bkgLabel.c_str(), normalize);
    }

    c->SaveAs(Form("final_bjets%s.png", outTag.c_str()));
    if (keepWindowsOpen) gCanvases.push_back(c);
    else delete c;
  }

  // =========================================
  // CANVAS 5: Final lepton + dRll (2x2 = 4)
  // =========================================
  {
    TCanvas* c = new TCanvas(Form("c_lep_dRll%s", outTag.c_str()), "Final lepton + dRll (Mode-2)", 1400, 900);
    c->Divide(2, 2);

    const char* hnames[4] = {
      "h_lep1_pt_final", "h_lep1_eta_final",
      "h_lep1_phi_final","h_dRll"
    };

    for (int i = 0; i < 4; ++i) {
      TH1* hs = CloneToROOT(fs, std::string("SIG_") + signalTag, hnames[i]);
      TH1* hb = CloneToROOT(fb, std::string("BKG_") + bkgTag,   hnames[i]);
      if (!hs || !hb) continue;
      DrawOverlayWithRatio(c->cd(i+1), hs, hb, sigLabel.c_str(), bkgLabel.c_str(), normalize);
    }

    c->SaveAs(Form("final_lepton_dRll%s.png", outTag.c_str()));
    if (keepWindowsOpen) gCanvases.push_back(c);
    else delete c;
  }

  // Close files safely (hists are cloned to gROOT)
  fb->Close(); delete fb;
  fs->Close(); delete fs;

  std::cout << "Done. Saved multi-pad overlays with tag: " << outTag << "\n";
  if (keepWindowsOpen) {
    std::cout << "Windows will stay open. To close/clear later: gCanvases.clear();\n";
  }
}
