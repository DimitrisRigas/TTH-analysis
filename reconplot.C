// reconplot.C
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TPad.h>
#include <TString.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TPaveStats.h>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>

// ============================================================================
// Sample selection helper
// ============================================================================
bool getSampleInfo(int sample,
                   const char* signalTag,
                   std::string& fname,
                   std::string& tag,
                   std::string& folderTag) {

  if (sample == 0) {
    fname     = std::string("output_signal_") + signalTag + ".root";
    tag       = std::string("_signal_") + signalTag;
    folderTag = std::string("signal_") + signalTag;
  }
  else if (sample == 1)  { fname = "output_ttbar.root";             tag = "_ttbar";             folderTag = "ttbar"; }
  else if (sample == 2)  { fname = "output_DYee.root";              tag = "_DYee";              folderTag = "DYee"; }
  else if (sample == 3)  { fname = "output_DYmumu.root";            tag = "_DYmumu";            folderTag = "DYmumu"; }

  // H->bb backgrounds
  else if (sample == 4)  { fname = "output_TTH_Hbb.root";           tag = "_TTH_Hbb";           folderTag = "TTH_Hbb"; }
  else if (sample == 5)  { fname = "output_VBF_Hbb.root";           tag = "_VBF_Hbb";           folderTag = "VBF_Hbb"; }
  else if (sample == 6)  { fname = "output_GGH_Hbb.root";           tag = "_GGH_Hbb";           folderTag = "GGH_Hbb"; }

  // Dibosons
  else if (sample == 7)  { fname = "output_WW.root";                tag = "_WW";                folderTag = "WW"; }
  else if (sample == 8)  { fname = "output_WZ.root";                tag = "_WZ";                folderTag = "WZ"; }
  else if (sample == 9)  { fname = "output_ZZ.root";                tag = "_ZZ";                folderTag = "ZZ"; }

  // Single-top / related
  else if (sample == 10) { fname = "output_TBbarQtoLNu.root";       tag = "_TBbarQtoLNu";       folderTag = "TBbarQtoLNu"; }
  else if (sample == 11) { fname = "output_TBbarQto2Q.root";        tag = "_TBbarQto2Q";        folderTag = "TBbarQto2Q"; }
  else if (sample == 12) { fname = "output_TTtoLNu2Q.root";         tag = "_TTtoLNu2Q";         folderTag = "TTtoLNu2Q"; }
  else if (sample == 13) { fname = "output_TbarWplusToNu2Q.root";   tag = "_TbarWplusToNu2Q";   folderTag = "TbarWplusToNu2Q"; }
  else if (sample == 14) { fname = "output_TbarWplusTo4Q.root";     tag = "_TbarWplusTo4Q";     folderTag = "TbarWplusTo4Q"; }

  // Extra backgrounds from your MyClass2 output logic, if present
  else if (sample == 15) { fname = "output_TTW.root";               tag = "_TTW";               folderTag = "TTW"; }
  else if (sample == 16) { fname = "output_TTZ.root";               tag = "_TTZ";               folderTag = "TTZ"; }

  else {
    std::cout << "ERROR: Unknown sample code: " << sample << std::endl;
    std::cout << "Valid sample codes:" << std::endl;
    std::cout << "  0  = signal" << std::endl;
    std::cout << "  1  = ttbar" << std::endl;
    std::cout << "  2  = DYee" << std::endl;
    std::cout << "  3  = DYmumu" << std::endl;
    std::cout << "  4  = TTH_Hbb" << std::endl;
    std::cout << "  5  = VBF_Hbb" << std::endl;
    std::cout << "  6  = GGH_Hbb" << std::endl;
    std::cout << "  7  = WW" << std::endl;
    std::cout << "  8  = WZ" << std::endl;
    std::cout << "  9  = ZZ" << std::endl;
    std::cout << "  10 = TBbarQtoLNu" << std::endl;
    std::cout << "  11 = TBbarQto2Q" << std::endl;
    std::cout << "  12 = TTtoLNu2Q" << std::endl;
    std::cout << "  13 = TbarWplusToNu2Q" << std::endl;
    std::cout << "  14 = TbarWplusTo4Q" << std::endl;
    std::cout << "  15 = TTW" << std::endl;
    std::cout << "  16 = TTZ" << std::endl;
    return false;
  }

  return true;
}

// ============================================================================
// Save canvas helper
// ============================================================================
void saveCanvasPDF(TCanvas* c, const std::string& outDir, const char* fileName) {
  if (!c) return;

  c->Modified();
  c->Update();

  gSystem->mkdir("reco_plots", kTRUE);
  gSystem->mkdir(outDir.c_str(), kTRUE);

  const std::string outPath = outDir + "/" + fileName;
  c->SaveAs(outPath.c_str());
}

// ============================================================================
// Main plotting function
// mode = 1 -> object-level reco plots
// mode = 2 -> final reconstruction plots
// sample = 0 signal, 1 ttbar, etc.
// ============================================================================
void reconplot(int mode = 1, int sample = 0, const char* signalTag = "tth12gev") {

  gStyle->SetOptStat(1111);

  // =======================================
  // Select input file by sample
  // =======================================
  std::string fname;
  std::string tag;
  std::string folderTag;

  if (!getSampleInfo(sample, signalTag, fname, tag, folderTag)) {
    return;
  }

  TFile *inFile = TFile::Open(fname.c_str(), "READ");

  if (!inFile || inFile->IsZombie()) {
    std::cout << "ERROR: Could not open " << fname << std::endl;
    return;
  }

  std::cout << "\n============================================" << std::endl;
  std::cout << "Loaded file: " << fname << std::endl;
  std::cout << "Mode: " << mode << std::endl;
  std::cout << "Output folder tag: " << folderTag << std::endl;
  std::cout << "============================================" << std::endl;

  // Output folder
  std::string modeFolder;

  if (mode == 1) {
    modeFolder = "mode1_object_level";
  }
  else if (mode == 2) {
    modeFolder = "mode2_final_reco";
  }
  else {
    std::cout << "Invalid mode. Use mode=1 or mode=2." << std::endl;
    inFile->Close();
    delete inFile;
    return;
  }

  std::string outDir = "reco_plots/" + modeFolder + "/" + folderTag;

  gSystem->mkdir("reco_plots", kTRUE);
  gSystem->mkdir(("reco_plots/" + modeFolder).c_str(), kTRUE);
  gSystem->mkdir(outDir.c_str(), kTRUE);

  // ============================================================
  // Persistent histogram handling
  // ============================================================
  static std::unordered_map<std::string, TH1*> hcache;
  static long long uniqueCounter = 0;

  auto GetH = [&](const char *hname) -> TH1* {

    // Cache key scoped by input file + mode + sample tag + histogram name
    const std::string key = fname + "::mode" + std::to_string(mode) + "::" + tag + "::" + hname;

    auto it = hcache.find(key);
    if (it != hcache.end() && it->second) return it->second;

    TH1 *h = (TH1*)inFile->Get(hname);

    if (!h) {
      std::cout << "WARNING: Missing histogram " << hname
                << " in " << fname << std::endl;
      return nullptr;
    }

    // Unique clone name
    const TString uniqName = Form("%s__uniq__%s__m%d__s%d__%lld",
                                  hname, tag.c_str(), mode, sample, uniqueCounter++);

    TH1 *hu = (TH1*)h->Clone(uniqName);
    hu->SetDirectory(gROOT);

    // Stable alias
    const TString aliasName = Form("%s__alias__%s__m%d",
                                   hname, tag.c_str(), mode);

    if (TObject *old = gROOT->FindObject(aliasName)) {
      old->Delete();
    }

    TH1 *ha = (TH1*)hu->Clone(aliasName);
    ha->SetDirectory(gROOT);

    hcache[key] = ha;
    return ha;
  };

  auto DrawHist = [&](int padNo, TCanvas *can, const char *hname) {
    can->cd(padNo);

    TH1 *hist = GetH(hname);

    if (hist) {
      hist->Draw("HIST");
      gPad->SetGrid();
    }
  };

  // =====================================================================
  // MODE 1 → OBJECT-LEVEL RECO PLOTS
  // =====================================================================
  if (mode == 1) {

    // ================= ELECTRONS =================
    TCanvas *c_ele = new TCanvas(Form("c_ele%s", tag.c_str()), "Electrons", 1600, 900);
    c_ele->Divide(3, 2);

    DrawHist(1, c_ele, "h_e1_pt");
    DrawHist(2, c_ele, "h_e1_eta");
    DrawHist(3, c_ele, "h_e1_phi");
    DrawHist(4, c_ele, "h_e2_pt");
    DrawHist(5, c_ele, "h_e2_eta");
    DrawHist(6, c_ele, "h_e2_phi");

    saveCanvasPDF(c_ele, outDir, "electrons.pdf");

    // ================= MUONS =================
    TCanvas *c_mu = new TCanvas(Form("c_mu%s", tag.c_str()), "Muons", 1600, 900);
    c_mu->Divide(3, 2);

    DrawHist(1, c_mu, "h_mu1_pt");
    DrawHist(2, c_mu, "h_mu1_eta");
    DrawHist(3, c_mu, "h_mu1_phi");
    DrawHist(4, c_mu, "h_mu2_pt");
    DrawHist(5, c_mu, "h_mu2_eta");
    DrawHist(6, c_mu, "h_mu2_phi");

    saveCanvasPDF(c_mu, outDir, "muons.pdf");

    // ================= JETS =================
    TCanvas *c_jet = new TCanvas(Form("c_jet%s", tag.c_str()), "Jets", 1600, 1200);
    c_jet->Divide(3, 4);

    for (int i = 1; i <= 4; ++i) {
      DrawHist(3*(i-1)+1, c_jet, Form("h_j%d_pt",  i));
      DrawHist(3*(i-1)+2, c_jet, Form("h_j%d_eta", i));
      DrawHist(3*(i-1)+3, c_jet, Form("h_j%d_phi", i));
    }

    saveCanvasPDF(c_jet, outDir, "jets.pdf");


    // ================= MULTIPLICITIES =================
    // Thesis-friendly 2x3 canvas
    {
      // Turn off stat boxes for thesis-quality plots
      gStyle->SetOptStat("emr");// entries, mean, RMS
      
      TCanvas *c_mult = new TCanvas(Form("c_mult%s", tag.c_str()),
				    "Object multiplicities", 1600, 1000);
      
      c_mult->Divide(3, 2, 0.01, 0.01);
      
      std::vector<std::pair<const char*, const char*> > multHists = {
	{"h_nEle",     "Electron multiplicity;N_{e};Events"},
	{"h_nMu",      "Muon multiplicity;N_{#mu};Events"},
	{"h_nJet",     "Jet multiplicity;N_{jets};Events"},
	{"h_nCJet",    "Clean jet multiplicity;N_{clean jets};Events"},
	{"h_nBjet",    "b-jet multiplicity;N_{b-jets};Events"},
	{"h_nDoubleB", "Double-b jet multiplicity;N_{double-b};Events"}
      };
      
      for (int i = 0; i < 6; ++i) {
	c_mult->cd(i + 1);
	
    gPad->SetGrid();
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.13);
    gPad->SetTopMargin(0.10);
    
    TH1 *h = GetH(multHists[i].first);
    
    if (!h) continue;
    
    h->SetTitle(multHists[i].second);
    h->SetLineWidth(2);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.6);
    
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetLabelSize(0.040);
    h->GetYaxis()->SetLabelSize(0.040);
    h->GetXaxis()->SetTitleOffset(1.10);
    h->GetYaxis()->SetTitleOffset(1.35);
    
    h->Draw("HIST");
    gPad->Update();

    TPaveStats *st = (TPaveStats*)h->FindObject("stats");
    if (st) {
      st->SetX1NDC(0.62);
      st->SetX2NDC(0.93);
      st->SetY1NDC(0.72);
      st->SetY2NDC(0.88);
      st->SetTextSize(0.030);
      st->SetFillColor(0);
      st->SetBorderSize(1);
    }
      }
      
      saveCanvasPDF(c_mult, outDir, "multiplicity_thesis.pdf");
 }
    
    // ================= ΔR CLEANING =================
    TCanvas *c_dR = new TCanvas(Form("c_dR%s", tag.c_str()), "Jet-Lepton dR", 1600, 1200);
    c_dR->Divide(2, 2);

    DrawHist(1, c_dR, "h_dR_jet_ele_before");
    DrawHist(2, c_dR, "h_dR_jet_mu_before");
    DrawHist(3, c_dR, "h_dR_jet_ele_after");
    DrawHist(4, c_dR, "h_dR_jet_mu_after");

    saveCanvasPDF(c_dR, outDir, "dR_jet_lepton.pdf");

    // ================= MET PRE-SELECTION =================
    TCanvas *c_met_pre = new TCanvas(Form("c_met_pre%s", tag.c_str()), "MET preselection", 1200, 600);
    c_met_pre->Divide(2, 1);

    DrawHist(1, c_met_pre, "h_pre_MET_pt");
    DrawHist(2, c_met_pre, "h_pre_MET_phi");

    saveCanvasPDF(c_met_pre, outDir, "MET_preselection.pdf");

    std::cout << "MODE 1 PDF plots saved in: " << outDir << std::endl;

    inFile->Close();
    delete inFile;
    return;
  }

  // =====================================================================
  // MODE 2 → FINAL RECONSTRUCTION
  // =====================================================================
  if (mode == 2) {

    // ================= HIGGS =================
    TCanvas *c_h = new TCanvas(Form("c_h%s", tag.c_str()), "Higgs", 1600, 900);
    c_h->Divide(3, 1);

    DrawHist(1, c_h, "h_Hdbb_mass");
    DrawHist(2, c_h, "h_Hdbb_pt");
    DrawHist(3, c_h, "h_Hdbb_eta");

    saveCanvasPDF(c_h, outDir, "Higgs_reco.pdf");

    // ================= MET PRE VS FINAL =================
    TCanvas *c_met = new TCanvas(Form("c_met%s", tag.c_str()), "MET", 1800, 600);
    c_met->Divide(3, 1);

    DrawHist(1, c_met, "h_pre_MET_pt");
    DrawHist(2, c_met, "h_MET_pt_final");
    DrawHist(3, c_met, "h_MET_phi_final");

    saveCanvasPDF(c_met, outDir, "MET.pdf");

    // ================= HT FINAL =================
    if (TH1 *hht = GetH("h_HT")) {
      TCanvas *c_ht = new TCanvas(Form("c_ht%s", tag.c_str()), "HT", 1000, 800);
      c_ht->cd();
      hht->Draw("HIST");
      gPad->SetGrid();

      saveCanvasPDF(c_ht, outDir, "HT.pdf");
    }

    // ================= DOUBLE-B JETS =================
    TCanvas *c_dbk = new TCanvas(Form("c_dbk%s", tag.c_str()), "DoubleB jets", 2000, 1000);
    c_dbk->Divide(4, 2);

    DrawHist(1, c_dbk, "h_dbj1_mass");
    DrawHist(2, c_dbk, "h_dbj1_pt");
    DrawHist(3, c_dbk, "h_dbj1_eta");
    DrawHist(4, c_dbk, "h_dbj1_phi_final");

    DrawHist(5, c_dbk, "h_dbj2_mass");
    DrawHist(6, c_dbk, "h_dbj2_pt");
    DrawHist(7, c_dbk, "h_dbj2_eta");
    DrawHist(8, c_dbk, "h_dbj2_phi_final");

    saveCanvasPDF(c_dbk, outDir, "doubleb_all_kinematics.pdf");

    // ================= ΔR(db1, db2) =================
    if (TH1 *h_dR = GetH("h_dR_dbj12_final")) {
      TCanvas *c_dRdb = new TCanvas(Form("c_dRdb%s", tag.c_str()), "DeltaR db1 db2", 1000, 800);
      c_dRdb->cd();
      h_dR->Draw("HIST");
      gPad->SetGrid();

      saveCanvasPDF(c_dRdb, outDir, "DeltaR_dbj1_dbj2.pdf");
    }

    // ================= Δm(db1, db2) =================
    if (TH1 *h_dM = GetH("h_dM_bj12_final")) {
      TCanvas *c_dm = new TCanvas(Form("c_dm%s", tag.c_str()), "DeltaM double-b jets", 1000, 800);
      c_dm->cd();
      h_dM->Draw("HIST");
      gPad->SetGrid();

      saveCanvasPDF(c_dm, outDir, "DeltaM_dbj1_dbj2.pdf");
    }

    // ================= FINAL B-JETS =================
    TCanvas *c_bj = new TCanvas(Form("c_bj%s", tag.c_str()), "Final b-jets", 1200, 800);
    c_bj->Divide(2, 2);

    DrawHist(1, c_bj, "h_bj1_pt_final");
    DrawHist(2, c_bj, "h_bj1_eta_final");
    DrawHist(3, c_bj, "h_bj2_pt_final");
    DrawHist(4, c_bj, "h_bj2_eta_final");

    saveCanvasPDF(c_bj, outDir, "final_bjets.pdf");

    // ================= FINAL LEPTONS =================
    TCanvas *c_lep = new TCanvas(Form("c_lep%s", tag.c_str()), "Final leptons", 1600, 900);
    c_lep->Divide(3, 2);

    DrawHist(1, c_lep, "h_lep1_pt_final");
    DrawHist(2, c_lep, "h_lep1_eta_final");
    DrawHist(3, c_lep, "h_lep1_phi_final");
    DrawHist(4, c_lep, "h_lep2_pt_final");
    DrawHist(5, c_lep, "h_lep2_eta_final");
    DrawHist(6, c_lep, "h_lep2_phi_final");

    saveCanvasPDF(c_lep, outDir, "final_leptons.pdf");

    // ================= DILEPTON =================
    TCanvas *c_ll = new TCanvas(Form("c_ll%s", tag.c_str()), "Dilepton", 1200, 600);
    c_ll->Divide(2, 1);

    DrawHist(1, c_ll, "h_mll");
    DrawHist(2, c_ll, "h_dRll");

    saveCanvasPDF(c_ll, outDir, "dilepton.pdf");

    std::cout << "MODE 2 PDF plots saved in: " << outDir << std::endl;

    inFile->Close();
    delete inFile;
    return;
  }
}

// ============================================================================
// Automatically find all output_signal_*.root files
// ============================================================================
std::vector<std::string> findSignalTags() {

  std::vector<std::string> signalTags;

  TSystemDirectory dir(".", ".");
  TList *files = dir.GetListOfFiles();

  if (files) {
    TIter next(files);
    TSystemFile *file;

    while ((file = (TSystemFile*)next())) {
      TString fname = file->GetName();

      if (file->IsDirectory()) continue;

      if (fname.BeginsWith("output_signal_") && fname.EndsWith(".root")) {
        fname.ReplaceAll("output_signal_", "");
        fname.ReplaceAll(".root", "");
        signalTags.push_back(std::string(fname.Data()));
      }
    }
  }

  std::sort(signalTags.begin(), signalTags.end());
  return signalTags;
}

// ============================================================================
// Run reco plots for all signal mass points for one mode
// Example:
//   reconplotAllSignals(1)
//   reconplotAllSignals(2)
// ============================================================================
void reconplotAllSignals(int mode = 1) {

  std::vector<std::string> signalTags = findSignalTags();

  if (signalTags.empty()) {
    std::cout << "ERROR: No signal files found." << std::endl;
    std::cout << "Expected files like output_signal_tth12gev.root" << std::endl;
    return;
  }

  std::cout << "Found " << signalTags.size() << " signal mass point(s)." << std::endl;

  for (const auto& sig : signalTags) {
    reconplot(mode, 0, sig.c_str());
  }
}

// ============================================================================
// Run reco plots for all available samples for one mode
// Signal mass point is selected by signalTag.
// Missing files are skipped safely.
// Example:
//   reconplotAllSamples(1, "tth12gev")
//   reconplotAllSamples(2, "tth12gev")
// ============================================================================
void reconplotAllSamples(int mode = 1, const char* signalTag = "tth12gev") {

  for (int sample = 0; sample <= 16; ++sample) {

    std::string fname;
    std::string tag;
    std::string folderTag;

    if (!getSampleInfo(sample, signalTag, fname, tag, folderTag)) continue;

    if (gSystem->AccessPathName(fname.c_str())) {
      std::cout << "Skipping missing file: " << fname << std::endl;
      continue;
    }

    reconplot(mode, sample, signalTag);
  }
}

// ============================================================================
// Convenience: run both modes for all signal mass points
// ============================================================================
void reconplotAllSignalsBothModes() {
  reconplotAllSignals(1);
  reconplotAllSignals(2);
}

// ============================================================================
// Convenience: run both modes for all available samples
// ============================================================================
void reconplotAllSamplesBothModes(const char* signalTag = "tth12gev") {
  reconplotAllSamples(1, signalTag);
  reconplotAllSamples(2, signalTag);
}
