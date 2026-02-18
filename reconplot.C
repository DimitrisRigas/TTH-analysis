// reconplot.C
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TPad.h>
#include <TString.h>
#include <iostream>
#include <string>
#include <unordered_map>

void reconplot(int mode = 1, int sample = 0, const char* signalTag = "tth12gev") {

  gStyle->SetOptStat(1111);

  // =======================================
  // Select input file by sample
  // =======================================
  std::string fname;
  std::string tag;

  // --- Signals (all masses handled via signalTag) ---
  // Examples:
  //   reconplot(1, 0, "tth12gev");
  //   reconplot(1, 0, "tth15gev");
  //   reconplot(1, 0, "tth20gev");
  //   reconplot(1, 0, "tth25gev");
  //   reconplot(1, 0, "tth30gev");
  //   reconplot(1, 0, "tth60gev");
  if (sample == 0) {
    fname = std::string("output_signal_") + signalTag + ".root";
    tag   = std::string("_signal_") + signalTag;
  }
  else if (sample == 1) { fname = "output_ttbar.root";  tag = "_ttbar"; }
  else if (sample == 2) { fname = "output_DYee.root";   tag = "_DYee"; }
  else if (sample == 3) { fname = "output_DYmumu.root"; tag = "_DYmumu"; }

  // --- H->bb backgrounds ---
  else if (sample == 4) { fname = "output_TTH_Hbb.root"; tag = "_TTH_Hbb"; }
  else if (sample == 5) { fname = "output_VBF_Hbb.root"; tag = "_VBF_Hbb"; }
  else if (sample == 6) { fname = "output_GGH_Hbb.root"; tag = "_GGH_Hbb"; }

  // --- Dibosons ---
  else if (sample == 7) { fname = "output_WW.root"; tag = "_WW"; }
  else if (sample == 8) { fname = "output_WZ.root"; tag = "_WZ"; }
  else if (sample == 9) { fname = "output_ZZ.root"; tag = "_ZZ"; }

  // --- Single-top / related ---
  else if (sample == 10) { fname = "output_TBbarQtoLNu.root";     tag = "_TBbarQtoLNu"; }
  else if (sample == 11) { fname = "output_TBbarQto2Q.root";      tag = "_TBbarQto2Q"; }
  else if (sample == 12) { fname = "output_TTtoLNu2Q.root";       tag = "_TTtoLNu2Q"; }
  else if (sample == 13) { fname = "output_TbarWplusToNu2Q.root"; tag = "_TbarWplusToNu2Q"; }
  else if (sample == 14) { fname = "output_TbarWplusTo4Q.root";   tag = "_TbarWplusTo4Q"; }

  else {
    std::cout << "ERROR: Unknown sample code: " << sample << std::endl;
    std::cout << "Valid: 0=signal, 1=ttbar, 2=DYee, 3=DYmumu, 4=TTH_Hbb, 5=VBF_Hbb, 6=GGH_Hbb, "
                 "7=WW, 8=WZ, 9=ZZ, 10=TBbarQtoLNu, 11=TBbarQto2Q, 12=TTtoLNu2Q, "
                 "13=TbarWplusToNu2Q, 14=TbarWplusTo4Q\n";
    return;
  }

  TFile *inFile = TFile::Open(fname.c_str(), "READ");
  if (!inFile || inFile->IsZombie()) {
    std::cout << "ERROR: Could not open " << fname << std::endl;
    return;
  }

  std::cout << "Loaded file: " << fname << std::endl;

  // ============================================================
  // Persistent histogram handling:
  //  - We CLONE histograms out of the file
  //  - We put the clones into gROOT so they survive after the macro ends
  //  - We keep a stable alias name so you can FindObject() later easily
  // ============================================================
  static std::unordered_map<std::string, TH1*> hcache;
  static long long uniqueCounter = 0;

  auto GetH = [&](const char *hname) -> TH1* {

    // Cache key scoped by sample tag + histogram name
    const std::string key = tag + "::" + hname;
    auto it = hcache.find(key);
    if (it != hcache.end() && it->second) return it->second;

    TH1 *h = (TH1*)inFile->Get(hname);
    if (!h) {
      std::cout << "WARNING: Missing histogram " << hname << " in " << fname << std::endl;
      return nullptr;
    }

    // Unique clone name (never collides)
    const TString uniqName = Form("%s__uniq__%s__m%d__s%d__%lld",
                                  hname, tag.c_str(), mode, sample, uniqueCounter++);

    // Clone and put into gROOT so it persists
    TH1 *hu = (TH1*)h->Clone(uniqName);
    hu->SetDirectory(gROOT);

    // Stable alias name (easy to retrieve later)
    const TString aliasName = Form("%s__alias__%s", hname, tag.c_str());

    // If alias exists from previous runs, delete it (avoid replacement warnings and stale pointers)
    if (TObject *old = gROOT->FindObject(aliasName)) {
      old->Delete(); // removes from ROOT lists
    }

    TH1 *ha = (TH1*)hu->Clone(aliasName);
    ha->SetDirectory(gROOT);

    // Cache the alias (stable name)
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

    c_ele->SaveAs(("electrons" + tag + ".png").c_str());

    // ================= MUONS =================
    TCanvas *c_mu = new TCanvas(Form("c_mu%s", tag.c_str()), "Muons", 1600, 900);
    c_mu->Divide(3, 2);

    DrawHist(1, c_mu, "h_mu1_pt");
    DrawHist(2, c_mu, "h_mu1_eta");
    DrawHist(3, c_mu, "h_mu1_phi");
    DrawHist(4, c_mu, "h_mu2_pt");
    DrawHist(5, c_mu, "h_mu2_eta");
    DrawHist(6, c_mu, "h_mu2_phi");

    c_mu->SaveAs(("muons" + tag + ".png").c_str());

    // ================= JETS =================
    TCanvas *c_jet = new TCanvas(Form("c_jet%s", tag.c_str()), "Jets", 1600, 1200);
    c_jet->Divide(3, 4);

    for (int i = 1; i <= 4; ++i) {
      DrawHist(3*(i-1)+1, c_jet, Form("h_j%d_pt",  i));
      DrawHist(3*(i-1)+2, c_jet, Form("h_j%d_eta", i));
      DrawHist(3*(i-1)+3, c_jet, Form("h_j%d_phi", i));
    }

    c_jet->SaveAs(("jets" + tag + ".png").c_str());

    // ================= MULTIPLICITIES =================
    TCanvas *c_mult = new TCanvas(Form("c_mult%s", tag.c_str()), "Multiplicities", 1800, 600);
    c_mult->Divide(6, 1);

    DrawHist(1, c_mult, "h_nEle");
    DrawHist(2, c_mult, "h_nMu");
    DrawHist(3, c_mult, "h_nJet");
    DrawHist(4, c_mult, "h_nCJet");
    DrawHist(5, c_mult, "h_nBjet");
    DrawHist(6, c_mult, "h_nDoubleB"); // NEW: was missing from the canvas

    c_mult->SaveAs(("multiplicity" + tag + ".png").c_str());

    // ================= ΔR CLEANING =================
    TCanvas *c_dR = new TCanvas(Form("c_dR%s", tag.c_str()), "Jet–Lepton dR", 1600, 1200);
    c_dR->Divide(2, 2);

    DrawHist(1, c_dR, "h_dR_jet_ele_before");
    DrawHist(2, c_dR, "h_dR_jet_mu_before");
    DrawHist(3, c_dR, "h_dR_jet_ele_after");
    DrawHist(4, c_dR, "h_dR_jet_mu_after");

    c_dR->SaveAs(("dR_jet_lepton" + tag + ".png").c_str());

    // ================= MET (PRE-SELECTION) =================
    // Filled in MyClass before the MET cut → shows MET distribution before cut
    TCanvas *c_met_pre = new TCanvas(Form("c_met_pre%s", tag.c_str()), "MET preselection", 1200, 600);
    c_met_pre->Divide(2, 1);

    DrawHist(1, c_met_pre, "h_pre_MET_pt");
    DrawHist(2, c_met_pre, "h_pre_MET_phi");

    c_met_pre->SaveAs(("MET_preselection" + tag + ".png").c_str());

    std::cout << "MODE 1 plots saved." << std::endl;

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

    c_h->SaveAs(("Higgs_reco" + tag + ".png").c_str());

    // ================= MET (PRE vs FINAL) =================
    // pre: filled before MET cut
    // final: filled after full selection
    TCanvas *c_met = new TCanvas(Form("c_met%s", tag.c_str()), "MET", 1800, 600);
    c_met->Divide(3, 1);

    DrawHist(1, c_met, "h_pre_MET_pt");
    DrawHist(2, c_met, "h_MET_pt_final");
    DrawHist(3, c_met, "h_MET_phi_final");

    c_met->SaveAs(("MET" + tag + ".png").c_str());

    // ================= HT (FINAL) =================
    if (TH1 *hht = GetH("h_HT")) {
      TCanvas *c_ht = new TCanvas(Form("c_ht%s", tag.c_str()), "HT", 1000, 800);
      c_ht->cd();
      hht->Draw("HIST");
      gPad->SetGrid();
      c_ht->SaveAs(("HT" + tag + ".png").c_str());
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

    c_dbk->SaveAs(("doubleb_all_kinematics" + tag + ".png").c_str());

    // ================= ΔR(db1, db2) =================
    if (TH1 *h_dR = GetH("h_dR_dbj12_final")) {
      TCanvas *c_dRdb = new TCanvas(Form("c_dRdb%s", tag.c_str()), "DeltaR db1 db2", 1000, 800);
      c_dRdb->cd();
      h_dR->Draw("HIST");
      gPad->SetGrid();
      c_dRdb->SaveAs(("DeltaR_dbj1_dbj2" + tag + ".png").c_str());
    }

    // ================= Δm(db1, db2) =================
    if (TH1 *h_dM = GetH("h_dM_bj12_final")) {
      TCanvas *c_dm = new TCanvas(Form("c_dm%s", tag.c_str()), "DeltaM double-b jets", 1000, 800);
      c_dm->cd();
      h_dM->Draw("HIST");
      gPad->SetGrid();
      c_dm->SaveAs(("DeltaM_dbj1_dbj2" + tag + ".png").c_str());
    }

    // ================= FINAL B-JETS =================
    TCanvas *c_bj = new TCanvas(Form("c_bj%s", tag.c_str()), "Final b-jets", 1200, 800);
    c_bj->Divide(2, 2);

    DrawHist(1, c_bj, "h_bj1_pt_final");
    DrawHist(2, c_bj, "h_bj1_eta_final");
    DrawHist(3, c_bj, "h_bj2_pt_final");
    DrawHist(4, c_bj, "h_bj2_eta_final");

    c_bj->SaveAs(("final_bjets" + tag + ".png").c_str());

    // ================= FINAL LEPTON =================
    TCanvas *c_lep = new TCanvas(Form("c_lep%s", tag.c_str()), "Final leptons", 1600, 900);
    c_lep->Divide(3, 2);

    DrawHist(1, c_lep, "h_lep1_pt_final");
    DrawHist(2, c_lep, "h_lep1_eta_final");
    DrawHist(3, c_lep, "h_lep1_phi_final");
    DrawHist(4, c_lep, "h_lep2_pt_final");
    DrawHist(5, c_lep, "h_lep2_eta_final");
    DrawHist(6, c_lep, "h_lep2_phi_final");

    c_lep->SaveAs(("final_leptons" + tag + ".png").c_str());

    // ================= DILEPTON =================
    TCanvas *c_ll = new TCanvas(Form("c_ll%s", tag.c_str()), "Dilepton", 1200, 600);
    c_ll->Divide(2, 1);

    DrawHist(1, c_ll, "h_mll");
    DrawHist(2, c_ll, "h_dRll");

    c_ll->SaveAs(("dilepton" + tag + ".png").c_str());

    inFile->Close();
    delete inFile;
    return;
  }

  std::cout << "Invalid mode." << std::endl;
  inFile->Close();
  delete inFile;
}
