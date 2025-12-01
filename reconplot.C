#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TPad.h>
#include <iostream>

// Simple RECO plotting macro for pre-selection histograms.
// It always reads the same histograms; whether they are raw or weighted
// depends on how you filled them in MyClass::Loop.

void reconplot() {

    gStyle->SetOptStat(1111);

    TFile *f = new TFile("ttH_aab_final_output.root");
    if (!f || f->IsZombie()) {
        std::cout << "ERROR: Could not open ttH_aab_final_output.root!" << std::endl;
        return;
    }

    // Helper to fetch histogram
    auto H = [&](std::string name) -> TH1F* {
        TH1F *h = (TH1F*)f->Get(name.c_str());
        if (!h) {
            std::cout << "WARNING: Missing histogram " << name << std::endl;
        }
        return h;
    };

    std::string tag = "_raw";  // just used for output file names

    // ============================================================
    // PRE-SELECTION CANVASES
    // ============================================================

    // --- Electrons (leading and subleading) ---
    TCanvas *c_ele = new TCanvas("c_ele", ("Electrons (pre-selection)" + tag).c_str(), 1600, 900);
    c_ele->Divide(3, 2);

    c_ele->cd(1); H("h_e1_pt") ->Draw("HIST"); gPad->SetGrid();
    c_ele->cd(2); H("h_e1_eta")->Draw("HIST"); gPad->SetGrid();
    c_ele->cd(3); H("h_e1_phi")->Draw("HIST"); gPad->SetGrid();

    c_ele->cd(4); H("h_e2_pt") ->Draw("HIST"); gPad->SetGrid();
    c_ele->cd(5); H("h_e2_eta")->Draw("HIST"); gPad->SetGrid();
    c_ele->cd(6); H("h_e2_phi")->Draw("HIST"); gPad->SetGrid();

    c_ele->SaveAs(("electrons_pre" + tag + ".png").c_str());

    // --- Muons (leading and subleading) ---
    TCanvas *c_mu = new TCanvas("c_mu", ("Muons (pre-selection)" + tag).c_str(), 1600, 900);
    c_mu->Divide(3, 2);

    c_mu->cd(1); H("h_mu1_pt") ->Draw("HIST"); gPad->SetGrid();
    c_mu->cd(2); H("h_mu1_eta")->Draw("HIST"); gPad->SetGrid();
    c_mu->cd(3); H("h_mu1_phi")->Draw("HIST"); gPad->SetGrid();

    c_mu->cd(4); H("h_mu2_pt") ->Draw("HIST"); gPad->SetGrid();
    c_mu->cd(5); H("h_mu2_eta")->Draw("HIST"); gPad->SetGrid();
    c_mu->cd(6); H("h_mu2_phi")->Draw("HIST"); gPad->SetGrid();

    c_mu->SaveAs(("muons_pre" + tag + ".png").c_str());

    // --- Jets j1..j6 ---
    TCanvas *c_jet = new TCanvas("c_jet", ("Jets (pre-selection)" + tag).c_str(), 1800, 2200);
    c_jet->Divide(3, 6);

    for (int i = 1; i <= 6; ++i) {
        c_jet->cd(3*(i-1) + 1);
        H(Form("h_j%d_pt",  i))->Draw("HIST"); gPad->SetGrid();

        c_jet->cd(3*(i-1) + 2);
        H(Form("h_j%d_eta", i))->Draw("HIST"); gPad->SetGrid();

        c_jet->cd(3*(i-1) + 3);
        H(Form("h_j%d_phi", i))->Draw("HIST"); gPad->SetGrid();
    }

    c_jet->SaveAs(("jets_pre" + tag + ".png").c_str());

    // --- Multiplicity histograms ---
    TCanvas *c_mult = new TCanvas("c_mult", ("Multiplicity (pre-selection)" + tag).c_str(), 1400, 600);
    c_mult->Divide(3, 1);

    c_mult->cd(1);
    if (H("h_nEle")) {
        H("h_nEle")->GetXaxis()->SetTitle("N electrons");
        H("h_nEle")->Draw("HIST");
    }
    gPad->SetGrid();

    c_mult->cd(2);
    if (H("h_nMu")) {
        H("h_nMu")->GetXaxis()->SetTitle("N muons");
        H("h_nMu")->Draw("HIST");
    }
    gPad->SetGrid();

    c_mult->cd(3);
    if (H("h_nJet")) {
        H("h_nJet")->GetXaxis()->SetTitle("N jets");
        H("h_nJet")->Draw("HIST");
    }
    gPad->SetGrid();

    c_mult->SaveAs(("multiplicity_pre" + tag + ".png").c_str());

    // ============================================================
    // NEW: B-JET PLOTS (RECO)
    // ============================================================

    // --- All b-jets (pt, eta, phi) ---
    TCanvas *c_bjet_all = new TCanvas("c_bjet_all", ("B-Jets (all, pre-selection)" + tag).c_str(), 1600, 900);
    c_bjet_all->Divide(3, 1);

    c_bjet_all->cd(1); H("h_bjet_pt") ->Draw("HIST");  gPad->SetGrid();
    c_bjet_all->cd(2); H("h_bjet_eta")->Draw("HIST");  gPad->SetGrid();
    c_bjet_all->cd(3); H("h_bjet_phi")->Draw("HIST");  gPad->SetGrid();

    c_bjet_all->SaveAs(("bjets_all_pre" + tag + ".png").c_str());

    // --- B-jet multiplicity ---
    TCanvas *c_bjet_mult = new TCanvas("c_bjet_mult", ("B-Jet Multiplicity" + tag).c_str(), 800, 600);
    if (H("h_nBjet")) {
      H("h_nBjet")->GetXaxis()->SetTitle("N b-jets");
      H("h_nBjet")->Draw("HIST");
    }
    gPad->SetGrid();
    c_bjet_mult->SaveAs(("bjet_multiplicity" + tag + ".png").c_str());

    // --- Ranked b-jets (leading → sixth) ---
    TCanvas *c_bjet_rank = new TCanvas("c_bjet_rank", ("Ranked B-Jets (pT-ordered)" + tag).c_str(), 1800, 2200);
    c_bjet_rank->Divide(3, 6);

    for (int i = 1; i <= 6; ++i) {

      c_bjet_rank->cd(3*(i-1) + 1);
      H(Form("h_bj%d_pt",  i))->Draw("HIST"); gPad->SetGrid();
      
      c_bjet_rank->cd(3*(i-1) + 2);
      H(Form("h_bj%d_eta", i))->Draw("HIST"); gPad->SetGrid();

      c_bjet_rank->cd(3*(i-1) + 3);
      H(Form("h_bj%d_phi", i))->Draw("HIST"); gPad->SetGrid();
    }

    c_bjet_rank->SaveAs(("bjets_ranked_pre" + tag + ".png").c_str());

    // ============================================================
    // MET
    // ============================================================

    TCanvas *c_met_pre = new TCanvas("c_met_pre", ("Puppi MET (pre-selection)" + tag).c_str(), 1200, 600);
    c_met_pre->Divide(2, 1);

    c_met_pre->cd(1); if (H("h_pre_MET_pt"))  H("h_pre_MET_pt")->Draw("HIST");
    c_met_pre->cd(2); if (H("h_pre_MET_phi")) H("h_pre_MET_phi")->Draw("HIST");

    gPad->SetGrid();
    c_met_pre->SaveAs(("met_pre" + tag + ".png").c_str());

    std::cout << "\nSaved PNGs (interpreting histograms as RAW; "
                 "if you filled them with w, these are weighted plots).\n";
}
