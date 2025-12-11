#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TPad.h>
#include <iostream>

void reconplot(int mode = 1, int sample = 0) {


    gStyle->SetOptStat(1111);

    // =======================================
    // Select input file by sample
    // sample = 0 → signal
    // sample = 1 → ttbar background
    // =======================================

    std::string fname;

    if (sample == 0)
      fname = "output_signal.root";
    else if (sample == 1)
      fname = "output_ttbar.root";
    else {
      std::cout << "ERROR: Unknown sample code! Use 0 (signal) or 1 (ttbar)." << std::endl;
      return;
    }

    TFile *f = TFile::Open(fname.c_str());
    if (!f || f->IsZombie()) {
      std::cout << "ERROR: Could not open " << fname << "!" << std::endl;
      return;
    }

    std::cout << "Loaded file: " << fname << std::endl;


    auto H = [&](std::string name) -> TH1F* {
        TH1F *h = (TH1F*)f->Get(name.c_str());
        if (!h) std::cout << "WARNING: Missing histogram " << name << std::endl;
        return h;
    };

    std::string tag = "_raw";

    // =====================================================================
    // MODE 1 → OBJECT-LEVEL RECO PLOTS (WITH ΔR CLEANING BACK)
    // =====================================================================
    if (mode == 1) {

        // ================= ELECTRONS =================
        TCanvas *c_ele = new TCanvas("c_ele", ("Electrons" + tag).c_str(), 1600, 900);
        c_ele->Divide(3, 2);

        c_ele->cd(1); H("h_e1_pt")->Draw("HIST"); gPad->SetGrid();
        c_ele->cd(2); H("h_e1_eta")->Draw("HIST"); gPad->SetGrid();
        c_ele->cd(3); H("h_e1_phi")->Draw("HIST"); gPad->SetGrid();

        c_ele->cd(4); H("h_e2_pt")->Draw("HIST"); gPad->SetGrid();
        c_ele->cd(5); H("h_e2_eta")->Draw("HIST"); gPad->SetGrid();
        c_ele->cd(6); H("h_e2_phi")->Draw("HIST"); gPad->SetGrid();

        c_ele->SaveAs(("electrons" + tag + ".png").c_str());

        // ================= MUONS =================
        TCanvas *c_mu = new TCanvas("c_mu", ("Muons" + tag).c_str(), 1600, 900);
        c_mu->Divide(3, 2);

        c_mu->cd(1); H("h_mu1_pt")->Draw("HIST"); gPad->SetGrid();
        c_mu->cd(2); H("h_mu1_eta")->Draw("HIST"); gPad->SetGrid();
        c_mu->cd(3); H("h_mu1_phi")->Draw("HIST"); gPad->SetGrid();

        c_mu->cd(4); H("h_mu2_pt")->Draw("HIST"); gPad->SetGrid();
        c_mu->cd(5); H("h_mu2_eta")->Draw("HIST"); gPad->SetGrid();
        c_mu->cd(6); H("h_mu2_phi")->Draw("HIST"); gPad->SetGrid();

        c_mu->SaveAs(("muons" + tag + ".png").c_str());

        // ================= JETS (J1–J4) =================
        TCanvas *c_jet = new TCanvas("c_jet", ("Jets" + tag).c_str(), 1600, 1200);
        c_jet->Divide(3, 4);

        for (int i = 1; i <= 4; ++i) {
            c_jet->cd(3*(i-1) + 1); H(Form("h_j%d_pt",  i))->Draw("HIST"); gPad->SetGrid();
            c_jet->cd(3*(i-1) + 2); H(Form("h_j%d_eta", i))->Draw("HIST"); gPad->SetGrid();
            c_jet->cd(3*(i-1) + 3); H(Form("h_j%d_phi", i))->Draw("HIST"); gPad->SetGrid();
        }

        c_jet->SaveAs(("jets" + tag + ".png").c_str());

        // ================= MULTIPLICITIES =================
        TCanvas *c_mult = new TCanvas("c_mult", ("Multiplicities" + tag).c_str(), 1800, 600);
        c_mult->Divide(5, 1);

        c_mult->cd(1); H("h_nEle")->Draw("HIST");   gPad->SetGrid();
        c_mult->cd(2); H("h_nMu")->Draw("HIST");    gPad->SetGrid();
        c_mult->cd(3); H("h_nJet")->Draw("HIST");   gPad->SetGrid();
        c_mult->cd(4); H("h_nCJet")->Draw("HIST");  gPad->SetGrid();
        c_mult->cd(5); H("h_nBjet")->Draw("HIST");  gPad->SetGrid();

        c_mult->SaveAs(("multiplicity" + tag + ".png").c_str());

        // ================= RANKED B-JETS =================
        TCanvas *c_bjets2 = new TCanvas("c_bjets2", ("2 leading b-jets" + tag).c_str(), 1600, 900);
        c_bjets2->Divide(3, 2);

        c_bjets2->cd(1); H("h_bj1_pt")->Draw("HIST");  gPad->SetGrid();
        c_bjets2->cd(2); H("h_bj1_eta")->Draw("HIST"); gPad->SetGrid();
        c_bjets2->cd(3); H("h_bj1_phi")->Draw("HIST"); gPad->SetGrid();

        c_bjets2->cd(4); H("h_bj2_pt")->Draw("HIST");  gPad->SetGrid();
        c_bjets2->cd(5); H("h_bj2_eta")->Draw("HIST"); gPad->SetGrid();
        c_bjets2->cd(6); H("h_bj2_phi")->Draw("HIST"); gPad->SetGrid();

        c_bjets2->SaveAs(("bjets_ranked2" + tag + ".png").c_str());

        // ================= RANKED DOUBLE-B JETS =================
        TCanvas *c_dbj2 = new TCanvas("c_dbj2", ("2 leading double-b-tagged jets" + tag).c_str(), 1600, 900);
        c_dbj2->Divide(3, 2);

        c_dbj2->cd(1); H("h_dbj1_pt")->Draw("HIST");  gPad->SetGrid();
        c_dbj2->cd(2); H("h_dbj1_eta")->Draw("HIST"); gPad->SetGrid();
        c_dbj2->cd(3); H("h_dbj1_phi")->Draw("HIST"); gPad->SetGrid();

        c_dbj2->cd(4); H("h_dbj2_pt")->Draw("HIST");  gPad->SetGrid();
        c_dbj2->cd(5); H("h_dbj2_eta")->Draw("HIST"); gPad->SetGrid();
        c_dbj2->cd(6); H("h_dbj2_phi")->Draw("HIST"); gPad->SetGrid();

        c_dbj2->SaveAs(("doublebjets_ranked2" + tag + ".png").c_str());

        // ============================================================
        // ✅ ΔR(jet, lepton) BEFORE & AFTER CLEANING (RESTORED)
        // ============================================================
        TCanvas *c_dRclean = new TCanvas("c_dRclean", "Jet–Lepton ΔR (Before/After Cleaning)", 1600, 1200);
        c_dRclean->Divide(2, 2);

        c_dRclean->cd(1); H("h_dR_jet_ele_before")->Draw("HIST"); gPad->SetGrid();
        c_dRclean->cd(2); H("h_dR_jet_mu_before")->Draw("HIST");  gPad->SetGrid();
        c_dRclean->cd(3); H("h_dR_jet_ele_after")->Draw("HIST");  gPad->SetGrid();
        c_dRclean->cd(4); H("h_dR_jet_mu_after")->Draw("HIST");   gPad->SetGrid();

        c_dRclean->SaveAs("dR_jet_lepton_before_after.png");

        std::cout << "\nAll PNGs saved for MODE 1.\n";
        return;
    }

    // =====================================================================
    // MODE 2 → FULL HIGGS + DOUBLE-B RECONSTRUCTION (FINAL)
    // =====================================================================
    if (mode == 2) {

        // ================= HIGGS MASS, pT, η =================
        TCanvas *c_h = new TCanvas("c_h", "Higgs_reco", 1600, 900);
        c_h->Divide(3, 1);

        c_h->cd(1); H("h_Hdbb_mass")->Draw("HIST"); gPad->SetGrid();
        c_h->cd(2); H("h_Hdbb_pt")->Draw("HIST");   gPad->SetGrid();
        c_h->cd(3); H("h_Hdbb_eta")->Draw("HIST");  gPad->SetGrid();

        c_h->SaveAs("Higgs_reco_mass_pt_eta.png");

        // ================= MET =================
        TCanvas *c_met = new TCanvas("c_met", "MET_final", 1200, 600);
        c_met->Divide(2, 1);

        c_met->cd(1); H("h_MET_pt_final")->Draw("HIST");  gPad->SetGrid();
        c_met->cd(2); H("h_MET_phi_final")->Draw("HIST"); gPad->SetGrid();

        c_met->SaveAs("MET_final.png");

        // ============================================================
        // ✅ DOUBLE-B MASS + pT + η (ALL IN ONE CANVAS)
        // ============================================================
        TCanvas *c_dbk = new TCanvas("c_dbk", "DoubleB_kinematics_and_mass", 1800, 1000);
        c_dbk->Divide(3, 2);

        c_dbk->cd(1); H("h_dbj1_mass")->Draw("HIST"); gPad->SetGrid();
        c_dbk->cd(2); H("h_dbj1_pt")->Draw("HIST");   gPad->SetGrid();
        c_dbk->cd(3); H("h_dbj1_eta")->Draw("HIST");  gPad->SetGrid();

        c_dbk->cd(4); H("h_dbj2_mass")->Draw("HIST"); gPad->SetGrid();
        c_dbk->cd(5); H("h_dbj2_pt")->Draw("HIST");   gPad->SetGrid();
        c_dbk->cd(6); H("h_dbj2_eta")->Draw("HIST");  gPad->SetGrid();

        c_dbk->SaveAs("doubleb_mass_pt_eta.png");

        // ================= HT =================
        TCanvas *c_ht = new TCanvas("c_ht", "HT_scalar", 1200, 800);
        H("h_HT")->Draw("HIST");
        gPad->SetGrid();

        c_ht->SaveAs("HT_scalar.png");

        std::cout << "\nAll PNGs saved for MODE 2 (Higgs, MET, double-b mass+pT+eta, HT).\n";
        return;
    }

    std::cout << "\nInvalid mode! Use reconplot(1) or reconplot(2).\n";
}
