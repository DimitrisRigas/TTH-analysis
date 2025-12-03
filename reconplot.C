#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TPad.h>
#include <iostream>

void reconplot() {

    gStyle->SetOptStat(1111);

    TFile *f = new TFile("ttH_aab_final_output.root");
    if (!f || f->IsZombie()) {
        std::cout << "ERROR: Could not open ttH_aab_final_output.root!" << std::endl;
        return;
    }

    // Fetch helper
    auto H = [&](std::string name) -> TH1F* {
        TH1F *h = (TH1F*)f->Get(name.c_str());
        if (!h) std::cout << "WARNING: Missing histogram " << name << std::endl;
        return h;
    };

    std::string tag = "_raw";

    // ============================================================
    // ELECTRONS
    // ============================================================

    TCanvas *c_ele = new TCanvas("c_ele", ("Electrons" + tag).c_str(), 1600, 900);
    c_ele->Divide(3, 2);

    c_ele->cd(1); H("h_e1_pt")->Draw("HIST"); gPad->SetGrid();
    c_ele->cd(2); H("h_e1_eta")->Draw("HIST"); gPad->SetGrid();
    c_ele->cd(3); H("h_e1_phi")->Draw("HIST"); gPad->SetGrid();

    c_ele->cd(4); H("h_e2_pt")->Draw("HIST"); gPad->SetGrid();
    c_ele->cd(5); H("h_e2_eta")->Draw("HIST"); gPad->SetGrid();
    c_ele->cd(6); H("h_e2_phi")->Draw("HIST"); gPad->SetGrid();

    c_ele->SaveAs(("electrons" + tag + ".png").c_str());

    // ============================================================
    // MUONS
    // ============================================================

    TCanvas *c_mu = new TCanvas("c_mu", ("Muons" + tag).c_str(), 1600, 900);
    c_mu->Divide(3, 2);

    c_mu->cd(1); H("h_mu1_pt")->Draw("HIST"); gPad->SetGrid();
    c_mu->cd(2); H("h_mu1_eta")->Draw("HIST"); gPad->SetGrid();
    c_mu->cd(3); H("h_mu1_phi")->Draw("HIST"); gPad->SetGrid();

    c_mu->cd(4); H("h_mu2_pt")->Draw("HIST"); gPad->SetGrid();
    c_mu->cd(5); H("h_mu2_eta")->Draw("HIST"); gPad->SetGrid();
    c_mu->cd(6); H("h_mu2_phi")->Draw("HIST"); gPad->SetGrid();

    c_mu->SaveAs(("muons" + tag + ".png").c_str());

    // ============================================================
    // JETS (ONLY J1–J4)
    // ============================================================

    TCanvas *c_jet = new TCanvas("c_jet", ("Jets" + tag).c_str(), 1600, 1200);
    c_jet->Divide(3, 4);

    for (int i = 1; i <= 4; ++i) {
        c_jet->cd(3*(i-1) + 1); H(Form("h_j%d_pt",  i))->Draw("HIST"); gPad->SetGrid();
        c_jet->cd(3*(i-1) + 2); H(Form("h_j%d_eta", i))->Draw("HIST"); gPad->SetGrid();
        c_jet->cd(3*(i-1) + 3); H(Form("h_j%d_phi", i))->Draw("HIST"); gPad->SetGrid();
    }

    c_jet->SaveAs(("jets" + tag + ".png").c_str());

    // ============================================================
    // MULTIPLICITIES
    // ============================================================

    TCanvas *c_mult = new TCanvas("c_mult", ("Multiplicities" + tag).c_str(), 1600, 600);
    c_mult->Divide(4, 1);

    c_mult->cd(1); if (H("h_nEle"))  H("h_nEle")->Draw("HIST");  gPad->SetGrid();
    c_mult->cd(2); if (H("h_nMu"))   H("h_nMu")->Draw("HIST");   gPad->SetGrid();
    c_mult->cd(3); if (H("h_nJet"))  H("h_nJet")->Draw("HIST");  gPad->SetGrid();
    c_mult->cd(4); if (H("h_nBjet")) H("h_nBjet")->Draw("HIST"); gPad->SetGrid();

    c_mult->SaveAs(("multiplicity" + tag + ".png").c_str());

    // ============================================================
    // ALL B-JETS (pt,eta,phi)
    // ============================================================

    TCanvas *c_bjet_all = new TCanvas("c_bjet_all", ("All b-jets" + tag).c_str(), 1600, 900);
    c_bjet_all->Divide(3, 1);

    c_bjet_all->cd(1); H("h_bjet_pt")->Draw("HIST");  gPad->SetGrid();
    c_bjet_all->cd(2); H("h_bjet_eta")->Draw("HIST"); gPad->SetGrid();
    c_bjet_all->cd(3); H("h_bjet_phi")->Draw("HIST"); gPad->SetGrid();

    c_bjet_all->SaveAs(("bjets_all" + tag + ".png").c_str());

    // ============================================================
    // RANKED B-JETS (bj1, bj2)
    // ============================================================

    TCanvas *c_bjets2 = new TCanvas("c_bjets2", ("2 leading b-jets" + tag).c_str(), 1600, 900);
    c_bjets2->Divide(3, 2);

    c_bjets2->cd(1); H("h_bj1_pt")->Draw("HIST");  gPad->SetGrid();
    c_bjets2->cd(2); H("h_bj1_eta")->Draw("HIST"); gPad->SetGrid();
    c_bjets2->cd(3); H("h_bj1_phi")->Draw("HIST"); gPad->SetGrid();

    c_bjets2->cd(4); H("h_bj2_pt")->Draw("HIST");  gPad->SetGrid();
    c_bjets2->cd(5); H("h_bj2_eta")->Draw("HIST"); gPad->SetGrid();
    c_bjets2->cd(6); H("h_bj2_phi")->Draw("HIST"); gPad->SetGrid();

    c_bjets2->SaveAs(("bjets_ranked2" + tag + ".png").c_str());

    // ============================================================
    // NEW: ALL DOUBLE B-TAGGED JETS (pt,eta,phi)
    // ============================================================

    TCanvas *c_dbj_all = new TCanvas("c_dbj_all", ("All double-b-tagged jets" + tag).c_str(), 1600, 900);
    c_dbj_all->Divide(3, 1);

    c_dbj_all->cd(1); H("h_dbjet_pt")->Draw("HIST");  gPad->SetGrid();
    c_dbj_all->cd(2); H("h_dbjet_eta")->Draw("HIST"); gPad->SetGrid();
    c_dbj_all->cd(3); H("h_dbjet_phi")->Draw("HIST"); gPad->SetGrid();

    c_dbj_all->SaveAs(("doublebjets_all" + tag + ".png").c_str());

    // ============================================================
    // NEW: RANKED DOUBLE-B-JETS (dbj1, dbj2)
    // ============================================================

    TCanvas *c_dbj2 = new TCanvas("c_dbj2", ("2 leading double-b-tagged jets" + tag).c_str(), 1600, 900);
    c_dbj2->Divide(3, 2);

    c_dbj2->cd(1); H("h_dbj1_pt")->Draw("HIST");  gPad->SetGrid();
    c_dbj2->cd(2); H("h_dbj1_eta")->Draw("HIST"); gPad->SetGrid();
    c_dbj2->cd(3); H("h_dbj1_phi")->Draw("HIST"); gPad->SetGrid();

    c_dbj2->cd(4); H("h_dbj2_pt")->Draw("HIST");  gPad->SetGrid();
    c_dbj2->cd(5); H("h_dbj2_eta")->Draw("HIST"); gPad->SetGrid();
    c_dbj2->cd(6); H("h_dbj2_phi")->Draw("HIST"); gPad->SetGrid();

    c_dbj2->SaveAs(("doublebjets_ranked2" + tag + ".png").c_str());

    std::cout << "\nAll PNGs saved.\n";
}
