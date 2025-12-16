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
    // =======================================
    std::string fname;

    if (sample == 0)
        fname = "output_signal.root";
    else if (sample == 1)
        fname = "output_ttbar.root";
    else {
        std::cout << "ERROR: Unknown sample code!" << std::endl;
        return;
    }

    TFile *f = TFile::Open(fname.c_str());
    if (!f || f->IsZombie()) {
        std::cout << "ERROR: Could not open " << fname << std::endl;
        return;
    }

    std::cout << "Loaded file: " << fname << std::endl;

    auto H = [&](const char *name) -> TH1F* {
        TH1F *h = (TH1F*)f->Get(name);
        if (!h) std::cout << "WARNING: Missing histogram " << name << std::endl;
        return h;
    };

    std::string tag = "_raw";

    // =====================================================================
    // MODE 1 → OBJECT-LEVEL RECO
    // =====================================================================
    if (mode == 1) {

        // ================= ELECTRONS =================
        TCanvas *c_ele = new TCanvas("c_ele", "Electrons", 1600, 900);
        c_ele->Divide(3, 2);

        c_ele->cd(1); if (auto h = H("h_e1_pt"))  { h->Draw("HIST"); gPad->SetGrid(); }
        c_ele->cd(2); if (auto h = H("h_e1_eta")) { h->Draw("HIST"); gPad->SetGrid(); }
        c_ele->cd(3); if (auto h = H("h_e1_phi")) { h->Draw("HIST"); gPad->SetGrid(); }
        c_ele->cd(4); if (auto h = H("h_e2_pt"))  { h->Draw("HIST"); gPad->SetGrid(); }
        c_ele->cd(5); if (auto h = H("h_e2_eta")) { h->Draw("HIST"); gPad->SetGrid(); }
        c_ele->cd(6); if (auto h = H("h_e2_phi")) { h->Draw("HIST"); gPad->SetGrid(); }

        c_ele->SaveAs(("electrons" + tag + ".png").c_str());

        // ================= MUONS =================
        TCanvas *c_mu = new TCanvas("c_mu", "Muons", 1600, 900);
        c_mu->Divide(3, 2);

        c_mu->cd(1); if (auto h = H("h_mu1_pt"))  { h->Draw("HIST"); gPad->SetGrid(); }
        c_mu->cd(2); if (auto h = H("h_mu1_eta")) { h->Draw("HIST"); gPad->SetGrid(); }
        c_mu->cd(3); if (auto h = H("h_mu1_phi")) { h->Draw("HIST"); gPad->SetGrid(); }
        c_mu->cd(4); if (auto h = H("h_mu2_pt"))  { h->Draw("HIST"); gPad->SetGrid(); }
        c_mu->cd(5); if (auto h = H("h_mu2_eta")) { h->Draw("HIST"); gPad->SetGrid(); }
        c_mu->cd(6); if (auto h = H("h_mu2_phi")) { h->Draw("HIST"); gPad->SetGrid(); }

        c_mu->SaveAs(("muons" + tag + ".png").c_str());

        // ================= JETS =================
        TCanvas *c_jet = new TCanvas("c_jet", "Jets", 1600, 1200);
        c_jet->Divide(3, 4);

        for (int i = 1; i <= 4; ++i) {
            c_jet->cd(3*(i-1)+1); if (auto h = H(Form("h_j%d_pt",i)))  { h->Draw("HIST"); gPad->SetGrid(); }
            c_jet->cd(3*(i-1)+2); if (auto h = H(Form("h_j%d_eta",i))) { h->Draw("HIST"); gPad->SetGrid(); }
            c_jet->cd(3*(i-1)+3); if (auto h = H(Form("h_j%d_phi",i))) { h->Draw("HIST"); gPad->SetGrid(); }
        }

        c_jet->SaveAs(("jets" + tag + ".png").c_str());

        // ================= ΔR CLEANING =================
        TCanvas *c_dR = new TCanvas("c_dRclean", "Jet–Lepton ΔR", 1600, 1200);
        c_dR->Divide(2, 2);

        c_dR->cd(1); if (auto h = H("h_dR_jet_ele_before")) { h->Draw("HIST"); gPad->SetGrid(); }
        c_dR->cd(2); if (auto h = H("h_dR_jet_mu_before"))  { h->Draw("HIST"); gPad->SetGrid(); }
        c_dR->cd(3); if (auto h = H("h_dR_jet_ele_after"))  { h->Draw("HIST"); gPad->SetGrid(); }
        c_dR->cd(4); if (auto h = H("h_dR_jet_mu_after"))   { h->Draw("HIST"); gPad->SetGrid(); }

        c_dR->SaveAs("dR_jet_lepton_before_after.png");

        return;
    }

    // =====================================================================
    // MODE 2 → FINAL RECONSTRUCTION
    // =====================================================================
    if (mode == 2) {

        // ================= HIGGS =================
        TCanvas *c_h = new TCanvas("c_h", "Higgs", 1600, 900);
        c_h->Divide(3,1);

        c_h->cd(1); if (auto h = H("h_Hdbb_mass")) { h->Draw("HIST"); gPad->SetGrid(); }
        c_h->cd(2); if (auto h = H("h_Hdbb_pt"))   { h->Draw("HIST"); gPad->SetGrid(); }
        c_h->cd(3); if (auto h = H("h_Hdbb_eta"))  { h->Draw("HIST"); gPad->SetGrid(); }

        c_h->SaveAs("Higgs_reco.png");
	// ================= DOUBLE-B JETS (ALL KINEMATICS, 4x2) =================
	TCanvas *c_dbk = new TCanvas("c_dbk", "DoubleB jets", 2000, 1000);
	c_dbk->Divide(4, 2);

	// --- double-b jet 1 ---
	c_dbk->cd(1); if (auto h = H("h_dbj1_mass"))      { h->Draw("HIST"); gPad->SetGrid(); }
	c_dbk->cd(2); if (auto h = H("h_dbj1_pt"))        { h->Draw("HIST"); gPad->SetGrid(); }
	c_dbk->cd(3); if (auto h = H("h_dbj1_eta"))       { h->Draw("HIST"); gPad->SetGrid(); }
	c_dbk->cd(4); if (auto h = H("h_dbj1_phi_final")) { h->Draw("HIST"); gPad->SetGrid(); }

	// --- double-b jet 2 ---
	c_dbk->cd(5); if (auto h = H("h_dbj2_mass"))      { h->Draw("HIST"); gPad->SetGrid(); }
	c_dbk->cd(6); if (auto h = H("h_dbj2_pt"))        { h->Draw("HIST"); gPad->SetGrid(); }
	c_dbk->cd(7); if (auto h = H("h_dbj2_eta"))       { h->Draw("HIST"); gPad->SetGrid(); }
	c_dbk->cd(8); if (auto h = H("h_dbj2_phi_final")) { h->Draw("HIST"); gPad->SetGrid(); }

	c_dbk->SaveAs("doubleb_all_kinematics.png");


        // ================= ΔR(db1,db2) =================
        if (auto h = H("h_dR_dbj12_final")) {
            TCanvas *c_dRdb = new TCanvas("c_dRdb", "DeltaR db1 db2", 1000, 800);
            h->Draw("HIST");
            gPad->SetGrid();
            c_dRdb->SaveAs("DeltaR_dbj1_dbj2.png");
        }

        // ================= FINAL LEPTONS =================
        TCanvas *c_lep = new TCanvas("c_lep", "Final leptons", 1600, 900);
        c_lep->Divide(3,2);

        c_lep->cd(1); if (auto h = H("h_lep1_pt_final"))  { h->Draw("HIST"); gPad->SetGrid(); }
        c_lep->cd(2); if (auto h = H("h_lep1_eta_final")) { h->Draw("HIST"); gPad->SetGrid(); }
        c_lep->cd(3); if (auto h = H("h_lep1_phi_final")) { h->Draw("HIST"); gPad->SetGrid(); }
        c_lep->cd(4); if (auto h = H("h_lep2_pt_final"))  { h->Draw("HIST"); gPad->SetGrid(); }
        c_lep->cd(5); if (auto h = H("h_lep2_eta_final")) { h->Draw("HIST"); gPad->SetGrid(); }
        c_lep->cd(6); if (auto h = H("h_lep2_phi_final")) { h->Draw("HIST"); gPad->SetGrid(); }

        c_lep->SaveAs("final_leptons.png");

        return;
    }

    std::cout << "Invalid mode." << std::endl;
}
