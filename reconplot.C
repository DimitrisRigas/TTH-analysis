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
        // (UNCHANGED — omitted for brevity)
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

        // ================= DOUBLE-B JETS (4x2) =================
        TCanvas *c_dbk = new TCanvas("c_dbk", "DoubleB jets", 2000, 1000);
        c_dbk->Divide(4, 2);

        c_dbk->cd(1); if (auto h = H("h_dbj1_mass"))      { h->Draw("HIST"); gPad->SetGrid(); }
        c_dbk->cd(2); if (auto h = H("h_dbj1_pt"))        { h->Draw("HIST"); gPad->SetGrid(); }
        c_dbk->cd(3); if (auto h = H("h_dbj1_eta"))       { h->Draw("HIST"); gPad->SetGrid(); }
        c_dbk->cd(4); if (auto h = H("h_dbj1_phi_final")) { h->Draw("HIST"); gPad->SetGrid(); }

        c_dbk->cd(5); if (auto h = H("h_dbj2_mass"))      { h->Draw("HIST"); gPad->SetGrid(); }
        c_dbk->cd(6); if (auto h = H("h_dbj2_pt"))        { h->Draw("HIST"); gPad->SetGrid(); }
        c_dbk->cd(7); if (auto h = H("h_dbj2_eta"))       { h->Draw("HIST"); gPad->SetGrid(); }
        c_dbk->cd(8); if (auto h = H("h_dbj2_phi_final")) { h->Draw("HIST"); gPad->SetGrid(); }

        c_dbk->SaveAs("doubleb_all_kinematics.png");

        // ================= ΔR(db1, db2) =================
        if (auto h = H("h_dR_dbj12_final")) {
            TCanvas *c_dRdb = new TCanvas("c_dRdb", "DeltaR db1 db2", 1000, 800);
            h->Draw("HIST");
            gPad->SetGrid();
            c_dRdb->SaveAs("DeltaR_dbj1_dbj2.png");
        }

        // ================= Δm = |m(db1) − m(db2)| =================
        if (auto h = H("h_dM_bj12_final")) {
            TCanvas *c_dm = new TCanvas("c_dm", "DeltaM double-b jets", 1000, 800);
            h->Draw("HIST");
            gPad->SetGrid();
	    c_dm->SaveAs("DeltaM_dbj1_dbj2.png");
        }

        // ================= FINAL B-JETS (2x2) =================
        TCanvas *c_bj = new TCanvas("c_bj", "Final b-jets", 1200, 800);
        c_bj->Divide(2,2);

        c_bj->cd(1); if (auto h = H("h_bj1_pt_final"))  { h->Draw("HIST"); gPad->SetGrid(); }
        c_bj->cd(2); if (auto h = H("h_bj1_eta_final")) { h->Draw("HIST"); gPad->SetGrid(); }
        c_bj->cd(3); if (auto h = H("h_bj2_pt_final"))  { h->Draw("HIST"); gPad->SetGrid(); }
        c_bj->cd(4); if (auto h = H("h_bj2_eta_final")) { h->Draw("HIST"); gPad->SetGrid(); }

        c_bj->SaveAs("final_bjets.png");

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
