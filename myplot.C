#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <iostream>

// Helper to style histograms
void styleHist(TH1F* h, int color) {
    if(!h) return;
    h->SetLineWidth(2);
    h->SetLineColor(color);
    h->GetYaxis()->SetTitleOffset(1.4);
    h->SetStats(0);
}

// Helper to make stat box
TPaveText* makeStat(TH1F* h) {
    if(!h) return nullptr;
    TPaveText *pt = new TPaveText(0.62,0.75,0.88,0.88,"NDC");
    pt->SetFillColor(kWhite);
    pt->SetBorderSize(1);
    pt->SetTextAlign(12);
    pt->AddText(Form("Entries = %.0f", h->GetEntries()));
    pt->AddText(Form("Mean    = %.3f", h->GetMean()));
    pt->AddText(Form("RMS     = %.3f", h->GetRMS()));
    return pt;
}

// Helper to Draw Histo + Stat
void drawWithStat(TH1F* h) {
    if(h) {
        h->Draw();
        TPaveText *st = makeStat(h);
        if(st) st->Draw();
    }
}

void myplot() {

    TFile *f = TFile::Open("output_signal.root");
    if (!f || f->IsZombie()) {
        std::cout << "Error: Cannot open ROOT file." << std::endl;
        return;
    }

    gStyle->SetOptStat(0);

    // --- Retrieve histograms ---

    // Higgs
    TH1F *h_H_pt   = (TH1F*)f->Get("h_H_pt");
    TH1F *h_H_eta  = (TH1F*)f->Get("h_H_eta");
    TH1F *h_H_phi  = (TH1F*)f->Get("h_H_phi");
    TH1F *h_H_mass = (TH1F*)f->Get("h_H_mass");

    // a bosons
    TH1F *h_a1_pt = (TH1F*)f->Get("h_a1_pt");
    TH1F *h_a1_eta = (TH1F*)f->Get("h_a1_eta");
    TH1F *h_a1_phi = (TH1F*)f->Get("h_a1_phi");

    TH1F *h_a2_pt = (TH1F*)f->Get("h_a2_pt");
    TH1F *h_a2_eta = (TH1F*)f->Get("h_a2_eta");
    TH1F *h_a2_phi = (TH1F*)f->Get("h_a2_phi");

    TH1F *h_dR_aa = (TH1F*)f->Get("h_dR_aa");
    TH1F *h_a_mass = (TH1F*)f->Get("h_a_mass");

    // b from a
    TH1F *h_b1_a_pt = (TH1F*)f->Get("h_b1_from_a_pt");
    TH1F *h_b1_a_eta = (TH1F*)f->Get("h_b1_from_a_eta");
    TH1F *h_b1_a_phi = (TH1F*)f->Get("h_b1_from_a_phi");

    TH1F *h_b2_a_pt = (TH1F*)f->Get("h_b2_from_a_pt");
    TH1F *h_b2_a_eta = (TH1F*)f->Get("h_b2_from_a_eta");
    TH1F *h_b2_a_phi = (TH1F*)f->Get("h_b2_from_a_phi");

    TH1F *h_b3_a_pt = (TH1F*)f->Get("h_b3_from_a_pt");
    TH1F *h_b3_a_eta = (TH1F*)f->Get("h_b3_from_a_eta");
    TH1F *h_b3_a_phi = (TH1F*)f->Get("h_b3_from_a_phi");

    TH1F *h_b4_a_pt = (TH1F*)f->Get("h_b4_from_a_pt");
    TH1F *h_b4_a_eta = (TH1F*)f->Get("h_b4_from_a_eta");
    TH1F *h_b4_a_phi = (TH1F*)f->Get("h_b4_from_a_phi");

    // dR(b,b)
    TH1F *h_dR_bb_a1 = (TH1F*)f->Get("h_dR_bb_a1");
    TH1F *h_dR_bb_a2 = (TH1F*)f->Get("h_dR_bb_a2");

    // Tops
    TH1F *h_t1_pt = (TH1F*)f->Get("h_t1_pt");
    TH1F *h_t1_eta = (TH1F*)f->Get("h_t1_eta");
    TH1F *h_t1_phi = (TH1F*)f->Get("h_t1_phi");
    TH1F *h_t1_mass = (TH1F*)f->Get("h_t1_mass");

    TH1F *h_t2_pt = (TH1F*)f->Get("h_t2_pt");
    TH1F *h_t2_eta = (TH1F*)f->Get("h_t2_eta");
    TH1F *h_t2_phi = (TH1F*)f->Get("h_t2_phi");
    TH1F *h_t2_mass = (TH1F*)f->Get("h_t2_mass");

    // b from top
    TH1F *h_b1_t_pt = (TH1F*)f->Get("h_b1_top_pt");
    TH1F *h_b1_t_eta = (TH1F*)f->Get("h_b1_top_eta");
    TH1F *h_b1_t_phi = (TH1F*)f->Get("h_b1_top_phi");

    TH1F *h_b2_t_pt = (TH1F*)f->Get("h_b2_top_pt");
    TH1F *h_b2_t_eta = (TH1F*)f->Get("h_b2_top_eta");
    TH1F *h_b2_t_phi = (TH1F*)f->Get("h_b2_top_phi");

    TH1F *h_dR_Wb = (TH1F*)f->Get("h_dR_Wb");

    // W
    TH1F *h_W_had_pt  = (TH1F*)f->Get("h_W_had_pt");
    TH1F *h_W_had_eta = (TH1F*)f->Get("h_W_had_eta");
    TH1F *h_W_lep_pt  = (TH1F*)f->Get("h_W_lep_pt");
    TH1F *h_W_lep_eta = (TH1F*)f->Get("h_W_lep_eta");

    TH1F *h_q_W_pt = (TH1F*)f->Get("h_q_from_W_pt");
    TH1F *h_lep_W_pt = (TH1F*)f->Get("h_lep_from_W_pt");
    TH1F *h_nu_W_pt = (TH1F*)f->Get("h_nu_from_W_pt");


    // -----------------------------
    //       STYLE HISTOGRAMS
    // -----------------------------
    auto styleAllPT = [&](TH1F* h){
        if(h) h->GetXaxis()->SetRangeUser(0,1000);
    };

    // Apply uniform pT range
    styleAllPT(h_H_pt);
    styleAllPT(h_a1_pt); styleAllPT(h_a2_pt);
    styleAllPT(h_b1_a_pt); styleAllPT(h_b2_a_pt);
    styleAllPT(h_b3_a_pt); styleAllPT(h_b4_a_pt);
    styleAllPT(h_t1_pt);   styleAllPT(h_t2_pt);
    styleAllPT(h_b1_t_pt); styleAllPT(h_b2_t_pt);
    styleAllPT(h_W_had_pt); styleAllPT(h_W_lep_pt);
    styleAllPT(h_q_W_pt);  styleAllPT(h_lep_W_pt); styleAllPT(h_nu_W_pt);

    styleHist(h_H_pt, kBlue); styleHist(h_H_eta, kBlue);
    styleHist(h_H_phi, kBlue); styleHist(h_H_mass, kBlue);

    styleHist(h_a1_pt, kRed); styleHist(h_a1_eta, kRed); styleHist(h_a1_phi, kRed);
    styleHist(h_a2_pt, kRed+2); styleHist(h_a2_eta, kRed+2); styleHist(h_a2_phi, kRed+2);

    styleHist(h_dR_aa, kMagenta+3);
    styleHist(h_dR_bb_a1, kAzure+7); styleHist(h_dR_bb_a2, kAzure+9);

    styleHist(h_t1_pt, kGreen+2); styleHist(h_t1_eta, kGreen+2); styleHist(h_t1_phi, kGreen+2); styleHist(h_t1_mass, kGreen+2);
    styleHist(h_t2_pt, kGreen+3); styleHist(h_t2_eta, kGreen+3); styleHist(h_t2_phi, kGreen+3); styleHist(h_t2_mass, kGreen+3);

    styleHist(h_b1_t_pt, kGreen+1); styleHist(h_b2_t_pt, kGreen+1);
    styleHist(h_dR_Wb, kGreen+4);

    styleHist(h_W_had_pt, kOrange+7); styleHist(h_W_had_eta, kOrange+7);
    styleHist(h_W_lep_pt, kMagenta);  styleHist(h_W_lep_eta, kMagenta);


    // -----------------------------
    //            CANVASES
    // -----------------------------

    // 1. Higgs
    TCanvas *cH = new TCanvas("cH", "Higgs Analysis", 1000, 800);
    cH->Divide(2,2);
    cH->cd(1); drawWithStat(h_H_pt);
    cH->cd(2); drawWithStat(h_H_eta);
    cH->cd(3); drawWithStat(h_H_phi);
    cH->cd(4); drawWithStat(h_H_mass);
    cH->SaveAs("Higgs_Analysis.png");

    // 2. a bosons
    TCanvas *cA = new TCanvas("cA", "a Boson Analysis", 1200, 800);
    cA->Divide(3,2);
    cA->cd(1); drawWithStat(h_a1_pt);
    cA->cd(2); drawWithStat(h_a1_eta);
    cA->cd(3); drawWithStat(h_a1_phi);
    cA->cd(4); drawWithStat(h_a2_pt);
    cA->cd(5); drawWithStat(h_a2_eta);
    cA->cd(6); drawWithStat(h_a2_phi);
    cA->SaveAs("a_Boson_Analysis.png");

    // 3. ΔR Diagrams
    TCanvas *cDR = new TCanvas("cDR", "ΔR Diagrams", 1200, 600);
    cDR->Divide(3,1);
    cDR->cd(1); drawWithStat(h_dR_aa);
    cDR->cd(2); drawWithStat(h_dR_bb_a1);
    cDR->cd(3); drawWithStat(h_dR_bb_a2);
    cDR->SaveAs("DeltaR_Diagrams.png");

    // 4. b-from-a
    TCanvas *cBa = new TCanvas("cBa", "b from a Analysis", 1200, 1000);
    cBa->Divide(3,4);
   
    cBa->cd(1); drawWithStat(h_b1_a_pt);gPad->SetLogy();
    cBa->cd(2); drawWithStat(h_b1_a_eta);
    cBa->cd(3); drawWithStat(h_b1_a_phi);

    cBa->cd(4); drawWithStat(h_b2_a_pt);
    cBa->cd(5); drawWithStat(h_b2_a_eta);
    cBa->cd(6); drawWithStat(h_b2_a_phi);

    cBa->cd(7); drawWithStat(h_b3_a_pt);
    cBa->cd(8); drawWithStat(h_b3_a_eta);
    cBa->cd(9); drawWithStat(h_b3_a_phi);

    cBa->cd(10); drawWithStat(h_b4_a_pt);
    cBa->cd(11); drawWithStat(h_b4_a_eta);
    cBa->cd(12); drawWithStat(h_b4_a_phi);
    cBa->SaveAs("b_from_a_Analysis.png");

    // 5. Tops
    TCanvas *cTop = new TCanvas("cTop", "Top Analysis", 1200, 800);
    cTop->Divide(3,2);
    cTop->cd(1); drawWithStat(h_t1_pt);
    cTop->cd(2); drawWithStat(h_t1_eta);
    cTop->cd(3); drawWithStat(h_t1_phi);
    cTop->cd(4); drawWithStat(h_t2_pt);
    cTop->cd(5); drawWithStat(h_t2_eta);
    cTop->cd(6); drawWithStat(h_t2_phi);
    cTop->SaveAs("Top_Analysis.png");

    // 6. Top masses
    TCanvas *cTm = new TCanvas("cTm", "Top Masses", 800, 600);
    cTm->Divide(2,1);
    cTm->cd(1); drawWithStat(h_t1_mass);
    cTm->cd(2); drawWithStat(h_t2_mass);
    cTm->SaveAs("Top_Masses.png");

    // 7. b-from-top
    TCanvas *cBt = new TCanvas("cBt", "b from Top Analysis", 1200, 600);
    cBt->Divide(3,2);
    cBt->cd(1); drawWithStat(h_b1_t_pt);
    cBt->cd(2); drawWithStat(h_b1_t_eta);
    cBt->cd(3); drawWithStat(h_b1_t_phi);
    cBt->cd(4); drawWithStat(h_b2_t_pt);
    cBt->cd(5); drawWithStat(h_b2_t_eta);
    cBt->cd(6); drawWithStat(h_dR_Wb);
    cBt->SaveAs("b_from_top_Analysis.png");

    // 8. NEW W-from-top canvas
    TCanvas *cWTop = new TCanvas("cWTop", "W from Top Analysis", 1200, 800);
    cWTop->Divide(2,2);
    cWTop->cd(1); drawWithStat(h_W_had_pt);
    cWTop->cd(2); drawWithStat(h_W_had_eta);
    cWTop->cd(3); drawWithStat(h_W_lep_pt);
    cWTop->cd(4); drawWithStat(h_W_lep_eta);
    cWTop->SaveAs("W_from_Top_Analysis.png");

    std::cout << "All canvases saved!" << std::endl;
}
