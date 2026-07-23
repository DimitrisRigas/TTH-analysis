// myplot.C
//
// Generator-level plotting macro.
//
// This version does TWO things:
//
// 1) Keeps ALL generator-level plots from your previous macro.
//    These are saved in:
//
//      gen_level_plots/<signalTag>/
//      gen_level_plots/overlay_tth12_vs_tth30/
//
// 2) Also creates a NEW reduced teacher-required set only with:
//
//      pT, eta of Higgs
//      pT, eta of one A-boson
//      pT, eta of the 4 b-quarks from H -> aa -> 4b
//      pT, eta of one representative top quark
//      pT, eta of one representative b-quark from top decay
//      DeltaR(A,A)
//      DeltaR(b,b) from one representative A-boson
//
//    These are saved in:
//
//      gen_level_plots/teacher_required/<signalTag>/
//      gen_level_plots/teacher_required/overlay_tth12_vs_tth30/
//
// Why 18 plots?
// - a1 and a2 are physically identical pseudoscalar bosons, so only a1 is shown.
// - t1 and t2 are top/anti-top from the same ttbar system, so one representative top is shown.
// - b1_top and b2_top are the two b-quarks from top decays, so one representative b-from-top is shown.
// - dR_bb_a1 and dR_bb_a2 describe the same a -> bb decay topology, so only one is shown.
// - All four b-quarks from H -> aa -> 4b are kept because the teacher explicitly asked for 4-b quarks.
//
// Run full old + teacher-required:
//   root -l
//   .L myplot.C
//   myplot()
//
// Run only teacher-required plots:
//   root -l
//   .L myplot.C
//   myplotTeacherRequired()
//
// Run only one signal:
//   myplotOne("tth12gev")
//
// Run only overlays:
//   myplotOverlay12vs30()

#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPad.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TString.h>
#include <TRandom.h>

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

// ============================================================================
// Histogram configuration
// ============================================================================
struct HistConfig {
    std::string histName;
    std::string outName;
    std::string xTitle;
    std::string yTitle;
    double xmin;
    double xmax;
    bool useRange;
};

// ============================================================================
// Configuration for all generator-level histograms
// ============================================================================
HistConfig getConfig(const std::string& hname) {
    HistConfig c;
    c.histName = hname;
    c.outName = hname;
    c.xTitle = "";
    c.yTitle = "Entries";
    c.xmin = 0.0;
    c.xmax = 0.0;
    c.useRange = false;

    // ------------------------------------------------------------------------
    // Higgs
    // ------------------------------------------------------------------------
    if (hname == "h_H_pt") {
        c.outName = "H_pt";
        c.xTitle = "p_{T}(H) [GeV]";
        c.xmin = 0.0;
        c.xmax = 500.0;
        c.useRange = true;
    }
    else if (hname == "h_H_eta") {
        c.outName = "H_eta";
        c.xTitle = "#eta(H)";
        c.xmin = -5.0;
        c.xmax = 5.0;
        c.useRange = true;
    }
    else if (hname == "h_H_phi") {
        c.outName = "H_phi";
        c.xTitle = "#phi(H)";
        c.xmin = -3.2;
        c.xmax = 3.2;
        c.useRange = true;
    }
    else if (hname == "h_H_mass") {
        c.outName = "H_mass";
        c.xTitle = "m(H) [GeV]";
        c.xmin = 100.0;
        c.xmax = 150.0;
        c.useRange = true;
    }

    // ------------------------------------------------------------------------
    // A bosons
    // ------------------------------------------------------------------------
    else if (hname == "h_a1_pt") {
        c.outName = "a1_pt";
        c.xTitle = "p_{T}(a) [GeV]";
        c.xmin = 0.0;
        c.xmax = 500.0;
        c.useRange = true;
    }
    else if (hname == "h_a1_eta") {
        c.outName = "a1_eta";
        c.xTitle = "#eta(a)";
        c.xmin = -5.0;
        c.xmax = 5.0;
        c.useRange = true;
    }
    else if (hname == "h_a1_phi") {
        c.outName = "a1_phi";
        c.xTitle = "#phi(a_{1})";
        c.xmin = -3.2;
        c.xmax = 3.2;
        c.useRange = true;
    }
    else if (hname == "h_a2_pt") {
        c.outName = "a2_pt";
        c.xTitle = "p_{T}(a_{2}) [GeV]";
        c.xmin = 0.0;
        c.xmax = 500.0;
        c.useRange = true;
    }
    else if (hname == "h_a2_eta") {
        c.outName = "a2_eta";
        c.xTitle = "#eta(a_{2})";
        c.xmin = -5.0;
        c.xmax = 5.0;
        c.useRange = true;
    }
    else if (hname == "h_a2_phi") {
        c.outName = "a2_phi";
        c.xTitle = "#phi(a_{2})";
        c.xmin = -3.2;
        c.xmax = 3.2;
        c.useRange = true;
    }
    else if (hname == "h_a_mass") {
        c.outName = "a_mass";
        c.xTitle = "m(a) [GeV]";
        c.xmin = 0.0;
        c.xmax = 70.0;
        c.useRange = true;
    }
    else if (hname == "h_dR_aa") {
        c.outName = "dR_aa";
        c.xTitle = "#DeltaR(a_{1}, a_{2})";
        c.xmin = 0.0;
        c.xmax = 6.0;
        c.useRange = true;
    }

    // ------------------------------------------------------------------------
    // b quarks from A bosons
    // ------------------------------------------------------------------------
    else if (hname == "h_b1_from_a_pt") {
        c.outName = "b1_from_a_pt";
        c.xTitle = "p_{T}(b_{1} from a) [GeV]";
        c.xmin = 0.0;
        c.xmax = 400.0;
        c.useRange = true;
    }
    else if (hname == "h_b2_from_a_pt") {
        c.outName = "b2_from_a_pt";
        c.xTitle = "p_{T}(b_{2} from a) [GeV]";
        c.xmin = 0.0;
        c.xmax = 400.0;
        c.useRange = true;
    }
    else if (hname == "h_b3_from_a_pt") {
        c.outName = "b3_from_a_pt";
        c.xTitle = "p_{T}(b_{3} from a) [GeV]";
        c.xmin = 0.0;
        c.xmax = 400.0;
        c.useRange = true;
    }
    else if (hname == "h_b4_from_a_pt") {
        c.outName = "b4_from_a_pt";
        c.xTitle = "p_{T}(b_{4} from a) [GeV]";
        c.xmin = 0.0;
        c.xmax = 400.0;
        c.useRange = true;
    }
    else if (hname == "h_b1_from_a_eta") {
        c.outName = "b1_from_a_eta";
        c.xTitle = "#eta(b_{1} from a)";
        c.xmin = -5.0;
        c.xmax = 5.0;
        c.useRange = true;
    }
    else if (hname == "h_b2_from_a_eta") {
        c.outName = "b2_from_a_eta";
        c.xTitle = "#eta(b_{2} from a)";
        c.xmin = -5.0;
        c.xmax = 5.0;
        c.useRange = true;
    }
    else if (hname == "h_b3_from_a_eta") {
        c.outName = "b3_from_a_eta";
        c.xTitle = "#eta(b_{3} from a)";
        c.xmin = -5.0;
        c.xmax = 5.0;
        c.useRange = true;
    }
    else if (hname == "h_b4_from_a_eta") {
        c.outName = "b4_from_a_eta";
        c.xTitle = "#eta(b_{4} from a)";
        c.xmin = -5.0;
        c.xmax = 5.0;
        c.useRange = true;
    }
    else if (hname == "h_b1_from_a_phi") {
        c.outName = "b1_from_a_phi";
        c.xTitle = "#phi(b_{1} from a)";
        c.xmin = -3.2;
        c.xmax = 3.2;
        c.useRange = true;
    }
    else if (hname == "h_b2_from_a_phi") {
        c.outName = "b2_from_a_phi";
        c.xTitle = "#phi(b_{2} from a)";
        c.xmin = -3.2;
        c.xmax = 3.2;
        c.useRange = true;
    }
    else if (hname == "h_b3_from_a_phi") {
        c.outName = "b3_from_a_phi";
        c.xTitle = "#phi(b_{3} from a)";
        c.xmin = -3.2;
        c.xmax = 3.2;
        c.useRange = true;
    }
    else if (hname == "h_b4_from_a_phi") {
        c.outName = "b4_from_a_phi";
        c.xTitle = "#phi(b_{4} from a)";
        c.xmin = -3.2;
        c.xmax = 3.2;
        c.useRange = true;
    }
    else if (hname == "h_dR_bb_a1") {
        c.outName = "dR_bb_a1";
        c.xTitle = "#DeltaR(b,b) from a";
        c.xmin = 0.0;
        c.xmax = 5.0;
        c.useRange = true;
    }
    else if (hname == "h_dR_bb_a2") {
        c.outName = "dR_bb_a2";
        c.xTitle = "#DeltaR(b,b) from a_{2}";
        c.xmin = 0.0;
        c.xmax = 5.0;
        c.useRange = true;
    }

    // ------------------------------------------------------------------------
    // Top quarks
    // ------------------------------------------------------------------------
    else if (hname == "h_t1_pt") {
        c.outName = "t1_pt";
        c.xTitle = "p_{T}(t) [GeV]";
        c.xmin = 0.0;
        c.xmax = 500.0;
        c.useRange = true;
    }
    else if (hname == "h_t1_eta") {
        c.outName = "t1_eta";
        c.xTitle = "#eta(t)";
        c.xmin = -5.0;
        c.xmax = 5.0;
        c.useRange = true;
    }
    else if (hname == "h_t1_phi") {
        c.outName = "t1_phi";
        c.xTitle = "#phi(t_{1})";
        c.xmin = -3.2;
        c.xmax = 3.2;
        c.useRange = true;
    }
    else if (hname == "h_t1_mass") {
        c.outName = "t1_mass";
        c.xTitle = "m(t_{1}) [GeV]";
        c.xmin = 140.0;
        c.xmax = 220.0;
        c.useRange = true;
    }
    else if (hname == "h_t2_pt") {
        c.outName = "t2_pt";
        c.xTitle = "p_{T}(t_{2}) [GeV]";
        c.xmin = 0.0;
        c.xmax = 500.0;
        c.useRange = true;
    }
    else if (hname == "h_t2_eta") {
        c.outName = "t2_eta";
        c.xTitle = "#eta(t_{2})";
        c.xmin = -5.0;
        c.xmax = 5.0;
        c.useRange = true;
    }
    else if (hname == "h_t2_phi") {
        c.outName = "t2_phi";
        c.xTitle = "#phi(t_{2})";
        c.xmin = -3.2;
        c.xmax = 3.2;
        c.useRange = true;
    }
    else if (hname == "h_t2_mass") {
        c.outName = "t2_mass";
        c.xTitle = "m(t_{2}) [GeV]";
        c.xmin = 140.0;
        c.xmax = 220.0;
        c.useRange = true;
    }

    // ------------------------------------------------------------------------
    // b quarks from top
    // ------------------------------------------------------------------------
    else if (hname == "h_b1_top_pt") {
        c.outName = "b1_top_pt";
        c.xTitle = "p_{T}(b from top) [GeV]";
        c.xmin = 0.0;
        c.xmax = 400.0;
        c.useRange = true;
    }
    else if (hname == "h_b1_top_eta") {
        c.outName = "b1_top_eta";
        c.xTitle = "#eta(b from top)";
        c.xmin = -5.0;
        c.xmax = 5.0;
        c.useRange = true;
    }
    else if (hname == "h_b1_top_phi") {
        c.outName = "b1_top_phi";
        c.xTitle = "#phi(b_{1} from top)";
        c.xmin = -3.2;
        c.xmax = 3.2;
        c.useRange = true;
    }
    else if (hname == "h_b2_top_pt") {
        c.outName = "b2_top_pt";
        c.xTitle = "p_{T}(b_{2} from top) [GeV]";
        c.xmin = 0.0;
        c.xmax = 400.0;
        c.useRange = true;
    }
    else if (hname == "h_b2_top_eta") {
        c.outName = "b2_top_eta";
        c.xTitle = "#eta(b_{2} from top)";
        c.xmin = -5.0;
        c.xmax = 5.0;
        c.useRange = true;
    }
    else if (hname == "h_b2_top_phi") {
        c.outName = "b2_top_phi";
        c.xTitle = "#phi(b_{2} from top)";
        c.xmin = -3.2;
        c.xmax = 3.2;
        c.useRange = true;
    }
    else if (hname == "h_dR_Wb") {
        c.outName = "dR_Wb";
        c.xTitle = "#DeltaR(W,b)";
        c.xmin = 0.0;
        c.xmax = 6.0;
        c.useRange = true;
    }

    // ------------------------------------------------------------------------
    // W bosons
    // ------------------------------------------------------------------------
    else if (hname == "h_W_had_pt") {
        c.outName = "W_had_pt";
        c.xTitle = "p_{T}(W_{had}) [GeV]";
        c.xmin = 0.0;
        c.xmax = 500.0;
        c.useRange = true;
    }
    else if (hname == "h_W_had_eta") {
        c.outName = "W_had_eta";
        c.xTitle = "#eta(W_{had})";
        c.xmin = -5.0;
        c.xmax = 5.0;
        c.useRange = true;
    }
    else if (hname == "h_W_lep_pt") {
        c.outName = "W_lep_pt";
        c.xTitle = "p_{T}(W_{lep}) [GeV]";
        c.xmin = 0.0;
        c.xmax = 500.0;
        c.useRange = true;
    }
    else if (hname == "h_W_lep_eta") {
        c.outName = "W_lep_eta";
        c.xTitle = "#eta(W_{lep})";
        c.xmin = -5.0;
        c.xmax = 5.0;
        c.useRange = true;
    }

    // ------------------------------------------------------------------------
    // W decay products
    // ------------------------------------------------------------------------
    else if (hname == "h_q_from_W_pt") {
        c.outName = "q_from_W_pt";
        c.xTitle = "p_{T}(q from W) [GeV]";
        c.xmin = 0.0;
        c.xmax = 400.0;
        c.useRange = true;
    }
    else if (hname == "h_lep_from_W_pt") {
        c.outName = "lep_from_W_pt";
        c.xTitle = "p_{T}(lepton from W) [GeV]";
        c.xmin = 0.0;
        c.xmax = 400.0;
        c.useRange = true;
    }
    else if (hname == "h_nu_from_W_pt") {
        c.outName = "nu_from_W_pt";
        c.xTitle = "p_{T}(#nu from W) [GeV]";
        c.xmin = 0.0;
        c.xmax = 400.0;
        c.useRange = true;
    }

    return c;
}

// ============================================================================
// All generator-level histograms from your previous macro
// ============================================================================
std::vector<std::string> allGeneratorHists() {
    return {
        // Higgs
        "h_H_pt",
        "h_H_eta",
        "h_H_phi",
        "h_H_mass",

        // A bosons
        "h_a1_pt",
        "h_a1_eta",
        "h_a1_phi",
        "h_a2_pt",
        "h_a2_eta",
        "h_a2_phi",
        "h_a_mass",
        "h_dR_aa",

        // b quarks from A bosons
        "h_b1_from_a_pt",
        "h_b1_from_a_eta",
        "h_b1_from_a_phi",
        "h_b2_from_a_pt",
        "h_b2_from_a_eta",
        "h_b2_from_a_phi",
        "h_b3_from_a_pt",
        "h_b3_from_a_eta",
        "h_b3_from_a_phi",
        "h_b4_from_a_pt",
        "h_b4_from_a_eta",
        "h_b4_from_a_phi",

        // DeltaR bb
        "h_dR_bb_a1",
        "h_dR_bb_a2",

        // Tops
        "h_t1_pt",
        "h_t1_eta",
        "h_t1_phi",
        "h_t1_mass",
        "h_t2_pt",
        "h_t2_eta",
        "h_t2_phi",
        "h_t2_mass",

        // b from top
        "h_b1_top_pt",
        "h_b1_top_eta",
        "h_b1_top_phi",
        "h_b2_top_pt",
        "h_b2_top_eta",
        "h_b2_top_phi",
        "h_dR_Wb",

        // W bosons
        "h_W_had_pt",
        "h_W_had_eta",
        "h_W_lep_pt",
        "h_W_lep_eta",

        // W decay products
        "h_q_from_W_pt",
        "h_lep_from_W_pt",
        "h_nu_from_W_pt"
    };
}

// ============================================================================
// Reduced teacher-required generator-level histograms
// Total: 18 plots
// ============================================================================
std::vector<std::string> teacherRequiredGeneratorHists() {
    return {
        // Higgs: pT and eta
        "h_H_pt",
        "h_H_eta",

        // A-boson: pT and eta
        // Only a1 is kept because a1 and a2 are identical particles.
        "h_a1_pt",
        "h_a1_eta",

        // Four b-quarks from H -> aa -> 4b: pT and eta
        // All four are kept because the teacher explicitly asked for 4-b quarks.
        "h_b1_from_a_pt",
        "h_b1_from_a_eta",
        "h_b2_from_a_pt",
        "h_b2_from_a_eta",
        "h_b3_from_a_pt",
        "h_b3_from_a_eta",
        "h_b4_from_a_pt",
        "h_b4_from_a_eta",

        // Top quark: pT and eta
        // Only one representative top is kept.
        "h_t1_pt",
        "h_t1_eta",

        // b-quark from top decay: pT and eta
        // Only one representative b-from-top is kept.
        "h_b1_top_pt",
        "h_b1_top_eta",

        // Angular separations
        // dR_aa is unique. Only one dR_bb is kept because a1 -> bb and a2 -> bb
        // are identical decay topologies.
        "h_dR_aa",
        "h_dR_bb_a1"
    };
}

// ============================================================================
// Clone histogram safely
// ============================================================================
TH1* getHistClone(TFile* f, const char* histName, const char* cloneName) {
    if (!f) return nullptr;

    TH1* h = (TH1*)f->Get(histName);

    if (!h) {
        std::cout << "WARNING: Histogram not found: " << histName << std::endl;
        return nullptr;
    }

    TH1* hClone = (TH1*)h->Clone(cloneName);
    hClone->SetDirectory(nullptr);
    hClone->SetStats(0);

    return hClone;
}

// ============================================================================
// Style histogram
// ============================================================================
void styleHist(TH1* h, const HistConfig& cfg, int color, int lineStyle = 1) {
    if (!h) return;

    h->SetTitle("");
    h->SetStats(0);
    h->SetLineColor(color);
    h->SetLineStyle(lineStyle);
    h->SetLineWidth(3);
    h->SetFillStyle(0);

    h->GetXaxis()->SetTitle(cfg.xTitle.c_str());
    h->GetYaxis()->SetTitle(cfg.yTitle.c_str());

    h->GetXaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetLabelSize(0.040);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetLabelSize(0.040);
    h->GetYaxis()->SetTitleOffset(1.25);

    if (cfg.useRange) {
        h->GetXaxis()->SetRangeUser(cfg.xmin, cfg.xmax);
    }
}

// ============================================================================
// Canvas style
// ============================================================================
void styleCanvas(TCanvas* c) {
    if (!c) return;

    c->SetTicks(1, 1);
    c->SetGridx(0);
    c->SetGridy(0);
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.06);
}

// ============================================================================
// Save one single-signal histogram
// ============================================================================
void saveSinglePlot(TFile* f,
                    const char* outDir,
                    const std::string& histName,
                    const char* signalTag) {
    HistConfig cfg = getConfig(histName);
    cfg.yTitle = "Entries";

    TH1* h = getHistClone(
        f,
        histName.c_str(),
        Form("%s_%s_single_clone_%u",
             histName.c_str(),
             signalTag,
             gRandom->Integer(1000000000))
    );

    if (!h) return;

    styleHist(h, cfg, kBlue + 1, 1);

    double ymax = h->GetMaximum();

    if (ymax > 0) {
        h->SetMaximum(1.25 * ymax);
        h->SetMinimum(0.0);
    }

    TCanvas* c = new TCanvas(
        Form("c_%s_%s_%u",
             cfg.outName.c_str(),
             signalTag,
             gRandom->Integer(1000000000)),
        cfg.outName.c_str(),
        900,
        700
    );

    styleCanvas(c);

    h->Draw("HIST");

    c->Modified();
    c->Update();

    std::string outFile =
        std::string(outDir) + "/" + cfg.outName + ".pdf";

    c->SaveAs(outFile.c_str());

    delete c;
    delete h;
}

// ============================================================================
// Save one overlay plot
// ============================================================================
void saveOverlayPlot(TFile* f12,
                     TFile* f30,
                     const char* outDir,
                     const std::string& histName,
                     bool teacherStyle = false) {
    HistConfig cfg = getConfig(histName);
    cfg.yTitle = "A.U.";

    TH1* h12 = getHistClone(
        f12,
        histName.c_str(),
        Form("%s_tth12_overlay_clone_%u",
             histName.c_str(),
             gRandom->Integer(1000000000))
    );

    TH1* h30 = getHistClone(
        f30,
        histName.c_str(),
        Form("%s_tth30_overlay_clone_%u",
             histName.c_str(),
             gRandom->Integer(1000000000))
    );

    if (!h12 || !h30) {
        delete h12;
        delete h30;
        return;
    }

    if (h12->Integral() > 0) h12->Scale(1.0 / h12->Integral());
    if (h30->Integral() > 0) h30->Scale(1.0 / h30->Integral());

    if (teacherStyle) {
        // Teacher-requested style:
        // mA = 12: black solid
        // mA = 30: black dashed
        styleHist(h12, cfg, kBlack, 1);
        styleHist(h30, cfg, kBlack, 2);
    }
    else {
        // Old style kept for the full old overlay set
        styleHist(h12, cfg, kBlue + 1, 1);
        styleHist(h30, cfg, kRed + 1, 1);
    }

    double ymax = std::max(h12->GetMaximum(), h30->GetMaximum());

    if (ymax > 0) {
        h12->SetMaximum(1.30 * ymax);
        h30->SetMaximum(1.30 * ymax);
        h12->SetMinimum(0.0);
        h30->SetMinimum(0.0);
    }

    TCanvas* c = new TCanvas(
        Form("c_overlay_%s_%u",
             cfg.outName.c_str(),
             gRandom->Integer(1000000000)),
        cfg.outName.c_str(),
        900,
        700
    );

    styleCanvas(c);

    h12->Draw("HIST");
    h30->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.58, 0.72, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);

    if (teacherStyle) {
        leg->AddEntry(h12, "Signal ttH (mA=12)", "l");
        leg->AddEntry(h30, "Signal ttH (mA=30)", "l");
    }
    else {
        leg->AddEntry(h12, Form("m_{a}=12 GeV, mean=%.2f", h12->GetMean()), "l");
        leg->AddEntry(h30, Form("m_{a}=30 GeV, mean=%.2f", h30->GetMean()), "l");
    }

    leg->Draw();

    c->Modified();
    c->Update();

    std::string outFile =
        std::string(outDir) + "/" + cfg.outName + "_Overlay_tth12_vs_tth30.pdf";

    c->SaveAs(outFile.c_str());

    delete leg;
    delete c;
    delete h12;
    delete h30;
}

// ============================================================================
// Draw all generator-level plots for one signal file
// Saves into gen_level_plots/<signalTag>/
// ============================================================================
void drawGenLevelForFile(const char* inputFileName, const char* signalTag) {
    std::cout << "\n============================================" << std::endl;
    std::cout << "Processing ALL plots: " << inputFileName << std::endl;
    std::cout << "Signal tag: " << signalTag << std::endl;
    std::cout << "============================================" << std::endl;

    TFile* f = TFile::Open(inputFileName);

    if (!f || f->IsZombie()) {
        std::cout << "ERROR: Cannot open ROOT file: " << inputFileName << std::endl;
        return;
    }

    gStyle->SetOptStat(0);

    std::string outDir = std::string("gen_level_plots/") + signalTag;

    gSystem->mkdir("gen_level_plots", kTRUE);
    gSystem->mkdir(outDir.c_str(), kTRUE);

    for (const auto& hname : allGeneratorHists()) {
        saveSinglePlot(f, outDir.c_str(), hname, signalTag);
    }

    f->Close();
    delete f;

    std::cout << "Saved ALL generator-level PDFs in:" << std::endl;
    std::cout << "  " << outDir << std::endl;
}

// ============================================================================
// Draw teacher-required generator-level plots for one signal file
// Saves into gen_level_plots/teacher_required/<signalTag>/
// ============================================================================
void drawTeacherRequiredForFile(const char* inputFileName, const char* signalTag) {
    std::cout << "\n============================================" << std::endl;
    std::cout << "Processing TEACHER-REQUIRED plots: " << inputFileName << std::endl;
    std::cout << "Signal tag: " << signalTag << std::endl;
    std::cout << "============================================" << std::endl;

    TFile* f = TFile::Open(inputFileName);

    if (!f || f->IsZombie()) {
        std::cout << "ERROR: Cannot open ROOT file: " << inputFileName << std::endl;
        return;
    }

    gStyle->SetOptStat(0);

    std::string baseDir = "gen_level_plots/teacher_required";
    std::string outDir = baseDir + "/" + signalTag;

    gSystem->mkdir("gen_level_plots", kTRUE);
    gSystem->mkdir(baseDir.c_str(), kTRUE);
    gSystem->mkdir(outDir.c_str(), kTRUE);

    std::vector<std::string> required = teacherRequiredGeneratorHists();

    for (const auto& hname : required) {
        saveSinglePlot(f, outDir.c_str(), hname, signalTag);
    }

    f->Close();
    delete f;

    std::cout << "Saved TEACHER-REQUIRED generator-level PDFs in:" << std::endl;
    std::cout << "  " << outDir << std::endl;
    std::cout << "Number of teacher-required plots: "
              << required.size()
              << std::endl;
}

// ============================================================================
// Draw all tth12gev vs tth30gev overlays
// Saves into gen_level_plots/overlay_tth12_vs_tth30/
// ============================================================================
void drawOverlay12vs30() {
    const char* file12Name = "output_signal_tth12gev.root";
    const char* file30Name = "output_signal_tth30gev.root";

    if (gSystem->AccessPathName(file12Name)) {
        std::cout << "ERROR: File not found: " << file12Name << std::endl;
        return;
    }

    if (gSystem->AccessPathName(file30Name)) {
        std::cout << "ERROR: File not found: " << file30Name << std::endl;
        return;
    }

    TFile* f12 = TFile::Open(file12Name);
    TFile* f30 = TFile::Open(file30Name);

    if (!f12 || f12->IsZombie()) {
        std::cout << "ERROR: Cannot open ROOT file: " << file12Name << std::endl;
        return;
    }

    if (!f30 || f30->IsZombie()) {
        std::cout << "ERROR: Cannot open ROOT file: " << file30Name << std::endl;

        if (f12) {
            f12->Close();
            delete f12;
        }

        return;
    }

    gStyle->SetOptStat(0);

    std::string outDir = "gen_level_plots/overlay_tth12_vs_tth30";

    gSystem->mkdir("gen_level_plots", kTRUE);
    gSystem->mkdir(outDir.c_str(), kTRUE);

    std::cout << "\n============================================" << std::endl;
    std::cout << "Creating ALL overlay plots directly in old folder" << std::endl;
    std::cout << "Output folder: " << outDir << std::endl;
    std::cout << "============================================" << std::endl;

    for (const auto& hname : allGeneratorHists()) {
        saveOverlayPlot(f12, f30, outDir.c_str(), hname, false);
    }

    f12->Close();
    f30->Close();

    delete f12;
    delete f30;

    std::cout << "ALL overlay PDFs are ready in:" << std::endl;
    std::cout << "  " << outDir << std::endl;
}

// ============================================================================
// Draw teacher-required tth12gev vs tth30gev overlays
// Saves into gen_level_plots/teacher_required/overlay_tth12_vs_tth30/
// ============================================================================
void drawTeacherRequiredOverlay12vs30() {
    const char* file12Name = "output_signal_tth12gev.root";
    const char* file30Name = "output_signal_tth30gev.root";

    if (gSystem->AccessPathName(file12Name)) {
        std::cout << "ERROR: File not found: " << file12Name << std::endl;
        return;
    }

    if (gSystem->AccessPathName(file30Name)) {
        std::cout << "ERROR: File not found: " << file30Name << std::endl;
        return;
    }

    TFile* f12 = TFile::Open(file12Name);
    TFile* f30 = TFile::Open(file30Name);

    if (!f12 || f12->IsZombie()) {
        std::cout << "ERROR: Cannot open ROOT file: " << file12Name << std::endl;
        return;
    }

    if (!f30 || f30->IsZombie()) {
        std::cout << "ERROR: Cannot open ROOT file: " << file30Name << std::endl;

        if (f12) {
            f12->Close();
            delete f12;
        }

        return;
    }

    gStyle->SetOptStat(0);

    std::string baseDir = "gen_level_plots/teacher_required";
    std::string outDir = baseDir + "/overlay_tth12_vs_tth30";

    gSystem->mkdir("gen_level_plots", kTRUE);
    gSystem->mkdir(baseDir.c_str(), kTRUE);
    gSystem->mkdir(outDir.c_str(), kTRUE);

    std::cout << "\n============================================" << std::endl;
    std::cout << "Creating TEACHER-REQUIRED overlay plots" << std::endl;
    std::cout << "Output folder: " << outDir << std::endl;
    std::cout << "============================================" << std::endl;

    std::vector<std::string> required = teacherRequiredGeneratorHists();

    for (const auto& hname : required) {
        saveOverlayPlot(f12, f30, outDir.c_str(), hname, true);
    }

    f12->Close();
    f30->Close();

    delete f12;
    delete f30;

    std::cout << "TEACHER-REQUIRED overlay PDFs are ready in:" << std::endl;
    std::cout << "  " << outDir << std::endl;
    std::cout << "Number of teacher-required overlay plots: "
              << required.size()
              << std::endl;
}

// ============================================================================
// Helper: find all signal ROOT files
// ============================================================================
std::vector<std::string> findSignalFiles() {
    std::vector<std::string> signalFiles;

    TSystemDirectory dir(".", ".");
    TList* files = dir.GetListOfFiles();

    if (files) {
        TIter next(files);
        TSystemFile* file;

        while ((file = (TSystemFile*)next())) {
            TString fname = file->GetName();

            if (file->IsDirectory()) continue;

            if (fname.BeginsWith("output_signal_") && fname.EndsWith(".root")) {
                signalFiles.push_back(std::string(fname.Data()));
            }
        }
    }

    std::sort(signalFiles.begin(), signalFiles.end());
    return signalFiles;
}

// ============================================================================
// Main function: full old set + teacher-required set
// ============================================================================
void myplot() {
    std::vector<std::string> signalFiles = findSignalFiles();

    if (signalFiles.empty()) {
        std::cout << "ERROR: No signal files found." << std::endl;
        std::cout << "Expected files like output_signal_tth12gev.root" << std::endl;
        return;
    }

    std::cout << "Found " << signalFiles.size() << " signal file(s)." << std::endl;

    for (const auto& fileName : signalFiles) {
        TString tag = fileName.c_str();
        tag.ReplaceAll("output_signal_", "");
        tag.ReplaceAll(".root", "");

        // Old full set
        drawGenLevelForFile(fileName.c_str(), tag.Data());

        // New teacher-required reduced set
        drawTeacherRequiredForFile(fileName.c_str(), tag.Data());
    }

    // Old full overlay set
    drawOverlay12vs30();

    // New teacher-required overlay set
    drawTeacherRequiredOverlay12vs30();

    std::cout << "\nAll generator-level PDFs are ready." << std::endl;
    std::cout << "Full set:" << std::endl;
    std::cout << "  gen_level_plots/<signalTag>/" << std::endl;
    std::cout << "  gen_level_plots/overlay_tth12_vs_tth30/" << std::endl;
    std::cout << "Teacher-required reduced set:" << std::endl;
    std::cout << "  gen_level_plots/teacher_required/<signalTag>/" << std::endl;
    std::cout << "  gen_level_plots/teacher_required/overlay_tth12_vs_tth30/" << std::endl;
    std::cout << "Teacher-required number of plots: "
              << teacherRequiredGeneratorHists().size()
              << " per mass point and "
              << teacherRequiredGeneratorHists().size()
              << " overlays." << std::endl;
}

// ============================================================================
// Optional helper: one signal mass point only
// ============================================================================
void myplotOne(const char* signalTag = "tth12gev") {
    TString fileName = Form("output_signal_%s.root", signalTag);

    if (gSystem->AccessPathName(fileName.Data())) {
        std::cout << "ERROR: File not found: " << fileName << std::endl;
        return;
    }

    drawGenLevelForFile(fileName.Data(), signalTag);
    drawTeacherRequiredForFile(fileName.Data(), signalTag);
}

// ============================================================================
// Optional helper: overlay only
// ============================================================================
void myplotOverlay12vs30() {
    drawOverlay12vs30();
    drawTeacherRequiredOverlay12vs30();
}

// ============================================================================
// Optional helper: only teacher-required plots
// ============================================================================
void myplotTeacherRequired() {
    std::vector<std::string> signalFiles = findSignalFiles();

    if (signalFiles.empty()) {
        std::cout << "ERROR: No signal files found." << std::endl;
        std::cout << "Expected files like output_signal_tth12gev.root" << std::endl;
        return;
    }

    std::cout << "Found " << signalFiles.size() << " signal file(s)." << std::endl;

    for (const auto& fileName : signalFiles) {
        TString tag = fileName.c_str();
        tag.ReplaceAll("output_signal_", "");
        tag.ReplaceAll(".root", "");

        drawTeacherRequiredForFile(fileName.c_str(), tag.Data());
    }

    drawTeacherRequiredOverlay12vs30();

    std::cout << "\nTeacher-required generator-level PDFs are ready in:" << std::endl;
    std::cout << "  gen_level_plots/teacher_required/" << std::endl;
    std::cout << "Teacher-required number of plots: "
              << teacherRequiredGeneratorHists().size()
              << " per mass point and "
              << teacherRequiredGeneratorHists().size()
              << " overlays." << std::endl;
}
