// redraw_brazilian_upto25.C
//
// Fast redraw of the Bayesian Brazilian plot up to m_a = 25 GeV.
// No MCMC. No toys. No log plot.
//
// It saves:
//
//   fits_and_limit_plots/step5_expected_limits_bayesian_upto25.pdf
//
// Run:
//
//   root -l
//   .L redraw_brazilian_upto25.C
//   redraw_brazilian_upto25()

#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TH1F.h>
#include <TAxis.h>
#include <TColor.h>

#include <vector>
#include <algorithm>
#include <iostream>

void redraw_brazilian_upto25() {

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    // Grid style
    gStyle->SetGridStyle(3);
    gStyle->SetGridColor(kGray + 1);
    gStyle->SetGridWidth(1);

    gSystem->mkdir(
        "fits_and_limit_plots",
        kTRUE
    );

    // Only keep mass points up to 25 GeV.
    const int n = 4;

    double mass[n] = {
        12.0,
        15.0,
        20.0,
        25.0
    };

    double expected[n] = {
        4.39218,
        4.2696,
        8.56049,
        20.5018
    };

    double sigma1_down[n] = {
        2.09082,
        1.68117,
        3.6413,
        8.77225
    };

    double sigma1_up[n] = {
        7.58048,
        7.67774,
        15.3496,
        33.2533
    };

    double sigma2_down[n] = {
        -0.475065,
        -0.580349,
        -0.987691,
        -2.77622
    };

    double sigma2_up[n] = {
        11.2525,
        10.886,
        23.6389,
        41.9963
    };

    double pseudo_observed[n] = {
        4.32045,
        7.23413,
        6.70276,
        21.1643
    };

    double exl[n] = {
        0.0,
        0.0,
        0.0,
        0.0
    };

    double exh[n] = {
        0.0,
        0.0,
        0.0,
        0.0
    };

    double eyl1[n];
    double eyh1[n];
    double eyl2[n];
    double eyh2[n];

    double ymax = 0.0;

    for (int i = 0; i < n; ++i) {
        eyl1[i] = expected[i] - sigma1_down[i];
        eyh1[i] = sigma1_up[i] - expected[i];

        eyl2[i] = expected[i] - sigma2_down[i];
        eyh2[i] = sigma2_up[i] - expected[i];

        ymax = std::max(
            ymax,
            sigma2_up[i]
        );

        ymax = std::max(
            ymax,
            pseudo_observed[i]
        );

        ymax = std::max(
            ymax,
            expected[i]
        );
    }

    ymax *= 1.25;

    TCanvas* c = new TCanvas(
        "c_brazilian_upto25",
        "Bayesian 95% upper limits up to 25 GeV",
        900,
        700
    );

    c->SetTicks(1, 1);

    c->SetGridx(1);
    c->SetGridy(1);
    c->SetGrid();

    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);
    c->SetLogy(0);

    // Exact x-axis range: 12–25 GeV.
    // There is no empty space before 12 or after 25.
    TH1F* frame = new TH1F(
        "frame_brazilian_upto25",
        ";m_{a} [GeV];Bayesian 95% upper limit on BR",
        130,
        12.0,
        25.0
    );

    frame->SetMinimum(0.0);
    frame->SetMaximum(ymax);

    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.040);

    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetTitleOffset(1.25);

    frame->Draw();

    TGraphAsymmErrors* g2 = new TGraphAsymmErrors(
        n,
        mass,
        expected,
        exl,
        exh,
        eyl2,
        eyh2
    );

    TGraphAsymmErrors* g1 = new TGraphAsymmErrors(
        n,
        mass,
        expected,
        exl,
        exh,
        eyl1,
        eyh1
    );

    TGraph* gmed = new TGraph(
        n,
        mass,
        expected
    );

    TGraph* gpseudo = new TGraph(
        n,
        mass,
        pseudo_observed
    );

    // Brazilian bands
    g2->SetFillColor(kYellow);
    g2->SetLineColor(kYellow);

    g1->SetFillColor(kGreen + 1);
    g1->SetLineColor(kGreen + 1);

    // Median expected
    gmed->SetLineColor(kBlack);
    gmed->SetLineWidth(2);
    gmed->SetLineStyle(2);
    gmed->SetMarkerStyle(20);
    gmed->SetMarkerSize(1.1);
    gmed->SetMarkerColor(kBlack);

    // Pseudo-observed
    gpseudo->SetLineColor(kRed + 1);
    gpseudo->SetLineWidth(2);
    gpseudo->SetMarkerStyle(24);
    gpseudo->SetMarkerSize(1.1);
    gpseudo->SetMarkerColor(kRed + 1);

    g2->Draw("3 SAME");
    g1->Draw("3 SAME");
    gmed->Draw("LP SAME");
    gpseudo->Draw("LP SAME");

    TLegend* leg = new TLegend(
        0.55,
        0.64,
        0.88,
        0.88
    );

    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);

    leg->AddEntry(
        gmed,
        "Median expected",
        "lp"
    );

    leg->AddEntry(
        gpseudo,
        "Pseudo-observed",
        "lp"
    );

    leg->AddEntry(
        g1,
        "#pm1#sigma",
        "f"
    );

    leg->AddEntry(
        g2,
        "#pm2#sigma",
        "f"
    );

    leg->Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.035);

    lat.DrawLatex(
        0.15,
        0.92,
        "Bayesian 95% upper limits"
    );

    // Redraw the grid and axes after the coloured bands.
    c->RedrawAxis("g");
    c->RedrawAxis();

    c->Modified();
    c->Update();

    c->SaveAs(
        "fits_and_limit_plots/"
        "step5_expected_limits_bayesian_upto25.pdf"
    );

    std::cout << "Saved:" << std::endl;
    std::cout
        << "  fits_and_limit_plots/"
        << "step5_expected_limits_bayesian_upto25.pdf"
        << std::endl;

    delete leg;
    delete gpseudo;
    delete gmed;
    delete g1;
    delete g2;
    delete frame;
    delete c;
}
