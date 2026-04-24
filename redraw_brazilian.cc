#include <vector>
#include <iostream>
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLatex.h"

void redraw_brazilian() {
  std::vector<double> x    = {12, 15, 20, 25, 30};

  std::vector<double> y    = {4.39218, 4.2696, 8.56049, 20.5018, 30.8669};
  std::vector<double> y1dn = {2.09082, 1.68117, 3.6413, 8.77225, 13.9704};
  std::vector<double> y1up = {7.58048, 7.67774, 15.3496, 33.2533, 55.7041};
  std::vector<double> y2dn = {-0.475065, -0.580349, -0.987691, -2.77622, -21.6915};
  std::vector<double> y2up = {11.2525, 10.886, 23.6389, 41.9963, 75.5984};
  std::vector<double> yobs = {4.32045, 7.23413, 6.70276, 21.1643, 21.9045};

  const int n = (int)x.size();

  std::vector<double> exl(n, 0.0), exh(n, 0.0);
  std::vector<double> eyl1(n), eyh1(n), eyl2(n), eyh2(n);

  for (int i = 0; i < n; ++i) {
    eyl1[i] = y[i] - y1dn[i];
    eyh1[i] = y1up[i] - y[i];
    eyl2[i] = y[i] - y2dn[i];
    eyh2[i] = y2up[i] - y[i];
  }

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.035);

  // ---------------- linear canvas ----------------
  TCanvas* c_lin = new TCanvas("c_brazil_lin", "Brazilian plot linear", 900, 700);
  c_lin->SetLogy(0);

  TGraphAsymmErrors* g2_lin = new TGraphAsymmErrors(
    n, x.data(), y.data(), exl.data(), exh.data(), eyl2.data(), eyh2.data()
  );
  TGraphAsymmErrors* g1_lin = new TGraphAsymmErrors(
    n, x.data(), y.data(), exl.data(), exh.data(), eyl1.data(), eyh1.data()
  );
  TGraph* gmed_lin = new TGraph(n, x.data(), y.data());
  TGraph* gobs_lin = new TGraph(n, x.data(), yobs.data());

  g2_lin->SetFillColor(kYellow);
  g2_lin->SetLineColor(kYellow);
  g2_lin->SetTitle(";m_{a} [GeV];Bayesian 95% upper limit on BR");
  g2_lin->SetMinimum(-10);

  g1_lin->SetFillColor(kGreen + 1);
  g1_lin->SetLineColor(kGreen + 1);

  gmed_lin->SetLineColor(kBlack);
  gmed_lin->SetLineWidth(2);
  gmed_lin->SetLineStyle(2);
  gmed_lin->SetMarkerStyle(20);

  gobs_lin->SetLineColor(kRed + 1);
  gobs_lin->SetLineWidth(2);
  gobs_lin->SetMarkerStyle(24);

  g2_lin->Draw("A3");
  g1_lin->Draw("3 SAME");
  gmed_lin->Draw("LP SAME");
  gobs_lin->Draw("LP SAME");

  TLegend* leg_lin = new TLegend(0.55, 0.64, 0.88, 0.88);
  leg_lin->SetBorderSize(0);
  leg_lin->SetFillStyle(0);
  leg_lin->AddEntry(gmed_lin, "Median expected", "lp");
  leg_lin->AddEntry(gobs_lin, "Pseudo-observed", "lp");
  leg_lin->AddEntry(g1_lin, "#pm1#sigma", "f");
  leg_lin->AddEntry(g2_lin, "#pm2#sigma", "f");
  leg_lin->Draw();

  lat.DrawLatex(0.15, 0.92, "Bayesian 95% upper limits");
  c_lin->SaveAs("brazilian_redraw_linear.pdf");

  // ---------------- log canvas ----------------
  TCanvas* c_log = new TCanvas("c_brazil_log", "Brazilian plot log", 900, 700);
  c_log->SetLogy(1);

  TGraphAsymmErrors* g2_log = new TGraphAsymmErrors(
    n, x.data(), y.data(), exl.data(), exh.data(), eyl2.data(), eyh2.data()
  );
  TGraphAsymmErrors* g1_log = new TGraphAsymmErrors(
    n, x.data(), y.data(), exl.data(), exh.data(), eyl1.data(), eyh1.data()
  );
  TGraph* gmed_log = new TGraph(n, x.data(), y.data());
  TGraph* gobs_log = new TGraph(n, x.data(), yobs.data());

  g2_log->SetFillColor(kYellow);
  g2_log->SetLineColor(kYellow);
  g2_log->SetTitle(";m_{a} [GeV];Bayesian 95% upper limit on BR");
  g2_log->SetMinimum(1e-1);

  g1_log->SetFillColor(kGreen + 1);
  g1_log->SetLineColor(kGreen + 1);

  gmed_log->SetLineColor(kBlack);
  gmed_log->SetLineWidth(2);
  gmed_log->SetLineStyle(2);
  gmed_log->SetMarkerStyle(20);

  gobs_log->SetLineColor(kRed + 1);
  gobs_log->SetLineWidth(2);
  gobs_log->SetMarkerStyle(24);

  g2_log->Draw("A3");
  g1_log->Draw("3 SAME");
  gmed_log->Draw("LP SAME");
  gobs_log->Draw("LP SAME");

  TLegend* leg_log = new TLegend(0.55, 0.64, 0.88, 0.88);
  leg_log->SetBorderSize(0);
  leg_log->SetFillStyle(0);
  leg_log->AddEntry(gmed_log, "Median expected", "lp");
  leg_log->AddEntry(gobs_log, "Pseudo-observed", "lp");
  leg_log->AddEntry(g1_log, "#pm1#sigma", "f");
  leg_log->AddEntry(g2_log, "#pm2#sigma", "f");
  leg_log->Draw();

  lat.DrawLatex(0.15, 0.92, "Bayesian 95% upper limits");
  c_log->SaveAs("brazilian_redraw_log.pdf");
}
