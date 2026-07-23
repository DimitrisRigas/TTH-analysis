// overlay_signal_allbkgs.C
// Mode-2 overlays: Signal vs selected dominant backgrounds using THStack,
// plus ratio pad S / (sum backgrounds).
//
// UPDATED AFTER MENTOR COMMENTS:
// - Signal ttH (mA=12): black solid line, SetLineStyle(1)
// - Signal ttH (mA=30): black dashed line, SetLineStyle(2)
// - Legend labels fixed:
//      Signal ttH (mA=12)
//      Signal ttH (mA=30)
// - Top pad has no grid.
// - Ratio pad has GridX and GridY.
// - Reconstructed phi plots are NOT produced.
// - Only three backgrounds are drawn:
//      1) TTH_Hbb
//      2) TT dileptonic
//      3) TT semileptonic
// - Output PDF names are kept in the same format as before.
// - Uses custom rebinning that preserves clean x-axis ranges.
//
// Run:
//   .L overlay_signal_allbkgs.C
//   overlay_signal_allbkgs("tth12gev", true,  true);   // shapes, keep windows
//   overlay_signal_allbkgs("tth12gev", false, false);  // raw entries, close windows
//
// Notes:
// - If normalize=true: each component is unit-area normalized and the y-axis is A.U.
// - If normalize=false: plots raw entries from each file.
// - Histograms are cloned into gROOT so files can be closed safely.

#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TString.h>
#include <TRandom.h>
#include <TSystem.h>

#include <cmath>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>

// Keep canvases alive
static std::vector<TCanvas*> gCanvases;

struct BkgInfo {
  std::string file;
  std::string tag;
  std::string label;
};

struct AxisRange {
  bool use;
  double xmin;
  double xmax;
};

// ----------------------------------------------------------------------------
// Thesis-friendly fixed x-axis ranges.
// ----------------------------------------------------------------------------
static AxisRange GetNiceXRange(const std::string& hname) {
  // Phi variables are kept here only for safety,
  // but reconstructed phi plots are removed from the final histogram list.
  if (hname.find("phi") != std::string::npos) {
    return {true, -3.2, 3.2};
  }

  // Eta variables
  if (hname == "h_Hdbb_eta") {
    return {true, -5.0, 5.0};
  }

  if (hname.find("eta") != std::string::npos) {
    return {true, -2.5, 2.5};
  }

  // Higgs candidate variables
  if (hname == "h_Hdbb_mass") {
    return {true, 0.0, 400.0};
  }

  if (hname == "h_Hdbb_pt") {
    return {true, 0.0, 500.0};
  }

  // Double-b jet masses
  if (hname == "h_dbj1_mass" || hname == "h_dbj2_mass") {
    return {true, 0.0, 80.0};
  }

  // Double-b jet pT
  if (hname == "h_dbj1_pt" || hname == "h_dbj2_pt") {
    return {true, 0.0, 400.0};
  }

  // Final b-jet pT
  if (hname == "h_bj1_pt_final" || hname == "h_bj2_pt_final") {
    return {true, 0.0, 400.0};
  }

  // Leading lepton pT
  if (hname == "h_lep1_pt_final") {
    return {true, 0.0, 400.0};
  }

  // DeltaR variables
  if (hname == "h_dR_dbj12_final" || hname == "h_dRll") {
    return {true, 0.0, 5.0};
  }

  // Delta mass variable
  if (hname == "h_dM_bj12_final") {
    return {true, 0.0, 150.0};
  }

  // MET
  if (hname == "h_MET_pt_final") {
    return {true, 0.0, 300.0};
  }

  // HT
  if (hname == "h_HT") {
    return {true, 0.0, 1200.0};
  }

  // Dilepton mass
  if (hname == "h_mll") {
    return {true, 0.0, 200.0};
  }

  return {false, 0.0, 0.0};
}

// ----------------------------------------------------------------------------
// Apply fixed x-axis range to an already existing histogram.
// ----------------------------------------------------------------------------
static void ApplyNiceXRange(TH1* h, const std::string& hname) {
  if (!h) return;

  AxisRange r = GetNiceXRange(hname);

  if (r.use) {
    h->GetXaxis()->SetRangeUser(r.xmin, r.xmax);
  }
}

// ----------------------------------------------------------------------------
// Clone histogram out of file into gROOT and cache it.
// ----------------------------------------------------------------------------
static TH1* CloneToROOT(TFile* f, const std::string& keyPrefix, const char* hname) {
  static std::unordered_map<std::string, TH1*> cache;
  static long long ctr = 0;

  const std::string key = keyPrefix + "::" + hname;
  auto it = cache.find(key);

  if (it != cache.end()) return it->second;

  TH1* h = (TH1*)f->Get(hname);

  if (!h) {
    std::cout << "WARNING: Missing histogram " << hname
              << " in file " << f->GetName() << std::endl;
    return nullptr;
  }

  const TString uniqName = Form("%s__keep__%s__%lld",
                                hname,
                                keyPrefix.c_str(),
                                ctr++);

  TH1* hc = (TH1*)h->Clone(uniqName);
  hc->SetDirectory(gROOT);
  hc->SetStats(0);

  cache[key] = hc;
  return hc;
}

// ----------------------------------------------------------------------------
// Custom rebinning that keeps the intended smoothing while using a clean range.
// ----------------------------------------------------------------------------
static TH1* RebinPreserveCleanRange(TH1* h,
                                    const std::string& hname,
                                    const std::string& tag,
                                    int rebinFactor) {
  if (!h) return nullptr;

  const int oldN = h->GetNbinsX();

  const double oldXmin = h->GetXaxis()->GetXmin();
  const double oldXmax = h->GetXaxis()->GetXmax();

  AxisRange niceRange = GetNiceXRange(hname);

  double xmin = oldXmin;
  double xmax = oldXmax;

  if (niceRange.use) {
    xmin = niceRange.xmin;
    xmax = niceRange.xmax;
  }

  double oldBinWidth = (oldXmax - oldXmin) / oldN;

  if (oldBinWidth <= 0.0) {
    oldBinWidth = (xmax - xmin) / oldN;
  }

  const int safeRebin = std::max(1, rebinFactor);
  const double targetBinWidth = oldBinWidth * safeRebin;

  int newN = (int)std::ceil((xmax - xmin) / targetBinWidth);

  if (newN < 1) newN = 1;

  TH1D* hnew = new TH1D(
    Form("%s__plotclone_%s_%u",
         h->GetName(),
         tag.c_str(),
         gRandom->Integer(1000000000)),
    "",
    newN,
    xmin,
    xmax
  );

  hnew->SetDirectory(gROOT);
  hnew->SetStats(0);
  hnew->Sumw2();

  hnew->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
  hnew->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());

  std::vector<double> err2(newN + 2, 0.0);

  for (int ib = 1; ib <= oldN; ++ib) {
    const double x = h->GetBinCenter(ib);

    if (x < xmin || x >= xmax) continue;

    const int jb = hnew->FindBin(x);

    if (jb < 1 || jb > newN) continue;

    const double oldContent = h->GetBinContent(ib);
    const double oldError   = h->GetBinError(ib);

    hnew->SetBinContent(jb, hnew->GetBinContent(jb) + oldContent);
    err2[jb] += oldError * oldError;
  }

  for (int jb = 1; jb <= newN; ++jb) {
    hnew->SetBinError(jb, std::sqrt(err2[jb]));
  }

  ApplyNiceXRange(hnew, hname);

  return hnew;
}

// ----------------------------------------------------------------------------
// Clone a histogram for one specific plot and rebin it safely.
// ----------------------------------------------------------------------------
static TH1* MakePlotClone(TH1* h,
                          const std::string& hname,
                          const std::string& tag,
                          int ibin) {
  if (!h) return nullptr;

  TH1* hc = RebinPreserveCleanRange(h, hname, tag, ibin);

  if (!hc) return nullptr;

  hc->SetDirectory(gROOT);
  hc->SetStats(0);

  ApplyNiceXRange(hc, hname);

  return hc;
}

// ----------------------------------------------------------------------------
// Normalize if requested.
// ----------------------------------------------------------------------------
static void Prep(TH1* h, bool normalize) {
  if (!h) return;

  h->SetStats(0);

  if (normalize) {
    const double I = h->Integral();

    if (I > 0) {
      h->Scale(1.0 / I);
    }
  }
}

// ----------------------------------------------------------------------------
// Nice axis labels/titles for each histogram key.
// ----------------------------------------------------------------------------
static void SetNiceLabels(TH1* h, const std::string& hname, bool normalize) {
  if (!h) return;

  h->SetTitle("");
  h->GetYaxis()->SetTitle(normalize ? "A.U." : "Entries");

  if (hname == "h_Hdbb_mass") {
    h->GetXaxis()->SetTitle("m_{bb} (Higgs candidate) [GeV]");
  }
  else if (hname == "h_Hdbb_pt") {
    h->GetXaxis()->SetTitle("p_{T}(Higgs candidate) [GeV]");
  }
  else if (hname == "h_Hdbb_eta") {
    h->GetXaxis()->SetTitle("#eta(Higgs candidate)");
  }

  else if (hname == "h_dbj1_mass") {
    h->GetXaxis()->SetTitle("m(double-b jet 1) [GeV]");
  }
  else if (hname == "h_dbj1_pt") {
    h->GetXaxis()->SetTitle("p_{T}(double-b jet 1) [GeV]");
  }
  else if (hname == "h_dbj1_eta") {
    h->GetXaxis()->SetTitle("#eta(double-b jet 1)");
  }

  else if (hname == "h_dbj2_mass") {
    h->GetXaxis()->SetTitle("m(double-b jet 2) [GeV]");
  }
  else if (hname == "h_dbj2_pt") {
    h->GetXaxis()->SetTitle("p_{T}(double-b jet 2) [GeV]");
  }
  else if (hname == "h_dbj2_eta") {
    h->GetXaxis()->SetTitle("#eta(double-b jet 2)");
  }

  else if (hname == "h_dR_dbj12_final") {
    h->GetXaxis()->SetTitle("#DeltaR(double-b jet 1, double-b jet 2)");
  }
  else if (hname == "h_dM_bj12_final") {
    h->GetXaxis()->SetTitle("|m(b_{1}) - m(b_{2})| [GeV]");
  }
  else if (hname == "h_MET_pt_final") {
    h->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  }
  else if (hname == "h_HT") {
    h->GetXaxis()->SetTitle("H_{T} (scalar #Sigma p_{T}^{jets}) [GeV]");
  }
  else if (hname == "h_mll") {
    h->GetXaxis()->SetTitle("m_{ll} [GeV]");
  }

  else if (hname == "h_bj1_pt_final") {
    h->GetXaxis()->SetTitle("p_{T}(b-jet 1) [GeV]");
  }
  else if (hname == "h_bj1_eta_final") {
    h->GetXaxis()->SetTitle("#eta(b-jet 1)");
  }
  else if (hname == "h_bj2_pt_final") {
    h->GetXaxis()->SetTitle("p_{T}(b-jet 2) [GeV]");
  }
  else if (hname == "h_bj2_eta_final") {
    h->GetXaxis()->SetTitle("#eta(b-jet 2)");
  }

  else if (hname == "h_lep1_pt_final") {
    h->GetXaxis()->SetTitle("p_{T}(leading lepton) [GeV]");
  }
  else if (hname == "h_lep1_eta_final") {
    h->GetXaxis()->SetTitle("#eta(leading lepton)");
  }
  else if (hname == "h_dRll") {
    h->GetXaxis()->SetTitle("#DeltaR(l_{1}, l_{2})");
  }

  h->GetXaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetLabelSize(0.040);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetLabelSize(0.040);
  h->GetYaxis()->SetTitleOffset(1.25);

  ApplyNiceXRange(h, hname);
}

// ----------------------------------------------------------------------------
// Prettier file naming.
// IMPORTANT: kept unchanged so existing Overleaf filenames still work.
// ----------------------------------------------------------------------------
static std::string PrettyKey(const std::string& hname) {
  if (hname == "h_Hdbb_mass") return "Higgs_mbb";
  if (hname == "h_Hdbb_pt")   return "Higgs_pt";
  if (hname == "h_Hdbb_eta")  return "Higgs_eta";

  if (hname == "h_dbj1_mass") return "dbjet1_mass";
  if (hname == "h_dbj1_pt")   return "dbjet1_pt";
  if (hname == "h_dbj1_eta")  return "dbjet1_eta";

  if (hname == "h_dbj2_mass") return "dbjet2_mass";
  if (hname == "h_dbj2_pt")   return "dbjet2_pt";
  if (hname == "h_dbj2_eta")  return "dbjet2_eta";

  if (hname == "h_dR_dbj12_final") return "dR_dbjets";
  if (hname == "h_dM_bj12_final")  return "dM_bjets";
  if (hname == "h_MET_pt_final")   return "MET_pt";
  if (hname == "h_HT")             return "HT";
  if (hname == "h_mll")            return "mll";

  if (hname == "h_bj1_pt_final")  return "bjet1_pt";
  if (hname == "h_bj1_eta_final") return "bjet1_eta";
  if (hname == "h_bj2_pt_final")  return "bjet2_pt";
  if (hname == "h_bj2_eta_final") return "bjet2_eta";

  if (hname == "h_lep1_pt_final")  return "lep1_pt";
  if (hname == "h_lep1_eta_final") return "lep1_eta";
  if (hname == "h_dRll")           return "dR_ll";

  return hname;
}

// ----------------------------------------------------------------------------
// Draw stacked backgrounds + two signal overlays with ratio pad.
// Ratio shows both: S(mA=12)/ΣB and S(mA=30)/ΣB.
// ----------------------------------------------------------------------------
static void DrawStackWithRatio2Sig(TCanvas* c,
                                   const char* hnameKey,
                                   TH1* hsA_in,
                                   TH1* hsB_in,
                                   const std::vector<TH1*>& hb_list_in,
                                   const std::vector<std::string>& bkgLabels,
                                   const char* sigALabel,
                                   const char* sigBLabel,
                                   bool normalize) {
  if (!c || !hsA_in || !hsB_in || hb_list_in.empty() || !hnameKey) return;

  const std::string hname = hnameKey;

  c->cd();

  c->SetGridx(0);
  c->SetGridy(0);

  // --------------------------------------------------------------------------
  // Top pad: no grid
  // --------------------------------------------------------------------------
  TPad* pTop = new TPad(Form("%s_top", c->GetName()), "", 0.0, 0.30, 1.0, 1.0);
  pTop->SetBottomMargin(0.02);
  pTop->SetLeftMargin(0.12);
  pTop->SetRightMargin(0.05);
  pTop->SetTopMargin(0.06);
  pTop->SetGridx(0);
  pTop->SetGridy(0);
  pTop->Draw();

  // --------------------------------------------------------------------------
  // Ratio pad: grid ON, as requested
  // --------------------------------------------------------------------------
  TPad* pBot = new TPad(Form("%s_bot", c->GetName()), "", 0.0, 0.00, 1.0, 0.30);
  pBot->SetTopMargin(0.02);
  pBot->SetBottomMargin(0.35);
  pBot->SetLeftMargin(0.12);
  pBot->SetRightMargin(0.05);
  pBot->SetGridx(1);
  pBot->SetGridy(1);
  pBot->Draw();

  TH1* hsA = (TH1*)hsA_in->Clone(Form("%s__draw_sigA__%u",
                                      hsA_in->GetName(),
                                      gRandom->Integer(1000000000)));

  TH1* hsB = (TH1*)hsB_in->Clone(Form("%s__draw_sigB__%u",
                                      hsB_in->GetName(),
                                      gRandom->Integer(1000000000)));

  hsA->SetDirectory(gROOT);
  hsB->SetDirectory(gROOT);

  Prep(hsA, normalize);
  Prep(hsB, normalize);

  SetNiceLabels(hsA, hname, normalize);
  SetNiceLabels(hsB, hname, normalize);

  THStack* st = new THStack(Form("st_%s_%u",
                                 hnameKey,
                                 gRandom->Integer(1000000000)), "");

  TH1* hbSum = nullptr;

  const int colors[] = {
    kGreen - 7,
    kAzure - 9,
    kOrange - 2
  };

  const int ncol = sizeof(colors) / sizeof(colors[0]);

  std::vector<TH1*> hb_drawn;
  hb_drawn.reserve(hb_list_in.size());

  for (size_t i = 0; i < hb_list_in.size(); ++i) {
    TH1* hb0 = hb_list_in[i];

    if (!hb0) continue;

    TH1* hb = (TH1*)hb0->Clone(Form("%s__draw_bkg_%zu__%u",
                                    hb0->GetName(),
                                    i,
                                    gRandom->Integer(1000000000)));

    hb->SetDirectory(gROOT);

    Prep(hb, normalize);
    SetNiceLabels(hb, hname, normalize);

    const int col = colors[i % ncol];

    hb->SetFillColor(col);
    hb->SetLineColor(col);
    hb->SetLineWidth(1);
    hb->SetFillStyle(1001);

    st->Add(hb, "HIST");
    hb_drawn.push_back(hb);

    if (!hbSum) {
      hbSum = (TH1*)hb->Clone(Form("hbSum_%s_%u",
                                   hb->GetName(),
                                   gRandom->Integer(1000000000)));
      hbSum->SetDirectory(gROOT);
    }
    else {
      hbSum->Add(hb);
    }
  }

  if (!hbSum) return;

  SetNiceLabels(hbSum, hname, normalize);

  // --------------------------------------------------------------------------
  // Signal style requested by mentor
  // --------------------------------------------------------------------------
  hsA->SetLineWidth(3);
  hsA->SetLineColor(kBlack);
  hsA->SetLineStyle(1);
  hsA->SetFillStyle(0);

  hsB->SetLineWidth(3);
  hsB->SetLineColor(kBlack);
  hsB->SetLineStyle(2);
  hsB->SetFillStyle(0);

  // --------------------------------------------------------------------------
  // Top pad drawing
  // --------------------------------------------------------------------------
  pTop->cd();
  pTop->SetGridx(0);
  pTop->SetGridy(0);

  hbSum->GetXaxis()->SetLabelSize(0.0);
  hbSum->GetXaxis()->SetTitleSize(0.0);

  const double ymax = std::max({
    hsA->GetMaximum(),
    hsB->GetMaximum(),
    hbSum->GetMaximum()
  });

  hbSum->SetMaximum(1.35 * ymax);
  hbSum->SetMinimum(0.0);

  ApplyNiceXRange(hbSum, hname);

  hbSum->Draw("HIST");
  st->Draw("HIST SAME");
  hsA->Draw("HIST SAME");
  hsB->Draw("HIST SAME");

  hbSum->Draw("AXIS SAME");

  TLegend* leg = new TLegend(0.62, 0.60, 0.92, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.032);

  leg->AddEntry(hsA, sigALabel, "l");
  leg->AddEntry(hsB, sigBLabel, "l");

  for (size_t i = 0; i < hb_drawn.size(); ++i) {
    const std::string& lab =
      (i < bkgLabels.size()) ? bkgLabels[i] : std::string(Form("Bkg%zu", i));

    leg->AddEntry(hb_drawn[i], lab.c_str(), "f");
  }

  leg->Draw();

  // --------------------------------------------------------------------------
  // Ratio pad drawing
  // --------------------------------------------------------------------------
  pBot->cd();
  pBot->SetGridx(1);
  pBot->SetGridy(1);

  auto makeRatio = [&](TH1* hs, const char* tag) {
    TH1* hr = (TH1*)hs->Clone(Form("%s__ratio_%s__%u",
                                   hs->GetName(),
                                   tag,
                                   gRandom->Integer(1000000000)));

    hr->SetDirectory(gROOT);

    for (int ib = 1; ib <= hr->GetNbinsX(); ++ib) {
      const double b  = hbSum->GetBinContent(ib);
      const double s  = hs->GetBinContent(ib);
      const double es = hs->GetBinError(ib);

      if (b > 0) {
        hr->SetBinContent(ib, s / b);
        hr->SetBinError(ib, es / b);
      }
      else {
        hr->SetBinContent(ib, 0.0);
        hr->SetBinError(ib, 0.0);
      }
    }

    hr->SetTitle("");
    hr->GetYaxis()->SetTitle("S/#SigmaB");
    hr->GetYaxis()->CenterTitle(true);

    hr->GetYaxis()->SetNdivisions(505);
    hr->GetYaxis()->SetTitleSize(0.12);
    hr->GetYaxis()->SetLabelSize(0.10);
    hr->GetYaxis()->SetTitleOffset(0.50);

    hr->GetXaxis()->SetTitleSize(0.14);
    hr->GetXaxis()->SetLabelSize(0.12);
    hr->GetXaxis()->SetTitleOffset(1.10);
    hr->GetXaxis()->SetTitle(hbSum->GetXaxis()->GetTitle());

    ApplyNiceXRange(hr, hname);

    hr->SetMinimum(0.0);
    hr->SetMaximum(2.0);

    return hr;
  };

  TH1* rA = makeRatio(hsA, "A");
  TH1* rB = makeRatio(hsB, "B");

  // Ratio styles match the signal styles
  rA->SetLineColor(kBlack);
  rA->SetLineStyle(1);
  rA->SetLineWidth(2);
  rA->SetMarkerStyle(20);
  rA->SetMarkerSize(0.6);
  rA->SetMarkerColor(kBlack);

  rB->SetLineColor(kBlack);
  rB->SetLineStyle(2);
  rB->SetLineWidth(2);
  rB->SetMarkerStyle(24);
  rB->SetMarkerSize(0.6);
  rB->SetMarkerColor(kBlack);

  rA->Draw("EPL");
  rB->Draw("EPL SAME");

  rA->Draw("AXIS SAME");

  c->Modified();
  c->Update();
}

// ----------------------------------------------------------------------------
// Helper to build and save one histogram plot.
// ----------------------------------------------------------------------------
static void MakeOnePlotPDF(const std::string& hname,
                           TFile* fs,
                           const std::string& signalTag,
                           TFile* fs30,
                           const std::vector<std::pair<BkgInfo,TFile*>>& obkgs,
                           bool normalize,
                           int ibin) {
  const std::string sigLabel   = "Signal ttH (mA=12)";
  const std::string sig30Tag   = "tth30gev";
  const std::string sig30Label = "Signal ttH (mA=30)";

  TH1* hs_original =
    CloneToROOT(fs, std::string("SIG_") + signalTag, hname.c_str());

  TH1* hs30_original =
    CloneToROOT(fs30, std::string("SIG_") + sig30Tag, hname.c_str());

  if (!hs_original || !hs30_original) return;

  TH1* hs =
    MakePlotClone(hs_original, hname, std::string("SIG_") + signalTag, ibin);

  TH1* hs30 =
    MakePlotClone(hs30_original, hname, std::string("SIG_") + sig30Tag, ibin);

  if (!hs || !hs30) return;

  std::vector<TH1*> hb_list;
  std::vector<std::string> labels;

  hb_list.reserve(obkgs.size());
  labels.reserve(obkgs.size());

  for (const auto& ob : obkgs) {
    const BkgInfo& info = ob.first;
    TFile* fb = ob.second;

    TH1* hb_original =
      CloneToROOT(fb, std::string("BKG_") + info.tag, hname.c_str());

    if (!hb_original) continue;

    TH1* hb =
      MakePlotClone(hb_original, hname, std::string("BKG_") + info.tag, ibin);

    if (!hb) continue;

    hb_list.push_back(hb);
    labels.push_back(info.label);
  }

  if (hb_list.empty()) {
    std::cout << "WARNING: no backgrounds found for " << hname << " -> skipping.\n";
    return;
  }

  const std::string mode = normalize ? "norm" : "raw";

  const std::string out =
    "plots/" + PrettyKey(hname) + "_" + signalTag + "_plus_tth30gev_" + mode + ".pdf";

  TCanvas* c = new TCanvas(Form("c_%s_%u",
                                hname.c_str(),
                                gRandom->Integer(1000000000)),
                           hname.c_str(),
                           900,
                           800);

  c->SetGridx(0);
  c->SetGridy(0);

  DrawStackWithRatio2Sig(c,
                         hname.c_str(),
                         hs,
                         hs30,
                         hb_list,
                         labels,
                         sigLabel.c_str(),
                         sig30Label.c_str(),
                         normalize);

  c->Modified();
  c->Update();
  c->Print(out.c_str());

  gCanvases.push_back(c);
}

// ----------------------------------------------------------------------------
// Main.
// ----------------------------------------------------------------------------
void overlay_signal_allbkgs(const char* signalTag = "tth12gev",
                            bool normalize = true,
                            bool keepWindowsOpen = true) {
  gStyle->SetOptStat(0);
  gStyle->SetCanvasDefW(900);
  gStyle->SetCanvasDefH(800);

  // Global pads have no grid by default.
  // Ratio pad grid is set explicitly inside DrawStackWithRatio2Sig.
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);

  gSystem->mkdir("plots", kTRUE);

  // Rebin factor.
  // ibin = 6 is the best balance between noisy and over-blocky.
  Int_t ibin = 6;

  const std::string fSig = std::string("output_signal_") + signalTag + ".root";

  TFile* fs = TFile::Open(fSig.c_str(), "READ");

  if (!fs || fs->IsZombie()) {
    std::cout << "ERROR: cannot open signal file " << fSig << "\n";
    return;
  }

  const std::string fSig30 = "output_signal_tth30gev.root";

  TFile* fs30 = TFile::Open(fSig30.c_str(), "READ");

  if (!fs30 || fs30->IsZombie()) {
    std::cout << "ERROR: cannot open signal30 file " << fSig30 << "\n";
    fs->Close();
    delete fs;
    return;
  }

  // Only the three backgrounds requested by the mentor.
  std::vector<BkgInfo> bkgs = {
    {"output_TTH_Hbb.root",   "TTH_Hbb",    "t#bar{t}H, H#rightarrow b#bar{b}"},
    {"output_ttbar.root",     "TT_dilep",   "t#bar{t} dileptonic"},
    {"output_TTtoLNu2Q.root", "TT_semilep", "t#bar{t} semileptonic"}
  };

  std::vector<std::pair<BkgInfo,TFile*>> obkgs;
  obkgs.reserve(bkgs.size());

  for (const auto& b : bkgs) {
    TFile* f = TFile::Open(b.file.c_str(), "READ");

    if (!f || f->IsZombie()) {
      std::cout << "WARNING: cannot open " << b.file << " (skipping)\n";

      if (f) {
        f->Close();
        delete f;
      }

      continue;
    }

    obkgs.push_back({b, f});
  }

  if (obkgs.empty()) {
    std::cout << "ERROR: no background files could be opened.\n";
    fs->Close();
    delete fs;
    fs30->Close();
    delete fs30;
    return;
  }

  // Histogram list.
  // Reconstructed phi plots are removed, as requested by the mentor.
  std::vector<std::string> allH = {
    // Higgs candidate
    "h_Hdbb_mass",
    "h_Hdbb_pt",
    "h_Hdbb_eta",

    // Double-b jets
    "h_dbj1_mass",
    "h_dbj1_pt",
    "h_dbj1_eta",
    "h_dbj2_mass",
    "h_dbj2_pt",
    "h_dbj2_eta",

    // Pair/MET/HT/dilepton variables
    "h_dR_dbj12_final",
    "h_dM_bj12_final",
    "h_MET_pt_final",
    "h_HT",
    "h_mll",

    // Final b-jets
    "h_bj1_pt_final",
    "h_bj1_eta_final",
    "h_bj2_pt_final",
    "h_bj2_eta_final",

    // Lepton and dilepton angular distance
    "h_lep1_pt_final",
    "h_lep1_eta_final",
    "h_dRll"
  };

  for (const auto& hname : allH) {
    MakeOnePlotPDF(hname, fs, signalTag, fs30, obkgs, normalize, ibin);
  }

  for (auto& ob : obkgs) {
    ob.second->Close();
    delete ob.second;
  }

  fs->Close();
  delete fs;

  fs30->Close();
  delete fs30;

  std::cout << "Done. PDFs saved in ./plots/ ("
            << (normalize ? "normalized" : "raw")
            << ", clean-range rebin="
            << ibin
            << ")\n";

  std::cout << "Note: reconstructed phi plots were not regenerated.\n";
  std::cout << "If old phi PDFs still exist in ./plots/, remove them manually before uploading to Overleaf.\n";

  if (!keepWindowsOpen) {
    for (auto* c : gCanvases) {
      if (c) c->Close();
    }

    gCanvases.clear();
  }
}
