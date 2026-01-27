// overlay_signal_allbkgs.C
// Mode-2 overlays: Signal vs ALL backgrounds at once using THStack,
// plus ratio pad S / (sum backgrounds).
//
// One histogram per canvas, saved as PDF into ./plots/
//
// Run:
//   .L overlay_signal_allbkgs.C
//   overlay_signal_allbkgs("tth12gev", true,  true);   // shapes (unit area), keep windows
//   overlay_signal_allbkgs("tth60gev", false, false);  // raw entries, close windows
//
// Notes:
// - If normalize=true: each component (signal + each bkg) is unit-area normalized (shape compare).
// - If normalize=false: plots raw entries from each file (only meaningful if your analysis filled weighted hists).
// - Histograms are cloned into gROOT so files can be closed safely.

#include <TFile.h>
#include <TH1.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TString.h>
#include <TRandom.h>
#include <TSystem.h>

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

// ----------------------------------------------------------------------------
// Clone histogram out of file into gROOT and cache it
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

  const TString uniqName = Form("%s__keep__%s__%lld", hname, keyPrefix.c_str(), ctr++);
  TH1* hc = (TH1*)h->Clone(uniqName);
  hc->SetDirectory(gROOT);

  cache[key] = hc;
  return hc;
}

static void Prep(TH1* h, bool normalize) {
  if (!h) return;
  h->SetStats(0);
  if (normalize) {
    const double I = h->Integral();
    if (I > 0) h->Scale(1.0 / I);
  }
}

// ----------------------------------------------------------------------------
// Nice axis labels/titles for each histogram key (hname)
// ----------------------------------------------------------------------------
static void SetNiceLabels(TH1* h, const std::string& hname, bool normalize) {
  if (!h) return;

  h->SetTitle("");
  h->GetYaxis()->SetTitle(normalize ? "A.U." : "Entries");

  if (hname == "h_Hdbb_mass")      h->GetXaxis()->SetTitle("m_{bb} (Higgs candidate) [GeV]");
  else if (hname == "h_Hdbb_pt")   h->GetXaxis()->SetTitle("p_{T}(Higgs candidate) [GeV]");
  else if (hname == "h_Hdbb_eta")  h->GetXaxis()->SetTitle("#eta(Higgs candidate)");

  else if (hname == "h_dbj1_mass")      h->GetXaxis()->SetTitle("m(double-b jet 1) [GeV]");
  else if (hname == "h_dbj1_pt")        h->GetXaxis()->SetTitle("p_{T}(double-b jet 1) [GeV]");
  else if (hname == "h_dbj1_eta")       h->GetXaxis()->SetTitle("#eta(double-b jet 1)");
  else if (hname == "h_dbj1_phi_final") h->GetXaxis()->SetTitle("#phi(double-b jet 1)");

  else if (hname == "h_dbj2_mass")      h->GetXaxis()->SetTitle("m(double-b jet 2) [GeV]");
  else if (hname == "h_dbj2_pt")        h->GetXaxis()->SetTitle("p_{T}(double-b jet 2) [GeV]");
  else if (hname == "h_dbj2_eta")       h->GetXaxis()->SetTitle("#eta(double-b jet 2)");
  else if (hname == "h_dbj2_phi_final") h->GetXaxis()->SetTitle("#phi(double-b jet 2)");

  else if (hname == "h_dR_dbj12_final") h->GetXaxis()->SetTitle("#DeltaR(double-b jet 1, double-b jet 2)");
  else if (hname == "h_dM_bj12_final")  h->GetXaxis()->SetTitle("|m(b_{1}) - m(b_{2})| [GeV]");
  else if (hname == "h_MET_pt_final")   h->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  else if (hname == "h_MET_phi_final")  h->GetXaxis()->SetTitle("#phi(E_{T}^{miss})");
  else if (hname == "h_HT")             h->GetXaxis()->SetTitle("H_{T} (scalar #Sigma p_{T}^{jets}) [GeV]");
  else if (hname == "h_mll")            h->GetXaxis()->SetTitle("m_{ll} [GeV]");

  else if (hname == "h_bj1_pt_final")   h->GetXaxis()->SetTitle("p_{T}(b-jet 1) [GeV]");
  else if (hname == "h_bj1_eta_final")  h->GetXaxis()->SetTitle("#eta(b-jet 1)");
  else if (hname == "h_bj2_pt_final")   h->GetXaxis()->SetTitle("p_{T}(b-jet 2) [GeV]");
  else if (hname == "h_bj2_eta_final")  h->GetXaxis()->SetTitle("#eta(b-jet 2)");

  else if (hname == "h_lep1_pt_final")  h->GetXaxis()->SetTitle("p_{T}(leading lepton) [GeV]");
  else if (hname == "h_lep1_eta_final") h->GetXaxis()->SetTitle("#eta(leading lepton)");
  else if (hname == "h_lep1_phi_final") h->GetXaxis()->SetTitle("#phi(leading lepton)");
  else if (hname == "h_dRll")           h->GetXaxis()->SetTitle("#DeltaR(l_{1}, l_{2})");
}

// ----------------------------------------------------------------------------
// Prettier file naming (short “physics meaning” keys)
// ----------------------------------------------------------------------------
static std::string PrettyKey(const std::string& hname) {
  if (hname == "h_Hdbb_mass") return "Higgs_mbb";
  if (hname == "h_Hdbb_pt")   return "Higgs_pt";
  if (hname == "h_Hdbb_eta")  return "Higgs_eta";

  if (hname == "h_dbj1_mass")      return "dbjet1_mass";
  if (hname == "h_dbj1_pt")        return "dbjet1_pt";
  if (hname == "h_dbj1_eta")       return "dbjet1_eta";
  if (hname == "h_dbj1_phi_final") return "dbjet1_phi";

  if (hname == "h_dbj2_mass")      return "dbjet2_mass";
  if (hname == "h_dbj2_pt")        return "dbjet2_pt";
  if (hname == "h_dbj2_eta")       return "dbjet2_eta";
  if (hname == "h_dbj2_phi_final") return "dbjet2_phi";

  if (hname == "h_dR_dbj12_final") return "dR_dbjets";
  if (hname == "h_dM_bj12_final")  return "dM_bjets";
  if (hname == "h_MET_pt_final")   return "MET_pt";
  if (hname == "h_MET_phi_final")  return "MET_phi";
  if (hname == "h_HT")             return "HT";
  if (hname == "h_mll")            return "mll";

  if (hname == "h_bj1_pt_final")   return "bjet1_pt";
  if (hname == "h_bj1_eta_final")  return "bjet1_eta";
  if (hname == "h_bj2_pt_final")   return "bjet2_pt";
  if (hname == "h_bj2_eta_final")  return "bjet2_eta";

  if (hname == "h_lep1_pt_final")  return "lep1_pt";
  if (hname == "h_lep1_eta_final") return "lep1_eta";
  if (hname == "h_lep1_phi_final") return "lep1_phi";
  if (hname == "h_dRll")           return "dR_ll";

  return hname; // fallback
}

// ----------------------------------------------------------------------------
// Draw stacked backgrounds + signal overlay with ratio S/sumB
// (Draws inside the CURRENT canvas; creates its own top/bottom pads)
// ----------------------------------------------------------------------------
static void DrawStackWithRatio(TCanvas* c,
                               const char* hnameKey,
                               TH1* hs_in,
                               const std::vector<TH1*>& hb_list_in,
                               const std::vector<std::string>& bkgLabels,
                               const char* sigLabel,
                               bool normalize) {
  if (!c || !hs_in || hb_list_in.empty() || !hnameKey) return;

  c->cd();

  // Create split pads on the canvas
  TPad* pTop = new TPad(Form("%s_top", c->GetName()), "", 0.0, 0.30, 1.0, 1.0);
  pTop->SetBottomMargin(0.02);
  pTop->SetLeftMargin(0.12);
  pTop->SetRightMargin(0.05);
  pTop->SetTopMargin(0.06);
  pTop->Draw();

  TPad* pBot = new TPad(Form("%s_bot", c->GetName()), "", 0.0, 0.00, 1.0, 0.30);
  pBot->SetTopMargin(0.02);
  pBot->SetBottomMargin(0.35);
  pBot->SetLeftMargin(0.12);
  pBot->SetRightMargin(0.05);
  pBot->Draw();

  // Clone signal for per-draw styling
  TH1* hs = (TH1*)hs_in->Clone(Form("%s__draw_sig__%u", hs_in->GetName(), gRandom->Integer(1e9)));
  hs->SetDirectory(gROOT);
  Prep(hs, normalize);
  SetNiceLabels(hs, hnameKey, normalize);

  // Build stack + sum
  THStack* st = new THStack(Form("st_%s_%u", hs_in->GetName(), gRandom->Integer(1e9)), "");
  TH1* hbSum = nullptr;

  const int colors[] = {kAzure-9, kGreen-7, kOrange-2, kMagenta-7, kCyan-7, kYellow-7,
                        kSpring-7, kViolet-7, kTeal-7, kPink-7, kGray+1};
  const int ncol = sizeof(colors)/sizeof(colors[0]);

  std::vector<TH1*> hb_drawn;
  hb_drawn.reserve(hb_list_in.size());

  for (size_t i = 0; i < hb_list_in.size(); ++i) {
    TH1* hb0 = hb_list_in[i];
    if (!hb0) continue;

    TH1* hb = (TH1*)hb0->Clone(Form("%s__draw_bkg_%zu__%u", hb0->GetName(), i, gRandom->Integer(1e9)));
    hb->SetDirectory(gROOT);
    Prep(hb, normalize);

    const int col = colors[i % ncol];
    hb->SetFillColor(col);
    hb->SetLineColor(col);
    hb->SetLineWidth(1);

    st->Add(hb, "HIST");
    hb_drawn.push_back(hb);

    if (!hbSum) {
      hbSum = (TH1*)hb->Clone(Form("hbSum_%s_%u", hb->GetName(), gRandom->Integer(1e9)));
      hbSum->SetDirectory(gROOT);
    } else {
      hbSum->Add(hb);
    }
  }

  if (!hbSum) return;
  SetNiceLabels(hbSum, hnameKey, normalize);

  // Signal style
  hs->SetLineWidth(3);
  hs->SetLineColor(kRed+1);
  hs->SetFillStyle(0);

  // TOP
  pTop->cd();
  pTop->SetGrid();

  hbSum->GetXaxis()->SetLabelSize(0.0);
  hbSum->GetXaxis()->SetTitleSize(0.0);

  const double ymax = std::max(hs->GetMaximum(), hbSum->GetMaximum());
  hbSum->SetMaximum(1.30 * ymax);
  hbSum->SetMinimum(0.0);

  hbSum->Draw("HIST");
  st->Draw("HIST SAME");
  hs->Draw("HIST SAME");

  // Legend
  TLegend* leg = new TLegend(0.55, 0.50, 0.90, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  leg->AddEntry(hs, sigLabel, "l");
  for (size_t i = 0; i < hb_drawn.size(); ++i) {
    const std::string& lab = (i < bkgLabels.size()) ? bkgLabels[i] : std::string(Form("Bkg%zu", i));
    leg->AddEntry(hb_drawn[i], lab.c_str(), "f");
  }
  leg->Draw();

  // RATIO: S / sumB
  pBot->cd();
  pBot->SetGrid();

  TH1* hr = (TH1*)hs->Clone(Form("%s__ratio__%u", hs->GetName(), gRandom->Integer(1e9)));
  hr->SetDirectory(gROOT);
  hr->Divide(hbSum);

  hr->SetLineColor(kBlack);
  hr->SetLineWidth(2);
  hr->SetMarkerStyle(20);
  hr->SetMarkerSize(0.6);

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

  hr->SetMinimum(0.0);
  hr->SetMaximum(2.0);

  hr->Draw("EP");

  // IMPORTANT: force paint for PDF
  c->Modified();
  c->Update();
}

// ----------------------------------------------------------------------------
// Helper to build & save ONE histogram plot (one canvas per histogram)
// ----------------------------------------------------------------------------
static void MakeOnePlotPDF(const std::string& hname,
                           TFile* fs,
                           const std::string& signalTag,
                           const std::vector<std::pair<BkgInfo,TFile*>>& obkgs,
                           bool normalize,
                           int ibin) {
  const std::string sigLabel = std::string("Signal ") + signalTag;

  TH1* hs = CloneToROOT(fs, std::string("SIG_") + signalTag, hname.c_str());
  if (!hs) return;
  hs->Rebin(ibin);

  std::vector<TH1*> hb_list;
  std::vector<std::string> labels;
  hb_list.reserve(obkgs.size());
  labels.reserve(obkgs.size());

  for (const auto& ob : obkgs) {
    const BkgInfo& info = ob.first;
    TFile* fb = ob.second;

    TH1* hb = CloneToROOT(fb, std::string("BKG_") + info.tag, hname.c_str());
    if (!hb) continue;
    hb->Rebin(ibin);

    hb_list.push_back(hb);
    labels.push_back(info.label);
  }

  if (hb_list.empty()) {
    std::cout << "WARNING: no backgrounds found for " << hname << " -> skipping.\n";
    return;
  }

  const std::string mode = normalize ? "norm" : "raw";
  const std::string out = "plots/" + PrettyKey(hname) + "_" + signalTag + "_" + mode + ".pdf";

  // One canvas per histogram
  TCanvas* c = new TCanvas(Form("c_%s_%u", hname.c_str(), gRandom->Integer(1e9)),
                           hname.c_str(), 900, 800);

  DrawStackWithRatio(c, hname.c_str(), hs, hb_list, labels, sigLabel.c_str(), normalize);

  // Save via the CANVAS (NOT gPad) -> avoids blank PDFs
  c->Modified();
  c->Update();
  c->Print(out.c_str());

  if (gROOT->IsBatch() == kFALSE) {
    // keepWindowsOpen handled by caller (push into gCanvases or delete)
  }
}

// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------
void overlay_signal_allbkgs(const char* signalTag = "tth12gev",
                           bool normalize = true,
                           bool keepWindowsOpen = true) {
  gStyle->SetOptStat(0);

  // Make sure plots/ exists (you said it already exists, but this is harmless)
  gSystem->mkdir("plots", kTRUE);

  Int_t ibin = 4;  // rebin factor for ALL histograms

  // Signal file
  const std::string fSig = std::string("output_signal_") + signalTag + ".root";
  TFile* fs = TFile::Open(fSig.c_str(), "READ");
  if (!fs || fs->IsZombie()) {
    std::cout << "ERROR: cannot open signal file " << fSig << "\n";
    return;
  }

  // Background list (edit freely)
  std::vector<BkgInfo> bkgs = {
    {"output_ttbar.root",            "TT_dilep",        "t#bar{t} dilepton"},
    // {"output_TTtoLNu2Q.root",        "TT_semilep",      "t#bar{t} semileptonic"},
    // {"output_DYee.root",             "DYee",            "DY#rightarrow ee"},
    {"output_DYmumu.root",           "DYmumu",          "DY#rightarrow #mu#mu"},
    {"output_WW.root",               "WW",              "WW"},
    {"output_WZ.root",               "WZ",              "WZ"},
    {"output_ZZ.root",               "ZZ",              "ZZ"},
    {"output_TTH_Hbb.root",          "TTH_Hbb",         "ttH#rightarrow bb"},
    {"output_VBF_Hbb.root",          "VBF_Hbb",         "VBF H#rightarrow bb"},
    {"output_GGH_Hbb.root",          "GGH_Hbb",         "ggH H#rightarrow bb"},
    {"output_TBbarQtoLNu.root",      "TBbarQtoLNu",     "t-channel (LNu)"},
    {"output_TBbarQto2Q.root",       "TBbarQto2Q",      "t-channel (2Q)"},
    {"output_TbarWplusToNu2Q.root",  "TbarWNu2Q",       "#bar{t}W (Nu2Q)"},
    {"output_TbarWplusTo4Q.root",    "TbarW4Q",         "#bar{t}W (4Q)"}
  };

  // Open backgrounds
  std::vector<std::pair<BkgInfo,TFile*>> obkgs;
  obkgs.reserve(bkgs.size());

  for (const auto& b : bkgs) {
    TFile* f = TFile::Open(b.file.c_str(), "READ");
    if (!f || f->IsZombie()) {
      std::cout << "WARNING: cannot open " << b.file << " (skipping)\n";
      if (f) { f->Close(); delete f; }
      continue;
    }
    obkgs.push_back({b, f});
  }

  if (obkgs.empty()) {
    std::cout << "ERROR: no background files could be opened.\n";
    fs->Close(); delete fs;
    return;
  }

  // List ALL histograms you want to plot (one per canvas / one PDF each)
  std::vector<std::string> allH = {
    // Higgs
    "h_Hdbb_mass", "h_Hdbb_pt", "h_Hdbb_eta",

    // double-b jets
    "h_dbj1_mass", "h_dbj1_pt", "h_dbj1_eta", "h_dbj1_phi_final",
    "h_dbj2_mass", "h_dbj2_pt", "h_dbj2_eta", "h_dbj2_phi_final",

    // pair/MET/HT
    "h_dR_dbj12_final", "h_dM_bj12_final", "h_MET_pt_final",
    "h_MET_phi_final", "h_HT", "h_mll",

    // final b-jets
    "h_bj1_pt_final", "h_bj1_eta_final", "h_bj2_pt_final", "h_bj2_eta_final",

    // lepton / dRll
    "h_lep1_pt_final", "h_lep1_eta_final", "h_lep1_phi_final", "h_dRll"
  };

  for (const auto& hname : allH) {
    MakeOnePlotPDF(hname, fs, signalTag, obkgs, normalize, ibin);
  }

  // Close files safely (histograms are cloned to gROOT)
  for (auto& ob : obkgs) {
    ob.second->Close();
    delete ob.second;
  }
  fs->Close();
  delete fs;

  std::cout << "Done. PDFs saved in ./plots/ ("
            << (normalize ? "normalized" : "raw") << ", rebin=" << ibin << ")\n";
}
