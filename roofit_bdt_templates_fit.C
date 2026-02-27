// File: roofit_bdt_templates_fit.C
#include <iostream>
#include <vector>
#include <memory>
#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace std;

struct Sample {
  std::string name;
  std::string file;
  bool isSignal;
};

static TH1* GetAndRebinBDT(TFile& f, const char* hname, int rebinFactor)
{
  TH1* h = dynamic_cast<TH1*>(f.Get(hname));
  if (!h) return nullptr;

  TH1* hc = dynamic_cast<TH1*>(h->Clone((std::string(h->GetName()) + "_clone").c_str()));
  hc->SetDirectory(nullptr);

  if (rebinFactor > 1) hc->Rebin(rebinFactor);
  return hc;
}

void roofit_bdt_templates_fit()
{
  const char* hname = "h_BDT";
  const int rebinFactor = 5; // 50 bins -> 10 bins

  // backgrounds + signals
  std::vector<Sample> samples = {
    // backgrounds
    {"TTW",       "output_TTW.root",        false},
    {"TTZ",       "output_TTZ.root",        false},
    {"TTH_Hbb",   "output_TTH_Hbb.root",    false},
    {"TTtoLNu2Q", "output_TTtoLNu2Q.root",  false},
    {"ttbar",     "output_ttbar.root",      false},

    // signals
    {"tth12", "output_signal_tth12gev.root", true},
    {"tth15", "output_signal_tth15gev.root", true},
    {"tth20", "output_signal_tth20gev.root", true},
    {"tth25", "output_signal_tth25gev.root", true},
    {"tth30", "output_signal_tth30gev.root", true},
  };

  std::vector<std::unique_ptr<TFile>> files;

  std::unique_ptr<TH1> hSigSum;
  std::unique_ptr<TH1> hBkgSum;

  for (const auto& s : samples) {
    auto tf = std::make_unique<TFile>(s.file.c_str(), "READ");
    if (!tf || tf->IsZombie()) {
      std::cerr << "[ERROR] Cannot open file: " << s.file << "\n";
      return;
    }

    TH1* h = GetAndRebinBDT(*tf, hname, rebinFactor);
    if (!h) {
      std::cerr << "[ERROR] Histogram '" << hname << "' not found in " << s.file << "\n";
      return;
    }

    // protect from negative bins (can happen with weights)
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
      if (h->GetBinContent(b) < 0) h->SetBinContent(b, 0.0);
    }

    if (s.isSignal) {
      if (!hSigSum) {
        hSigSum.reset(dynamic_cast<TH1*>(h->Clone("hSigSum")));
        hSigSum->SetDirectory(nullptr);
        hSigSum->Reset("ICES");
      }
      hSigSum->Add(h);
    } else {
      if (!hBkgSum) {
        hBkgSum.reset(dynamic_cast<TH1*>(h->Clone("hBkgSum")));
        hBkgSum->SetDirectory(nullptr);
        hBkgSum->Reset("ICES");
      }
      hBkgSum->Add(h);
    }

    delete h;
    files.push_back(std::move(tf));
  }

  if (!hSigSum || !hBkgSum) {
    std::cerr << "[ERROR] Could not build summed signal/background histograms.\n";
    return;
  }

  // Optional: total = S + B (dashed)
  std::unique_ptr<TH1> hTot(dynamic_cast<TH1*>(hBkgSum->Clone("hTot")));
  hTot->SetDirectory(nullptr);
  hTot->Add(hSigSum.get());

  // style
  hBkgSum->SetLineColor(kBlue);
  hBkgSum->SetLineWidth(3);

  hSigSum->SetLineColor(kRed);
  hSigSum->SetLineWidth(3);

  hTot->SetLineColor(kBlack);
  hTot->SetLineWidth(2);
  hTot->SetLineStyle(2);

  hBkgSum->GetXaxis()->SetTitle("BDT score");
  hBkgSum->GetYaxis()->SetTitle("Events (weighted)");

  double ymax = std::max(hBkgSum->GetMaximum(), hSigSum->GetMaximum());
  hBkgSum->SetMaximum(1.25 * ymax);

  TCanvas c("c", "c", 1000, 800);

  hBkgSum->Draw("HIST");
  hSigSum->Draw("HIST SAME");
  hTot->Draw("HIST SAME");

  TLegend leg(0.55, 0.70, 0.88, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.AddEntry(hSigSum.get(), "Signal sum (tth12+15+20+25+30)", "l");
  leg.AddEntry(hBkgSum.get(), "Background sum (TTW+TTZ+TTH_Hbb+TTtoLNu2Q+ttbar)", "l");
  leg.AddEntry(hTot.get(),    "Total (S+B)", "l");
  leg.Draw();

  c.SaveAs("roofit_bdt_templates_fit.pdf");

  std::cout << "[OK] Saved: roofit_bdt_templates_fit.pdf\n";
  std::cout << "     Signal integral: " << hSigSum->Integral() << "\n";
  std::cout << "     Bkg integral   : " << hBkgSum->Integral() << "\n";
  std::cout << "     Total integral : " << hTot->Integral() << "\n";
}
