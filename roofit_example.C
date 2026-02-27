// roofit_example.C
// Run with: root -l -q roofit_example.C+

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooAbsData.h"

using namespace RooFit;

void roofit_example()  // <-- REQUIRED: name must match the filename
{
  // --- Change these two lines to your file/histogram name ---
const char* fname = "TTH12Gev.root";
const char* hname = "hJetPt_jec";
  // ---------------------------------------------------------

  TFile f(fname, "READ");
  TH1* h = dynamic_cast<TH1*>(f.Get(hname));
  if (!h) {
    Error("roofit_example", "Histogram '%s' not found in file '%s'", hname, fname);
    return;
  }

  // Observable with histogram range
  const double xmin = h->GetXaxis()->GetXmin();
  const double xmax = h->GetXaxis()->GetXmax();
  RooRealVar x("x", "observable", xmin, xmax);
  x.setBins(h->GetNbinsX());
  // TH1 -> RooDataHist (binned dataset)
  RooDataHist dataHist("dataHist", "binned data from TH1", RooArgList(x), Import(*h));

  // Model (RooAbsPdf) from the histogram
  RooHistPdf model("model", "hist-based pdf", RooArgSet(x), dataHist);

  // Generate data from model (p.13)
  std::unique_ptr<RooDataSet> toyUnbinned(model.generate(RooArgSet(x), 20000));
  std::unique_ptr<RooDataHist> toyBinned(model.generateBinned(RooArgSet(x), 20000));

  // Plot with RooFit
  RooPlot* frame = x.frame(Title("TH1 -> RooDataHist, RooHistPdf, generate, plot"));

  dataHist.plotOn(frame, DataError(RooAbsData::Poisson));
  model.plotOn(frame);
  toyUnbinned->plotOn(frame, MarkerStyle(kOpenCircle));
  toyBinned->plotOn(frame, MarkerStyle(kFullTriangleUp));

  TCanvas c("c", "c", 900, 700);
  frame->Draw();
  c.SaveAs("roofit_example.pdf");

  Info("roofit_example", "Saved plot to roofit_example.pdf");
}
