// fit_and_limits3.cc 1 background
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <memory>

#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMatrixDSym.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"

#include "RooConstVar.h"
#include "RooRealVar.h"
#include "RooBinning.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooGaussian.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooMinimizer.h"
#include "RooWorkspace.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooAbsReal.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/MCMCInterval.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/ProposalHelper.h"
#include "RooStats/ProposalFunction.h"
#include "RooRandom.h"
using namespace RooFit;
using namespace RooStats;

// ============================================================
// Options
// ============================================================
struct AnalysisOptions {
  std::vector<int> mass_points = {12, 15, 20, 25, 30};

  int Ntoys_pulls  = 10000;
  int Ntoys_limits = 1001;

  double CL_95 = 0.95;
  int numIters = 10000;
  int numBurnInSteps = 1000;
  int nbins = 5;
  double var_bins[6] = {-1.0, -0.70, -0.4, -0.08, 0.52, 1.0};

  double rel_tt     = 0.15;
  double rel_ttsemi = 0.30;

  bool runToyStudy = false;   // true = run pull study, false = skip it
  int toyStudyMode = 1;       // 0 = existing case: toy from model_0, fit with model_0
                              // 1 = new case: toy from model_1, fit with model_1
  int toyStudySignalMass = 12; // only used for mode 1
};

AnalysisOptions opt;

// ============================================================
// Small helpers
// ============================================================
static TH1* getHistOrDie(TFile* f, const char* hname, const char* tag) {
  if (!f || f->IsZombie()) {
    std::cerr << "[FATAL] File open failed for: " << tag << std::endl;
    return nullptr;
  }
  TH1* h = dynamic_cast<TH1*>(f->Get(hname));
  if (!h) {
    std::cerr << "[FATAL] Histogram '" << hname << "' not found in " << tag << std::endl;
    return nullptr;
  }
  return h;
}

static double safeSigma(double nominal, double rel) {
  const double s = std::fabs(nominal) * rel;
  return (s > 1e-6 ? s : 1.0);
}

static TString signalFileForMass(int mass) {
  return TString::Format("output_signal_tth%dgev.root", mass);
}

static TMatrixDSym scaleCovMatrix(const TMatrixDSym& cov, double scale) {
  TMatrixDSym c(cov);
  c *= (scale * scale);
  return c;
}

// ============================================================
// File container
// ============================================================
struct InputFiles {
  TFile *fsig   = nullptr;
  TFile *ftt    = nullptr;
  TFile *ftt2q  = nullptr;
  TFile *fTTHbb = nullptr;
  TFile *fTTW   = nullptr;
  TFile *fTTZ   = nullptr;

  void open(int mass) {
    fsig   = TFile::Open(signalFileForMass(mass), "READ");
    ftt    = TFile::Open("output_ttbar.root", "READ");
    ftt2q  = TFile::Open("output_TTtoLNu2Q.root", "READ");
    fTTHbb = TFile::Open("output_TTH_Hbb.root", "READ");
    fTTW   = TFile::Open("output_TTW.root", "READ");
    fTTZ   = TFile::Open("output_TTZ.root", "READ");
  }

  void openBackgroundOnly() {
    ftt    = TFile::Open("output_ttbar.root", "READ");
    ftt2q  = TFile::Open("output_TTtoLNu2Q.root", "READ");
    fTTHbb = TFile::Open("output_TTH_Hbb.root", "READ");
    fTTW   = TFile::Open("output_TTW.root", "READ");
    fTTZ   = TFile::Open("output_TTZ.root", "READ");
  }

  void close() {
    if (fsig)   { fsig->Close();   delete fsig;   fsig   = nullptr; }
    if (ftt)    { ftt->Close();    delete ftt;    ftt    = nullptr; }
    if (ftt2q)  { ftt2q->Close();  delete ftt2q;  ftt2q  = nullptr; }
    if (fTTHbb) { fTTHbb->Close(); delete fTTHbb; fTTHbb = nullptr; }
    if (fTTW)   { fTTW->Close();   delete fTTW;   fTTW   = nullptr; }
    if (fTTZ)   { fTTZ->Close();   delete fTTZ;   fTTZ   = nullptr; }
  }
};

// ============================================================
// Histogram container
// ============================================================
struct Histograms {
  TH1* h_sig      = nullptr;
  TH1* h_tt       = nullptr;
  TH1* h_tt2q     = nullptr;
  TH1* h_tthbb    = nullptr;
  TH1* h_ttw      = nullptr;
  TH1* h_ttz      = nullptr;

  TH1* h_sig_r    = nullptr;
  TH1* h_tt_r     = nullptr;
  TH1* h_tt2q_r   = nullptr;
  TH1* h_tthbb_r  = nullptr;
  TH1* h_ttw_r    = nullptr;
  TH1* h_ttz_r    = nullptr;
  TH1* h_ttsemi_r = nullptr;
  TH1* h_bkg_r    = nullptr;

  double Nsig        = 0.0;
  double Ntt         = 0.0;
  double Nttsemi     = 0.0;
  double Nbkg        = 0.0;
  double denominator = 1.0;
};

// ============================================================
// Statistical model container
// ============================================================
struct AnalysisModel {
  RooRealVar* output_BDT = nullptr;
  RooBinning* customBinning = nullptr;

  RooDataHist* dh_sig    = nullptr;
  RooDataHist* dh_tt     = nullptr;
  RooDataHist* dh_ttsemi = nullptr;
  RooDataHist* dh_bkg    = nullptr;

  RooHistPdf* pdf_sig    = nullptr;
  RooHistPdf* pdf_tt     = nullptr;
  RooHistPdf* pdf_ttsemi = nullptr;
  RooHistPdf* pdf_bkg    = nullptr;

  RooRealVar* Nexp_sig    = nullptr;
  RooRealVar* Nexp_tt     = nullptr;
  RooRealVar* Nexp_ttsemi = nullptr;
  RooRealVar* Nexp_bkg    = nullptr;

  RooConstVar* tt_nom     = nullptr;
  RooConstVar* tt_sig     = nullptr;
  RooConstVar* ttsemi_nom = nullptr;
  RooConstVar* ttsemi_sig = nullptr;
  RooConstVar* bkg_nom    = nullptr;
  RooConstVar* bkg_sig    = nullptr;

  RooGaussian* c_tt     = nullptr;
  RooGaussian* c_ttsemi = nullptr;
  RooGaussian* c_bkg    = nullptr;

  RooAddPdf* model_0 = nullptr;
  RooAddPdf* model_1 = nullptr;

  RooProdPdf* total_model_0 = nullptr;
  RooProdPdf* total_model_1 = nullptr;

  RooWorkspace* wr = nullptr;
  ModelConfig* mc  = nullptr;
  ProposalFunction* pf = nullptr;

  void cleanup() {
    delete pf;

    delete mc;
    delete wr;

    delete total_model_1;
    delete total_model_0;
    delete model_1;
    delete model_0;

    delete c_bkg;
    delete c_ttsemi;
    delete c_tt;

    delete bkg_sig;
    delete bkg_nom;
    delete ttsemi_sig;
    delete ttsemi_nom;
    delete tt_sig;
    delete tt_nom;

    delete Nexp_bkg;
    delete Nexp_ttsemi;
    delete Nexp_tt;
    delete Nexp_sig;

    delete pdf_bkg;
    delete pdf_ttsemi;
    delete pdf_tt;
    delete pdf_sig;

    delete dh_bkg;
    delete dh_ttsemi;
    delete dh_tt;
    delete dh_sig;

    delete customBinning;
    delete output_BDT;

    pf = nullptr;
    mc = nullptr;
    wr = nullptr;
    total_model_1 = nullptr;
    total_model_0 = nullptr;
    model_1 = nullptr;
    model_0 = nullptr;
    c_bkg = nullptr;
    c_ttsemi = nullptr;
    c_tt = nullptr;
    bkg_sig = nullptr;
    bkg_nom = nullptr;
    ttsemi_sig = nullptr;
    ttsemi_nom = nullptr;
    tt_sig = nullptr;
    tt_nom = nullptr;
    Nexp_bkg = nullptr;
    Nexp_ttsemi = nullptr;
    Nexp_tt = nullptr;
    Nexp_sig = nullptr;
    pdf_bkg = nullptr;
    pdf_ttsemi = nullptr;
    pdf_tt = nullptr;
    pdf_sig = nullptr;
    dh_bkg = nullptr;
    dh_ttsemi = nullptr;
    dh_tt = nullptr;
    dh_sig = nullptr;
    customBinning = nullptr;
    output_BDT = nullptr;
  }
};

// ============================================================
// Limit result container
// ============================================================
struct LimitSummary {
  std::vector<double> mass_points_d;
  std::vector<double> expected_br_limits;
  std::vector<double> sigma_1_br_down;
  std::vector<double> sigma_1_br_up;
  std::vector<double> sigma_2_br_down;
  std::vector<double> sigma_2_br_up;
  std::vector<double> pseudo_observed_br_limits;
};

// ============================================================
// 1) Construct input histograms and rebin
// ============================================================
bool buildInputHistograms(int mass, InputFiles& files, Histograms& hs) {
  files.open(mass);

  const char* HNAME = "h_BDT";

  hs.h_sig   = getHistOrDie(files.fsig,   HNAME, files.fsig   ? files.fsig->GetName()   : "sig");
  hs.h_tt    = getHistOrDie(files.ftt,    HNAME, files.ftt    ? files.ftt->GetName()    : "ttbar");
  hs.h_tt2q  = getHistOrDie(files.ftt2q,  HNAME, files.ftt2q  ? files.ftt2q->GetName()  : "tt2q");
  hs.h_tthbb = getHistOrDie(files.fTTHbb, HNAME, files.fTTHbb ? files.fTTHbb->GetName() : "tthbb");
  hs.h_ttw   = getHistOrDie(files.fTTW,   HNAME, files.fTTW   ? files.fTTW->GetName()   : "ttw");
  hs.h_ttz   = getHistOrDie(files.fTTZ,   HNAME, files.fTTZ   ? files.fTTZ->GetName()   : "ttz");

  if (!hs.h_sig || !hs.h_tt || !hs.h_tt2q || !hs.h_tthbb || !hs.h_ttw || !hs.h_ttz) {
    std::cerr << "[FATAL] Missing input histogram(s).\n";
    return false;
  }

  hs.h_sig_r   = hs.h_sig  ->Rebin(opt.nbins, "h_sig_r",   opt.var_bins);
  hs.h_tt_r    = hs.h_tt   ->Rebin(opt.nbins, "h_tt_r",    opt.var_bins);
  hs.h_tt2q_r  = hs.h_tt2q ->Rebin(opt.nbins, "h_tt2q_r",  opt.var_bins);
  hs.h_tthbb_r = hs.h_tthbb->Rebin(opt.nbins, "h_tthbb_r", opt.var_bins);
  hs.h_ttw_r   = hs.h_ttw  ->Rebin(opt.nbins, "h_ttw_r",   opt.var_bins);
  hs.h_ttz_r   = hs.h_ttz  ->Rebin(opt.nbins, "h_ttz_r",   opt.var_bins);

  hs.h_ttsemi_r = (TH1*)hs.h_tt2q_r->Clone("h_ttsemi_r");
  hs.h_ttsemi_r->Reset();
  hs.h_ttsemi_r->Add(hs.h_tt2q_r);
  hs.h_ttsemi_r->Add(hs.h_tthbb_r);
  hs.h_ttsemi_r->Add(hs.h_ttw_r);
  hs.h_ttsemi_r->Add(hs.h_ttz_r);

  hs.h_bkg_r = (TH1*)hs.h_tt_r->Clone("h_bkg_r");
  hs.h_bkg_r->Reset();
  hs.h_bkg_r->Add(hs.h_tt_r);
  hs.h_bkg_r->Add(hs.h_tt2q_r);
  hs.h_bkg_r->Add(hs.h_tthbb_r);
  hs.h_bkg_r->Add(hs.h_ttw_r);
  hs.h_bkg_r->Add(hs.h_ttz_r);

  hs.Nsig    = hs.h_sig_r->Integral();
  hs.Ntt     = hs.h_tt_r->Integral();
  hs.Nttsemi = hs.h_ttsemi_r->Integral();

  hs.Nbkg = hs.h_bkg_r->Integral();
  hs.denominator = (hs.Nsig > 0.0 ? hs.Nsig : 1.0);

  return true;
}

// ============================================================
// 1b) Construct input histograms and rebin (background only)
//     COMBINED BACKGROUND = tt + tt2q + tthbb + ttw + ttz
// ============================================================
bool buildInputHistogramsBackgroundOnly(InputFiles& files, Histograms& hs) {
  files.openBackgroundOnly();

  const char* HNAME = "h_BDT";

  hs.h_tt    = getHistOrDie(files.ftt,    HNAME, files.ftt    ? files.ftt->GetName()    : "ttbar");
  hs.h_tt2q  = getHistOrDie(files.ftt2q,  HNAME, files.ftt2q  ? files.ftt2q->GetName()  : "tt2q");
  hs.h_tthbb = getHistOrDie(files.fTTHbb, HNAME, files.fTTHbb ? files.fTTHbb->GetName() : "tthbb");
  hs.h_ttw   = getHistOrDie(files.fTTW,   HNAME, files.fTTW   ? files.fTTW->GetName()   : "ttw");
  hs.h_ttz   = getHistOrDie(files.fTTZ,   HNAME, files.fTTZ   ? files.fTTZ->GetName()   : "ttz");

  if (!hs.h_tt || !hs.h_tt2q || !hs.h_tthbb || !hs.h_ttw || !hs.h_ttz) {
    std::cerr << "[FATAL] Missing background input histogram(s).\n";
    return false;
  }

  hs.h_tt_r    = hs.h_tt   ->Rebin(opt.nbins, "h_tt_r",    opt.var_bins);
  hs.h_tt2q_r  = hs.h_tt2q ->Rebin(opt.nbins, "h_tt2q_r",  opt.var_bins);
  hs.h_tthbb_r = hs.h_tthbb->Rebin(opt.nbins, "h_tthbb_r", opt.var_bins);
  hs.h_ttw_r   = hs.h_ttw  ->Rebin(opt.nbins, "h_ttw_r",   opt.var_bins);
  hs.h_ttz_r   = hs.h_ttz  ->Rebin(opt.nbins, "h_ttz_r",   opt.var_bins);

  hs.h_ttsemi_r = (TH1*)hs.h_tt2q_r->Clone("h_ttsemi_r");
  hs.h_ttsemi_r->Reset();
  hs.h_ttsemi_r->Add(hs.h_tt2q_r);
  hs.h_ttsemi_r->Add(hs.h_tthbb_r);
  hs.h_ttsemi_r->Add(hs.h_ttw_r);
  hs.h_ttsemi_r->Add(hs.h_ttz_r);

  hs.h_bkg_r = (TH1*)hs.h_tt_r->Clone("h_bkg_r");
  hs.h_bkg_r->Reset();
  hs.h_bkg_r->Add(hs.h_tt_r);
  hs.h_bkg_r->Add(hs.h_tt2q_r);
  hs.h_bkg_r->Add(hs.h_tthbb_r);
  hs.h_bkg_r->Add(hs.h_ttw_r);
  hs.h_bkg_r->Add(hs.h_ttz_r);

  hs.Ntt     = hs.h_tt_r->Integral();
  hs.Nttsemi = hs.h_ttsemi_r->Integral();
  hs.Nbkg    = hs.h_bkg_r->Integral();

  return true;
}

// ============================================================
// 1) Construct RooFit model
// ============================================================
bool constructModel(int mass, const Histograms& hs, AnalysisModel& fm) {
  fm.output_BDT = new RooRealVar("output_BDT", "BDT score", -1.0, 1.0);
  fm.customBinning = new RooBinning(opt.nbins, opt.var_bins);
  fm.output_BDT->setBinning(*fm.customBinning, "customBinning");

  fm.dh_sig = new RooDataHist("sig", "sig", *fm.output_BDT, hs.h_sig_r);
  fm.dh_bkg = new RooDataHist("bkg", "bkg", *fm.output_BDT, hs.h_bkg_r);

  fm.pdf_sig = new RooHistPdf("sig_pdf", "sig_pdf", *fm.output_BDT, *fm.dh_sig);
  fm.pdf_bkg = new RooHistPdf("bkg_pdf", "bkg_pdf", *fm.output_BDT, *fm.dh_bkg);

 const double sigRangeScale = (mass == 30 ? 100.0 : 30.0);

 fm.Nexp_sig = new RooRealVar("Nexp_sig", "Expected signal events",
                             hs.Nsig,
                             -sigRangeScale * std::max(1.0, hs.Nsig),
                              sigRangeScale * std::max(1.0, hs.Nsig));
 
  fm.Nexp_bkg = new RooRealVar("Nexp_bkg", "Expected combined background",
                             hs.Nbkg, -0.5 * hs.Nbkg, 2.0 * hs.Nbkg);
  fm.bkg_nom = new RooConstVar("bkg_nom", "bkg_nom", hs.Nbkg);
  fm.bkg_sig = new RooConstVar("bkg_sig", "bkg_sig", safeSigma(hs.Nbkg, opt.rel_tt));
  fm.c_bkg   = new RooGaussian("c_bkg", "constraint bkg", *fm.Nexp_bkg, *fm.bkg_nom, *fm.bkg_sig);

  fm.model_0 = new RooAddPdf("model_0", "Background-only (extended)",
                             RooArgList(*fm.pdf_bkg),
                             RooArgList(*fm.Nexp_bkg));

  fm.model_1 = new RooAddPdf("model_1", "Signal+Background (extended)",
                             RooArgList(*fm.pdf_sig, *fm.pdf_bkg),
                             RooArgList(*fm.Nexp_sig, *fm.Nexp_bkg));

  // Optional constrained versions:
  // fm.total_model_0 = new RooProdPdf("total_model_0", "b-only with constraints",
  //                                   RooArgList(*fm.model_0, *fm.c_bkg));
  // fm.total_model_1 = new RooProdPdf("total_model_1", "s+b with constraints",
  //                                   RooArgList(*fm.model_1, *fm.c_bkg));

  return true;
}

// ============================================================
// 1b) Construct RooFit model (background only)
// ============================================================
bool constructBackgroundOnlyModel(const Histograms& hs, AnalysisModel& fm) {
  fm.output_BDT = new RooRealVar("output_BDT", "BDT score", -1.0, 1.0);
  fm.customBinning = new RooBinning(opt.nbins, opt.var_bins);
  fm.output_BDT->setBinning(*fm.customBinning, "customBinning");

  fm.dh_bkg = new RooDataHist("bkg", "bkg", *fm.output_BDT, hs.h_bkg_r);
  fm.pdf_bkg = new RooHistPdf("bkg_pdf", "bkg_pdf", *fm.output_BDT, *fm.dh_bkg);

  fm.Nexp_bkg = new RooRealVar("Nexp_bkg", "Expected combined background",
			       hs.Nbkg, -0.5 * hs.Nbkg, 2.0 * hs.Nbkg);
  fm.bkg_nom = new RooConstVar("bkg_nom", "bkg_nom", hs.Nbkg);
  fm.bkg_sig = new RooConstVar("bkg_sig", "bkg_sig", safeSigma(hs.Nbkg, opt.rel_tt));
  fm.c_bkg   = new RooGaussian("c_bkg", "constraint bkg", *fm.Nexp_bkg, *fm.bkg_nom, *fm.bkg_sig);

  fm.model_0 = new RooAddPdf("model_0", "Background-only (extended)",
                             RooArgList(*fm.pdf_bkg),
                             RooArgList(*fm.Nexp_bkg));

  return true;
}

// ============================================================
// 2) Generate pseudo-data
// ============================================================
RooDataHist* generateAsimovB(AnalysisModel& fm) {
  return fm.model_0->generateBinned(
    *fm.output_BDT,
    RooFit::ExpectedData(true),
    RooFit::Name("asimov_B")
  );
}

RooDataHist* generateToyB(AnalysisModel& fm, const char* name = "toy_B") {
  return fm.model_0->generateBinned(
    *fm.output_BDT,
    RooFit::Extended(true),
    RooFit::Name(name)
  );
}

//Observed

RooDataHist* generateToyBWithObservedCount(AnalysisModel& fm,
                                           int mass,
                                           double Nexp_nominal,
                                           const char* name = "toy_Bobs") {
  const int Nobs = RooRandom::randomGenerator()->Poisson(Nexp_nominal);

  std::cout << "[OBSERVED-LIKE TOY] m_a=" << mass
            << "  Nexp_nominal=" << Nexp_nominal
            << "  Nobs_poisson=" << Nobs << "\n";

  return fm.model_0->generateBinned(
    *fm.output_BDT,
    RooFit::NumEvents(Nobs),
    RooFit::Name(name)
  );
}

RooDataHist* generateAsimovSB(AnalysisModel& fm) {
  return fm.model_1->generateBinned(
    *fm.output_BDT,
    RooFit::ExpectedData(true),
    RooFit::Name("asimov_SB")
  );
}

RooDataHist* generateToySB(AnalysisModel& fm, const char* name = "toy_SB") {
  return fm.model_1->generateBinned(
    *fm.output_BDT,
    RooFit::Extended(true),
    RooFit::Name(name)
  );
}

// ============================================================
// Fit helpers
// ============================================================
RooFitResult* fitSBModel(AnalysisModel& fm, RooDataHist& data) {
  return fm.model_1->fitTo(
    data,
    Save(),
    Extended(kTRUE),
    SumW2Error(kFALSE),
    RooFit::MaxCalls(20000),
    Strategy(2),
    RooFit::Optimize(kTRUE)
  );
}
RooFitResult* fitBOnlyModel(AnalysisModel& fm, RooDataHist& data) {
  return fm.model_0->fitTo(
    data,
    Save(),
    Extended(kTRUE),
    SumW2Error(kFALSE),
    RooFit::MaxCalls(10000),
    Strategy(1),
    RooFit::Optimize(kTRUE)
  );
}

// ============================================================
// Correlation matrix plot helper
// ============================================================
void saveCorrelationPlot(RooFitResult* fitres, const TString& outname, const char* title) {
  if (!fitres) return;

  TH2* hcorr = fitres->correlationHist();
  if (!hcorr) return;

  TCanvas* c = new TCanvas(TString::Format("c_%s", outname.Data()), title, 900, 800);
  c->SetRightMargin(0.15);
  c->SetLeftMargin(0.15);
  c->SetBottomMargin(0.15);

  hcorr->SetTitle(title);
  hcorr->GetZaxis()->SetRangeUser(-1.0, 1.0);
  hcorr->LabelsOption("v", "X");
  hcorr->SetMarkerSize(1.5);
  gStyle->SetPaintTextFormat(".2f");
  hcorr->Draw("COLZ TEXT");

  c->SaveAs(outname);

  delete c;
}

// ============================================================
// 3) Print one-fit outputs
// ============================================================
void printOneFitSummary(int mass, const Histograms& hs, AnalysisModel& fm, RooFitResult* fitres) {
  std::cout << "\n============================================================\n";
  std::cout << " ONE-FIT SUMMARY  (m_a = " << mass << " GeV)\n";
  std::cout << "============================================================\n";

  if (!fitres) {
    std::cout << "[ERROR] Null fit result.\n";
    return;
  }

  std::cout << "status   = " << fitres->status()  << "\n";
  std::cout << "covQual  = " << fitres->covQual() << "\n";
  std::cout << "EDM      = " << fitres->edm()     << "\n";
  std::cout << "minNll   = " << fitres->minNll()  << "\n";
  std::cout << "numInvalidNLL = " << fitres->numInvalidNLL() << "\n";

  std::cout << "\n--- Parameter of interest in this fit ---\n";
  std::cout << "Nexp_bkg = " << fm.Nexp_bkg->getVal() << " +/- " << fm.Nexp_bkg->getError() << "\n";

  std::cout << "\n--- Nominal value ---\n";
  std::cout << "true combined background nominal = " << hs.Nbkg << "\n";

  const double pull_bkg =
    (fm.Nexp_bkg->getError() > 0.0) ? (hs.Nbkg - fm.Nexp_bkg->getVal()) / fm.Nexp_bkg->getError() : 0.0;

  std::cout << "\n--- Pull ---\n";
  std::cout << "pull(combined background) = " << pull_bkg << "\n";

  std::cout << "\n--- Final floating parameters ---\n";
  fitres->floatParsFinal().Print("v");

  std::cout << "\n--- Correlation matrix ---\n";
  fitres->correlationMatrix().Print();

  std::cout << "============================================================\n";
}

void printOneFitSummarySB(int mass, const Histograms& hs, AnalysisModel& fm, RooFitResult* fitres) {
  std::cout << "\n============================================================\n";
  std::cout << " ONE-FIT SUMMARY S+B  (m_a = " << mass << " GeV)\n";
  std::cout << "============================================================\n";

  if (!fitres) {
    std::cout << "[ERROR] Null fit result.\n";
    return;
  }

  std::cout << "status   = " << fitres->status()  << "\n";
  std::cout << "covQual  = " << fitres->covQual() << "\n";
  std::cout << "EDM      = " << fitres->edm()     << "\n";
  std::cout << "minNll   = " << fitres->minNll()  << "\n";
  std::cout << "numInvalidNLL = " << fitres->numInvalidNLL() << "\n";

  std::cout << "\n--- Parameters in this fit ---\n";
  std::cout << "Nexp_sig = " << fm.Nexp_sig->getVal() << " +/- " << fm.Nexp_sig->getError() << "\n";
  std::cout << "Nexp_bkg = " << fm.Nexp_bkg->getVal() << " +/- " << fm.Nexp_bkg->getError() << "\n";

  std::cout << "\n--- Nominal values ---\n";
  std::cout << "true signal nominal           = " << hs.Nsig << "\n";
  std::cout << "true total background nominal = " << hs.Nbkg << "\n";

  const double pull_sig =
    (fm.Nexp_sig->getError() > 0.0) ? (hs.Nsig - fm.Nexp_sig->getVal()) / fm.Nexp_sig->getError() : 0.0;

  const double pull_bkg =
    (fm.Nexp_bkg->getError() > 0.0) ? (hs.Nbkg - fm.Nexp_bkg->getVal()) / fm.Nexp_bkg->getError() : 0.0;

  std::cout << "\n--- Pulls ---\n";
  std::cout << "pull(signal)           = " << pull_sig << "\n";
  std::cout << "pull(total background) = " << pull_bkg << "\n";

  std::cout << "\n--- Final floating parameters ---\n";
  fitres->floatParsFinal().Print("v");

  std::cout << "\n--- Correlation matrix ---\n";
  fitres->correlationMatrix().Print();

  std::cout << "============================================================\n";
}

// ============================================================
// 3) Plot one toy fit
// ============================================================
void plotOneToyFit(int mass, AnalysisModel& fm, RooDataHist& toyData) {
  TCanvas* c_bdt = new TCanvas(TString::Format("c_bdt_ma%d", mass), "Pseudodata + fit", 900, 800);
  RooPlot* fr = fm.output_BDT->frame();

  toyData.plotOn(fr, Name("toydata"));
  fm.model_0->plotOn(fr, Name("fullfit"), LineColor(kBlue));

  fm.model_0->plotOn(fr, Components("bkg_pdf"),
                     LineColor(kRed + 1), LineStyle(kSolid), Name("bkgcomp"));

  fr->SetTitle(TString::Format("One pseudodata fit (B-only), m_{a}=%d GeV", mass));
  fr->GetXaxis()->SetTitle("BDT score");
  fr->GetYaxis()->SetTitle("Events");

  c_bdt->cd();
  c_bdt->SetLogy(0);
  fr->SetMinimum(0.0);
  fr->Draw();

  TLegend* leg = new TLegend(0.56, 0.60, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(fr->findObject("toydata"), "Pseudodata", "pe");
  leg->AddEntry(fr->findObject("fullfit"), "B-only fit", "l");
  leg->AddEntry(fr->findObject("bkgcomp"), "Combined background", "l");
  leg->Draw();

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.035);
  lat.DrawLatex(0.15, 0.92, TString::Format("m_{a}=%d GeV", mass));

  c_bdt->SaveAs(TString::Format("fits_and_limit_plots/step3_onefit_bonly_bdt_ma%d.pdf", mass));

  c_bdt->cd();
  c_bdt->SetLogy(1);
  fr->SetMinimum(0.1);
  fr->SetMaximum(1000.0);
  fr->Draw();
  leg->Draw();

  TLatex lat2;
  lat2.SetNDC();
  lat2.SetTextSize(0.035);
  lat2.DrawLatex(0.15, 0.92, TString::Format("m_{a}=%d GeV", mass));

  c_bdt->SaveAs(TString::Format("fits_and_limit_plots/step3_onefit_bonly_bdt_ma%d_log.pdf", mass));
  delete leg;
  delete fr;
  delete c_bdt;
}

void plotOneToyFitSB(int mass, AnalysisModel& fm, RooDataHist& toyData) {
  TCanvas* c_bdt = new TCanvas(TString::Format("c_bdt_sb_ma%d", mass), "Pseudodata + fit S+B", 900, 800);
  RooPlot* fr = fm.output_BDT->frame();

  toyData.plotOn(fr, Name("toydata"));
  fm.model_1->plotOn(fr, Name("fullfit"), LineColor(kBlue));

  fm.model_1->plotOn(fr, Components("sig_pdf"),
                     LineColor(kRed + 1), LineStyle(kSolid), Name("sigcomp"));
  fm.model_1->plotOn(fr, Components("bkg_pdf"),
                     LineColor(kGreen + 2), LineStyle(kSolid), Name("bkgcomp"));

  fr->SetTitle(TString::Format("One pseudodata fit (S+B), m_{a}=%d GeV", mass));
  fr->GetXaxis()->SetTitle("BDT score");
  fr->GetYaxis()->SetTitle("Events");

  c_bdt->cd();
  c_bdt->SetLogy(0);
  fr->SetMinimum(0.0);
  fr->Draw();

  TLegend* leg = new TLegend(0.54, 0.60, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(fr->findObject("toydata"), "Pseudodata", "pe");
  leg->AddEntry(fr->findObject("fullfit"), "S+B fit", "l");
  leg->AddEntry(fr->findObject("sigcomp"), "Signal", "l");
  leg->AddEntry(fr->findObject("bkgcomp"), "Combined background", "l");
  leg->Draw();

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.035);
  lat.DrawLatex(0.15, 0.92, TString::Format("m_{a}=%d GeV", mass));

  c_bdt->SaveAs(TString::Format("fits_and_limit_plots/step3_onefit_splusb_bdt_ma%d.pdf", mass));

  c_bdt->cd();
  c_bdt->SetLogy(1);
  fr->SetMinimum(0.1);
  fr->SetMaximum(1000.0);
  fr->Draw();
  leg->Draw();

  TLatex lat2;
  lat2.SetNDC();
  lat2.SetTextSize(0.035);
  lat2.DrawLatex(0.15, 0.92, TString::Format("m_{a}=%d GeV", mass));

  c_bdt->SaveAs(TString::Format("fits_and_limit_plots/step3_onefit_splusb_bdt_ma%d_log.pdf", mass));

  delete leg;
  delete fr;
  delete c_bdt;
}

// ============================================================
// 4) Toy MC study for pulls
// ============================================================
void runToyPullStudy(const Histograms& hs, AnalysisModel& fm) {
  if (opt.toyStudyMode == 0) {
    TH1F* h_pullDist_bkg = new TH1F(
      "h_pullDist_bkg_bonly",
      "Combined background pull distribution;Pull;Toys",
      60, -5, 5
    );

    TH1F* h_fitval_bkg = new TH1F(
      "h_fitval_bkg_bonly",
      "Fitted Nexp_bkg;Nexp_bkg;Toys",
      60, 0.0, 2.0 * std::max(1.0, hs.Nbkg)
    );

    TH1F* h_err_bkg = new TH1F(
      "h_err_bkg_bonly",
      "Fit error on Nexp_bkg;#sigma(Nexp_bkg);Toys",
      60, 0.0, std::max(1.0, hs.Nbkg)
    );

    int nFitFail = 0;
    const double Nbkg_nominal = hs.Nbkg;

    std::cout << "\n============================================================\n";
    std::cout << " TOY MC PULL STUDY FROM model_0 FIT WITH model_0\n";
    std::cout << " nToys = " << opt.Ntoys_pulls << "\n";
    std::cout << "============================================================\n";

    for (int it = 0; it < opt.Ntoys_pulls; ++it) {
      fm.Nexp_bkg->setVal(Nbkg_nominal);
      RooDataHist* toyData = generateToyB(fm, TString::Format("toy_pull_%d", it));

      RooFitResult* fitres = fitBOnlyModel(fm, *toyData);

      if (!fitres || fitres->status() != 0 || fitres->covQual() < 2) {
        ++nFitFail;
        delete fitres;
        delete toyData;
        continue;
      }

      const double val_bkg = fm.Nexp_bkg->getVal();
      const double err_bkg = fm.Nexp_bkg->getError();

      const double pull_bkg =
        (err_bkg > 0.0) ? (hs.Nbkg - val_bkg) / err_bkg : 999.0;

      if (std::fabs(pull_bkg) > 3.0) {
        std::cout << "\n[TAIL TOY] it = " << it
                  << "  status = " << fitres->status()
                  << "  covQual = " << fitres->covQual()
                  << "  edm = " << fitres->edm()
                  << "\n"
                  << "  BKG: true = " << hs.Nbkg
                  << "  fit = " << val_bkg
                  << "  err = " << err_bkg
                  << "  pull = " << pull_bkg
                  << "\n";
      }

      h_fitval_bkg->Fill(val_bkg);
      h_err_bkg->Fill(err_bkg);

      if (err_bkg > 0.0) {
        h_pullDist_bkg->Fill((hs.Nbkg - val_bkg) / err_bkg);
      }

      delete fitres;
      delete toyData;
    }

    std::cout << "Number of problematic fits = " << nFitFail << " / " << opt.Ntoys_pulls << "\n";
    std::cout << "Combined background pull: mean = " << h_pullDist_bkg->GetMean()
              << "   RMS = " << h_pullDist_bkg->GetRMS() << "\n";

    {
      TCanvas* c1 = new TCanvas("c_pullDist_bkg_bonly", "Combined background pull", 800, 600);
      h_pullDist_bkg->Fit("gaus", "Q", "");
      h_pullDist_bkg->Draw();
      c1->SaveAs("fits_and_limit_plots/step4_pull_distribution_bkg_bonly.pdf");
      delete c1;
    }

    {
      TCanvas* c2 = new TCanvas("c_fitval_bkg_bonly", "Fitted combined background yield", 800, 600);
      h_fitval_bkg->Draw();
      c2->SaveAs("fits_and_limit_plots/step4_fitval_bkg_bonly.pdf");
      delete c2;
    }

    {
      TCanvas* c3 = new TCanvas("c_err_bkg_bonly", "Error combined background yield", 800, 600);
      h_err_bkg->Draw();
      c3->SaveAs("fits_and_limit_plots/step4_fiterr_bkg_bonly.pdf");
      delete c3;
    }

    delete h_pullDist_bkg;
    delete h_fitval_bkg;
    delete h_err_bkg;
  }

  else if (opt.toyStudyMode == 1) {
    TH1F* h_pullDist_sig = new TH1F(
      "h_pullDist_sig_splusb",
      "Signal pull distribution;Pull;Toys",
      60, -5, 5
    );

    TH1F* h_fitval_sig = new TH1F(
      "h_fitval_sig_splusb",
      "Fitted Nexp_sig;Nexp_sig;Toys",
      60, 0.0, 2.0 * std::max(1.0, hs.Nsig)
    );

    TH1F* h_err_sig = new TH1F(
      "h_err_sig_splusb",
      "Fit error on Nexp_sig;#sigma(Nexp_sig);Toys",
      120, 0.0, 15.0
    );

    TH1F* h_pullDist_bkg = new TH1F(
      "h_pullDist_bkg_splusb",
      "Total background pull distribution;Pull;Toys",
      60, -5, 5
    );

    TH1F* h_fitval_bkg = new TH1F(
      "h_fitval_bkg_splusb",
      "Fitted total background yield;Nexp_bkg;Toys",
      60, 0.0, 2.0 * std::max(1.0, hs.Nbkg)
    );

    TH1F* h_err_bkg = new TH1F(
      "h_err_bkg_splusb",
      "Fit error on total background;#sigma(Nexp_bkg);Toys",
      60, 0.0, std::max(1.0, hs.Nbkg)
    );

    int nFitFail = 0;

    const double Nsig_nominal = hs.Nsig;
    const double Nbkg_nominal = hs.Nbkg;

    std::cout << "\n============================================================\n";
    std::cout << " TOY MC PULL STUDY FROM model_1 FIT WITH model_1\n";
    std::cout << " nToys = " << opt.Ntoys_pulls << "\n";
    std::cout << " signal mass = " << opt.toyStudySignalMass << " GeV\n";
    std::cout << "============================================================\n";

    for (int it = 0; it < opt.Ntoys_pulls; ++it) {
      fm.Nexp_sig->setVal(Nsig_nominal);
      fm.Nexp_bkg->setVal(Nbkg_nominal);

      RooDataHist* toyData = generateToySB(fm, TString::Format("toy_pull_sb_%d", it));
      RooFitResult* fitres = fitSBModel(fm, *toyData);

      if (!fitres || fitres->status() != 0 || fitres->covQual() < 2) {
        ++nFitFail;
        delete fitres;
        delete toyData;
        continue;
      }

      const double val_sig = fm.Nexp_sig->getVal();
      const double err_sig = fm.Nexp_sig->getError();

      const double val_bkg = fm.Nexp_bkg->getVal();
      const double err_bkg = fm.Nexp_bkg->getError();

      const double pull_sig =
        (err_sig > 0.0) ? (hs.Nsig - val_sig) / err_sig : 999.0;

      const double pull_bkg =
        (err_bkg > 0.0) ? (hs.Nbkg - val_bkg) / err_bkg : 999.0;

      if (std::fabs(pull_sig) > 3.0 || std::fabs(pull_bkg) > 3.0) {
        std::cout << "\n[TAIL TOY] it = " << it
                  << "  status = " << fitres->status()
                  << "  covQual = " << fitres->covQual()
                  << "  edm = " << fitres->edm()
                  << "\n"
                  << "  SIG: true = " << hs.Nsig
                  << "  fit = " << val_sig
                  << "  err = " << err_sig
                  << "  pull = " << pull_sig
                  << "\n"
                  << "  BKG: true = " << hs.Nbkg
                  << "  fit = " << val_bkg
                  << "  err = " << err_bkg
                  << "  pull = " << pull_bkg
                  << "\n";
      }

      h_fitval_sig->Fill(val_sig);
      h_err_sig->Fill(err_sig);

      if (err_sig > 0.0) {
        h_pullDist_sig->Fill((hs.Nsig - val_sig) / err_sig);
      }

      h_fitval_bkg->Fill(val_bkg);
      h_err_bkg->Fill(err_bkg);

      if (err_bkg > 0.0) {
        h_pullDist_bkg->Fill((hs.Nbkg - val_bkg) / err_bkg);
      }

      delete fitres;
      delete toyData;
    }

    std::cout << "Number of problematic fits = " << nFitFail << " / " << opt.Ntoys_pulls << "\n";
    std::cout << "Signal pull: mean = " << h_pullDist_sig->GetMean()
              << "   RMS = " << h_pullDist_sig->GetRMS() << "\n";
    std::cout << "Total background pull: mean = " << h_pullDist_bkg->GetMean()
              << "   RMS = " << h_pullDist_bkg->GetRMS() << "\n";

    {
      TCanvas* c1 = new TCanvas("c_pullDist_sig_splusb", "Signal pull", 800, 600);
      h_pullDist_sig->Fit("gaus", "Q", "");
      h_pullDist_sig->Draw();
      c1->SaveAs(TString::Format("fits_and_limit_plots/step4_pull_distribution_sig_splusb_ma%d.pdf",
                                 opt.toyStudySignalMass));
      delete c1;
    }

    {
      TCanvas* c2 = new TCanvas("c_fitval_sig_splusb", "Fitted signal yield", 800, 600);
      h_fitval_sig->Draw();
      c2->SaveAs(TString::Format("fits_and_limit_plots/step4_fitval_sig_splusb_ma%d.pdf",
                                 opt.toyStudySignalMass));
      delete c2;
    }

    {
      TCanvas* c3 = new TCanvas("c_err_sig_splusb", "Error signal yield", 800, 600);
      h_err_sig->Draw();
      c3->SaveAs(TString::Format("fits_and_limit_plots/step4_fiterr_sig_splusb_ma%d.pdf",
                                 opt.toyStudySignalMass));
      delete c3;
    }

    {
      TCanvas* c4 = new TCanvas("c_pullDist_bkg_splusb", "Total background pull", 800, 600);
      h_pullDist_bkg->Fit("gaus", "Q", "");
      h_pullDist_bkg->Draw();
      c4->SaveAs(TString::Format("fits_and_limit_plots/step4_pull_distribution_bkg_splusb_ma%d.pdf",
                                 opt.toyStudySignalMass));
      delete c4;
    }

    {
      TCanvas* c5 = new TCanvas("c_fitval_bkg_splusb", "Fitted total background yield", 800, 600);
      h_fitval_bkg->Draw();
      c5->SaveAs(TString::Format("fits_and_limit_plots/step4_fitval_bkg_splusb_ma%d.pdf",
                                 opt.toyStudySignalMass));
      delete c5;
    }

    {
      TCanvas* c6 = new TCanvas("c_err_bkg_splusb", "Error total background yield", 800, 600);
      h_err_bkg->Draw();
      c6->SaveAs(TString::Format("fits_and_limit_plots/step4_fiterr_bkg_splusb_ma%d.pdf",
                                 opt.toyStudySignalMass));
      delete c6;
    }

    delete h_pullDist_sig;
    delete h_fitval_sig;
    delete h_err_sig;
    delete h_pullDist_bkg;
    delete h_fitval_bkg;
    delete h_err_bkg;
  }
}

// ============================================================
// 5) Build workspace and proposal for limits
// FIXED for one-background model
// ============================================================
void buildWorkspaceAndProposal(AnalysisModel& fm, RooDataHist& data_asimov_B) {
  RooFitResult* fit_asimov = fitSBModel(fm, data_asimov_B);
  if (!fit_asimov) {
    std::cerr << "[ERROR] buildWorkspaceAndProposal: null fit result\n";
    return;
  }

  fm.wr = new RooWorkspace("wr");

  // import the unconstrained S+B model you are actually using
  fm.wr->import(*fm.model_1);

  fm.mc = new ModelConfig("ModelConfig", fm.wr);
  fm.mc->SetPdf(*fm.wr->pdf("model_1"));
  fm.mc->SetParametersOfInterest(RooArgSet(*fm.Nexp_sig));

  RooArgSet nuis;
  nuis.add(*fm.Nexp_bkg);
  fm.mc->SetNuisanceParameters(nuis);

  fm.wr->import(*fm.mc);

  TMatrixDSym cov = fit_asimov->covarianceMatrix();
  TMatrixDSym covSmall = scaleCovMatrix(cov, 0.20);

  RooArgSet proposalVars(*fm.Nexp_sig, *fm.Nexp_bkg);

  ProposalHelper ph;
  ph.SetVariables(proposalVars);
  ph.SetCovMatrix(covSmall);
  ph.SetUpdateProposalParameters(kTRUE);
  ph.SetCacheSize(100);
  fm.pf = ph.GetProposalFunction();

  delete fit_asimov;
}

void makePosteriorPlot(int mass, AnalysisModel& fm, RooDataHist& data_asimov_B) {

  MCMCCalculator post_mcmc(data_asimov_B, *fm.mc);
  post_mcmc.SetProposalFunction(*fm.pf);
  post_mcmc.SetConfidenceLevel(opt.CL_95);
  post_mcmc.SetNumIters(opt.numIters);
  post_mcmc.SetNumBurnInSteps(opt.numBurnInSteps);
  post_mcmc.SetLeftSideTailFraction(0.0);

  MCMCInterval* post_interval = post_mcmc.GetInterval();
  if (!post_interval) {
    std::cerr << "[ERROR] makePosteriorPlot: null MCMC interval for mass " << mass << "\n";
    return;
  }

  TCanvas* c_post = new TCanvas(TString::Format("c_post_ma%d", mass), "Posterior", 800, 700);
  MCMCIntervalPlot postPlot(*post_interval);
  postPlot.SetLineColor(kBlue + 1);
  postPlot.Draw();

  c_post->Update();

  RooPlot* frame = fm.Nexp_sig->frame();
  frame->SetTitle("");
  frame->GetXaxis()->SetTitle("Nexp_sig");
  frame->GetYaxis()->SetTitle("Posterior probability");
  frame->SetMinimum(0.0);
  frame->Draw();

  postPlot.Draw("same");

  TLatex latp;
  latp.SetNDC();
  latp.SetTextSize(0.035);
  latp.DrawLatex(0.16, 0.92, TString::Format("Posterior for m_{a}=%d GeV", mass));
  latp.DrawLatex(0.16, 0.87,
                 TString::Format("95%% upper limit = %.4g",
                                 post_interval->UpperLimit(*fm.Nexp_sig)));

  c_post->Update();
  c_post->SaveAs(TString::Format("fits_and_limit_plots/step5_posterior_ma%d.pdf", mass));
  c_post->SaveAs(TString::Format("fits_and_limit_plots/step5_posterior_ma%d.root", mass));
  delete frame;
  delete c_post;
  delete post_interval;
}
// ============================================================
// 5) Expected limits from toys
// ============================================================
void runExpectedLimitsForMass(int mass, const Histograms& hs, AnalysisModel& fm, LimitSummary& summary) {
  std::vector<double> toy_br_limits;
  toy_br_limits.reserve(opt.Ntoys_limits);

  std::cout << "\n====================================================\n";
  std::cout << " LIMIT STUDY  (m_a = " << mass << " GeV)\n";
  std::cout << " denominator (Nsig) = " << hs.denominator << "\n";
  std::cout << " Nbkg = " << hs.Nbkg << "\n";
  std::cout << "====================================================\n";

  for (int it = 0; it < opt.Ntoys_limits; ++it) {
    RooDataHist* toyData = generateToyB(fm, TString::Format("toy_limit_%d", it));

    MCMCCalculator mcmc(*toyData, *fm.mc);
    mcmc.SetProposalFunction(*fm.pf);
    mcmc.SetConfidenceLevel(opt.CL_95);
    mcmc.SetNumIters(opt.numIters);
    mcmc.SetNumBurnInSteps(opt.numBurnInSteps);
    mcmc.SetLeftSideTailFraction(0.0);

    MCMCInterval* interval = mcmc.GetInterval();

    const double Nsig_up = interval->UpperLimit(*fm.Nexp_sig);
    const double BR_up   = Nsig_up / hs.denominator;

    toy_br_limits.push_back(BR_up);

    if ((it + 1) % 50 == 0) {
      std::cout << "Toy " << (it + 1) << "/" << opt.Ntoys_limits
                << "  Nsig^95 = " << Nsig_up
                << "  BR^95 = " << BR_up << "\n";
    }

    delete interval;
    delete toyData;
  }

  std::sort(toy_br_limits.begin(), toy_br_limits.end());

  auto q = [&](double p) -> double {
    const int idx = std::max(0, std::min(int(std::floor(p * (opt.Ntoys_limits - 1))), opt.Ntoys_limits - 1));
    return toy_br_limits[idx];
  };

  const double median = q(0.50);

  summary.mass_points_d.push_back((double)mass);
  summary.expected_br_limits.push_back(median);
  summary.sigma_1_br_down.push_back(q(0.16));
  summary.sigma_1_br_up.push_back(q(0.84));
  summary.sigma_2_br_down.push_back(q(0.025));
  summary.sigma_2_br_up.push_back(q(0.975));

  std::cout << "\n[RESULT] m_a=" << mass
            << "  expected BR^95 median=" << median
            << "  (-1σ=" << q(0.16) << ", +1σ=" << q(0.84) << ")"
            << "  (-2σ=" << q(0.025) << ", +2σ=" << q(0.975) << ")\n";
}

void runPseudoObservedLimitForMass(int mass, const Histograms& hs, AnalysisModel& fm, LimitSummary& summary) {
 RooDataHist* pseudoData =
  generateToyBWithObservedCount(fm, mass, hs.Nbkg, TString::Format("pseudo_obs_%d", mass));
  
  MCMCCalculator mcmc(*pseudoData, *fm.mc);
  mcmc.SetProposalFunction(*fm.pf);
  mcmc.SetConfidenceLevel(opt.CL_95);
  mcmc.SetNumIters(opt.numIters);
  mcmc.SetNumBurnInSteps(opt.numBurnInSteps);
  mcmc.SetLeftSideTailFraction(0.0);

  MCMCInterval* interval = mcmc.GetInterval();
  if (!interval) {
    std::cerr << "[ERROR] pseudo-observed interval failed for mass " << mass << "\n";
    summary.pseudo_observed_br_limits.push_back(-999.0);
    delete pseudoData;
    return;
  }

  const double Nsig_up = interval->UpperLimit(*fm.Nexp_sig);
  const double BR_up   = Nsig_up / hs.denominator;

  summary.pseudo_observed_br_limits.push_back(BR_up);

std::cout << "[PSEUDO-OBSERVED] m_a=" << mass
          << "  Nbkg_nominal=" << hs.Nbkg
          << "  Nsig^95=" << Nsig_up
          << "  BR^95=" << BR_up << "\n";
  delete interval;
  delete pseudoData;
}

// ============================================================
// 5) Final expected-limit plot
// ============================================================
void makeFinalLimitPlot(const LimitSummary& summary) {
  const int n = (int)summary.mass_points_d.size();
  if (n <= 0) return;

  std::vector<double> x(n), y(n), yobs(n);
  std::vector<double> exl(n, 0.0), exh(n, 0.0);
  std::vector<double> eyl1(n), eyh1(n), eyl2(n), eyh2(n);

  for (int i = 0; i < n; ++i) {
    x[i] = summary.mass_points_d[i];
    y[i] = summary.expected_br_limits[i];
    yobs[i] = summary.pseudo_observed_br_limits[i];

    eyl1[i] = summary.expected_br_limits[i] - summary.sigma_1_br_down[i];
    eyh1[i] = summary.sigma_1_br_up[i]      - summary.expected_br_limits[i];

    eyl2[i] = summary.expected_br_limits[i] - summary.sigma_2_br_down[i];
    eyh2[i] = summary.sigma_2_br_up[i]      - summary.expected_br_limits[i];
  }

  TCanvas* c_lim = new TCanvas("c_lim", "Expected Bayesian limits", 900, 700);

  TGraphAsymmErrors* g2 = new TGraphAsymmErrors(
    n, &x[0], &y[0], &exl[0], &exh[0], &eyl2[0], &eyh2[0]
  );
  TGraphAsymmErrors* g1 = new TGraphAsymmErrors(
    n, &x[0], &y[0], &exl[0], &exh[0], &eyl1[0], &eyh1[0]
  );
  TGraph* gmed = new TGraph(n, &x[0], &y[0]);
  TGraph* gpseudo = new TGraph(n, &x[0], &yobs[0]);

  g2->SetFillColor(kYellow);
  g2->SetLineColor(kYellow);
  g2->SetTitle(";m_{a} [GeV];Bayesian 95% upper limit on BR");

  g1->SetFillColor(kGreen + 1);
  g1->SetLineColor(kGreen + 1);

  gmed->SetLineColor(kBlack);
  gmed->SetLineWidth(2);
  gmed->SetLineStyle(2);
  gmed->SetMarkerStyle(20);
  gmed->SetMarkerSize(1.1);

  gpseudo->SetLineColor(kRed + 1);
  gpseudo->SetLineWidth(2);
  gpseudo->SetMarkerStyle(24);
  gpseudo->SetMarkerSize(1.1);

g2->SetMinimum(-10.0);

g2->Draw("A3");
g1->Draw("3 SAME");
gmed->Draw("LP SAME");
gpseudo->Draw("LP SAME");

TLegend* leglim = new TLegend(0.55, 0.64, 0.88, 0.88);
leglim->SetBorderSize(0);
leglim->SetFillStyle(0);
leglim->AddEntry(gmed, "Median expected", "lp");
leglim->AddEntry(gpseudo, "Pseudo-observed", "lp");
leglim->AddEntry(g1, "#pm1#sigma", "f");
leglim->AddEntry(g2, "#pm2#sigma", "f");
leglim->Draw();

TLatex latlim;
latlim.SetNDC();
latlim.SetTextSize(0.035);
latlim.DrawLatex(0.15, 0.92, "Bayesian 95% upper limits");

c_lim->SaveAs("fits_and_limit_plots/step5_expected_limits_bayesian.pdf");

// log plot
c_lim->Clear();
c_lim->SetLogy();

g2->SetMinimum(1e-1);

g2->Draw("A3");
g1->Draw("3 SAME");
gmed->Draw("LP SAME");
gpseudo->Draw("LP SAME");
leglim->Draw();
latlim.DrawLatex(0.15, 0.92, "Bayesian 95% upper limits");

c_lim->SaveAs("fits_and_limit_plots/step5_expected_limits_bayesian_log.pdf");
  
  delete leglim;
  delete gpseudo;
  delete gmed;
  delete g1;
  delete g2;
  delete c_lim;
}
// ============================================================
// Main steering function
// ============================================================
void fit_and_limits() {
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(0);
  gSystem->mkdir("fits_and_limit_plots", kTRUE);

  LimitSummary summary;
std::vector<double> Nexp_sig_nominal_per_mass;
  // ========================================================
  // STEP 1-4: keep your current study flow unchanged
  // ========================================================
  {
    InputFiles files;
    Histograms hs;
    AnalysisModel fm;

    std::cout << "\n\n############################################################\n";
    std::cout << " Background-only / signal+background study\n";
    std::cout << "############################################################\n";

    std::cout << "\n[STEP 1] Construct input histograms and model\n";

    if (opt.toyStudyMode == 1) {
      if (!buildInputHistograms(opt.toyStudySignalMass, files, hs)) {
        std::cerr << "[FATAL] Could not build full input histograms for mass "
                  << opt.toyStudySignalMass << "\n";
        files.close();
        return;
      }
      constructModel(opt.toyStudySignalMass, hs, fm);
    } else {
      if (!buildInputHistogramsBackgroundOnly(files, hs)) {
        std::cerr << "[FATAL] Could not build background-only histograms\n";
        files.close();
        return;
      }
      constructBackgroundOnlyModel(hs, fm);
    }

    std::cout << "[STEP 2] Generate pseudo-data\n";

    RooDataHist* data_asimov_B = nullptr;
    RooDataHist* data_toy_B    = nullptr;

    if (opt.toyStudyMode == 1) {
      data_asimov_B = generateAsimovSB(fm);
      data_toy_B    = generateToySB(fm, "toy_SB");
    } else {
      data_asimov_B = generateAsimovB(fm);
      data_toy_B    = generateToyB(fm, "toy_B");
    }

    std::cout << "[STEP 3] One fit and detailed output\n";
    RooFitResult* fit_one = nullptr;

    if (opt.toyStudyMode == 1) {
      fit_one = fitSBModel(fm, *data_toy_B);
      printOneFitSummarySB(opt.toyStudySignalMass, hs, fm, fit_one);
      plotOneToyFitSB(opt.toyStudySignalMass, fm, *data_toy_B);
      saveCorrelationPlot(
        fit_one,
        TString::Format("fits_and_limit_plots/step3_correlation_matrix_splusb_ma%d.pdf",
                        opt.toyStudySignalMass),
        TString::Format("Correlation matrix (S+B fit), m_{a}=%d GeV",
                        opt.toyStudySignalMass)
      );
    } else {
      fit_one = fitBOnlyModel(fm, *data_toy_B);
      printOneFitSummary(0, hs, fm, fit_one);
      plotOneToyFit(0, fm, *data_toy_B);
      saveCorrelationPlot(
        fit_one,
        "fits_and_limit_plots/step3_correlation_matrix_bonly.pdf",
        "Correlation matrix (B-only fit)"
      );
    }

    if (opt.runToyStudy) {
      std::cout << "[STEP 4] Toy MC pull study\n";
      runToyPullStudy(hs, fm);
    } else {
      std::cout << "[STEP 4] Toy MC pull study skipped\n";
    }

    delete fit_one;
    delete data_asimov_B;
    delete data_toy_B;

    delete hs.h_bkg_r;
    delete hs.h_ttsemi_r;

    fm.cleanup();
    files.close();
  }

  // ========================================================
  // STEP 5: Limits over all masses
  // ========================================================
  std::cout << "\n[STEP 5] Bayesian limits\n";

  for (int mass : opt.mass_points) {
    InputFiles files_lim;
    Histograms hs_lim;
    AnalysisModel fm_lim;

    if (!buildInputHistograms(mass, files_lim, hs_lim)) {
      std::cerr << "[FATAL] Could not build input histograms for limit mass " << mass << "\n";
      files_lim.close();
      continue;
    }

    constructModel(mass, hs_lim, fm_lim);
    Nexp_sig_nominal_per_mass.push_back(hs_lim.Nsig);
    RooDataHist* data_asimov_B = generateAsimovB(fm_lim);

    buildWorkspaceAndProposal(fm_lim, *data_asimov_B);

   if (fm_lim.mc && fm_lim.pf) {
  makePosteriorPlot(mass, fm_lim, *data_asimov_B);
  runExpectedLimitsForMass(mass, hs_lim, fm_lim, summary);
  runPseudoObservedLimitForMass(mass, hs_lim, fm_lim, summary);
} else {
  std::cerr << "[ERROR] Workspace/proposal not built for mass " << mass << "\n";
}

    delete data_asimov_B;
    delete hs_lim.h_bkg_r;
    delete hs_lim.h_ttsemi_r;

    fm_lim.cleanup();
    files_lim.close();
  }

  std::cout << "\nexpected_br_limits {";
  for (size_t i = 0; i < summary.expected_br_limits.size(); ++i) {
    std::cout << summary.expected_br_limits[i]
              << (i + 1 < summary.expected_br_limits.size() ? ", " : "");
  }
  std::cout << "}\n";

  std::cout << "sigma_1_br_down {";
  for (size_t i = 0; i < summary.sigma_1_br_down.size(); ++i) {
    std::cout << summary.sigma_1_br_down[i]
              << (i + 1 < summary.sigma_1_br_down.size() ? ", " : "");
  }
  std::cout << "}\n";

  std::cout << "sigma_1_br_up {";
  for (size_t i = 0; i < summary.sigma_1_br_up.size(); ++i) {
    std::cout << summary.sigma_1_br_up[i]
              << (i + 1 < summary.sigma_1_br_up.size() ? ", " : "");
  }
  std::cout << "}\n";

  std::cout << "sigma_2_br_down {";
  for (size_t i = 0; i < summary.sigma_2_br_down.size(); ++i) {
    std::cout << summary.sigma_2_br_down[i]
              << (i + 1 < summary.sigma_2_br_down.size() ? ", " : "");
  }
  std::cout << "}\n";

  std::cout << "sigma_2_br_up {";
  for (size_t i = 0; i < summary.sigma_2_br_up.size(); ++i) {
    std::cout << summary.sigma_2_br_up[i]
              << (i + 1 < summary.sigma_2_br_up.size() ? ", " : "");
  }
  std::cout << "}\n";
std::cout << "\nNexp_sig_nominal {";
for (size_t i = 0; i < Nexp_sig_nominal_per_mass.size(); ++i) {
  std::cout << Nexp_sig_nominal_per_mass[i]
            << (i + 1 < Nexp_sig_nominal_per_mass.size() ? ", " : "");
}
std::cout << "}\n";
 
 std::cout << "pseudo_observed_br_limits {";
for (size_t i = 0; i < summary.pseudo_observed_br_limits.size(); ++i) {
  std::cout << summary.pseudo_observed_br_limits[i]
            << (i + 1 < summary.pseudo_observed_br_limits.size() ? ", " : "");
}
std::cout << "}\n";
  makeFinalLimitPlot(summary);
}
