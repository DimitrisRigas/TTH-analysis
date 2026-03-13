// fit_and_limits_stepwise.cc
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <memory>

#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
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

using namespace RooFit;
using namespace RooStats;

// ============================================================
// Options
// ============================================================
struct AnalysisOptions {
  std::vector<int> mass_points = {12, 15, 20, 25, 30};

  int Ntoys_pulls  = 1000;
  int Ntoys_limits = 1001;

  double CL_95 = 0.95;
  int numIters = 10000;
  int numBurnInSteps = 0;

  int nbins = 4;
  double var_bins[5] = {-1.0, -0.4, -0.08, 0.52, 1.0};

  double rel_tt     = 0.15;
  double rel_ttsemi = 0.30;
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

  RooHistPdf* pdf_sig    = nullptr;
  RooHistPdf* pdf_tt     = nullptr;
  RooHistPdf* pdf_ttsemi = nullptr;

  RooRealVar* Nexp_sig    = nullptr;
  RooRealVar* Nexp_tt     = nullptr;
  RooRealVar* Nexp_ttsemi = nullptr;

  RooConstVar* tt_nom     = nullptr;
  RooConstVar* tt_sig     = nullptr;
  RooConstVar* ttsemi_nom = nullptr;
  RooConstVar* ttsemi_sig = nullptr;

  RooGaussian* c_tt     = nullptr;
  RooGaussian* c_ttsemi = nullptr;

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

    delete c_ttsemi;
    delete c_tt;

    delete ttsemi_sig;
    delete ttsemi_nom;
    delete tt_sig;
    delete tt_nom;

    delete Nexp_ttsemi;
    delete Nexp_tt;
    delete Nexp_sig;

    delete pdf_ttsemi;
    delete pdf_tt;
    delete pdf_sig;

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

  hs.Nsig    = hs.h_sig_r->Integral();
  hs.Ntt     = hs.h_tt_r->Integral();
  hs.Nttsemi = hs.h_ttsemi_r->Integral();

  hs.Nbkg = hs.Ntt + hs.Nttsemi;
  hs.denominator = (hs.Nsig > 0.0 ? hs.Nsig : 1.0);

  return true;
}

// ============================================================
// 1) Construct RooFit model
// ============================================================
bool constructModel(const Histograms& hs, AnalysisModel& fm) {
  fm.output_BDT = new RooRealVar("output_BDT", "BDT score", -1.0, 1.0);
  fm.customBinning = new RooBinning(opt.nbins, opt.var_bins);
  fm.output_BDT->setBinning(*fm.customBinning, "customBinning");

  fm.dh_sig    = new RooDataHist("sig",    "sig",    *fm.output_BDT, hs.h_sig_r);
  fm.dh_tt     = new RooDataHist("tt",     "tt",     *fm.output_BDT, hs.h_tt_r);
  fm.dh_ttsemi = new RooDataHist("ttsemi", "ttsemi", *fm.output_BDT, hs.h_ttsemi_r);

  fm.pdf_sig    = new RooHistPdf("sig_pdf",    "sig_pdf",    *fm.output_BDT, *fm.dh_sig);
  fm.pdf_tt     = new RooHistPdf("tt_pdf",     "tt_pdf",     *fm.output_BDT, *fm.dh_tt);
  fm.pdf_ttsemi = new RooHistPdf("ttsemi_pdf", "ttsemi_pdf", *fm.output_BDT, *fm.dh_ttsemi);

  fm.Nexp_sig = new RooRealVar("Nexp_sig", "Expected signal events",
                               hs.Nsig, 0.0, 1000.0 * std::max(1.0, hs.Nsig));

  fm.Nexp_tt = new RooRealVar("Nexp_tt", "Expected TT dileptonic",
                              hs.Ntt, 0.0, 10.0 * std::max(1.0, hs.Ntt));

  fm.Nexp_ttsemi = new RooRealVar("Nexp_ttsemi", "Expected TT semileptonic",
                                  hs.Nttsemi, 0.0, 10.0 * std::max(1.0, hs.Nttsemi));

  fm.tt_nom     = new RooConstVar("tt_nom",     "tt_nom",     hs.Ntt);
  fm.tt_sig     = new RooConstVar("tt_sig",     "tt_sig",     safeSigma(hs.Ntt,     opt.rel_tt));
  fm.ttsemi_nom = new RooConstVar("ttsemi_nom", "ttsemi_nom", hs.Nttsemi);
  fm.ttsemi_sig = new RooConstVar("ttsemi_sig", "ttsemi_sig", safeSigma(hs.Nttsemi, opt.rel_ttsemi));

  fm.c_tt     = new RooGaussian("c_tt",     "constraint tt",     *fm.Nexp_tt,     *fm.tt_nom,     *fm.tt_sig);
  fm.c_ttsemi = new RooGaussian("c_ttsemi", "constraint ttsemi", *fm.Nexp_ttsemi, *fm.ttsemi_nom, *fm.ttsemi_sig);

  fm.model_0 = new RooAddPdf("model_0", "Background-only (extended)",
                             RooArgList(*fm.pdf_tt, *fm.pdf_ttsemi),
                             RooArgList(*fm.Nexp_tt, *fm.Nexp_ttsemi));

  fm.model_1 = new RooAddPdf("model_1", "Signal+Background (extended)",
                             RooArgList(*fm.pdf_sig, *fm.pdf_tt, *fm.pdf_ttsemi),
                             RooArgList(*fm.Nexp_sig, *fm.Nexp_tt, *fm.Nexp_ttsemi));

  fm.total_model_0 = new RooProdPdf("total_model_0", "b-only with constraints",
                                    RooArgList(*fm.model_0, *fm.c_tt, *fm.c_ttsemi));

  fm.total_model_1 = new RooProdPdf("total_model_1", "s+b with constraints",
                                    RooArgList(*fm.model_1, *fm.c_tt, *fm.c_ttsemi));

  return true;
}

// ============================================================
// 2) Generate pseudo-data
// ============================================================
RooDataHist* generateAsimovB(AnalysisModel& fm) {
  return fm.total_model_0->generateBinned(
    *fm.output_BDT,
    RooFit::ExpectedData(true),
    RooFit::Name("asimov_B")
  );
}

RooDataHist* generateToyB(AnalysisModel& fm, const char* name = "toy_B") {
  return fm.total_model_0->generateBinned(
    *fm.output_BDT,
    RooFit::Extended(true),
    RooFit::Name(name)
  );
}

// ============================================================
// Fit helpers
// ============================================================
RooFitResult* fitSBModel(AnalysisModel& fm, RooDataHist& data) {
  return fm.total_model_1->fitTo(
    data,
    Save(),
    Extended(kTRUE),
    RooFit::MaxCalls(20000),
    Strategy(2),
    RooFit::Optimize(kTRUE),
    RooFit::SumW2Error(kTRUE)
  );
}

RooFitResult* fitBOnlyModel(AnalysisModel& fm, RooDataHist& data) {
  return fm.total_model_0->fitTo(
    data,
    Save(),
    Extended(kTRUE),
    RooFit::MaxCalls(10000),
    Strategy(1),
    RooFit::Optimize(kTRUE),
    RooFit::SumW2Error(kTRUE)
  );
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

  std::cout << "\n--- Parameters of interest / nuisance parameters ---\n";
  //  std::cout << "Nexp_sig    = " << fm.Nexp_sig->getVal()    << " +/- " << fm.Nexp_sig->getError()    << "\n";
  std::cout << "Nexp_tt     = " << fm.Nexp_tt->getVal()     << " +/- " << fm.Nexp_tt->getError()     << "\n";
  std::cout << "Nexp_ttsemi = " << fm.Nexp_ttsemi->getVal() << " +/- " << fm.Nexp_ttsemi->getError() << "\n";

  std::cout << "\n--- Nominal values ---\n";
  std::cout << "true TTbar nominal           = " << hs.Ntt << "\n";
  std::cout << "true TT semileptonic nominal = " << hs.Nttsemi << "\n";
  std::cout << "signal nominal               = " << hs.Nsig << "\n";

  const double pull_tt =
    (fm.Nexp_tt->getError() > 0.0) ? (hs.Ntt - fm.Nexp_tt->getVal()) / fm.Nexp_tt->getError() : 0.0;
  const double pull_ttsemi =
    (fm.Nexp_ttsemi->getError() > 0.0) ? (hs.Nttsemi - fm.Nexp_ttsemi->getVal()) / fm.Nexp_ttsemi->getError() : 0.0;

  std::cout << "\n--- Pulls ---\n";
  std::cout << "pull(TTbar)        = " << pull_tt << "\n";
  std::cout << "pull(TT semilept.) = " << pull_ttsemi << "\n";

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
  fm.total_model_0->plotOn(fr, Name("fullfit"), LineColor(kBlue));

  fm.total_model_0->plotOn(fr, Components("tt_pdf"),
                           LineColor(kRed + 1), LineStyle(kSolid), Name("ttcomp"));

  fm.total_model_0->plotOn(fr, Components("ttsemi_pdf"),
                           LineColor(kMagenta + 1), LineStyle(kSolid), Name("ttsemicomp"));

  fr->SetTitle(TString::Format("One pseudodata fit (B-only), m_{a}=%d GeV", mass));
  fr->GetXaxis()->SetTitle("BDT score");
  fr->GetYaxis()->SetTitle("Events");
  fr->Draw();

  TLegend* leg = new TLegend(0.56, 0.60, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(fr->findObject("toydata"),    "Pseudodata", "pe");
  leg->AddEntry(fr->findObject("fullfit"),    "B-only fit", "l");
  leg->AddEntry(fr->findObject("ttcomp"),     "t \bar{t} (dil)", "l");
  leg->AddEntry(fr->findObject("ttsemicomp"), "t \bar{t} (semi)", "l");
  leg->Draw();

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.035);
  lat.DrawLatex(0.15, 0.92, TString::Format("m_{a}=%d GeV", mass));

  c_bdt->SaveAs(TString::Format("fits_and_limit_plots/step3_onefit_bonly_bdt_ma%d.pdf", mass));

  delete leg;
  delete fr;
  delete c_bdt;
}

// ============================================================
// 4) Toy MC study for pulls
// ============================================================
void runToyPullStudy(int mass, const Histograms& hs, AnalysisModel& fm) {
  TH1F* h_pullDist_tt = new TH1F(
    TString::Format("h_pullDist_tt_ma%d", mass),
    TString::Format("TTbar pull distribution, m_{a}=%d;Pull;Toys", mass),
    60, -5, 5
  );

  TH1F* h_pullDist_ttsemi = new TH1F(
    TString::Format("h_pullDist_ttsemi_ma%d", mass),
    TString::Format("TT semileptonic pull distribution, m_{a}=%d;Pull;Toys", mass),
    60, -5, 5
  );

  TH1F* h_fitval_tt = new TH1F(
    TString::Format("h_fitval_tt_ma%d", mass),
    TString::Format("Fitted Nexp_tt, m_{a}=%d;Nexp_tt;Toys", mass),
    60, 0.0, 2.0 * std::max(1.0, hs.Ntt)
  );

  TH1F* h_fitval_ttsemi = new TH1F(
    TString::Format("h_fitval_ttsemi_ma%d", mass),
    TString::Format("Fitted Nexp_ttsemi, m_{a}=%d;Nexp_ttsemi;Toys", mass),
    60, 0.0, 2.0 * std::max(1.0, hs.Nttsemi)
  );

  TH1F* h_err_tt = new TH1F(
    TString::Format("h_err_tt_ma%d", mass),
    TString::Format("Fit error on Nexp_tt, m_{a}=%d;#sigma(Nexp_tt);Toys", mass),
    60, 0.0, std::max(1.0, hs.Ntt)
  );

  TH1F* h_err_ttsemi = new TH1F(
    TString::Format("h_err_ttsemi_ma%d", mass),
    TString::Format("Fit error on Nexp_ttsemi, m_{a}=%d;#sigma(Nexp_ttsemi);Toys", mass),
    60, 0.0, std::max(1.0, hs.Nttsemi)
  );

  int nFitFail = 0;

  std::cout << "\n============================================================\n";
  std::cout << " TOY MC PULL STUDY  (m_a = " << mass << " GeV)\n";
  std::cout << " nToys = " << opt.Ntoys_pulls << "\n";
  std::cout << "============================================================\n";

  for (int it = 0; it < opt.Ntoys_pulls; ++it) {
    RooDataHist* toyData = generateToyB(fm, TString::Format("toy_pull_%d", it));

    RooFitResult* fitres = fitBOnlyModel(fm, *toyData);

    if (!fitres || fitres->status() != 0 || fitres->covQual() < 2) {
      ++nFitFail;
    }

    const double val_tt     = fm.Nexp_tt->getVal();
    const double err_tt     = fm.Nexp_tt->getError();
    const double val_ttsemi = fm.Nexp_ttsemi->getVal();
    const double err_ttsemi = fm.Nexp_ttsemi->getError();

    h_fitval_tt->Fill(val_tt);
    h_fitval_ttsemi->Fill(val_ttsemi);
    h_err_tt->Fill(err_tt);
    h_err_ttsemi->Fill(err_ttsemi);

    if (err_tt > 0.0) {
      h_pullDist_tt->Fill((hs.Ntt - val_tt) / err_tt);
    }
    if (err_ttsemi > 0.0) {
      h_pullDist_ttsemi->Fill((hs.Nttsemi - val_ttsemi) / err_ttsemi);
    }

    if ((it + 1) % 100 == 0) {
      std::cout << "Toy " << (it + 1) << "/" << opt.Ntoys_pulls << "\n";
    }

    delete fitres;
    delete toyData;
  }

  std::cout << "Number of problematic fits = " << nFitFail << " / " << opt.Ntoys_pulls << "\n";

  {
    TCanvas* c1 = new TCanvas(TString::Format("c_pullDist_tt_ma%d", mass), "TTbar pull", 800, 600);
    h_pullDist_tt->Fit("gaus", "Q");
    h_pullDist_tt->Draw();
    c1->SaveAs(TString::Format("fits_and_limit_plots/step4_pull_distribution_ttbar_ma%d.pdf", mass));
    delete c1;
  }

  {
    TCanvas* c2 = new TCanvas(TString::Format("c_pullDist_ttsemi_ma%d", mass), "TTsemi pull", 800, 600);
    h_pullDist_ttsemi->Fit("gaus", "Q");
    h_pullDist_ttsemi->Draw();
    c2->SaveAs(TString::Format("fits_and_limit_plots/step4_pull_distribution_ttsemi_ma%d.pdf", mass));
    delete c2;
  }

  {
    TCanvas* c3 = new TCanvas(TString::Format("c_fitval_tt_ma%d", mass), "Fitted TTbar yield", 800, 600);
    h_fitval_tt->Draw();
    c3->SaveAs(TString::Format("fits_and_limit_plots/step4_fitval_ttbar_ma%d.pdf", mass));
    delete c3;
  }

  {
    TCanvas* c4 = new TCanvas(TString::Format("c_fitval_ttsemi_ma%d", mass), "Fitted TTsemi yield", 800, 600);
    h_fitval_ttsemi->Draw();
    c4->SaveAs(TString::Format("fits_and_limit_plots/step4_fitval_ttsemi_ma%d.pdf", mass));
    delete c4;
  }

  {
    TCanvas* c5 = new TCanvas(TString::Format("c_err_tt_ma%d", mass), "Error TTbar yield", 800, 600);
    h_err_tt->Draw();
    c5->SaveAs(TString::Format("fits_and_limit_plots/step4_fiterr_ttbar_ma%d.pdf", mass));
    delete c5;
  }

  {
    TCanvas* c6 = new TCanvas(TString::Format("c_err_ttsemi_ma%d", mass), "Error TTsemi yield", 800, 600);
    h_err_ttsemi->Draw();
    c6->SaveAs(TString::Format("fits_and_limit_plots/step4_fiterr_ttsemi_ma%d.pdf", mass));
    delete c6;
  }

  delete h_pullDist_tt;
  delete h_pullDist_ttsemi;
  delete h_fitval_tt;
  delete h_fitval_ttsemi;
  delete h_err_tt;
  delete h_err_ttsemi;
}

// ============================================================
// 5) Build workspace and proposal for limits
// ============================================================
void buildWorkspaceAndProposal(AnalysisModel& fm, RooDataHist& data_asimov_B) {
  RooFitResult* fit_asimov = fitSBModel(fm, data_asimov_B);

  fm.wr = new RooWorkspace("wr");
  fm.wr->import(*fm.total_model_1);

  fm.mc = new ModelConfig("ModelConfig", fm.wr);
  fm.mc->SetPdf(*fm.wr->pdf("total_model_1"));
  fm.mc->SetParametersOfInterest(RooArgSet(*fm.Nexp_sig));

  RooArgSet nuis(*fm.Nexp_tt, *fm.Nexp_ttsemi);
  fm.mc->SetNuisanceParameters(nuis);
  fm.wr->import(*fm.mc);

  TMatrixDSym cov = fit_asimov->covarianceMatrix();
  TMatrixDSym covSmall = scaleCovMatrix(cov, 0.20);

  ProposalHelper ph;
  ph.SetVariables((RooArgSet&)fit_asimov->floatParsFinal());
  ph.SetCovMatrix(covSmall);
  ph.SetUpdateProposalParameters(kTRUE);
  ph.SetCacheSize(100);
  fm.pf = ph.GetProposalFunction();

  delete fit_asimov;
}

// ============================================================
// 5) Posterior plot from Asimov
// ============================================================
void makePosteriorPlot(int mass, AnalysisModel& fm, RooDataHist& data_asimov_B) {
  MCMCCalculator post_mcmc(data_asimov_B, *fm.mc);
  post_mcmc.SetProposalFunction(*fm.pf);
  post_mcmc.SetConfidenceLevel(opt.CL_95);
  post_mcmc.SetNumIters(opt.numIters);
  post_mcmc.SetNumBurnInSteps(opt.numBurnInSteps);
  post_mcmc.SetLeftSideTailFraction(0.0);

  MCMCInterval* post_interval = post_mcmc.GetInterval();

  TCanvas* c_post = new TCanvas(TString::Format("c_post_ma%d", mass), "Posterior", 800, 700);
  MCMCIntervalPlot postPlot(*post_interval);
  postPlot.SetLineColor(kBlue + 1);
  postPlot.Draw();

  TLatex latp;
  latp.SetNDC();
  latp.SetTextSize(0.035);
  latp.DrawLatex(0.16, 0.92, TString::Format("Posterior for m_{a}=%d GeV", mass));
  latp.DrawLatex(0.16, 0.87, TString::Format("95%% upper limit = %.4g", post_interval->UpperLimit(*fm.Nexp_sig)));

  c_post->SaveAs(TString::Format("fits_and_limit_plots/step5_posterior_ma%d.pdf", mass));

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

// ============================================================
// 5) Final expected-limit plot
// ============================================================
void makeFinalLimitPlot(const LimitSummary& summary) {
  const int n = (int)summary.mass_points_d.size();
  if (n <= 0) return;

  std::vector<double> x(n), y(n);
  std::vector<double> exl(n, 0.0), exh(n, 0.0);
  std::vector<double> eyl1(n), eyh1(n), eyl2(n), eyh2(n);

  for (int i = 0; i < n; ++i) {
    x[i] = summary.mass_points_d[i];
    y[i] = summary.expected_br_limits[i];

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

  g2->SetFillColor(kYellow);
  g2->SetLineColor(kYellow);
  g2->SetTitle(";m_{a} [GeV];Expected Bayesian 95% upper limit");

  g1->SetFillColor(kGreen + 1);
  g1->SetLineColor(kGreen + 1);

  gmed->SetLineColor(kBlack);
  gmed->SetLineWidth(2);
  gmed->SetMarkerStyle(20);
  gmed->SetMarkerSize(1.1);

  g2->Draw("A3");
  g1->Draw("3 SAME");
  gmed->Draw("LP SAME");

  TLegend* leglim = new TLegend(0.58, 0.68, 0.88, 0.88);
  leglim->SetBorderSize(0);
  leglim->SetFillStyle(0);
  leglim->AddEntry(gmed, "Median expected", "lp");
  leglim->AddEntry(g1, "#pm1#sigma", "f");
  leglim->AddEntry(g2, "#pm2#sigma", "f");
  leglim->Draw();

  TLatex latlim;
  latlim.SetNDC();
  latlim.SetTextSize(0.035);
  latlim.DrawLatex(0.15, 0.92, "Expected Bayesian 95% upper limits");

  c_lim->SaveAs("fits_and_limit_plots/step5_expected_limits_bayesian.pdf");

  delete leglim;
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
  gSystem->mkdir("fits_and_limit_plots", kTRUE);

  LimitSummary summary;

  for (size_t im = 0; im < opt.mass_points.size(); ++im) {
    const int mass = opt.mass_points[im];

    std::cout << "\n\n############################################################\n";
    std::cout << " Processing mass point m_a = " << mass << " GeV\n";
    std::cout << "############################################################\n";

    InputFiles files;
    Histograms hs;
    AnalysisModel fm;

    // ========================================================
    // 1) Construct models
    // ========================================================
    std::cout << "\n[STEP 1] Construct input histograms and statistical model\n";
    if (!buildInputHistograms(mass, files, hs)) {
      std::cerr << "[FATAL] Could not build histograms for mass " << mass << "\n";
      files.close();
      continue;
    }
    constructModel(hs, fm);

    // ========================================================
    // 2) Generate pseudo-data
    // ========================================================
    std::cout << "[STEP 2] Generate pseudo-data\n";
    RooDataHist* data_asimov_B = generateAsimovB(fm);
    RooDataHist* data_toy_B    = generateToyB(fm, "toy_B");

    // ========================================================
    // 3) Perform one fit and inspect outputs
    // ========================================================
    std::cout << "[STEP 3] One fit and detailed output\n";
    RooFitResult* fit_one = fitBOnlyModel(fm, *data_toy_B);
    printOneFitSummary(mass, hs, fm, fit_one);
    plotOneToyFit(mass, fm, *data_toy_B);

    delete fit_one;
    delete data_asimov_B;
    delete data_toy_B;

    delete hs.h_ttsemi_r;

    fm.cleanup();
    files.close();

    return;

    // ========================================================
    // 4) Toy MC study for pulls
    // ========================================================
    std::cout << "[STEP 4] Toy MC pull study\n";
    runToyPullStudy(mass, hs, fm);

    // ========================================================
    // 5) Limits
    // ========================================================
    std::cout << "[STEP 5] Bayesian limits\n";
    buildWorkspaceAndProposal(fm, *data_asimov_B);
    makePosteriorPlot(mass, fm, *data_asimov_B);
    runExpectedLimitsForMass(mass, hs, fm, summary);

    delete fit_one;
    delete data_asimov_B;
    delete data_toy_B;

    delete hs.h_ttsemi_r;

    fm.cleanup();
    files.close();
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

  makeFinalLimitPlot(summary);
}
