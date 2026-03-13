// fit_and_limits.cc
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

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

#include "RooStats/ModelConfig.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/MCMCInterval.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/ProposalHelper.h"
#include "RooStats/ProposalFunction.h"

using namespace RooFit;
using namespace RooStats;

// ------------------------------
// MCMC options
// ------------------------------
struct BayesianMCMCOptions {
  double CL_95 = 0.95;
  int numIters = 10000;
  int numBurnInSteps = 0; // safest for your ROOT version
};
BayesianMCMCOptions optMCMC;

// ------------------------------
// Helpers
// ------------------------------
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

static TMatrixDSym scaleCovMatrix(const TMatrixDSym& cov, double scale)
{
  TMatrixDSym c(cov);
  c *= (scale * scale);
  return c;
}

// ------------------------------
// Main
// ------------------------------
void fit_and_limits()
{
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(1111);

  // create output folder
  gSystem->mkdir("fits_and_limit_plots", kTRUE);

  // mass points
  std::vector<int> mass_points = {12,15,20,25,30};

  // toys
  const int Ntoys = 1001;

  // DO NOT TOUCH BINS
  const int b = 4;
  double var_bins[b+1];
  var_bins[0] = -1;
  var_bins[1] = -0.4;
  var_bins[2] = -0.08;
  var_bins[3] = 0.52;
  var_bins[4] = 1;

  // results
  std::vector<double> expected_br_limits;
  std::vector<double> sigma_1_br_up, sigma_1_br_down;
  std::vector<double> sigma_2_br_up, sigma_2_br_down;

  // keep masses as doubles for final limit plot
  std::vector<double> mass_points_d;

  for (size_t im = 0; im < mass_points.size(); ++im) {
    const int mass = mass_points[im];
    mass_points_d.push_back((double)mass);

    // =========================================================
    // 1) Input files
    // =========================================================
    TFile *fsig   = TFile::Open(signalFileForMass(mass), "READ");
    TFile *ftt    = TFile::Open("output_ttbar.root", "READ");
    TFile *ftt2q  = TFile::Open("output_TTtoLNu2Q.root", "READ");
    TFile *fTTHbb = TFile::Open("output_TTH_Hbb.root", "READ");
    TFile *fTTW   = TFile::Open("output_TTW.root", "READ");
    TFile *fTTZ   = TFile::Open("output_TTZ.root", "READ");

    // =========================================================
    // 2) Histograms
    // =========================================================
    const char* HNAME = "h_BDT";

    TH1* h_sig   = getHistOrDie(fsig,   HNAME, fsig   ? fsig->GetName()   : "sig");
    TH1* h_tt    = getHistOrDie(ftt,    HNAME, ftt    ? ftt->GetName()    : "ttbar");
    TH1* h_tt2q  = getHistOrDie(ftt2q,  HNAME, ftt2q  ? ftt2q->GetName()  : "tt2q");
    TH1* h_tthbb = getHistOrDie(fTTHbb, HNAME, fTTHbb ? fTTHbb->GetName() : "tthbb");
    TH1* h_ttw   = getHistOrDie(fTTW,   HNAME, fTTW   ? fTTW->GetName()   : "ttw");
    TH1* h_ttz   = getHistOrDie(fTTZ,   HNAME, fTTZ   ? fTTZ->GetName()   : "ttz");

    if (!h_sig || !h_tt || !h_tt2q || !h_tthbb || !h_ttw || !h_ttz) {
      std::cerr << "[FATAL] Missing input histogram(s). Fix file names / hist names.\n";
      return;
    }

    // =========================================================
    // 3) Rebin
    // =========================================================
    TH1* h_sig_r   = h_sig  ->Rebin(b, "h_sig_r",   var_bins);
    TH1* h_tt_r    = h_tt   ->Rebin(b, "h_tt_r",    var_bins);
    TH1* h_tt2q_r  = h_tt2q ->Rebin(b, "h_tt2q_r",  var_bins);
    TH1* h_tthbb_r = h_tthbb->Rebin(b, "h_tthbb_r", var_bins);
    TH1* h_ttw_r   = h_ttw  ->Rebin(b, "h_ttw_r",   var_bins);
    TH1* h_ttz_r   = h_ttz  ->Rebin(b, "h_ttz_r",   var_bins);

    // =========================================================
    // 3b) Combined TT semileptonic histogram:
    //     TTtoLNu2Q + TTH(H->bb) + TTW + TTZ
    // =========================================================
    TH1* h_ttsemi_r = (TH1*)h_tt2q_r->Clone("h_ttsemi_r");
    h_ttsemi_r->Reset();
    h_ttsemi_r->Add(h_tt2q_r);
    h_ttsemi_r->Add(h_tthbb_r);
    h_ttsemi_r->Add(h_ttw_r);
    h_ttsemi_r->Add(h_ttz_r);

    // =========================================================
    // 4) Observable
    // =========================================================
    RooRealVar output_BDT("output_BDT", "BDT score", -1, 1);
    RooBinning customBinning(b, var_bins);
    output_BDT.setBinning(customBinning, "customBinning");

    // =========================================================
    // 5) RooDataHists + RooHistPdf
    // =========================================================
    RooDataHist dh_sig    ("sig",    "sig",    output_BDT, h_sig_r);
    RooDataHist dh_tt     ("tt",     "tt",     output_BDT, h_tt_r);
    RooDataHist dh_ttsemi ("ttsemi", "ttsemi", output_BDT, h_ttsemi_r);

    RooHistPdf pdf_sig    ("sig_pdf",    "sig_pdf",    output_BDT, dh_sig);
    RooHistPdf pdf_tt     ("tt_pdf",     "tt_pdf",     output_BDT, dh_tt);
    RooHistPdf pdf_ttsemi ("ttsemi_pdf", "ttsemi_pdf", output_BDT, dh_ttsemi);

    // =========================================================
    // 6) Yields
    // =========================================================
    const double Nsig    = h_sig_r->Integral();
    const double Ntt     = h_tt_r->Integral();
    const double Nttsemi = h_ttsemi_r->Integral();

    const double Nbkg = Ntt + Nttsemi;
    const double denominator = (Nsig > 0 ? Nsig : 1.0);

    // =========================================================
    // 7) Parameters + constraints
    // =========================================================
    RooRealVar Nexp_sig   ("Nexp_sig",   "Expected signal events",
                           Nsig, 0.0, 1000.0 * std::max(1.0, Nsig));
    RooRealVar Nexp_tt    ("Nexp_tt",    "Expected TT dileptonic",
                           Ntt, 0.0, 10.0 * std::max(1.0, Ntt));
    RooRealVar Nexp_ttsemi("Nexp_ttsemi","Expected TT semileptonic",
                           Nttsemi, 0.0, 10.0 * std::max(1.0, Nttsemi));

    const double rel_tt     = 0.15;
    const double rel_ttsemi = 0.30;

    RooConstVar tt_nom     ("tt_nom",     "tt_nom",     Ntt);
    RooConstVar tt_sig     ("tt_sig",     "tt_sig",     safeSigma(Ntt, rel_tt));
    RooConstVar ttsemi_nom ("ttsemi_nom", "ttsemi_nom", Nttsemi);
    RooConstVar ttsemi_sig ("ttsemi_sig", "ttsemi_sig", safeSigma(Nttsemi, rel_ttsemi));

    RooGaussian c_tt    ("c_tt",    "constraint tt",     Nexp_tt,     tt_nom,     tt_sig);
    RooGaussian c_ttsemi("c_ttsemi","constraint ttsemi", Nexp_ttsemi, ttsemi_nom, ttsemi_sig);

    // =========================================================
    // 8) Models
    // =========================================================
    RooAddPdf model_0("model_0", "Background-only (extended)",
      RooArgList(pdf_tt, pdf_ttsemi),
      RooArgList(Nexp_tt, Nexp_ttsemi)
    );

    RooAddPdf model_1("model_1", "Signal+Background (extended)",
      RooArgList(pdf_sig, pdf_tt, pdf_ttsemi),
      RooArgList(Nexp_sig, Nexp_tt, Nexp_ttsemi)
    );

    RooProdPdf total_model_0("total_model_0", "b-only with constraints",
      RooArgList(model_0, c_tt, c_ttsemi)
    );

    RooProdPdf total_model_1("total_model_1", "s+b with constraints",
      RooArgList(model_1, c_tt, c_ttsemi)
    );

    // =========================================================
    // 9) Asimov dataset
    // =========================================================
    RooDataHist* data_asimov_B = total_model_0.generateBinned(
      output_BDT,
      RooFit::ExpectedData(true),
      RooFit::Name("asimov_B")
    );

    RooFitResult* fit_asimov = total_model_1.fitTo(
      *data_asimov_B,
      Save(),
      Extended(kTRUE),
      RooFit::MaxCalls(20000),
      Strategy(2),
      RooFit::Optimize(kTRUE),
      RooFit::SumW2Error(kTRUE)
    );

    // =========================================================
    // 10) One pseudodata fit
    // =========================================================
    RooDataHist* data_toy_B = total_model_0.generateBinned(
      output_BDT,
      RooFit::Extended(true),
      RooFit::Name("toy_B")
    );

    RooFitResult* fit_toy = total_model_1.fitTo(
      *data_toy_B,
      Save(),
      Extended(kTRUE),
      RooFit::MaxCalls(20000),
      Strategy(2),
      RooFit::Optimize(kTRUE),
      RooFit::SumW2Error(kTRUE)
    );

    // =========================================================
    // 10b) Plot 1 pseudodata + fit + normalized signal overlay
    //      backgrounds shown as TTbar and TT semileptonic
    // =========================================================
    {
      TCanvas* c_bdt = new TCanvas("c_bdt","Pseudodata + fit",900,800);
      RooPlot* fr = output_BDT.frame();

      data_toy_B->plotOn(fr, Name("toydata"));
      total_model_1.plotOn(fr, Name("fullfit"), LineColor(kBlue));

      total_model_1.plotOn(fr, Components("tt_pdf"),
                           LineColor(kRed+1), LineStyle(kSolid), Name("ttcomp"));
      total_model_1.plotOn(fr, Components("ttsemi_pdf"),
                           LineColor(kMagenta+1), LineStyle(kSolid), Name("ttsemicomp"));

      // signal superimposed with fitted normalization
      total_model_1.plotOn(fr, Components("sig_pdf"),
                           LineColor(kBlack), LineStyle(kDashed), Name("sigonly"));

      fr->SetTitle(TString::Format("Pseudodata + fit, m_{a}=%d GeV", mass));
      fr->GetXaxis()->SetTitle("BDT score");
      fr->GetYaxis()->SetTitle("Events");
      fr->Draw();

      TLegend* leg = new TLegend(0.56,0.60,0.88,0.88);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->AddEntry(fr->findObject("toydata"),     "Pseudodata", "pe");
      leg->AddEntry(fr->findObject("fullfit"),     "S+B fit", "l");
      leg->AddEntry(fr->findObject("sigonly"),     "Signal", "l");
      leg->AddEntry(fr->findObject("ttcomp"),      "TTbar", "l");
      leg->AddEntry(fr->findObject("ttsemicomp"),  "TT semileptonic", "l");
      leg->Draw();

      TLatex lat;
      lat.SetNDC();
      lat.SetTextSize(0.035);
      lat.DrawLatex(0.15, 0.92, TString::Format("m_{a}=%d GeV", mass));

      c_bdt->SaveAs(TString::Format("fits_and_limit_plots/pseudodata_fit_bdt_ma%d.pdf", mass));
      delete c_bdt;
    }

    // =========================================================
    // 10c) Separate pull plot:
    //      Pull = (Nexp - Nfit) / dNfit
    //      shown for TTbar and TT semileptonic
    // =========================================================
    {
      const int npulls = 2;

      const char* labels[npulls] = {
        "TTbar",
        "TT semileptonic"
      };

      double nom[npulls] = {
        Ntt,
        Nttsemi
      };

      double fitv[npulls] = {
        Nexp_tt.getVal(),
        Nexp_ttsemi.getVal()
      };

      double fite[npulls] = {
        Nexp_tt.getError(),
        Nexp_ttsemi.getError()
      };

      double pulls[npulls];
      double x[npulls];
      double ex[npulls];
      double ey[npulls];

      for (int i = 0; i < npulls; ++i) {
        x[i]  = i + 1;
        ex[i] = 0.0;
        ey[i] = 0.0;
        pulls[i] = (fite[i] > 0.0) ? (nom[i] - fitv[i]) / fite[i] : 0.0;
      }

      TCanvas* c_pull = new TCanvas("c_pull","Parameter pulls",800,600);
      TH1F* hframe = new TH1F("hframe",";Process;Pull = (N_{exp} - N_{fit}) / #delta N_{fit}", npulls, 0.5, npulls + 0.5);

      for (int i = 1; i <= npulls; ++i) {
        hframe->GetXaxis()->SetBinLabel(i, labels[i-1]);
      }

      hframe->SetMinimum(-5.0);
      hframe->SetMaximum(5.0);
      hframe->Draw();

      TGraphErrors* grPull = new TGraphErrors(npulls, x, pulls, ex, ey);
      grPull->SetMarkerStyle(20);
      grPull->SetMarkerSize(1.2);
      grPull->Draw("P SAME");

      TLine* l0  = new TLine(0.5,  0.0, npulls + 0.5,  0.0);
      TLine* lp1 = new TLine(0.5,  1.0, npulls + 0.5,  1.0);
      TLine* lm1 = new TLine(0.5, -1.0, npulls + 0.5, -1.0);
      TLine* lp2 = new TLine(0.5,  2.0, npulls + 0.5,  2.0);
      TLine* lm2 = new TLine(0.5, -2.0, npulls + 0.5, -2.0);

      l0->SetLineColor(kBlack);
      lp1->SetLineColor(kBlue);  lm1->SetLineColor(kBlue);
      lp2->SetLineColor(kRed);   lm2->SetLineColor(kRed);

      lp1->SetLineStyle(2); lm1->SetLineStyle(2);
      lp2->SetLineStyle(2); lm2->SetLineStyle(2);

      l0->Draw("SAME");
      lp1->Draw("SAME");
      lm1->Draw("SAME");
      lp2->Draw("SAME");
      lm2->Draw("SAME");

      TLatex lat2;
      lat2.SetNDC();
      lat2.SetTextSize(0.035);
      lat2.DrawLatex(0.15, 0.92, TString::Format("Post-fit pulls, m_{a}=%d GeV", mass));

      c_pull->SaveAs(TString::Format("fits_and_limit_plots/pulls_ma%d.pdf", mass));

      delete grPull;
      delete hframe;
      delete l0;
      delete lp1;
      delete lm1;
      delete lp2;
      delete lm2;
      delete c_pull;
    }

    // =========================================================
    // 10d) NEW: bin-by-bin pull plot like your colleague had
    // =========================================================
    {
      RooPlot* frame_for_pull = output_BDT.frame();
      data_toy_B->plotOn(frame_for_pull, Name("toydata_pull"));
      total_model_1.plotOn(frame_for_pull, Name("fullfit_pull"));

      RooHist* hpull = frame_for_pull->pullHist("toydata_pull", "fullfit_pull");

      TCanvas* c_binpull = new TCanvas("c_binpull","Bin pulls",800,600);
      RooPlot* frame_pull = output_BDT.frame();
      frame_pull->addPlotable(hpull, "P");

      frame_pull->SetTitle(TString::Format("Bin-by-bin pulls, m_{a}=%d GeV", mass));
      frame_pull->GetXaxis()->SetTitle("BDT score");
      frame_pull->GetYaxis()->SetTitle("Pull");
      frame_pull->GetYaxis()->SetTitleOffset(1.2);
      frame_pull->SetMinimum(-5.0);
      frame_pull->SetMaximum(5.0);
      frame_pull->Draw();

      TLine* l0b  = new TLine(-1.0,  0.0, 1.0,  0.0);
      TLine* lp1b = new TLine(-1.0,  1.0, 1.0,  1.0);
      TLine* lm1b = new TLine(-1.0, -1.0, 1.0, -1.0);
      TLine* lp2b = new TLine(-1.0,  2.0, 1.0,  2.0);
      TLine* lm2b = new TLine(-1.0, -2.0, 1.0, -2.0);

      l0b->SetLineColor(kBlack);
      lp1b->SetLineColor(kBlue);  lm1b->SetLineColor(kBlue);
      lp2b->SetLineColor(kRed);   lm2b->SetLineColor(kRed);

      lp1b->SetLineStyle(2); lm1b->SetLineStyle(2);
      lp2b->SetLineStyle(2); lm2b->SetLineStyle(2);

      l0b->Draw("SAME");
      lp1b->Draw("SAME");
      lm1b->Draw("SAME");
      lp2b->Draw("SAME");
      lm2b->Draw("SAME");

      TLatex latb;
      latb.SetNDC();
      latb.SetTextSize(0.035);
      latb.DrawLatex(0.15, 0.92, TString::Format("Bin pulls, m_{a}=%d GeV", mass));

      c_binpull->SaveAs(TString::Format("fits_and_limit_plots/bin_pulls_ma%d.pdf", mass));

      delete l0b;
      delete lp1b;
      delete lm1b;
      delete lp2b;
      delete lm2b;
      delete frame_for_pull;
      delete frame_pull;
      delete c_binpull;
    }

    // =========================================================
    // 11) Workspace + ModelConfig
    // =========================================================
    RooWorkspace wr("wr");
    wr.import(total_model_1);
    ModelConfig mc("ModelConfig", &wr);
    mc.SetPdf(*wr.pdf("total_model_1"));
    mc.SetParametersOfInterest(RooArgSet(Nexp_sig));

    RooArgSet nuis(Nexp_tt, Nexp_ttsemi);
    mc.SetNuisanceParameters(nuis);
    wr.import(mc);

    // =========================================================
    // 12) Proposal
    // =========================================================
    TMatrixDSym cov = fit_asimov->covarianceMatrix();
    TMatrixDSym covSmall = scaleCovMatrix(cov, 0.20);

    ProposalHelper ph;
    ph.SetVariables((RooArgSet&)fit_asimov->floatParsFinal());
    ph.SetCovMatrix(covSmall);
    ph.SetUpdateProposalParameters(kTRUE);
    ph.SetCacheSize(100);
    ProposalFunction* pf = ph.GetProposalFunction();

    // =========================================================
    // A) Posterior plot for this mass
    // =========================================================
    {
      MCMCCalculator post_mcmc(*data_asimov_B, mc);
      post_mcmc.SetProposalFunction(*pf);
      post_mcmc.SetConfidenceLevel(optMCMC.CL_95);
      post_mcmc.SetNumIters(optMCMC.numIters);
      post_mcmc.SetNumBurnInSteps(optMCMC.numBurnInSteps);
      post_mcmc.SetLeftSideTailFraction(0.0);

      MCMCInterval* post_interval = post_mcmc.GetInterval();

      TCanvas* c_post = new TCanvas("c_post","Posterior",800,700);
      MCMCIntervalPlot postPlot(*post_interval);
      postPlot.SetLineColor(kBlue+1);
      postPlot.Draw();

      TLatex latp;
      latp.SetNDC();
      latp.SetTextSize(0.035);
      latp.DrawLatex(0.16, 0.92, TString::Format("Posterior for m_{a}=%d GeV", mass));
      latp.DrawLatex(0.16, 0.87, TString::Format("95%% upper limit = %.4g", post_interval->UpperLimit(Nexp_sig)));

      c_post->SaveAs(TString::Format("fits_and_limit_plots/posterior_ma%d.pdf", mass));

      delete c_post;
      delete post_interval;
    }

    // =========================================================
    // 13) Expected limits + pull distributions over toys
    // =========================================================
    std::vector<double> toy_br_limits;
    toy_br_limits.reserve(Ntoys);

    TH1F* h_pullDist_tt = new TH1F(
      TString::Format("h_pullDist_tt_ma%d",mass),
      TString::Format("TTbar pull distribution, m_{a}=%d;Pull;Toys",mass),
      60,-5,5
    );

    TH1F* h_pullDist_ttsemi = new TH1F(
      TString::Format("h_pullDist_ttsemi_ma%d",mass),
      TString::Format("TT semileptonic pull distribution, m_{a}=%d;Pull;Toys",mass),
      60,-5,5
    );

    std::cout << "\n====================================================\n";
    std::cout << " Expected limits for m_a = " << mass << " GeV\n";
    std::cout << " denominator (Nsig) = " << denominator << "\n";
    std::cout << " Nbkg = " << Nbkg << "\n";
    std::cout << "====================================================\n";

    for (int it = 0; it < Ntoys; ++it) {
      RooDataHist* toyData = total_model_0.generateBinned(
        output_BDT,
        RooFit::Extended(true)
      );

      // fit MODEL 0 for pull distributions
      RooFitResult* fit_pull_toy = total_model_0.fitTo(
        *toyData,
        Save(),
        Extended(kTRUE),
        RooFit::MaxCalls(10000),
        Strategy(1),
        RooFit::Optimize(kTRUE),
        RooFit::SumW2Error(kTRUE)
      );

      if (Nexp_tt.getError() > 0.0) {
        h_pullDist_tt->Fill((Ntt - Nexp_tt.getVal()) / Nexp_tt.getError());
      }

      if (Nexp_ttsemi.getError() > 0.0) {
        h_pullDist_ttsemi->Fill((Nttsemi - Nexp_ttsemi.getVal()) / Nexp_ttsemi.getError());
      }

      delete fit_pull_toy;

      MCMCCalculator mcmc(*toyData, mc);
      mcmc.SetProposalFunction(*pf);
      mcmc.SetConfidenceLevel(optMCMC.CL_95);
      mcmc.SetNumIters(optMCMC.numIters);
      mcmc.SetNumBurnInSteps(optMCMC.numBurnInSteps);
      mcmc.SetLeftSideTailFraction(0.0);

      MCMCInterval* interval = mcmc.GetInterval();
      const double Nsig_up = interval->UpperLimit(Nexp_sig);
      const double BR_up   = Nsig_up / denominator;

      toy_br_limits.push_back(BR_up);

      if ((it+1) % 50 == 0) {
        std::cout << "Toy " << (it+1) << "/" << Ntoys
                  << "  Nsig^95 = " << Nsig_up
                  << "  BR^95 = " << BR_up << "\n";
      }

      delete interval;
      delete toyData;
    }

    // =========================================================
    // B) Save pull distribution histograms for this mass
    //    each histogram separately, with stat + fit boxes
    // =========================================================
    {
      TCanvas* c1 = new TCanvas(TString::Format("c_pullDist_tt_ma%d",mass),
                                "TTbar pull distribution",800,600);
      h_pullDist_tt->Fit("gaus","Q");
      h_pullDist_tt->Draw();
      c1->Update();
      c1->SaveAs(TString::Format("fits_and_limit_plots/pull_distribution_ttbar_ma%d.pdf", mass));
      delete c1;
    }

    {
      TCanvas* c2 = new TCanvas(TString::Format("c_pullDist_ttsemi_ma%d",mass),
                                "TT semileptonic pull distribution",800,600);
      h_pullDist_ttsemi->Fit("gaus","Q");
      h_pullDist_ttsemi->Draw();
      c2->Update();
      c2->SaveAs(TString::Format("fits_and_limit_plots/pull_distribution_ttsemi_ma%d.pdf", mass));
      delete c2;
    }

    std::sort(toy_br_limits.begin(), toy_br_limits.end());

    auto q = [&](double p)->double {
      const int idx = std::max(0, std::min(int(std::floor(p * (Ntoys - 1))), Ntoys - 1));
      return toy_br_limits[idx];
    };

    const double median = q(0.50);
    expected_br_limits.push_back(median);
    sigma_1_br_down.push_back(q(0.16));
    sigma_1_br_up  .push_back(q(0.84));
    sigma_2_br_down.push_back(q(0.025));
    sigma_2_br_up  .push_back(q(0.975));

    std::cout << "\n[RESULT] m_a=" << mass
              << "  expected BR^95 median=" << median
              << "  (-1σ=" << q(0.16) << ", +1σ=" << q(0.84) << ")"
              << "  (-2σ=" << q(0.025) << ", +2σ=" << q(0.975) << ")"
              << "\n";

    delete h_pullDist_tt;
    delete h_pullDist_ttsemi;
    delete h_ttsemi_r;

    delete fit_asimov;
    delete fit_toy;
    delete data_asimov_B;
    delete data_toy_B;

    if (fsig)   fsig->Close();
    if (ftt)    ftt->Close();
    if (ftt2q)  ftt2q->Close();
    if (fTTHbb) fTTHbb->Close();
    if (fTTW)   fTTW->Close();
    if (fTTZ)   fTTZ->Close();
  }

  // final vectors
  std::cout << "\nexpected_br_limits {";
  for (size_t i = 0; i < expected_br_limits.size(); ++i) {
    std::cout << expected_br_limits[i] << (i+1<expected_br_limits.size()? ", ":"");
  }
  std::cout << "}\n";

  std::cout << "sigma_1_br_down {";
  for (size_t i = 0; i < sigma_1_br_down.size(); ++i) {
    std::cout << sigma_1_br_down[i] << (i+1<sigma_1_br_down.size()? ", ":"");
  }
  std::cout << "}\n";

  std::cout << "sigma_1_br_up {";
  for (size_t i = 0; i < sigma_1_br_up.size(); ++i) {
    std::cout << sigma_1_br_up[i] << (i+1<sigma_1_br_up.size()? ", ":"");
  }
  std::cout << "}\n";

  std::cout << "sigma_2_br_down {";
  for (size_t i = 0; i < sigma_2_br_down.size(); ++i) {
    std::cout << sigma_2_br_down[i] << (i+1<sigma_2_br_down.size()? ", ":"");
  }
  std::cout << "}\n";

  std::cout << "sigma_2_br_up {";
  for (size_t i = 0; i < sigma_2_br_up.size(); ++i) {
    std::cout << sigma_2_br_up[i] << (i+1<sigma_2_br_up.size()? ", ":"");
  }
  std::cout << "}\n";

  // =========================================================
  // C) Final expected limit plot using all masses
  // =========================================================
  {
    const int n = (int)mass_points_d.size();

    std::vector<double> x(n), y(n);
    std::vector<double> exl(n,0.0), exh(n,0.0);
    std::vector<double> eyl1(n), eyh1(n), eyl2(n), eyh2(n);

    for (int i = 0; i < n; ++i) {
      x[i] = mass_points_d[i];
      y[i] = expected_br_limits[i];

      eyl1[i] = expected_br_limits[i] - sigma_1_br_down[i];
      eyh1[i] = sigma_1_br_up[i]      - expected_br_limits[i];

      eyl2[i] = expected_br_limits[i] - sigma_2_br_down[i];
      eyh2[i] = sigma_2_br_up[i]      - expected_br_limits[i];
    }

    TCanvas* c_lim = new TCanvas("c_lim","Expected Bayesian limits",900,700);

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

    g1->SetFillColor(kGreen+1);
    g1->SetLineColor(kGreen+1);

    gmed->SetLineColor(kBlack);
    gmed->SetLineWidth(2);
    gmed->SetMarkerStyle(20);
    gmed->SetMarkerSize(1.1);

    g2->Draw("A3");
    g1->Draw("3 SAME");
    gmed->Draw("LP SAME");

    TLegend* leglim = new TLegend(0.58,0.68,0.88,0.88);
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

    c_lim->SaveAs("fits_and_limit_plots/expected_limits_bayesian.pdf");

    delete leglim;
    delete gmed;
    delete g1;
    delete g2;
    delete c_lim;
  }
}
