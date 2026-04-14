#include <iostream>
#include <vector>
#include <cmath>

#include <TH1F.h>
#include <TFile.h>
#include <TMath.h>

void binopt()
{
  bool verbose(true);
  bool printout(true);

  // Already weighted in MyClass
  float int_lumi = 1.0;

  double bg_syst_mult = 0.0;
  double bg_syst_add  = 0.0;

  // ==========================================================
  // OPEN FILES (YOUR OUTPUTS)
  // ==========================================================
  TFile *Signal12 = TFile::Open("output_signal_tth12gev.root");
  TFile *Signal15 = TFile::Open("output_signal_tth15gev.root");
  TFile *Signal20 = TFile::Open("output_signal_tth20gev.root");
  TFile *Signal25 = TFile::Open("output_signal_tth25gev.root");
  TFile *Signal30 = TFile::Open("output_signal_tth30gev.root");

  TFile *TTbar   = TFile::Open("output_ttbar.root");
  TFile *TTlnu2q = TFile::Open("output_TTtoLNu2Q.root");
  TFile *TTHHbb  = TFile::Open("output_TTH_Hbb.root");

  // ==========================================================
  // GET HISTOGRAMS
  // ==========================================================
  TH1F *h_Signal =
    (TH1F*)Signal12->Get("h_BDT")->Clone("h_Signal");

  h_Signal->Add((TH1F*)Signal15->Get("h_BDT"));
  h_Signal->Add((TH1F*)Signal20->Get("h_BDT"));
  h_Signal->Add((TH1F*)Signal25->Get("h_BDT"));
  h_Signal->Add((TH1F*)Signal30->Get("h_BDT"));

  h_Signal->Scale(int_lumi);

  TH1F *h_TTbar   = (TH1F*)TTbar->Get("h_BDT");
  TH1F *h_TTlnu2q = (TH1F*)TTlnu2q->Get("h_BDT");
  TH1F *h_TTHHbb  = (TH1F*)TTHHbb->Get("h_BDT");

  h_TTbar->Scale(int_lumi);
  h_TTlnu2q->Scale(int_lumi);
  h_TTHHbb->Scale(int_lumi);

  TH1F *h_Bkg = (TH1F*)h_TTbar->Clone("h_Bkg");
  h_Bkg->Add(h_TTlnu2q);
  h_Bkg->Add(h_TTHHbb);

  std::vector<TH1*> h_BkgProcesses =
    {h_TTbar, h_TTlnu2q, h_TTHHbb};

  // ==========================================================
  // BASIC INFO
  // ==========================================================
  std::cout << "Signal integral after scaling: "
            << h_Signal->Integral() << std::endl;

  std::cout << "Background integral after scaling: "
            << h_Bkg->Integral() << std::endl;

  double totalSignal     = h_Signal->Integral();
  double totalBackground = h_Bkg->Integral();

  // ==========================================================
  // BIN OPTIMIZATION
  // ==========================================================
  int nBins = h_Signal->GetNbinsX();

  double maxSOverSigmaBTotal = 0;
  int bestBin1 = 0;
  int bestBin2 = 0;

  for (int bin1 = nBins; bin1 >= 1; --bin1)
  {
    bool validBin1 = true;
    for (auto& hProcess : h_BkgProcesses)
    {
      if (hProcess->Integral(bin1, nBins) <= 0)
      {
        validBin1 = false;
        break;
      }
    }
    if (!validBin1) continue;

    double s1 = h_Signal->Integral(bin1, nBins);
    double b1 = h_Bkg->Integral(bin1, nBins);

    double bgStatError1;
    h_Bkg->IntegralAndError(bin1, nBins, bgStatError1);

    double sigmaB1 =
      sqrt(pow(bgStatError1,2)+pow(bg_syst_mult*b1,2));

    if (sigmaB1 == 0 || b1 <= 0) continue;

    double sOverSigmaB1 = s1 / sigmaB1;

    for (int bin2 = bin1 - 1; bin2 >= 1; --bin2)
    {
      bool validBin2 = true;
      for (auto& hProcess : h_BkgProcesses)
      {
        if (hProcess->Integral(bin2, bin1 - 1) <= 0.0)
        {
          validBin2 = false;
          break;
        }
      }
      if (!validBin2) continue;

      double s2 = h_Signal->Integral(bin2, bin1 - 1);
      double b2 = h_Bkg->Integral(bin2, bin1 - 1);

      double bgStatError2;
      h_Bkg->IntegralAndError(bin2, bin1 - 1, bgStatError2);

      double sigmaB2 =
        sqrt(pow(bgStatError2,2)+pow(bg_syst_mult*b2,2));

      if (sigmaB2 == 0 || b2 <= 0) continue;

      double sOverSigmaB2 = s2 / sigmaB2;

      double sOverSigmaBTotal =
        sqrt(pow(sOverSigmaB1,2) + pow(sOverSigmaB2,2));

      if (sOverSigmaBTotal > maxSOverSigmaBTotal)
      {
        maxSOverSigmaBTotal = sOverSigmaBTotal;
        bestBin1 = bin1;
        bestBin2 = bin2;
      }
    }
  }

  // ==========================================================
  // ===== RESTORED VERBOSE PRINTING BLOCK (UNCHANGED STYLE)
  // ==========================================================
  std::cout << "Best bin pair: "
            << bestBin1 << " and " << bestBin2
            << " with S/B total: "
            << maxSOverSigmaBTotal << std::endl;

  if (bestBin1 > 0 && bestBin2 > 0)
  {
    std::cout <<"tot sig: "<< h_Signal->Integral()<<std::endl;
    std::cout <<"tot bkg: "<< h_Bkg->Integral()<<std::endl;

    double s1 = h_Signal->Integral(bestBin1, nBins);
    double b1 = h_Bkg->Integral(bestBin1, nBins);

    double bgStatError1;
    h_Bkg->IntegralAndError(bestBin1, nBins, bgStatError1);

    double sigmaB1 =
      sqrt(pow(bgStatError1,2)+pow(bg_syst_mult*b1,2));

    double sOverSigmaB1 = s1 / sigmaB1;

    double s2 = h_Signal->Integral(bestBin2, bestBin1 - 1);
    double b2 = h_Bkg->Integral(bestBin2, bestBin1 - 1);

    double bgStatError2;
    h_Bkg->IntegralAndError(bestBin2, bestBin1 - 1, bgStatError2);

    double sigmaB2 =
      sqrt(pow(bgStatError2,2)+pow(bg_syst_mult*b2,2));

    double sOverSigmaB2 = s2 / sigmaB2;

    std::cout << "Optimal bins found:\n";

    std::cout << "Bin 1: ["
              << h_Signal->GetXaxis()->GetBinLowEdge(bestBin1)
              << ", "
              << h_Signal->GetXaxis()->GetBinUpEdge(nBins)
              << "]\n";

    std::cout << "S: " << s1
              << ", Bkg: " << b1
              << ", S/sigmaB: " << sOverSigmaB1
              << ", S/sqrt(B+S): "
              << s1 / sqrt(b1 + s1) << "\n";

    std::cout << "Bin 2: ["
              << h_Signal->GetXaxis()->GetBinLowEdge(bestBin2)
              << ", "
              << h_Signal->GetXaxis()->GetBinUpEdge(bestBin1 - 1)
              << "]\n";

    std::cout << "S: " << s2
              << ", Bkg: " << b2
              << ", S/sigmaB: " << sOverSigmaB2
              << ", S/sqrt(B+S): "
              << s2 / sqrt(b2 + s2) << "\n";

    for (size_t i = 0; i < h_BkgProcesses.size(); ++i)
    {
      double processB1 =
        h_BkgProcesses[i]->Integral(bestBin1, nBins);
      double processB2 =
        h_BkgProcesses[i]->Integral(bestBin2, bestBin1 - 1);

      std::cout << "Background process "
                << i + 1
                << " in Bin 1: " << processB1
                << ", in Bin 2: " << processB2 << "\n";
    }

    std::cout << "Total S/sigmaB: "
              << maxSOverSigmaBTotal << "\n";

    for (size_t i = 0; i < h_BkgProcesses.size(); ++i)
    {
      double process = h_BkgProcesses[i]->Integral();
      std::cout << "Background process "
                << i + 1
                << " ratio: "
                << process/h_Bkg->Integral()
                << std::endl;
    }
  }
  else
  {
    std::cout << "No valid bin combination found." << std::endl;
  }
}
