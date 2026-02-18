// File: TMVA_Classification_MyAnalysis_AllSignals.C
//
// Train: Signal = {tth12, tth15, tth20, tth25, tth30}
//        Background = {TTtoLNu2Q, TTto2L2Nu, TTHHtobb}
//
// IMPORTANT: This expects your REDUCED output files produced by MyClass,
// i.e. files that contain TTree "AnalysisTree" with branch "weight".
//
// Usage:
//   root -l
//   .L TMVA_Classification_MyAnalysis_AllSignals.C
//   TMVClassification("", "TMVA_allSignals.root");

#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TCut.h"

#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVAGui.h"

static TTree* GetTreeOrReport(TFile* f, const char* tname, const TString& label) {
  if (!f || f->IsZombie()) {
    std::cerr << "[ERROR] Could not open file for: " << label << std::endl;
    return nullptr;
  }
  TTree* t = dynamic_cast<TTree*>(f->Get(tname));
  if (!t) {
    std::cerr << "[ERROR] Tree '" << tname << "' not found in: " << label << std::endl;
    std::cerr << "        Listing file contents:\n";
    f->ls();
    return nullptr;
  }
  std::cout << "[OK] Found tree '" << tname << "' in: " << label
            << "  (entries=" << t->GetEntries() << ")\n";
  return t;
}

void TMVClassification(TString dir = "", TString outFileName = "TMVA_allSignals.root")
{
  TMVA::Tools::Instance();

  const TString treeName = "AnalysisTree";

  // Output file
  TFile* outputFile = TFile::Open(dir + outFileName, "RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    std::cerr << "[ERROR] Cannot create output file: " << (dir + outFileName) << std::endl;
    return;
  }

  // Factory
  TMVA::Factory* factory = new TMVA::Factory(
    "MVAnalysis",
    outputFile,
    "!V:!Silent:Color:DrawProgressBar:Transformations=I;G:AnalysisType=Classification"
  );

  // DataLoader
  TMVA::DataLoader* dataloader = new TMVA::DataLoader("dataset");

  // ---------------------------------------------------------
  // INPUT FILES (reduced outputs from your analysis)
  // ---------------------------------------------------------
  // Signals
  TFile* f_sig12 = TFile::Open(dir + "output_signal_tth12gev.root");
  TFile* f_sig15 = TFile::Open(dir + "output_signal_tth15gev.root");
  TFile* f_sig20 = TFile::Open(dir + "output_signal_tth20gev.root");
  TFile* f_sig25 = TFile::Open(dir + "output_signal_tth25gev.root");
  TFile* f_sig30 = TFile::Open(dir + "output_signal_tth30gev.root");

  TTree* t_sig12 = GetTreeOrReport(f_sig12, treeName, "signal tth12 (output_signal_tth12gev.root)");
  TTree* t_sig15 = GetTreeOrReport(f_sig15, treeName, "signal tth15 (output_signal_tth15gev.root)");
  TTree* t_sig20 = GetTreeOrReport(f_sig20, treeName, "signal tth20 (output_signal_tth20gev.root)");
  TTree* t_sig25 = GetTreeOrReport(f_sig25, treeName, "signal tth25 (output_signal_tth25gev.root)");
  TTree* t_sig30 = GetTreeOrReport(f_sig30, treeName, "signal tth30 (output_signal_tth30gev.root)");

  if (!t_sig12 || !t_sig15 || !t_sig20 || !t_sig25 || !t_sig30) {
    std::cerr << "[FATAL] One or more SIGNAL trees missing. Aborting.\n";
    return;
  }

  // Backgrounds
  TFile* f_bkg_tt2l2nu = TFile::Open(dir + "output_ttbar.root");        // from TTto2L2Nu.root
  TFile* f_bkg_ttlnu2q = TFile::Open(dir + "output_TTtoLNu2Q.root");     // from TTtoLNu2Q.root
  TFile* f_bkg_tthhbb  = TFile::Open(dir + "output_TTH_Hbb.root");       // from TTHHtobb.root

  TTree* t_bkg_tt2l2nu = GetTreeOrReport(f_bkg_tt2l2nu, treeName, "TTto2L2Nu (output_ttbar.root)");
  TTree* t_bkg_ttlnu2q = GetTreeOrReport(f_bkg_ttlnu2q, treeName, "TTtoLNu2Q (output_TTtoLNu2Q.root)");
  TTree* t_bkg_tthhbb  = GetTreeOrReport(f_bkg_tthhbb,  treeName, "TTHHtobb (output_TTH_Hbb.root)");

  if (!t_bkg_tt2l2nu || !t_bkg_ttlnu2q || !t_bkg_tthhbb) {
    std::cerr << "[FATAL] One or more BACKGROUND trees missing. Aborting.\n";
    return;
  }

  // ---------------------------------------------------------
  // Event weights from your reduced trees
  // ---------------------------------------------------------
  dataloader->SetSignalWeightExpression("weight");
  dataloader->SetBackgroundWeightExpression("weight");

  // ---------------------------------------------------------
  // Add Signal trees (all masses) and Background trees
  // ---------------------------------------------------------
  dataloader->AddSignalTree(t_sig12, 1.0);
  dataloader->AddSignalTree(t_sig15, 1.0);
  dataloader->AddSignalTree(t_sig20, 1.0);
  dataloader->AddSignalTree(t_sig25, 1.0);
  dataloader->AddSignalTree(t_sig30, 1.0);

  dataloader->AddBackgroundTree(t_bkg_tt2l2nu, 1.0);
  dataloader->AddBackgroundTree(t_bkg_ttlnu2q, 1.0);
  dataloader->AddBackgroundTree(t_bkg_tthhbb,  1.0);

  // ---------------------------------------------------------
  // Variables (match AnalysisTree branches; your floats are /D)
  // ---------------------------------------------------------
  dataloader->AddVariable("hbb_m",    'D');
  dataloader->AddVariable("hbb_pt",   'D');
  // dataloader->AddVariable("hbb_eta",  'D');

  //dataloader->AddVariable("dbj1_m",   'D');
  //dataloader->AddVariable("dbj2_m",   'D');
  dataloader->AddVariable("dbj1_pt",  'D');
  dataloader->AddVariable("dbj2_pt",  'D');

  dataloader->AddVariable("dR_dbj12", 'D');
  dataloader->AddVariable("dM_dbj12", 'D');

  dataloader->AddVariable("met_pt",   'D');
  dataloader->AddVariable("HT",       'D');

  dataloader->AddVariable("njet",     'I');
  dataloader->AddVariable("nb",       'I');
  //dataloader->AddVariable("ndb",      'I');

  dataloader->AddVariable("mll",      'D');
  dataloader->AddVariable("dRll",     'D');

  dataloader->AddVariable("lep1_pt",  'D');
  dataloader->AddVariable("lep2_pt",  'D');

  // ---------------------------------------------------------
  // Optional preselection (currently none)
  // ---------------------------------------------------------
  TCut preselectionCut = "";

  // NOTE:
  // NormMode=NumEvents is fine; your "weight" expression still acts per-event.
  // If you want equalized training (common), use NormMode=EqualNumEvents.
 dataloader->PrepareTrainingAndTestTree(
  preselectionCut,
  "nTrain_Signal=0:nTrain_Background=0:"
  "SplitMode=Random:SplitSeed=12345:"
  "TrainTestSplit_Signal=0.6:TrainTestSplit_Background=0.6:"
  "NormMode=NumEvents:!V"
);


  // ---------------------------------------------------------
  // BDT
  // ---------------------------------------------------------
  factory->BookMethod(
    dataloader,
    TMVA::Types::kBDT,
    "BDT",
    "!H:!V:"
    "NTrees=800:"
    "MinNodeSize=2.5%:"
    "MaxDepth=3:"
    "BoostType=AdaBoost:"
    "AdaBoostBeta=0.5:"
    "UseBaggedBoost:"
    "BaggedSampleFraction=0.5:"
    "SeparationType=GiniIndex:"
    "nCuts=40"
  );

  // Run
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();

  outputFile->Close();

  delete factory;
  delete dataloader;

  // Close inputs
  if (f_sig12) f_sig12->Close();
  if (f_sig15) f_sig15->Close();
  if (f_sig20) f_sig20->Close();
  if (f_sig25) f_sig25->Close();
  if (f_sig30) f_sig30->Close();

  if (f_bkg_tt2l2nu) f_bkg_tt2l2nu->Close();
  if (f_bkg_ttlnu2q) f_bkg_ttlnu2q->Close();
  if (f_bkg_tthhbb)  f_bkg_tthhbb->Close();

  std::cout << "\n[OK] TMVA training finished. Output written to: " << (dir + outFileName) << "\n";

  TMVA::TMVAGui(dir + outFileName);
}
