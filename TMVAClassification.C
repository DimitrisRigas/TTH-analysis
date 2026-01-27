// File: TMVA_Classification_MyAnalysis_tth12.C
//
// STEP-BY-STEP VERSION: only Signal (tth12gev) vs Background (TTH_Hbb)
// Other backgrounds are kept but COMMENTED OUT so you can re-enable later.
//
// Usage (ROOT):
//   root -l
//   .L TMVA_Classification.C
//   TMVClassification()
//
// Output: TMVA_tth12.root (default)

#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TCut.h"
#include "TSystem.h"
#include "TROOT.h"

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

// Helper: confirm tree name by trying a list of common names and printing contents
static void ConfirmTreeName(const TString& filePath) {
  std::cout << "\n============================================================\n";
  std::cout << "Checking tree name in file: " << filePath << "\n";
  std::cout << "============================================================\n";

  TFile* f = TFile::Open(filePath);
  if (!f || f->IsZombie()) {
    std::cerr << "[ERROR] Cannot open " << filePath << "\n";
    return;
  }

  std::cout << "[INFO] Top-level contents (f->ls()):\n";
  f->ls();

  const char* candidates[] = {"AnalysisTree", "my_tree", "tree", "Events"};
  for (const char* name : candidates) {
    TObject* obj = f->Get(name);
    if (obj && obj->InheritsFrom(TTree::Class())) {
      std::cout << "[FOUND] TTree name = '" << name << "'"
                << "  (entries=" << ((TTree*)obj)->GetEntries() << ")\n";
      std::cout << "        Branch list:\n";
      ((TTree*)obj)->Print();
      break;
    }
  }

  f->Close();
  delete f;
}

void TMVClassification(TString dir = "", TString outFileName = "TMVA_tth12.root")
{
  TMVA::Tools::Instance();

  const TString treeName = "AnalysisTree";

  // ----------------------------------------------------------------------
  // OPTIONAL: Quick tree-name confirmation (run once if needed)
  // ----------------------------------------------------------------------
  // ConfirmTreeName(dir + "output_signal_tth12gev.root");
  // ConfirmTreeName(dir + "output_TTH_Hbb.root");
  //
  // --- commented backgrounds you may add back later:
  // ConfirmTreeName(dir + "output_ttbar.root");
  // ConfirmTreeName(dir + "output_TTtoLNu2Q.root");
  // ----------------------------------------------------------------------

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

  // ----------------------------
  // Input files (ONLY what you want now)
  // ----------------------------
  TFile* f_sig = TFile::Open(dir + "output_signal_tth12gev.root");
  TTree* t_sig = GetTreeOrReport(f_sig, treeName, "signal tth12gev (output_signal_tth12gev.root)");
  if (!t_sig) {
    std::cerr << "[FATAL] Signal tree missing. Aborting.\n";
    return;
  }

  TFile* f_bkg_tthhbb = TFile::Open(dir + "output_TTH_Hbb.root");
  TTree* t_bkg_tthhbb = GetTreeOrReport(f_bkg_tthhbb, treeName, "ttH(H->bb) (output_TTH_Hbb.root)");
  if (!t_bkg_tthhbb) {
    std::cerr << "[FATAL] TTH_Hbb background tree missing. Aborting.\n";
    return;
  }

  // --- other backgrounds (KEEP COMMENTED FOR LATER) ---
   TFile* f_bkg_tt2l2nu = TFile::Open(dir + "output_ttbar.root");
   TFile* f_bkg_ttlnu2q = TFile::Open(dir + "output_TTtoLNu2Q.root");
   TTree* t_bkg_tt2l2nu = GetTreeOrReport(f_bkg_tt2l2nu, treeName, "ttbar dilepton (output_ttbar.root)");
   TTree* t_bkg_ttlnu2q = GetTreeOrReport(f_bkg_ttlnu2q, treeName, "ttbar semileptonic (output_TTtoLNu2Q.root)");
   if (!t_bkg_tt2l2nu || !t_bkg_ttlnu2q) { std::cerr << "[FATAL] ttbar background tree missing.\n"; return; }

  // ----------------------------
  // Use event-by-event weights from your reduced tree
  // ----------------------------
  dataloader->SetSignalWeightExpression("weight");
  dataloader->SetBackgroundWeightExpression("weight");

  // Add trees
  dataloader->AddSignalTree(t_sig, 1.0);

  // ONLY this background for now
  dataloader->AddBackgroundTree(t_bkg_tthhbb, 1.0);

  //--- add back later ---
   dataloader->AddBackgroundTree(t_bkg_tt2l2nu, 1.0);
   dataloader->AddBackgroundTree(t_bkg_ttlnu2q, 1.0);

  // ----------------------------
  // Variables: MUST exist in your AnalysisTree
  // ----------------------------
  dataloader->AddVariable("hbb_m",    'F');
  dataloader->AddVariable("hbb_pt",   'F');
  dataloader->AddVariable("hbb_eta",  'F');

  dataloader->AddVariable("dbj1_m",   'F');
  dataloader->AddVariable("dbj2_m",   'F');
  dataloader->AddVariable("dbj1_pt",  'F');
  dataloader->AddVariable("dbj2_pt",  'F');

  dataloader->AddVariable("dR_dbj12", 'F');
  dataloader->AddVariable("dM_dbj12", 'F');

  dataloader->AddVariable("met_pt",   'F');
  dataloader->AddVariable("HT",       'F');

  dataloader->AddVariable("njet",     'I');
  dataloader->AddVariable("nb",       'I');
  dataloader->AddVariable("ndb",      'I');

  dataloader->AddVariable("mll",      'F');
  dataloader->AddVariable("dRll",     'F');

  dataloader->AddVariable("lep1_pt",  'F');
  dataloader->AddVariable("lep2_pt",  'F');

  // ----------------------------
  // Preselection (start with NONE; add later if you want)
  // ----------------------------
  TCut preselectionCut = ""; // e.g. "weight!=0 && met_pt>0 && hbb_m>0";

  dataloader->PrepareTrainingAndTestTree(
    preselectionCut,
    "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V"
  );

  // ----------------------------
  // Methods (ONLY BDT for now)
  // ----------------------------
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

  // --- Likelihood (commented out for now; can re-enable later) ---
  // factory->BookMethod(
  //   dataloader,
  //   TMVA::Types::kLikelihood,
  //   "LikelihoodD",
  //   "!H:!V:TransformOutput:"
  //   "PDFInterpol=Spline2:"
  //   "NSmoothSig[0]=20:"
  //   "NSmoothBkg[0]=20:"
  //   "NSmooth=5:"
  //   "NAvEvtPerBin=50:"
  //   "VarTransform=Decorrelate"
  // );

  // Run
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();

  outputFile->Close();

  delete factory;
  delete dataloader;

  if (f_sig) f_sig->Close();
  if (f_bkg_tthhbb) f_bkg_tthhbb->Close();

  // --- close other backgrounds later ---
   if (f_bkg_tt2l2nu) f_bkg_tt2l2nu->Close();
   if (f_bkg_ttlnu2q) f_bkg_ttlnu2q->Close();

  std::cout << "\n[OK] TMVA training finished. Output written to: " << (dir + outFileName) << "\n";

  // If running in batch (-b), comment this line out.
  TMVA::TMVAGui(dir + outFileName);
}
