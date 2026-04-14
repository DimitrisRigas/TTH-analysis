// runAnalysis.C
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>

// -------------------------------------------------------------------
// Sample selector
// -------------------------------------------------------------------
enum SampleType {
  kSignal12,
  kSignal15,
  kSignal20,
  kSignal25,
  kSignal30,
  kSignal60,
  kTTbar,
  kDYee,
  kDYmumu,
  kVBF_Hbb,
  kGGH_Hbb,
  kTTH_Hbb,
  kWW,
  kWZ,
  kZZ,
  kTBbarQtoLNu,
  kTBbarQto2Q,
  kTTtoLNu2Q,
  kTbarWplustoNu2Q,
  kTbarWplusto4Q,
  kTTW,   // NEW
  kTTZ    // NEW
};

void runAnalysis(SampleType sample = kSignal12)
{
  TString infile;

  switch (sample) {
    case kSignal12: infile="TTH12Gev.root"; break;
    case kSignal15: infile="TTH15Gev.root"; break;
    case kSignal20: infile="TTH20Gev.root"; break;
    case kSignal25: infile="TTH25Gev.root"; break;
    case kSignal30: infile="TTH30Gev.root"; break;
    case kSignal60: infile="TTH60Gev.root"; break;

    case kTTbar:    infile="TTto2L2Nu.root"; break;
    case kDYee:     infile="DYto2E4JETS.root"; break;
    case kDYmumu:   infile="DYto2MU4JETS.root"; break;

    case kVBF_Hbb:  infile="VBFHto2b.root"; break;
    case kGGH_Hbb:  infile="glugluHtobb.root"; break;
    case kTTH_Hbb:  infile="TTHHtobb.root"; break;

    case kWW:       infile="WW.root"; break;
    case kWZ:       infile="WZ.root"; break;
    case kZZ:       infile="ZZ.root"; break;

    case kTBbarQtoLNu:     infile="TBbarQtoLNu.root"; break;
    case kTBbarQto2Q:      infile="TBbarQto2Q.root"; break;
    case kTTtoLNu2Q:       infile="TTtoLNu2Q.root"; break;
    case kTbarWplustoNu2Q: infile="TbarWplustoNu2Q.root"; break;
    case kTbarWplusto4Q:   infile="TbarWplusto4Q.root"; break;

    case kTTW:      infile="TTW.root"; break;   // NEW
    case kTTZ:      infile="TTZ.root"; break;   // NEW

    default:
      std::cerr << "ERROR: Unknown sample type!\n";
      return;
  }

  std::cout << "Running on: " << infile << std::endl;

  TFile *f = TFile::Open(infile, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "ERROR: Cannot open file " << infile << std::endl;
    return;
  }

  TTree *tree = (TTree*)f->Get("Events");
  if (!tree) {
    std::cerr << "ERROR: Tree 'Events' not found in file!\n";
    f->Close();
    return;
  }

  // IMPORTANT: Do NOT include MyClass.h here.
  // We call MyClass through the interpreter (it is provided by MyClass2.C+).
  gROOT->ProcessLine(Form("MyClass analysis((TTree*)0x%lx); analysis.Loop();",
                          (ULong_t)tree));

  std::cout << "Analysis finished for file: " << infile << std::endl;
  f->Close();
}
