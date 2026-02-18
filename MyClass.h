//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jan 28 11:35:43 2026 by ROOT version 6.36.04
// from TTree Events/
// found on file: TTH12Gev.root
//////////////////////////////////////////////////////////

#ifndef MyClass_h
#define MyClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class MyClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   ULong64_t       event;
   UInt_t          luminosityBlock;
   Int_t           trigger_type;
   Bool_t          has_trigger;
   Float_t         PuppiMET_pt;
   Float_t         PuppiMET_phi;
   Int_t           nMuon;
   Float_t         Muon_pt[13];   //[nMuon]
   Float_t         Muon_eta[13];   //[nMuon]
   Float_t         Muon_phi[13];   //[nMuon]
   Int_t           Muon_charge[13];   //[nMuon]
   Bool_t          Muon_tightId[13];   //[nMuon]
   Bool_t          Muon_looseId[13];   //[nMuon]
   Float_t         Muon_mass[13];   //[nMuon]
   Float_t         Muon_pfRelIso04_all[13];   //[nMuon]
   Int_t           nElectron;
   Float_t         Electron_pt[8];   //[nElectron]
   Float_t         Electron_eta[8];   //[nElectron]
   Float_t         Electron_phi[8];   //[nElectron]
   Int_t           Electron_charge[8];   //[nElectron]
   UChar_t         Electron_cutBased[8];   //[nElectron]
   Float_t         Electron_mass[8];   //[nElectron]
   Float_t         Electron_pfRelIso03_all[8];   //[nElectron]
   UChar_t         Electron_seedGain[8];   //[nElectron]
   Float_t         Electron_r9[8];   //[nElectron]
   Float_t         Electron_superclusterEta[8];   //[nElectron]
   Bool_t          Electron_mvaIso_WP90[8];   //[nElectron]
   Bool_t          Electron_mvaIso_WP80[8];   //[nElectron]
   Int_t           nJet;
   Double_t        Jet_pt[21];   //[nJet]
   Float_t         Jet_eta[21];   //[nJet]
   Float_t         Jet_phi[21];   //[nJet]
   Double_t        Jet_mass[21];   //[nJet]
   Float_t         Jet_btagUParTAK4probbb[21];   //[nJet]
   Float_t         Jet_btagUParTAK4B[21];   //[nJet]
   Bool_t          Jet_passJetIdTightLepVeto[21];   //[nJet]
   UChar_t         Jet_hadronFlavour[21];   //[nJet]
   Float_t         Pileup_nTrueInt;
   Int_t           Pileup_nPU;
   UChar_t         PV_npvsGood;
   UChar_t         PV_npvs;
   Int_t           nGenPart;
   Float_t         GenPart_pt[128];   //[nGenPart]
   Float_t         GenPart_eta[128];   //[nGenPart]
   Float_t         GenPart_phi[128];   //[nGenPart]
   Float_t         GenPart_mass[128];   //[nGenPart]
   Int_t           GenPart_pdgId[128];   //[nGenPart]
   UShort_t        GenPart_statusFlags[128];   //[nGenPart]
   Short_t         GenPart_genPartIdxMother[128];   //[nGenPart]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_trigger_type;   //!
   TBranch        *b_has_trigger;   //!
   TBranch        *b_PuppiMET_pt;   //!
   TBranch        *b_PuppiMET_phi;   //!
   TBranch        *b_nMuon;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_tightId;   //!
   TBranch        *b_Muon_looseId;   //!
   TBranch        *b_Muon_mass;   //!
   TBranch        *b_Muon_pfRelIso04_all;   //!
   TBranch        *b_nElectron;   //!
   TBranch        *b_Electron_pt;   //!
   TBranch        *b_Electron_eta;   //!
   TBranch        *b_Electron_phi;   //!
   TBranch        *b_Electron_charge;   //!
   TBranch        *b_Electron_cutBased;   //!
   TBranch        *b_Electron_mass;   //!
   TBranch        *b_Electron_pfRelIso03_all;   //!
   TBranch        *b_Electron_seedGain;   //!
   TBranch        *b_Electron_r9;   //!
   TBranch        *b_Electron_superclusterEta;   //!
   TBranch        *b_Electron_mvaIso_WP90;   //!
   TBranch        *b_Electron_mvaIso_WP80;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_btagUParTAK4probbb;   //!
   TBranch        *b_Jet_btagUParTAK4B;   //!
   TBranch        *b_Jet_passJetIdTightLepVeto;   //!
   TBranch        *b_Jet_hadronFlavour;   //!
   TBranch        *b_Pileup_nTrueInt;   //!
   TBranch        *b_Pileup_nPU;   //!
   TBranch        *b_PV_npvsGood;   //!
   TBranch        *b_PV_npvs;   //!
   TBranch        *b_nGenPart;   //!
   TBranch        *b_GenPart_pt;   //!
   TBranch        *b_GenPart_eta;   //!
   TBranch        *b_GenPart_phi;   //!
   TBranch        *b_GenPart_mass;   //!
   TBranch        *b_GenPart_pdgId;   //!
   TBranch        *b_GenPart_statusFlags;   //!
   TBranch        *b_GenPart_genPartIdxMother;   //!

   MyClass(TTree *tree=0);
   virtual ~MyClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MyClass_cxx
MyClass::MyClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("TTH12Gev.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("TTH12Gev.root");
      }
      f->GetObject("Events",tree);

   }
   Init(tree);
}

MyClass::~MyClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MyClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MyClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MyClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("trigger_type", &trigger_type, &b_trigger_type);
   fChain->SetBranchAddress("has_trigger", &has_trigger, &b_has_trigger);
   fChain->SetBranchAddress("PuppiMET_pt", &PuppiMET_pt, &b_PuppiMET_pt);
   fChain->SetBranchAddress("PuppiMET_phi", &PuppiMET_phi, &b_PuppiMET_phi);
   fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
   fChain->SetBranchAddress("Muon_pt", Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon_eta", Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon_phi", Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon_charge", Muon_charge, &b_Muon_charge);
   fChain->SetBranchAddress("Muon_tightId", Muon_tightId, &b_Muon_tightId);
   fChain->SetBranchAddress("Muon_looseId", Muon_looseId, &b_Muon_looseId);
   fChain->SetBranchAddress("Muon_mass", Muon_mass, &b_Muon_mass);
   fChain->SetBranchAddress("Muon_pfRelIso04_all", Muon_pfRelIso04_all, &b_Muon_pfRelIso04_all);
   fChain->SetBranchAddress("nElectron", &nElectron, &b_nElectron);
   fChain->SetBranchAddress("Electron_pt", Electron_pt, &b_Electron_pt);
   fChain->SetBranchAddress("Electron_eta", Electron_eta, &b_Electron_eta);
   fChain->SetBranchAddress("Electron_phi", Electron_phi, &b_Electron_phi);
   fChain->SetBranchAddress("Electron_charge", Electron_charge, &b_Electron_charge);
   fChain->SetBranchAddress("Electron_cutBased", Electron_cutBased, &b_Electron_cutBased);
   fChain->SetBranchAddress("Electron_mass", Electron_mass, &b_Electron_mass);
   fChain->SetBranchAddress("Electron_pfRelIso03_all", Electron_pfRelIso03_all, &b_Electron_pfRelIso03_all);
   fChain->SetBranchAddress("Electron_seedGain", Electron_seedGain, &b_Electron_seedGain);
   fChain->SetBranchAddress("Electron_r9", Electron_r9, &b_Electron_r9);
   fChain->SetBranchAddress("Electron_superclusterEta", Electron_superclusterEta, &b_Electron_superclusterEta);
   fChain->SetBranchAddress("Electron_mvaIso_WP90", Electron_mvaIso_WP90, &b_Electron_mvaIso_WP90);
   fChain->SetBranchAddress("Electron_mvaIso_WP80", Electron_mvaIso_WP80, &b_Electron_mvaIso_WP80);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_btagUParTAK4probbb", Jet_btagUParTAK4probbb, &b_Jet_btagUParTAK4probbb);
   fChain->SetBranchAddress("Jet_btagUParTAK4B", Jet_btagUParTAK4B, &b_Jet_btagUParTAK4B);
   fChain->SetBranchAddress("Jet_passJetIdTightLepVeto", Jet_passJetIdTightLepVeto, &b_Jet_passJetIdTightLepVeto);
   fChain->SetBranchAddress("Jet_hadronFlavour", Jet_hadronFlavour, &b_Jet_hadronFlavour);
   fChain->SetBranchAddress("Pileup_nTrueInt", &Pileup_nTrueInt, &b_Pileup_nTrueInt);
   fChain->SetBranchAddress("Pileup_nPU", &Pileup_nPU, &b_Pileup_nPU);
   fChain->SetBranchAddress("PV_npvsGood", &PV_npvsGood, &b_PV_npvsGood);
   fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   fChain->SetBranchAddress("nGenPart", &nGenPart, &b_nGenPart);
   fChain->SetBranchAddress("GenPart_pt", GenPart_pt, &b_GenPart_pt);
   fChain->SetBranchAddress("GenPart_eta", GenPart_eta, &b_GenPart_eta);
   fChain->SetBranchAddress("GenPart_phi", GenPart_phi, &b_GenPart_phi);
   fChain->SetBranchAddress("GenPart_mass", GenPart_mass, &b_GenPart_mass);
   fChain->SetBranchAddress("GenPart_pdgId", GenPart_pdgId, &b_GenPart_pdgId);
   fChain->SetBranchAddress("GenPart_statusFlags", GenPart_statusFlags, &b_GenPart_statusFlags);
   fChain->SetBranchAddress("GenPart_genPartIdxMother", GenPart_genPartIdxMother, &b_GenPart_genPartIdxMother);
   Notify();
}

bool MyClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void MyClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MyClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MyClass_cxx
