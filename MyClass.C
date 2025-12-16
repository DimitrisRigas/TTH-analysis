#define MyClass_cxx
#include "MyClass.h"

#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TMath.h>
#include <TFile.h>
#include <TLorentzVector.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdio>

using namespace std;

// ============================================================================
// Helper Functions (geometry & kinematics)
// ============================================================================

static inline float DeltaPhi(float a, float b) {
   float d = a - b;
   while (d >  TMath::Pi()) d -= 2.f * TMath::Pi();
   while (d <= -TMath::Pi()) d += 2.f * TMath::Pi();
   return d;
}

static inline float DeltaR(float eta1, float phi1, float eta2, float phi2) {
   const float dphi = DeltaPhi(phi1, phi2);
   const float deta = eta1 - eta2;
   return std::sqrt(deta * deta + dphi * dphi);
}

// (currently unused, but kept for completeness)
static inline double massFromPtEtaPhiM(double pt1, double eta1, double phi1, double m1,
                                       double pt2, double eta2, double phi2, double m2) {
   const double px1 = pt1 * std::cos(phi1);
   const double py1 = pt1 * std::sin(phi1);
   const double pz1 = pt1 * std::sinh(eta1);
   const double E1  = std::sqrt(m1 * m1 + px1 * px1 + py1 * py1 + pz1 * pz1);

   const double px2 = pt2 * std::cos(phi2);
   const double py2 = pt2 * std::sin(phi2);
   const double pz2 = pt2 * std::sinh(eta2);
   const double E2  = std::sqrt(m2 * m2 + px2 * px2 + py2 * py2 + pz2 * pz2);

   const double E  = E1 + E2;
   const double px = px1 + px2;
   const double py = py1 + py2;
   const double pz = pz1 + pz2;
   const double M2 = E * E - (px * px + py * py + pz * pz);

   return (M2 > 0.0) ? std::sqrt(M2) : 0.0;
}

// ============================================================================
// Main Analysis Loop
// ============================================================================

void MyClass::Loop()
{
   if (fChain == nullptr) return;

   const Long64_t nentries = fChain->GetEntriesFast();

   std::cout << "Processing " << nentries << " entries from Tree." << std::endl;

   // detect if current file is the signal (TTH12Gev.root) or ttbar (TTto2L2Nu)
   Bool_t isSignal = false;
   Bool_t isTTbar  = false;

   TFile *currentFile = fChain->GetCurrentFile();
   std::string fname = currentFile ? std::string(currentFile->GetName()) : "";
   //identify file names
   if (fname.find("TTH12Gev") != std::string::npos)
     {
    isSignal = true;
     }
   else if (fname.find("TTto2L2Nu") != std::string::npos)
     {
       isTTbar = true;
     }


   // Output file name
   std::string outFileName = "output_unknown.root";
   if (isSignal)
      outFileName = "output_signal.root";
   else if (isTTbar)
      outFileName = "output_ttbar.root";

   // Create output ROOT file
   TFile *out = new TFile(outFileName.c_str(), "RECREATE");

   //
   //  wgt_prime: use for RAW histograms      (weight = 1)
   //  wgt      : physics per-event weight    (Nexp / Nstat)
   //
   //  Inside the RECO loop we use:
   //      const double weight = wgt_prime;
   //
   //  To switch to WEIGHTED histograms:
   //      const double weight = wgt;
   //
   // ------------------------------------------------------------------------
   float wgt_prime = 1.0f;
   float wgt       = 1.0f;   // will be set to Nexp/Nstat
   float sigma_pb  = 0.0f;   // cross section in pb
   const float L_int   = 108.96f;   // fb^-1
   const float Br_Wlnu = 0.1086f;   // BR(W->lnu) per lepton flavour (used only for signal)

   // Assign cross-section depending on sample
   if (isSignal) {
      sigma_pb = 0.5071f;   // ttH(aa) signal
   } else if (isTTbar) {
      sigma_pb = 98.04f;    // ttbar background
   }

   // Convert pb → fb
   const float sigma_fb = sigma_pb * 1000.0f;

   // Expected events
   float sigma = sigma_fb;
   if (isSignal) {
      // For signal apply BR(W->lnu)^2
      sigma *= (Br_Wlnu * Br_Wlnu);
   }

   float Nexp  = sigma * L_int;                     // expected events
   float Nstat = static_cast<float>(nentries);      // MC statistics

   // Correct per-event event weight:
   //   wgt = Nexp / Nstat
   if (Nstat > 0.0f)
      wgt = Nexp / Nstat;
   else
      wgt = 1.0f;

   const long long Nexp_int = static_cast<long long>(std::llround(Nexp));
   const double w = static_cast<double>(wgt);       // alias to keep your printout using "w"

   // ------------------------------------------------------------------------
   // 1. HISTOGRAM DEFINITIONS (RECO + GEN)
   // ------------------------------------------------------------------------

   // --- HIGGS (GEN) ---
   TH1F *h_H_pt   = new TH1F("h_H_pt",  "Higgs p_{T}; p_{T} [GeV]; Entries", 100, 0, 1000);
   TH1F *h_H_eta  = new TH1F("h_H_eta", "Higgs #eta; #eta; Entries",         60, -5, 5);
   TH1F *h_H_phi  = new TH1F("h_H_phi", "Higgs #phi; #phi; Entries",         64, -3.2, 3.2);
   TH1F *h_H_mass = new TH1F("h_H_mass","Higgs mass; m [GeV]; Entries",      80, 100, 150);

   // --- a BOSONS (a1, a2) (GEN) ---
   TH1F *h_a1_pt   = new TH1F("h_a1_pt", "Leading a (a1) p_{T}; p_{T} [GeV]; Entries", 100, 0, 1000);
   TH1F *h_a1_eta  = new TH1F("h_a1_eta","Leading a (a1) #eta; #eta; Entries",         60, -5, 5);
   TH1F *h_a1_phi  = new TH1F("h_a1_phi","Leading a (a1) #phi; #phi; Entries",         64, -3.2, 3.2);

   TH1F *h_a2_pt   = new TH1F("h_a2_pt", "Subleading a (a2) p_{T}; p_{T} [GeV]; Entries", 100, 0, 1000);
   TH1F *h_a2_eta  = new TH1F("h_a2_eta","Subleading a (a2) #eta; #eta; Entries",         60, -5, 5);
   TH1F *h_a2_phi  = new TH1F("h_a2_phi","Subleading a (a2) #phi; #phi; Entries",         64, -3.2, 3.2);

   TH1F *h_a_mass  = new TH1F("h_a_mass","a Boson Mass (Combined); m [GeV]; Entries", 100, 0, 100);
   TH1F *h_dR_aa   = new TH1F("h_dR_aa", "#Delta R(a1, a2); #Delta R; Entries",       100, 0, 5);

   // --- b FROM a (GEN) ---
   TH1F *h_dR_bb_a1 = new TH1F("h_dR_bb_a1","#Delta R(b,b) from Leading a (a1); #Delta R; Entries",   100, 0, 4);
   TH1F *h_dR_bb_a2 = new TH1F("h_dR_bb_a2","#Delta R(b,b) from Subleading a (a2); #Delta R; Entries",100, 0, 4);

   TH1F *h_b1_from_a_pt  = new TH1F("h_b1_from_a_pt","b1 from a p_{T}; p_{T}; Entries",100,0,1000);
   TH1F *h_b1_from_a_eta = new TH1F("h_b1_from_a_eta","b1 from a #eta; #eta; Entries",60,-5,5);
   TH1F *h_b1_from_a_phi = new TH1F("h_b1_from_a_phi","b1 from a #phi; #phi; Entries",64,-3.2,3.2);

   TH1F *h_b2_from_a_pt  = new TH1F("h_b2_from_a_pt","b2 from a p_{T}; p_{T}; Entries",100,0,1000);
   TH1F *h_b2_from_a_eta = new TH1F("h_b2_from_a_eta","b2 from a #eta; #eta; Entries",60,-5,5);
   TH1F *h_b2_from_a_phi = new TH1F("h_b2_from_a_phi","b2 from a #phi; #phi; Entries",64,-3.2,3.2);

   TH1F *h_b3_from_a_pt  = new TH1F("h_b3_from_a_pt","b3 from a p_{T}; p_{T}; Entries",100,0,1000);
   TH1F *h_b3_from_a_eta = new TH1F("h_b3_from_a_eta","b3 from a #eta; #eta; Entries",60,-5,5);
   TH1F *h_b3_from_a_phi = new TH1F("h_b3_from_a_phi","b3 from a #phi; #phi; Entries",64,-3.2,3.2);

   TH1F *h_b4_from_a_pt  = new TH1F("h_b4_from_a_pt","b4 from a p_{T}; p_{T}; Entries",100,0,1000);
   TH1F *h_b4_from_a_eta = new TH1F("h_b4_from_a_eta","b4 from a #eta; #eta; Entries",60,-5,5);
   TH1F *h_b4_from_a_phi = new TH1F("h_b4_from_a_phi","b4 from a #phi; #phi; Entries",64,-3.2,3.2);

   // --- TOPS (GEN, t1, t2 by ±6 selection) ---
   TH1F *h_t1_pt   = new TH1F("h_t1_pt",  "Top 1 p_{T}; p_{T} [GeV]; Entries",100,0,1000);
   TH1F *h_t1_eta  = new TH1F("h_t1_eta", "Top 1 #eta; #eta; Entries",        60,-5,5);
   TH1F *h_t1_phi  = new TH1F("h_t1_phi", "Top 1 #phi; #phi; Entries",        64,-3.2,3.2);
   TH1F *h_t1_mass = new TH1F("h_t1_mass","Top 1 mass; m [GeV]; Entries",     100,100,250);

   TH1F *h_t2_pt   = new TH1F("h_t2_pt",  "Top 2 p_{T}; p_{T} [GeV]; Entries",100,0,1000);
   TH1F *h_t2_eta  = new TH1F("h_t2_eta", "Top 2 #eta; #eta; Entries",        60,-5,5);
   TH1F *h_t2_phi  = new TH1F("h_t2_phi", "Top 2 #phi; #phi; Entries",        64,-3.2,3.2);
   TH1F *h_t2_mass = new TH1F("h_t2_mass","Top 2 mass; m [GeV]; Entries",     100,100,250);

   // --- b FROM TOPS (GEN) ---
   TH1F *h_b1_top_pt   = new TH1F("h_b1_top_pt","b from Top 1 p_{T}; p_{T}; Entries",100,0,1000);
   TH1F *h_b1_top_eta  = new TH1F("h_b1_top_eta","b from Top 1 #eta; #eta; Entries",60,-5,5);
   TH1F *h_b1_top_phi  = new TH1F("h_b1_top_phi","b from Top 1 #phi; #phi; Entries",64,-3.2,3.2);

   TH1F *h_b2_top_pt   = new TH1F("h_b2_top_pt","b from Top 2 p_{T}; p_{T}; Entries",100,0,1000);
   TH1F *h_b2_top_eta  = new TH1F("h_b2_top_eta","b from Top 2 #eta; #eta; Entries",60,-5,5);
   TH1F *h_b2_top_phi  = new TH1F("h_b2_top_phi","b from Top 2 #phi; #phi; Entries",64,-3.2,3.2);

   TH1F *h_dR_Wb = new TH1F("h_dR_Wb","#Delta R(W, b) from Top; #Delta R; Entries",100,0,5);

   // --- W BOSONS (GEN) ---
   TH1F *h_W_had_pt  = new TH1F("h_W_had_pt", "Hadronic W p_{T}; p_{T} [GeV]; Entries",100,0,1000);
   TH1F *h_W_had_eta = new TH1F("h_W_had_eta","Hadronic W #eta; #eta; Entries",        60,-5,5);
   TH1F *h_W_lep_pt  = new TH1F("h_W_lep_pt", "Leptonic W p_{T}; p_{T} [GeV]; Entries",      100,0,1000);
   TH1F *h_W_lep_eta = new TH1F("h_W_lep_eta","Leptonic W #eta; #eta; Entries",        60,-5,5);

   TH1F *h_q_from_W_pt   = new TH1F("h_q_from_W_pt",  "Quarks from W p_{T}; p_{T}; Entries",  100,0,1000);
   TH1F *h_lep_from_W_pt = new TH1F("h_lep_from_W_pt","Leptons from W p_{T}; p_{T}; Entries",100,0,1000);
   TH1F *h_nu_from_W_pt  = new TH1F("h_nu_from_W_pt", "Neutrinos from W p_{T}; p_{T}; Entries",100,0,1000);

   // --- RECO: ELECTRONS, MUONS, JETS (PRE-SELECTION only pT/eta) ---

   // multiplicities
   TH1F *h_nEle = new TH1F("h_nEle","Electron multiplicity; N_{e}; Entries",   10, 0, 10);
   TH1F *h_nMu  = new TH1F("h_nMu", "Muon multiplicity; N_{#mu}; Entries",     10, 0, 10);
   TH1F *h_nJet = new TH1F("h_nJet","Jet multiplicity; N_{jets}; Entries",     30, 0, 30);
   TH1F *h_nCJet = new TH1F("h_nCJet","Clean jet multiplicity; N_{clean jets}; Entries",30,0,30);

   // all electrons
   TH1F *h_ele_pt  = new TH1F("h_ele_pt", "Electron p_{T}; p_{T} [GeV]; Entries",100,0,500);
   TH1F *h_ele_eta = new TH1F("h_ele_eta","Electron #eta; #eta; Entries",       50,-2.5,2.5);
   TH1F *h_ele_phi = new TH1F("h_ele_phi","Electron #phi; #phi; Entries",       64,-3.2,3.2);

   // all muons
   TH1F *h_mu_pt  = new TH1F("h_mu_pt", "Muon p_{T}; p_{T} [GeV]; Entries",100,0,500);
   TH1F *h_mu_eta = new TH1F("h_mu_eta","Muon #eta; #eta; Entries",       50,-2.5,2.5);
   TH1F *h_mu_phi = new TH1F("h_mu_phi","Muon #phi; #phi; Entries",       64,-3.2,3.2);

   // all jets
   TH1F *h_jet_pt  = new TH1F("h_jet_pt", "Jet p_{T}; p_{T} [GeV]; Entries",100,0,500);
   TH1F *h_jet_eta = new TH1F("h_jet_eta","Jet #eta; #eta; Entries",       50,-2.5,2.5);
   TH1F *h_jet_phi = new TH1F("h_jet_phi","Jet #phi; #phi; Entries",       64,-3.2,3.2);

   // NEW: ΔR(jet, lepton) histograms (before and after cleaning)
   TH1F *h_dR_jet_ele_before = new TH1F("h_dR_jet_ele_before",
                                        "#Delta R(jet, electron) BEFORE cleaning; #Delta R; Entries",
                                        100, 0, 5);
   TH1F *h_dR_jet_mu_before  = new TH1F("h_dR_jet_mu_before",
                                        "#Delta R(jet, muon) BEFORE cleaning; #Delta R; Entries",
                                        100, 0, 5);
   TH1F *h_dR_jet_ele_after  = new TH1F("h_dR_jet_ele_after",
                                        "#Delta R(jet, electron) AFTER cleaning; #Delta R; Entries",
                                        100, 0, 5);
   TH1F *h_dR_jet_mu_after   = new TH1F("h_dR_jet_mu_after",
                                        "#Delta R(jet, muon) AFTER cleaning; #Delta R; Entries",
                                        100, 0, 5);

   // --- Reconstructed Higgs kinematics ---
   TH1F *h_Hdbb_pt  = new TH1F("h_Hdbb_pt",
                               "Reconstructed Higgs p_{T}; p_{T} [GeV]; Entries",
                               100, 0, 1000);

   TH1F *h_Hdbb_eta = new TH1F("h_Hdbb_eta",
                               "Reconstructed Higgs #eta; #eta; Entries",
                               60, -5, 5);

   // --- Double-b jet kinematics (leading + subleading) ---
   TH1F *h_dbj1_pt  = new TH1F("h_dbj1_pt",
                               "Double-b jet 1 p_{T}; p_{T} [GeV]; Entries",
                               100, 0, 500);

   TH1F *h_dbj1_eta = new TH1F("h_dbj1_eta",
                               "Double-b jet 1 #eta; #eta; Entries",
                               60, -2.5, 2.5);

   TH1F *h_dbj2_pt  = new TH1F("h_dbj2_pt",
                               "Double-b jet 2 p_{T}; p_{T} [GeV]; Entries",
                               100, 0, 500);

   TH1F *h_dbj2_eta = new TH1F("h_dbj2_eta",
                               "Double-b jet 2 #eta; #eta; Entries",
                               60, -2.5, 2.5);

   // ranked objects: 4 jets, 2 electrons, 2 muons
   const int MAXJET   = 4;
   const int MAXELE   = 2;
   const int MAXMU    = 2;
   const int MAXBJET  = 2;  // ranked b-jets
   const int MAXDBJET = 2;  // ranked double-b-tagged jets

   TH1F *h_j_pt [MAXJET];
   TH1F *h_j_eta[MAXJET];
   TH1F *h_j_phi[MAXJET];

   for (int i = 0; i < MAXJET; ++i) {
      std::string idx = std::to_string(i + 1);
      h_j_pt[i]  = new TH1F(("h_j" + idx + "_pt" ).c_str(), ("Jet " + idx + " p_{T}; p_{T} [GeV]; Entries").c_str(),100,0,500);
      h_j_eta[i] = new TH1F(("h_j" + idx + "_eta").c_str(), ("Jet " + idx + " #eta; #eta; Entries").c_str(),        50,-2.5,2.5);
      h_j_phi[i] = new TH1F(("h_j" + idx + "_phi").c_str(), ("Jet " + idx + " #phi; #phi; Entries").c_str(),        64,-3.2,3.2);
   }

   // --- NEW: B-JET HISTOGRAMS (RECO) ---
   TH1F *h_nBjet    = new TH1F("h_nBjet",    "b-jet multiplicity; N_{bjets}; Entries", 15, 0, 15);
   TH1F *h_bjet_pt  = new TH1F("h_bjet_pt",  "All b-jets p_{T}; p_{T} [GeV]; Entries", 100, 0, 500);
   TH1F *h_bjet_eta = new TH1F("h_bjet_eta", "All b-jets #eta; #eta; Entries",         50, -2.5, 2.5);
   TH1F *h_bjet_phi = new TH1F("h_bjet_phi", "All b-jets #phi; #phi; Entries",         64, -3.2, 3.2);

   // Ranked b-jets (top 2)
   TH1F *h_bj_pt [MAXBJET];
   TH1F *h_bj_eta[MAXBJET];
   TH1F *h_bj_phi[MAXBJET];

   for (int i = 0; i < MAXBJET; ++i) {
      std::string idx = std::to_string(i + 1);
      h_bj_pt[i]  = new TH1F(("h_bj" + idx + "_pt" ).c_str(), ("b-jet " + idx + " p_{T}; p_{T} [GeV]; Entries").c_str(), 100, 0, 500);
      h_bj_eta[i] = new TH1F(("h_bj" + idx + "_eta").c_str(), ("b-jet " + idx + " #eta; #eta; Entries").c_str(),         50, -2.5, 2.5);
      h_bj_phi[i] = new TH1F(("h_bj" + idx + "_phi").c_str(), ("b-jet " + idx + " #phi; #phi; Entries").c_str(),         64, -3.2, 3.2);
   }

   // --- NEW: DOUBLE-B-TAGGED JET HISTOGRAMS (RECO) ---
   TH1F *h_nDoubleB    = new TH1F("h_nDoubleB",    "double-b jet multiplicity; N_{dbjets}; Entries", 10, 0, 10);
   TH1F *h_dbjet_pt    = new TH1F("h_dbjet_pt",    "All double-b jets p_{T}; p_{T} [GeV]; Entries",  100, 0, 500);
   TH1F *h_dbjet_eta   = new TH1F("h_dbjet_eta",   "All double-b jets #eta; #eta; Entries",          50, -2.5, 2.5);
   TH1F *h_dbjet_phi   = new TH1F("h_dbjet_phi",   "All double-b jets #phi; #phi; Entries",          64, -3.2, 3.2);

   TH1F *h_dbj_pt [MAXDBJET];
   TH1F *h_dbj_eta[MAXDBJET];
   TH1F *h_dbj_phi[MAXDBJET];

   for (int i = 0; i < MAXDBJET; ++i) {
      std::string idx = std::to_string(i + 1);
      h_dbj_pt[i]  = new TH1F(("h_dbj" + idx + "_pt" ).c_str(), ("double-b jet " + idx + " p_{T}; p_{T} [GeV]; Entries").c_str(), 100, 0, 500);
      h_dbj_eta[i] = new TH1F(("h_dbj" + idx + "_eta").c_str(), ("double-b jet " + idx + " #eta; #eta; Entries").c_str(),         50, -2.5, 2.5);
      h_dbj_phi[i] = new TH1F(("h_dbj" + idx + "_phi").c_str(), ("double-b jet " + idx + " #phi; #phi; Entries").c_str(),         64, -3.2, 3.2);
   }

   TH1F *h_e_pt_rank [MAXELE];
   TH1F *h_e_eta_rank[MAXELE];
   TH1F *h_e_phi_rank[MAXELE];

   for (int i = 0; i < MAXELE; ++i) {
      std::string idx = std::to_string(i + 1);
      h_e_pt_rank[i]  = new TH1F(("h_e" + idx + "_pt" ).c_str(), ("Electron " + idx + " p_{T}; p_{T} [GeV]; Entries").c_str(),100,0,500);
      h_e_eta_rank[i] = new TH1F(("h_e" + idx + "_eta").c_str(), ("Electron " + idx + " #eta; #eta; Entries").c_str(),        50,-2.5,2.5);
      h_e_phi_rank[i] = new TH1F(("h_e" + idx + "_phi").c_str(), ("Electron " + idx + " #phi; #phi; Entries").c_str(),        64,-3.2,3.2);
   }

   TH1F *h_mu_pt_rank [MAXMU];
   TH1F *h_mu_eta_rank[MAXMU];
   TH1F *h_mu_phi_rank[MAXMU];

   for (int i = 0; i < MAXMU; ++i) {
      std::string idx = std::to_string(i + 1);
      h_mu_pt_rank[i]  = new TH1F(("h_mu" + idx + "_pt" ).c_str(), ("Muon " + idx + " p_{T}; p_{T} [GeV]; Entries").c_str(),100,0,500);
      h_mu_eta_rank[i] = new TH1F(("h_mu" + idx + "_eta").c_str(), ("Muon " + idx + " #eta; #eta; Entries").c_str(),        50,-2.5,2.5);
      h_mu_phi_rank[i] = new TH1F(("h_mu" + idx + "_phi").c_str(), ("Muon " + idx + " #phi; #phi; Entries").c_str(),        64,-3.2,3.2);
   }

   // --- MET (RECO, pre-selection, before cuts) ---
   TH1F *h_pre_MET_pt  = new TH1F("h_pre_MET_pt", "Puppi MET E_{T} (pre-selection); MET [GeV]; Entries", 100, 0, 500);
   TH1F *h_pre_MET_phi = new TH1F("h_pre_MET_phi","Puppi MET #phi (pre-selection); #phi; Entries",        64,-3.2, 3.2);

   // --- NEW HISTOGRAMS FOR REAL ANALYSIS ---
   TH1F *h_Hdbb_mass = new TH1F("h_Hdbb_mass",
                                "Reconstructed Higgs mass from double-b jets; m_{bb} [GeV]; Entries",
                                80, 0, 500);
   TH1F *h_MET_pt_final = new TH1F("h_MET_pt_final",
                                   "MET (final selection); MET [GeV]; Entries",
                                   100, 0, 500);
   TH1F *h_MET_phi_final = new TH1F("h_MET_phi_final",
                                    "MET #phi (final selection); #phi; Entries",
                                    64, -3.2, 3.2);
   TH1F *h_dbj1_mass = new TH1F("h_dbj1_mass",
                                "Double-b jet 1 mass; m [GeV]; Entries",
                                80, 0, 80);
   TH1F *h_dbj2_mass = new TH1F("h_dbj2_mass",
                                "Double-b jet 2 mass; m [GeV]; Entries",
                                80, 0, 80);
   TH1F *h_HT = new TH1F("h_HT",
                         "H_{T} scalar sum of jet p_{T}; H_{T} [GeV]; Entries",
                         100, 0, 2000);

   // ==========================================================
   // NEW REAL-ANALYSIS HISTOGRAMS (after full selection)
   // ==========================================================
   //
   // IMPORTANT FIX: these ROOT names must be UNIQUE (not "h_dbj1_phi" etc),
   // because you already have ranked histograms named "h_dbj1_phi" from h_dbj_phi[].
   //
   TH1F *h_dR_dbj12 = new TH1F(
       "h_dR_dbj12_final",
       "#Delta R(double-b jet 1, double-b jet 2) (final); #Delta R; Entries",
       100, 0, 5
   );

   TH1F *h_dbj1_phi = new TH1F(
       "h_dbj1_phi_final",
       "Double-b jet 1 #phi (final); #phi; Entries",
       64, -3.2, 3.2
   );

   TH1F *h_dbj2_phi = new TH1F(
       "h_dbj2_phi_final",
       "Double-b jet 2 #phi (final); #phi; Entries",
       64, -3.2, 3.2
   );

   TH1F *h_lep1_pt = new TH1F(
       "h_lep1_pt_final",
       "Leading lepton p_{T} (final); p_{T} [GeV]; Entries",
       100, 0, 500
   );

   TH1F *h_lep1_eta = new TH1F(
       "h_lep1_eta_final",
       "Leading lepton #eta (final); #eta; Entries",
       60, -2.5, 2.5
   );

   TH1F *h_lep1_phi = new TH1F(
       "h_lep1_phi_final",
       "Leading lepton #phi (final); #phi; Entries",
       64, -3.2, 3.2
   );

   TH1F *h_lep2_pt = new TH1F(
       "h_lep2_pt_final",
       "Subleading lepton p_{T} (final); p_{T} [GeV]; Entries",
       100, 0, 500
   );

   TH1F *h_lep2_eta = new TH1F(
       "h_lep2_eta_final",
       "Subleading lepton #eta (final); #eta; Entries",
       60, -2.5, 2.5
   );

   TH1F *h_lep2_phi = new TH1F(
       "h_lep2_phi_final",
       "Subleading lepton #phi (final); #phi; Entries",
       64, -3.2, 3.2
   );

   // ------------------------------------------------------------------------
   // Cut-flow counters (for RECO only)
   // ------------------------------------------------------------------------

   Long64_t nCut1 = 0; // N leptons >= 2
   Long64_t nCut2 = 0; // OS lepton pair
   Long64_t nCut3 = 0; // N jets >= 4 (using CLEAN jets)
   Long64_t nCut4 = 0; // N b-jets >= 2
   Long64_t nCut5 = 0; // N double-b-jets >= 2
   Long64_t nCut6 = 0; // MET >= 40

   // Helper for sorting GEN particles by pT
   auto sortParticles = [&](std::vector<int> &indices) {
       std::sort(indices.begin(), indices.end(),
                 [&](int i, int j) { return GenPart_pt[i] > GenPart_pt[j]; });
   };

   // Struct to unify electrons/muons into a single lepton container
   struct RecoLepton {
       TLorentzVector p4;
       int charge;   // +/-1
       int flavour;  // 0 = electron, 1 = muon
   };

   // ========================================================================
   // 2. EVENT LOOP
   // ========================================================================
   for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      fChain->GetEntry(jentry);

      // Run GENERATOR LEVEL analysis , ONLY IF input file is the Signal:
      if (isSignal) {

         // ------------------------------------------------------------------
         // 2B. GEN-LEVEL ANALYSIS
         // ------------------------------------------------------------------

         std::vector<int> final_tops;
         std::vector<int> final_higgs;

         int t1 = -1;
         int t2 = -1;

         for (int i = 0; i < nGenPart; ++i) {
             if      (GenPart_pdgId[i] ==  6 && t1 == -1) t1 = i;
             else if (GenPart_pdgId[i] == -6 && t2 == -1) t2 = i;
             if (t1 != -1 && t2 != -1) break;
         }

         if (t1 != -1) final_tops.push_back(t1);
         if (t2 != -1) final_tops.push_back(t2);

         for (int i = 0; i < nGenPart; ++i) {
            const int pdg = std::abs(GenPart_pdgId[i]);
            if (pdg != 25) continue;

            bool hasHiggsDaughter = false;
            for (int j = 0; j < nGenPart; ++j) {
                if (GenPart_genPartIdxMother[j] == i &&
                    std::abs(GenPart_pdgId[j]) == 25) {
                    hasHiggsDaughter = true;
                    break;
                }
            }
            if (!hasHiggsDaughter) final_higgs.push_back(i);
         }

         for (int ih : final_higgs) {

             h_H_pt  ->Fill(GenPart_pt[ih]);
             h_H_eta ->Fill(GenPart_eta[ih]);
             h_H_phi ->Fill(GenPart_phi[ih]);
             h_H_mass->Fill(GenPart_mass[ih]);

             std::vector<int> a_indices;
             for (int j = 0; j < nGenPart; ++j) {
                if (GenPart_genPartIdxMother[j] == ih &&
                    std::abs(GenPart_pdgId[j]) == 36) {
                    a_indices.push_back(j);
                }
             }

             if (!a_indices.empty()) {
                 int ia1 = a_indices[0];

                 h_a1_pt  ->Fill(GenPart_pt[ia1]);
                 h_a1_eta ->Fill(GenPart_eta[ia1]);
                 h_a1_phi ->Fill(GenPart_phi[ia1]);
                 h_a_mass ->Fill(GenPart_mass[ia1]);

                 std::vector<int> b_from_a1;
                 for (int k = 0; k < nGenPart; ++k) {
                     if (GenPart_genPartIdxMother[k] == ia1 &&
                         std::abs(GenPart_pdgId[k]) == 5) {
                         b_from_a1.push_back(k);
                     }
                 }
                 if (b_from_a1.size() == 2) {
                     h_dR_bb_a1->Fill(
                         DeltaR(GenPart_eta[b_from_a1[0]], GenPart_phi[b_from_a1[0]],
                                GenPart_eta[b_from_a1[1]], GenPart_phi[b_from_a1[1]])
                     );
                 }
             }

             if (a_indices.size() > 1) {
                 int ia2 = a_indices[1];

                 h_a2_pt  ->Fill(GenPart_pt[ia2]);
                 h_a2_eta ->Fill(GenPart_eta[ia2]);
                 h_a2_phi ->Fill(GenPart_phi[ia2]);
                 h_a_mass ->Fill(GenPart_mass[ia2]);

                 const double dRaa = DeltaR(GenPart_eta[a_indices[0]], GenPart_phi[a_indices[0]],
                                            GenPart_eta[a_indices[1]], GenPart_phi[a_indices[1]]);
                 h_dR_aa->Fill(dRaa);

                 std::vector<int> b_from_a2;
                 for (int k = 0; k < nGenPart; ++k) {
                     if (GenPart_genPartIdxMother[k] == ia2 &&
                         std::abs(GenPart_pdgId[k]) == 5) {
                         b_from_a2.push_back(k);
                     }
                 }
                 if (b_from_a2.size() == 2) {
                     h_dR_bb_a2->Fill(
                         DeltaR(GenPart_eta[b_from_a2[0]], GenPart_phi[b_from_a2[0]],
                                GenPart_eta[b_from_a2[1]], GenPart_phi[b_from_a2[1]])
                     );
                 }
             }

             std::vector<int> all_b_from_a;
             for (int ia : a_indices) {
                 for (int k = 0; k < nGenPart; ++k) {
                     if (GenPart_genPartIdxMother[k] == ia &&
                         std::abs(GenPart_pdgId[k]) == 5) {
                         all_b_from_a.push_back(k);
                     }
                 }
             }

             sortParticles(all_b_from_a);

             if (all_b_from_a.size() > 0) {
                 h_b1_from_a_pt ->Fill(GenPart_pt[all_b_from_a[0]]);
                 h_b1_from_a_eta->Fill(GenPart_eta[all_b_from_a[0]]);
                 h_b1_from_a_phi->Fill(GenPart_phi[all_b_from_a[0]]);
             }
             if (all_b_from_a.size() > 1) {
                 h_b2_from_a_pt ->Fill(GenPart_pt[all_b_from_a[1]]);
                 h_b2_from_a_eta->Fill(GenPart_eta[all_b_from_a[1]]);
                 h_b2_from_a_phi->Fill(GenPart_phi[all_b_from_a[1]]);
             }
             if (all_b_from_a.size() > 2) {
                 h_b3_from_a_pt ->Fill(GenPart_pt[all_b_from_a[2]]);
                 h_b3_from_a_eta->Fill(GenPart_eta[all_b_from_a[2]]);
                 h_b3_from_a_phi->Fill(GenPart_phi[all_b_from_a[2]]);
             }
             if (all_b_from_a.size() > 3) {
                 h_b4_from_a_pt ->Fill(GenPart_pt[all_b_from_a[3]]);
                 h_b4_from_a_eta->Fill(GenPart_eta[all_b_from_a[3]]);
                 h_b4_from_a_phi->Fill(GenPart_phi[all_b_from_a[3]]);
             }
         }

         auto processTop = [&](int idx, int rank) {

             if (rank == 1) {
                 h_t1_pt  ->Fill(GenPart_pt[idx]);
                 h_t1_eta ->Fill(GenPart_eta[idx]);
                 h_t1_phi ->Fill(GenPart_phi[idx]);
                 h_t1_mass->Fill(GenPart_mass[idx]);
             } else {
                 h_t2_pt  ->Fill(GenPart_pt[idx]);
                 h_t2_eta ->Fill(GenPart_eta[idx]);
                 h_t2_phi ->Fill(GenPart_phi[idx]);
                 h_t2_mass->Fill(GenPart_mass[idx]);
             }

             int  final_top = idx;
             bool search    = true;
             while (search) {
                 search = false;
                 for (int j = 0; j < nGenPart; ++j) {
                     if (GenPart_genPartIdxMother[j] == final_top &&
                         std::abs(GenPart_pdgId[j]) == 6) {
                         final_top = j;
                         search    = true;
                         break;
                     }
                 }
             }

             int b_idx = -1;
             int w_idx = -1;

             for (int j = 0; j < nGenPart; ++j) {
                 if (GenPart_genPartIdxMother[j] != final_top) continue;
                 if (std::abs(GenPart_pdgId[j]) == 5)  b_idx = j;
                 if (std::abs(GenPart_pdgId[j]) == 24) w_idx = j;
             }

             if (b_idx != -1) {
                 if (rank == 1) {
                     h_b1_top_pt ->Fill(GenPart_pt[b_idx]);
                     h_b1_top_eta->Fill(GenPart_eta[b_idx]);
                     h_b1_top_phi->Fill(GenPart_phi[b_idx]);
                 } else {
                     h_b2_top_pt ->Fill(GenPart_pt[b_idx]);
                     h_b2_top_eta->Fill(GenPart_eta[b_idx]);
                     h_b2_top_phi->Fill(GenPart_phi[b_idx]);
                 }
             }

             if (w_idx != -1) {

                 int  final_w = w_idx;
                 bool isFinal = false;

                 while (!isFinal) {
                     bool hasDaughter = false;
                     for (int k = 0; k < nGenPart; ++k) {
                         if (GenPart_genPartIdxMother[k] == final_w &&
                             std::abs(GenPart_pdgId[k]) == 24) {
                             final_w     = k;
                             hasDaughter = true;
                             break;
                         }
                     }
                     if (!hasDaughter) isFinal = true;
                 }

                 bool isHad = false;
                 bool isLep = false;

                 for (int k = 0; k < nGenPart; ++k) {
                     if (GenPart_genPartIdxMother[k] != final_w) continue;

                     const int absId = std::abs(GenPart_pdgId[k]);

                     if (absId <= 5) {
                         isHad = true;
                         h_q_from_W_pt->Fill(GenPart_pt[k]);
                     }

                     if (absId >= 11 && absId <= 16) {
                         isLep = true;
                         if (absId % 2 == 1)
                             h_lep_from_W_pt->Fill(GenPart_pt[k]);
                         else
                             h_nu_from_W_pt ->Fill(GenPart_pt[k]);
                     }
                 }

                 if (isHad) {
                     h_W_had_pt ->Fill(GenPart_pt[final_w]);
                     h_W_had_eta->Fill(GenPart_eta[final_w]);
                 }
                 if (isLep) {
                     h_W_lep_pt ->Fill(GenPart_pt[final_w]);
                     h_W_lep_eta->Fill(GenPart_eta[final_w]);
                 }

                 if (b_idx != -1) {
                     h_dR_Wb->Fill(
                         DeltaR(GenPart_eta[final_w], GenPart_phi[final_w],
                                GenPart_eta[b_idx],    GenPart_phi[b_idx])
                     );
                 }
             }
         };

         if (final_tops.size() > 0) processTop(final_tops[0], 1);
         if (final_tops.size() > 1) processTop(final_tops[1], 2);

      } // end if (isSignal)

      // ---------------------------------------------------------------------
      // 2A. RECO-LEVEL ANALYSIS
      // ---------------------------------------------------------------------

      const double weight = static_cast<double>(wgt_prime);

      h_pre_MET_pt ->Fill(PuppiMET_pt,  weight);
      h_pre_MET_phi->Fill(PuppiMET_phi, weight);

      std::vector<TLorentzVector> vec_ele;
      std::vector<TLorentzVector> vec_muons;
      std::vector<TLorentzVector> vec_jet;
      std::vector<TLorentzVector> vec_cjet;
      std::vector<TLorentzVector> vec_bjets;
      std::vector<TLorentzVector> vec_doublebjets;
      std::vector<RecoLepton>     leptons;

      for (int i = 0; i < nElectron; ++i) {
         TLorentzVector p4;
         p4.SetPtEtaPhiM(Electron_pt[i], Electron_eta[i], Electron_phi[i], Electron_mass[i]);

         if (p4.Pt() <= 20.0) continue;
         if (std::fabs(p4.Eta()) >= 2.5) continue;

         vec_ele.push_back(p4);

         RecoLepton lep;
         lep.p4      = p4;
         lep.charge  = Electron_charge[i];
         lep.flavour = 0;
         leptons.push_back(lep);
      }

      for (int i = 0; i < nMuon; ++i) {
         TLorentzVector p4;
         p4.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);

         if (p4.Pt() <= 20.0) continue;
         if (std::fabs(p4.Eta()) >= 2.5) continue;

         vec_muons.push_back(p4);

         RecoLepton lep;
         lep.p4      = p4;
         lep.charge  = Muon_charge[i];
         lep.flavour = 1;
         leptons.push_back(lep);
      }

      for (int i = 0; i < nJet; ++i) {
         float pt  = Jet_pt[i];
         float eta = Jet_eta[i];
         float phi = Jet_phi[i];

         if (pt <= 20.0) continue;
         if (std::fabs(eta) >= 2.5) continue;

         TLorentzVector j;
         j.SetPtEtaPhiM(pt, eta, phi, Jet_mass[i]);
         vec_jet.push_back(j);

         bool isDoubleB = (Jet_btagUParTAK4probbb[i] > 0.12f);

         if (isDoubleB) {
           vec_doublebjets.push_back(j);
         }
         else if (Jet_btagUParTAK4B[i] > 0.4648f) {
           vec_bjets.push_back(j);
         }
      }

      auto cmpPt = [](const TLorentzVector &a, const TLorentzVector &b) { return a.Pt() > b.Pt(); };

      std::sort(vec_ele.begin(),        vec_ele.end(),        cmpPt);
      std::sort(vec_muons.begin(),      vec_muons.end(),      cmpPt);
      std::sort(vec_jet.begin(),        vec_jet.end(),        cmpPt);
      std::sort(vec_bjets.begin(),      vec_bjets.end(),      cmpPt);
      std::sort(vec_doublebjets.begin(),vec_doublebjets.end(),cmpPt);

      std::sort(leptons.begin(), leptons.end(),
                [](const RecoLepton &a, const RecoLepton &b) { return a.p4.Pt() > b.p4.Pt(); });

      // PRE-SELECTION HISTOGRAMS
      h_nEle->Fill(vec_ele.size(),   weight);
      h_nMu ->Fill(vec_muons.size(), weight);
      h_nJet->Fill(vec_jet.size(),   weight);

      h_nBjet->Fill(vec_bjets.size(), weight);
      h_nDoubleB->Fill(vec_doublebjets.size(), weight);

      for (const auto &el : vec_ele) {
         h_ele_pt ->Fill(el.Pt(),  weight);
         h_ele_eta->Fill(el.Eta(), weight);
         h_ele_phi->Fill(el.Phi(), weight);
      }

      for (const auto &mu : vec_muons) {
         h_mu_pt ->Fill(mu.Pt(),  weight);
         h_mu_eta->Fill(mu.Eta(), weight);
         h_mu_phi->Fill(mu.Phi(), weight);
      }

      for (const auto &jet : vec_jet) {
         h_jet_pt ->Fill(jet.Pt(),  weight);
         h_jet_eta->Fill(jet.Eta(), weight);
         h_jet_phi->Fill(jet.Phi(), weight);
      }

      for (const auto &bj : vec_bjets) {
         h_bjet_pt ->Fill(bj.Pt(),  weight);
         h_bjet_eta->Fill(bj.Eta(), weight);
         h_bjet_phi->Fill(bj.Phi(), weight);
      }

      for (const auto &dbj : vec_doublebjets) {
         h_dbjet_pt ->Fill(dbj.Pt(),  weight);
         h_dbjet_eta->Fill(dbj.Eta(), weight);
         h_dbjet_phi->Fill(dbj.Phi(), weight);
      }

      {
         const int maxEle = std::min<int>(vec_ele.size(), MAXELE);
         for (int i = 0; i < maxEle; ++i) {
           h_e_pt_rank[i] ->Fill(vec_ele[i].Pt(),  weight);
           h_e_eta_rank[i]->Fill(vec_ele[i].Eta(), weight);
           h_e_phi_rank[i]->Fill(vec_ele[i].Phi(), weight);
         }
      }

      {
         const int maxMu = std::min<int>(vec_muons.size(), MAXMU);
         for (int i = 0; i < maxMu; ++i) {
           h_mu_pt_rank[i] ->Fill(vec_muons[i].Pt(),  weight);
           h_mu_eta_rank[i]->Fill(vec_muons[i].Eta(), weight);
           h_mu_phi_rank[i]->Fill(vec_muons[i].Phi(), weight);
         }
      }

      {
         const int maxJet = std::min<int>(vec_jet.size(), MAXJET);
         for (int i = 0; i < maxJet; ++i) {
           h_j_pt[i] ->Fill(vec_jet[i].Pt(),  weight);
           h_j_eta[i]->Fill(vec_jet[i].Eta(), weight);
           h_j_phi[i]->Fill(vec_jet[i].Phi(), weight);
         }
      }

      {
         const int maxB = std::min<int>(vec_bjets.size(), MAXBJET);
         for (int i = 0; i < maxB; ++i) {
           h_bj_pt[i] ->Fill(vec_bjets[i].Pt(),  weight);
           h_bj_eta[i]->Fill(vec_bjets[i].Eta(), weight);
           h_bj_phi[i]->Fill(vec_bjets[i].Phi(), weight);
         }
      }

      {
         const int maxDB = std::min<int>(vec_doublebjets.size(), MAXDBJET);
         for (int i = 0; i < maxDB; ++i) {
           h_dbj_pt[i] ->Fill(vec_doublebjets[i].Pt(),  weight);
           h_dbj_eta[i]->Fill(vec_doublebjets[i].Eta(), weight);
           h_dbj_phi[i]->Fill(vec_doublebjets[i].Phi(), weight);
         }
      }

      // CUT FLOW
      const int Nleptons = leptons.size();
      const int Nbjets   = vec_bjets.size();
      const int Ndoubleb = vec_doublebjets.size();

      if (Nleptons < 2) continue;
      ++nCut1;

      for (const auto &jet : vec_jet) {
         for (const auto &el : vec_ele)  h_dR_jet_ele_before->Fill(jet.DeltaR(el), weight);
         for (const auto &mu : vec_muons) h_dR_jet_mu_before ->Fill(jet.DeltaR(mu), weight);
      }

      const double dRclean = 0.4;
      for (const auto &jet : vec_jet) {
         bool keep = true;

         for (const auto &el : vec_ele) {
            if (jet.DeltaR(el) < dRclean) { keep = false; break; }
         }
         if (keep) {
            for (const auto &mu : vec_muons) {
               if (jet.DeltaR(mu) < dRclean) { keep = false; break; }
            }
         }
         if (keep) vec_cjet.push_back(jet);
      }

      h_nCJet->Fill(vec_cjet.size(), weight);

      for (const auto &jet : vec_cjet) {
         for (const auto &el : vec_ele)   h_dR_jet_ele_after->Fill(jet.DeltaR(el), weight);
         for (const auto &mu : vec_muons) h_dR_jet_mu_after ->Fill(jet.DeltaR(mu), weight);
      }

      const int Njets = vec_cjet.size();

      bool os_flavour_ok = false;
      RecoLepton &l1 = leptons[0];
      RecoLepton &l2 = leptons[1];

      if (l1.charge * l2.charge < 0) os_flavour_ok = true;

      if (!os_flavour_ok) continue;
      ++nCut2;

      if (Njets < 4) continue;
      ++nCut3;

      if (Nbjets < 2) continue;
      ++nCut4;

      if (Ndoubleb < 2) continue;
      ++nCut5;

      if (PuppiMET_pt < 40.0) continue;
      ++nCut6;

      // =====================================================================
      // REAL ANALYSIS (FINAL SELECTION): FILL ONLY (NO DECLARATIONS)
      // =====================================================================

      const TLorentzVector &db1 = vec_doublebjets[0];
      const TLorentzVector &db2 = vec_doublebjets[1];

      TLorentzVector HiggsCand = db1 + db2;

      h_Hdbb_mass->Fill(HiggsCand.M(), weight);
      h_Hdbb_pt  ->Fill(HiggsCand.Pt(),  weight);
      h_Hdbb_eta ->Fill(HiggsCand.Eta(), weight);

      h_MET_pt_final ->Fill(PuppiMET_pt,  weight);
      h_MET_phi_final->Fill(PuppiMET_phi, weight);

      h_dbj1_mass->Fill(db1.M(), weight);
      h_dbj2_mass->Fill(db2.M(), weight);

      h_dbj1_pt ->Fill(db1.Pt(),  weight);
      h_dbj1_eta->Fill(db1.Eta(), weight);

      h_dbj2_pt ->Fill(db2.Pt(),  weight);
      h_dbj2_eta->Fill(db2.Eta(), weight);

      double HT = 0.0;
      for (const auto &j : vec_jet) HT += j.Pt();
      h_HT->Fill(HT, weight);

      // --- NEW FINAL histograms (fixed) ---
      h_dR_dbj12->Fill(db1.DeltaR(db2), weight);

      h_dbj1_phi->Fill(db1.Phi(), weight);
      h_dbj2_phi->Fill(db2.Phi(), weight);

      h_lep1_pt ->Fill(l1.p4.Pt(),  weight);
      h_lep1_eta->Fill(l1.p4.Eta(), weight);
      h_lep1_phi->Fill(l1.p4.Phi(), weight);

      h_lep2_pt ->Fill(l2.p4.Pt(),  weight);
      h_lep2_eta->Fill(l2.p4.Eta(), weight);
      h_lep2_phi->Fill(l2.p4.Phi(), weight);

   } // end event loop

   // ========================================================================
   // 3. EVENT YIELD CALCULATION & WEIGHT (RECO ONLY)
   // ========================================================================

   auto fmt3 = [&](double x) {
     char buf[32];
     std::snprintf(buf, sizeof(buf), "%.3f", x);
     return std::string(buf);
   };

   std::cout << "\n=============================================================\n";
   std::cout << "                   EXPECTED EVENT YIELDS (RECO)              \n";
   std::cout << "=============================================================\n";
   std::cout << " sigma(ttH)           = " << sigma_pb          << " pb\n";
   std::cout << " BR(W->lnu)^2         = " << fmt3(Br_Wlnu * Br_Wlnu) << "\n";
   std::cout << " L                    = " << L_int             << " fb^-1\n";
   std::cout << " -------------------------------------------\n";
   std::cout << " Nexp                 = " << Nexp_int          << " events\n";
   std::cout << " Nstat (MC)           = " << fmt3(Nstat)       << " events\n";
   std::cout << " per-event weight w   = Nexp/Nstat = " << fmt3(w) << "\n";
   std::cout << "=============================================================\n\n";

   auto eff_decimal = [&](Long64_t n) {
     if (nentries == 0) return std::string("0.000");
     const double e = static_cast<double>(n) / static_cast<double>(nentries);
     char buf[16];
     std::snprintf(buf, sizeof(buf), "%.3f", e);
     return std::string(buf);
   };

   std::cout << "=============================================================\n";
   std::cout << "                       CUT FLOW TABLE (RECO)                 \n";
   std::cout << "=============================================================\n";
   std::cout << std::left << std::setw(35) << "Step / Requirement"
             << std::setw(15) << "Events"
             << std::setw(15) << "Weighted"
             << std::setw(15) << "Eff" << std::endl;
   std::cout << "-------------------------------------------------------------\n";

   std::cout << std::left << std::setw(35) << "Step 0) Raw events"
             << std::setw(15) << nentries
             << std::setw(15) << fmt3(static_cast<double>(nentries) * w)
             << std::setw(15) << "1.000" << std::endl;

   std::cout << std::left << std::setw(35) << "Step 1) N leptons >= 2"
             << std::setw(15) << nCut1
             << std::setw(15) << fmt3(static_cast<double>(nCut1) * w)
             << std::setw(15) << eff_decimal(nCut1) << std::endl;

   std::cout << std::left << std::setw(35) << "Step 2) OS lepton pair (e/μ)"
             << std::setw(15) << nCut2
             << std::setw(15) << fmt3(static_cast<double>(nCut2) * w)
             << std::setw(15) << eff_decimal(nCut2) << std::endl;

   std::cout << std::left << std::setw(35) << "Step 3) Clean N  jets >= 4"
             << std::setw(15) << nCut3
             << std::setw(15) << fmt3(static_cast<double>(nCut3) * w)
             << std::setw(15) << eff_decimal(nCut3) << std::endl;

   std::cout << std::left << std::setw(35) << "Step 4) N b-jets >= 2"
             << std::setw(15) << nCut4
             << std::setw(15) << fmt3(static_cast<double>(nCut4) * w)
             << std::setw(15) << eff_decimal(nCut4) << std::endl;

   std::cout << std::left << std::setw(35) << "Step 5) N double-b-jets >= 2"
             << std::setw(15) << nCut5
             << std::setw(15) << fmt3(static_cast<double>(nCut5) * w)
             << std::setw(15) << eff_decimal(nCut5) << std::endl;

   std::cout << std::left << std::setw(35) << "Step 6) MET >= 40 GeV"
             << std::setw(15) << nCut6
             << std::setw(15) << fmt3(static_cast<double>(nCut6) * w)
             << std::setw(15) << eff_decimal(nCut6) << std::endl;

   std::cout << "=============================================================\n";

   out->Write();
   out->Close();

   std::cout << "Analysis complete. Histograms saved to " << outFileName << std::endl;

}
