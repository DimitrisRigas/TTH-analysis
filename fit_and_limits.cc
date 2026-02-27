#include <iostream>
#include "TCanvas.h"
#include "TH1F.h"
#include "TFile.h"
#include "TGraph.h"
#include "TSpline.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooArgList.h"
#include "RooHist.h"
#include "RooBinning.h"
#include "RooMCStudy.h"
#include "RooStats/SimpleInterval.h"
#include "RooStats/BayesianCalculator.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/ProposalHelper.h"
#include "TH2.h"
#include "TCutG.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TSpline.h"
#include "TMatrixDSym.h"


using namespace RooFit;
using namespace RooStats;

// Bayesian Markov Chain Monte Carlo (MCMC) analysis
struct BayesianMCMCOptions
{
  double CL_95 = 0.95;
  double CL_68 = 0.68;
 
  // MCMCInterval::IntervalType intervalType = MCMCInterval::Upper;
  // type of interval (0 is shortest, 1 central, 2 upper limit)
  double maxPOI = -999;      // force different values of POI (Parameters of Interest) for doing the scan (default is given value)
  double minPOI = -999;
  int numIters = 10000;    // number of iterations
  int numBurnInSteps = 1000; // number of burn in steps to be ignored
};



BayesianMCMCOptions optMCMC;

void fit_and_limits()
{

  
  //  bool obs("true");
   bool exp("true");
  TString lep("0");
  TString sr("2");
  TString boost("grad");
  std::vector<int> mass_points = {60};
  
  
  std::vector<double> observed_limits;
  std::vector<double> expected_limits;
  std::vector<double> sigma_1_up, sigma_1_down;
  std::vector<double> sigma_2_up, sigma_2_down;
  std::vector<double> observed_br_limits;
  std::vector<double> expected_br_limits;
  std::vector<double> toy_br_limits;
  std::vector<double> sigma_1_br_up, sigma_1_br_down;
  std::vector<double> sigma_2_br_up, sigma_2_br_down;
  
  const int Ntoys = 1001;

  for(size_t i = 0; i < mass_points.size(); ++i)
    {
      int mass = mass_points[i];
      TString amass = TString::Format("%d", mass);
      //input files
      TFile *sigm = TFile::Open("sig"+amass+"_"+lep+"lep.root");
      TFile *ftt1 = TFile::Open("hist_TT_filt1_"+lep+"lep.root");
      TFile *ftt4 = TFile::Open("hist_TT_filt4_"+lep+"lep.root");
      TFile *ftt5 = TFile::Open("hist_TT_filt5_"+lep+"lep.root");
      TFile *fznn = TFile::Open("hist_znn.root");
      TFile *fw = TFile::Open("hist_w.root");
      TFile *fqcd = TFile::Open("hist_QCDpl.root");
      TFile *foth = TFile::Open("hist_others_"+lep+"lep.root");

      //define bin edges for variable binning of bdt score
      int b=4;
      double var_bins[b+1] ;
   

      //bins

	
	var_bins[0] = -1;
	var_bins[1] = -0.4;
	var_bins[2] = -0.08;
	var_bins[3] = 0.52;
	var_bins[4] = 1;
 

      //input bdt histos
      TH1 *h_m =(TH1*)sigm->Get("sr"+sr+"_"+boost+"_bdt");
      TH1 *h_tt1 =(TH1*)ftt1->Get("sr"+sr+"_"+boost+"_bdt");
      TH1 *h_tt4 =(TH1*)ftt4->Get("sr"+sr+"_"+boost+"_bdt");
      TH1 *h_tt5 =(TH1*)ftt5->Get("sr"+sr+"_"+boost+"_bdt");
      TH1 *h_znn =(TH1*)fznn->Get("sr"+sr+"_"+boost+"_bdt");
      TH1 *h_w =(TH1*)fw->Get("sr"+sr+"_"+boost+"_bdt");
      TH1 *h_qcd =(TH1*)fqcd->Get("sr"+sr+"_"+boost+"_bdt");
      TH1 *h_oth =(TH1*)foth->Get("sr"+sr+"_"+boost+"_bdt");
      //rebin input histos using the variable binning
      TH1*h_m_r= h_m->Rebin(b, "h_m_r", var_bins);
      TH1*h_znn_r= h_znn->Rebin(b, "h_znn_r", var_bins);
      TH1*h_w_r= h_w->Rebin( b, "h_w_r", var_bins);
      TH1*h_qcd_r= h_qcd->Rebin(b, "h_qcd_r", var_bins);
      TH1*h_oth_r= h_oth->Rebin(b, "h_oth_r", var_bins);
      TH1*h_tt1_r= h_tt1->Rebin(b, "h_tt1_r", var_bins);
      TH1*h_tt4_r= h_tt4->Rebin(b, "h_tt4_r", var_bins);
      TH1*h_tt5_r= h_tt5->Rebin(b, "h_tt5_r", var_bins);

      //observable (rebinned bdt score)
      RooRealVar output_BDT("output_BDT","BDT score",-1,1);
      RooBinning customBinning(4, var_bins); // bins with specified edges
      output_BDT.setBinning(customBinning, "customBinning");
      
     //define roodatahists from hist
     RooDataHist sig("sig","sig",output_BDT,h_m_r);
     RooDataHist tt1("tt1","tt1",output_BDT,h_tt1_r);
     RooDataHist tt4("tt4","tt4",output_BDT,h_tt4_r);
     RooDataHist tt5("tt5","tt5",output_BDT,h_tt5_r);
     RooDataHist w("w","w",output_BDT,h_w_r);
     RooDataHist qcd("qcd","qcd",output_BDT,h_qcd_r);
     RooDataHist znn("znn","znn",output_BDT,h_znn_r);
     //RooDataHist dy("dy","dy",output_BDT,h_dy_r);
     RooDataHist oth("oth","oth",output_BDT,h_oth_r);
    //define pdfs 
     RooHistPdf znn_pdf("znn_pdf","znn_pdf",output_BDT,znn);
     RooHistPdf w_pdf("w_pdf","w_pdf",output_BDT,w);
     RooHistPdf qcd_pdf("qcd_pdf","qcd_pdf",output_BDT,qcd);
     RooHistPdf tt1_pdf("tt1_pdf","tt1_pdf",output_BDT,tt1);
     RooHistPdf tt4_pdf("tt4_pdf","tt4_pdf",output_BDT,tt4);
     RooHistPdf tt5_pdf("tt5_pdf","tt5_pdf",output_BDT,tt5);
     //RooHistPdf dy_pdf("dy_pdf","dy_pdf",output_BDT,dy);
     RooHistPdf oth_pdf("oth_pdf","oth_pdf",output_BDT,oth);
     RooHistPdf sig_pdf("sig_pdf","sig_pdf",output_BDT,sig);
    
     
     //define nominal yields for each process
     double Nsig=h_m->Integral();
     double Ntt1=h_tt1->Integral();
     double Ntt4=h_tt4->Integral();
     double Ntt5=h_tt5->Integral();
     double Nw=h_w->Integral();
     double Nqcd=h_qcd->Integral();
     double Nznn=h_znn->Integral();
     double Noth=h_oth->Integral();
     double Nbkg=Ntt1+Ntt4+Ntt5+Nznn+Nw+Nqcd+Noth;
     double Nsb=Nsig+Nbkg;

     //define variables and  gaussian constraints
     RooRealVar Nexp_sig("Nexp_sig","Expected number of signal events",Nsig,-5*Nsig,5*Nsig );
    RooRealVar Nexp_sig_th("Nexp_sig_th","Expected number of signal events",Nsig);
    //TTBAR
  
    //tt5
    RooRealVar Nexp_tt5("Nexp_tt5", "Expected number of TT5 bkg", Ntt5,-1.1*Ntt5,3.2*Ntt5);
    //RooRealVar Nexp_tt5("Nexp_tt5", "Expected number of TT5 bkg", Ntt5);
    RooRealVar Nexp_tt5_th("Nexp_tt5_th", "Expected number of TT  bkg", Ntt5);
    RooGaussian tt5_constraint("tt5_constraint", "Gaussian syst  constraint on ttbar yield",Nexp_tt5, RooConst(Ntt5), RooConst(0.7*Ntt5));
    
    //w
    RooRealVar Nexp_w("Nexp_w", "Expected number of W bkg", Nw,0.*Nw,1.9*Nw);
    RooRealVar Nexp_w_th("Nexp_w_th", "Expected number of  W bkg", Nw);
    RooGaussian w_constraint("w_constraint", "Gaussian   constraint on W yield", Nexp_w, RooConst(Nw), RooConst(0.3*Nw));
    //qcd
    RooRealVar Nexp_qcd("Nexp_qcd", "Expected number of QCD bkg", Nqcd,-0.5*Nqcd,2.5*Nqcd);
    RooRealVar Nexp_qcd_th("Nexp_qcd_th", "Expected number of  QCD bkg", Nqcd);
    RooGaussian qcd_constraint("qcd_constraint", "Gaussian constraint on qcd yield", Nexp_qcd, RooConst(Nqcd), RooConst(0.5*Nqcd));
    
    //z->vv
    RooRealVar Nexp_znn("Nexp_znn", "Expected number of Zvv bkg", Nznn,0.*Nznn,1.9*Nznn);
    RooRealVar Nexp_znn_th("Nexp_znn_th", "Expected number of Zvv bkg", Nznn);
    RooGaussian znn_constraint("znn_constraint", "Gaussian  constraint on Zvv yield", Nexp_znn, RooConst(Nznn), RooConst(0.3*Nznn));

    //other bkg
    RooRealVar Nexp_oth("Nexp_oth", "Expected number of other bkg", Noth,-0.5*Noth,2.5*Noth);
    RooRealVar Nexp_oth_th("Nexp_oth_th", "Expected number of oth bkg", Noth);
    RooGaussian oth_constraint("oth_constraint", "Gaussian  constraint on others yield", Nexp_oth, RooConst(Noth), RooConst(0.5*Noth));
    //init 
    Nexp_sig.setVal(Nexp_sig_th.getVal());
    Nexp_tt1.setVal(Nexp_tt1_th.getVal());
    Nexp_tt4.setVal(Nexp_tt4_th.getVal());
    Nexp_tt5.setVal(Nexp_tt5_th.getVal());
    Nexp_znn.setVal(Nexp_znn_th.getVal());
    Nexp_w.setVal(Nexp_w_th.getVal());
    Nexp_qcd.setVal(Nexp_qcd_th.getVal());
    Nexp_oth.setVal(Nexp_oth_th.getVal());
    //define model0 (bkg only)
    RooAddPdf model_0("model_0", "Background", RooArgList(tt1_pdf,tt4_pdf,tt5_pdf,w_pdf,znn_pdf,qcd_pdf,oth_pdf), RooArgList(Nexp_tt1,Nexp_tt4,Nexp_tt5,Nexp_w,Nexp_znn,Nexp_qcd,Nexp_oth));
    //define model0 (s+b)
    RooAddPdf model_1("model_1", "Signal + Background", RooArgList(sig_pdf,tt1_pdf,tt4_pdf,tt5_pdf,w_pdf,znn_pdf,qcd_pdf,oth_pdf), RooArgList(Nexp_sig,Nexp_tt1,Nexp_tt4,Nexp_tt5,Nexp_w,Nexp_znn,Nexp_qcd,Nexp_oth));
    //model0 and model1 with costraints  
    RooProdPdf total_model_0("total_model_0", " Background model with constraints", RooArgList(model_0,tt1_constraint,tt4_constraint,tt5_constraint,w_constraint,qcd_constraint,znn_constraint,oth_constraint));
    RooProdPdf total_model_1("total_model_1", "Signal + Background model with constraints",RooArgList(model_1,tt1_constraint,tt4_constraint,tt5_constraint,w_constraint,qcd_constraint,znn_constraint,oth_constraint));
    //convert model0 into roodatahist 
     TH1* modelHist = model_0.createHistogram("totalModelHist", output_BDT, RooFit::Binning(customBinning));
     for (int i = 1; i <= modelHist->GetNbinsX(); ++i) {
       double binWidth = modelHist->GetBinWidth(i);
       double binContent = modelHist->GetBinContent(i);
       modelHist->SetBinContent(i, binContent * binWidth); // Adjust by bin width
     }
     //dataHist_B has  number of total events=Nbkg and has exactly the same per bin yields with the initial bdt distribution  
     RooDataHist dataHist_B("data_hist", "Data Histogram from Model", output_BDT, modelHist );
     //fit model1 into dataB
     RooFitResult *fit_model_1_B_asimov = total_model_1.fitTo(dataHist_B,Save(),Extended(kTRUE),RooFit::MaxCalls(10000), Strategy(2),RooFit::Optimize(kTRUE),SumW2Error(true));
     //compute nll
    RooAbsReal* nll = total_model_1.createNLL(dataHist_B,RooFit::Extended(kTRUE));
    RooMinimizer minim(*nll);

    // Run the minimization
    minim.minimize("Minuit2", "migrad");
    minim.minos();
    fit_model_1_B_asimov->Print("v");
    //return to init values
     Nexp_sig.setVal(Nexp_sig_th.getVal());
    Nexp_tt1.setVal(Nexp_tt1_th.getVal());
    Nexp_tt4.setVal(Nexp_tt4_th.getVal());
    Nexp_tt5.setVal(Nexp_tt5_th.getVal());
    Nexp_znn.setVal(Nexp_znn_th.getVal());
    Nexp_w.setVal(Nexp_w_th.getVal());
    Nexp_qcd.setVal(Nexp_qcd_th.getVal());
    Nexp_oth.setVal(Nexp_oth_th.getVal());
    
    //data_B  has Nbkg events but has PER BIN statistical  fluctuations (so will be used in the MCMC chain)
    RooDataHist  *data_B =model_0.generateBinned(output_BDT,Nbkg);
    RooFitResult *fit_model_1_B = total_model_1.fitTo(*data_B,Save(),Extended(kTRUE),RooFit::MaxCalls(10000), Strategy(2),RooFit::Optimize(kTRUE));

    RooNLLVar* nlls = (RooNLLVar*) total_model_1.createNLL( *data_B, RooFit::Extended(kTRUE));

    RooPlot* frameLS = Nexp_sig.frame();
    //nll plot
    // Plot the likelihood surface (profiling)
    nlls->plotOn(frameLS, RooFit::ShiftToZero(), RooFit::PrintEvalErrors(-1));

    // Draw the plot
    TCanvas* canvas = new TCanvas("canvas", "Likelihood Surface", 600, 600);
    frameLS->SetTitle("Profile Likelihood for Nexp_signal");
    frameLS->Draw();
    TLatex *latex = new TLatex();
    latex->SetNDC();            // Use Normalized Device Coordinates
    latex->SetTextSize(0.03);   // Adjust text size
    latex->SetTextFont(42);    
    latex->DrawLatex(0.78, 0.96, "M_{a}="+amassamass+" GeV"); 
  
    TCanvas* cp = new TCanvas("cp", "Fit and pulls of model 1 on B-only data", 800, 800);
    TPad *t1 = new TPad("t1","t1", 0.0, 0.25, 1.0, 1.0);
    // gPad->SetLogy();
    t1->SetFillColor(0);
    t1->SetBorderMode(0);
    t1->SetBorderSize(2);
    //  t1->SetTickx(1);
    // t1->SetTicky(1);
    t1->SetLeftMargin(0.10);
    t1->SetRightMargin(0.05);
    t1->SetTopMargin(0.05);
    t1->SetBottomMargin(0.05); 

    t1->SetFrameFillStyle(0);
    t1->SetFrameBorderMode(0);
    t1->SetFrameFillStyle(0);
    t1->SetFrameBorderMode(0);
    t1->Draw();
    t1->cd();
    RooPlot *frame_B1 = output_BDT.frame();
    dataHist_B.plotOn(frame_B1);
    total_model_1.plotOn(frame_B1);
    total_model_1.plotOn(frame_B1, LineColor(kMagenta-10),Components("oth_pdf"));
  
    total_model_1.plotOn(frame_B1, LineColor(kGreen-9),Components("tt1_pdf"));
    total_model_1.plotOn(frame_B1, LineColor(kGreen),Components("tt4_pdf"));
    total_model_1.plotOn(frame_B1, LineColor(kGreen+2),Components("tt5_pdf"));
    total_model_1.plotOn(frame_B1, LineColor(wC),Components("dy_pdf"));
  
    total_model_1.plotOn(frame_B1, LineColor(kBlack),LineStyle(1),Components("sig_pdf"));
   
  
 
    RooCurve* curve_2 = dynamic_cast<RooCurve*>(frame_B1->getObject(1));

    frame_B1->GetXaxis()->SetTitleSize(0);
    frame_B1->GetXaxis()->SetTitle(0);
    frame_B1->GetXaxis()->SetTitleOffset(0.4);
    frame_B1->SetMaximum(1000000);
    frame_B1->SetMinimum(0.01);
    gPad->SetLogy();
    frame_B1->Draw();
    t1->Update();
    TLegend *leg =  new TLegend(0.30,0.77,0.93,0.96, "NDC");
    leg->SetHeader("");
    leg->SetNColumns(3);   
    leg->SetBorderSize(0);
    leg->SetTextFont(42);   leg->SetTextSize(0.03);
    leg->SetLineColor(0);   leg->SetLineStyle(1);   leg->SetLineWidth(1);
    leg->SetFillColor(0);   leg->SetFillStyle(0);

    leg->AddEntry(data_B, "Bgk only Data", "p"); 
    leg->AddEntry(curve_2, "Fit model_1", "l")->SetLineColor(kBlue); // "l" for line
    leg->AddEntry(data_B, "other Bkgs", "l")->SetLineColor(kMagenta-10);
    leg->AddEntry(data_B,"t #bar{t}+light","l")->SetLineColor(kGreen-9);
    leg->AddEntry(data_B,"t #bar{t} +c #bar{c}","l")->SetLineColor(kGreen);
    leg->AddEntry(data_B,"t #bar{t} +b #bar{b}","l")->SetLineColor(kGreen+2);
    leg->AddEntry(data_B,"DY","l") ->SetLineColor(wC);
    leg->AddEntry(data_B,"Zh( M_{a}="+amass+")","l") ->SetLineColor(kBlack);
  
   
    leg->Draw("SAME");
    t1->Update();
    cp->cd();
    TPad *t2 = new TPad("t2", "t2",0.0,0.0, 1.0,0.25);

    t2->SetFillColor(0);
    t2->SetBorderMode(0);
    t2->SetBorderSize(2);
    t2->SetLeftMargin(0.10);     
    t2->SetRightMargin(0.05);    
    t2->SetTopMargin(0.02);      
    t2->SetBottomMargin(0.40);   
    t2->SetFrameFillStyle(0);
    t2->SetFrameBorderMode(0);
    t2->Draw();
    t2->cd();

    RooPlot* frame01 =output_BDT.frame();
    dataHist_B.plotOn(frame01);
    total_model_1.plotOn(frame01);
    RooHist* pulls01 = frame01->pullHist(); 
    // Create a new frame for residuals
    RooPlot* residualFrame01 =output_BDT.frame();
    residualFrame01->addPlotable(pulls01, "P");  // Plot the pulls
    // Customize the residual frame
    residualFrame01->GetYaxis()->SetTitle("Pulls"); 
    residualFrame01->GetYaxis()->SetTitleSize(0.12); 
    residualFrame01->GetYaxis()->SetTitleOffset(0.4); 
    residualFrame01->GetYaxis()->SetLabelSize(0.10);  
    residualFrame01->GetXaxis()->SetLabelSize(0.10);  
    residualFrame01->GetXaxis()->SetTitleSize(0.15);
    residualFrame01->GetXaxis()->SetTitleOffset(0.7);
    // Draw the pulls frame
    residualFrame01->Draw();
    t2->Update();
    // Draw a horizontal line at y = 0
    TLine *l1 = new TLine(residualFrame01->GetXaxis()->GetXmin(), 0, residualFrame01->GetXaxis()->GetXmax(), 0);
    l1->SetLineColor(kRed);
    l1->SetLineWidth(2);
    l1->Draw("same");

    // Update and redraw the second pad
    t2->Update();
    t2->RedrawAxis();
    t2->Update();
    // Update the main canvas
    cp->cd();
    cp->RedrawAxis();
    cp->Update();
    cp->SaveAs("toy"+lep+"l"+sr+".pdf");

    
    Nexp_sig.setVal(Nexp_sig_th.getVal());
    Nexp_tt1.setVal(Nexp_tt1_th.getVal());
    Nexp_tt4.setVal(Nexp_tt4_th.getVal());
    Nexp_tt5.setVal(Nexp_tt5_th.getVal());
    Nexp_znn.setVal(Nexp_znn_th.getVal());
    Nexp_w.setVal(Nexp_w_th.getVal());
    Nexp_qcd.setVal(Nexp_qcd_th.getVal());
    Nexp_oth.setVal(Nexp_oth_th.getVal());
    
    //TOWARDS EXTRACTING LIMITS:
    
    // Create RooWorkspace to hold models and parameters and ModelConfig
    RooWorkspace wr("wr");
    // Import models into the workspace
    wr.import(total_model_1, RooFit::RenameConflictNodes("total_model_1"));

    // Set up ModelConfig after importing models
    ModelConfig mc("ModelConfig", &wr);
    mc.SetPdf(*wr.pdf("total_model_1"));
    mc.SetParametersOfInterest(RooArgSet(Nexp_sig));
    mc.SetNuisanceParameters(RooArgSet(Nexp_tt1,Nexp_tt4,Nexp_tt5,Nexp_w,Nexp_znn,Nexp_qcd,Nexp_oth));
    wr.import(mc);
    double denominator =Nsig;
     
      
        //===================================================//
      // E X P E C T E D  L I M I T  C A L C U L A T I O N //
      //===================================================//
       if (exp){
	cout << "\n>>>>>> START with expected limit calculation for mass " << mass << endl;
	cout << "\n>>>>>> Nominal value : NSignal = " << Nexp_sig.getVal() << "\n" << endl;
	
	// Return to the initial values of Nexpected for the next fit
	Nexp_sig.setVal(Nexp_sig_th.getVal());
	Nexp_tt1.setVal(Nexp_tt1_th.getVal());
	Nexp_tt4.setVal(Nexp_tt4_th.getVal());
	Nexp_tt5.setVal(Nexp_tt5_th.getVal());
	Nexp_znn.setVal(Nexp_znn_th.getVal());
	Nexp_w.setVal(Nexp_w_th.getVal());
	Nexp_qcd.setVal(Nexp_qcd_th.getVal());
	Nexp_oth.setVal(Nexp_oth_th.getVal());
	ProposalHelper ph_exp;
	ph_exp.SetVariables((RooArgSet&)fit_model_1_B_asimov->floatParsFinal());
	ph_exp.SetCovMatrix(fit_model_1_B_asimov->covarianceMatrix());
	ph_exp.SetUpdateProposalParameters(kTRUE); 
	ph_exp.SetCacheSize(100);
	ProposalFunction* pf_exp = ph_exp.GetProposalFunction();
	
	std::vector<double> toy_limits;
        
	for (int i = 0; i < Ntoys; ++i)
	  {
	    
	    RooDataHist* toyData =model_0.generateBinned(output_BDT, Nbkg);
	    MCMCCalculator toy_mcmc_calculator(*toyData, mc);
	    toy_mcmc_calculator.SetProposalFunction(*pf_exp);
	    toy_mcmc_calculator.SetConfidenceLevel(optMCMC.CL_95); 
	    toy_mcmc_calculator.SetNumIters(optMCMC.numIters);
	    toy_mcmc_calculator.SetNumBurnInSteps(optMCMC.numBurnInSteps);
	    toy_mcmc_calculator.SetLeftSideTailFraction(0.0);
	    
	    MCMCInterval* toy_interval = toy_mcmc_calculator.GetInterval();
	    double Nexp_Signal_exp_upper = toy_interval->UpperLimit(Nexp_sig);
	    toy_limits.push_back(Nexp_Signal_exp_upper);
	    
	    double BR_exp = Nexp_Signal_exp_upper / denominator;
	    toy_br_limits.push_back(BR_exp);
	    
	    std::cout << "\n>>>> RESULT : " << optMCMC.CL_95 * 100 << "% interval is [" << toy_interval->LowerLimit(Nexp_sig) <<  ", " << toy_interval->UpperLimit(Nexp_sig) << "] \n" << std::endl;
	    std::cout << "\n>>>> Toy " << i + 1 << " / " << Ntoys << " | Upper limit: " << Nexp_Signal_exp_upper << " | BR: " << BR_exp << "\n" << std::endl;
	    // toy_limits.clear();
	    // toy_br_limits.clear();
	  }    
	
	// Calculate median and sigma bands from toy limits
	std::sort(toy_br_limits.begin(), toy_br_limits.end());
	double expected_95_CL = toy_br_limits[Ntoys / 2];
	expected_br_limits.push_back(expected_95_CL);
	sigma_1_up.push_back(toy_br_limits[Ntoys*0.84] );
	sigma_1_down.push_back( toy_br_limits[Ntoys* 0.16]);
	sigma_2_up.push_back(toy_br_limits[Ntoys* 0.975]  );
	sigma_2_down.push_back( toy_br_limits[Ntoys*0.025]);  }	 
       										
    }
  
  // Print the vector
  if(exp){
    std::cout << "exp Upper limits are  expected_br_limits {";
    for (size_t i = 0; i <expected_br_limits.size(); ++i) {
      std::cout << expected_br_limits[i];
      if (i != expected_br_limits.size() - 1) {
	std::cout << ", "; 
      }
    }
    std::cout << "}" << std::endl;
    // Print the vector
    std::cout << "Upper limits 1S up {";
    for (size_t i = 0; i <sigma_1_up.size(); ++i) {
      std::cout << sigma_1_up[i];
	if (i !=sigma_1_up.size() - 1) {
	  std::cout << ", "; 
        }
    }
    std::cout << "}" << std::endl;
    
    
    // Print the vector
    std::cout << "Upper limits 1S down {";
    for (size_t i = 0; i <sigma_1_down.size(); ++i) {
      std::cout << sigma_1_down[i];
      if (i !=sigma_1_down.size() - 1) {
	std::cout << ", "; 
        }
    }
    std::cout << "}" << std::endl;
    // Print the vector
    std::cout << "Upper limits 2S up {";
    for (size_t i = 0; i <sigma_2_up.size(); ++i) {
      std::cout << sigma_2_up[i];
      if (i !=sigma_2_up.size() - 1) {
	std::cout << ", "; 
      }
    }
    std::cout << "}" << std::endl;

    
    // Print the vector
    std::cout << "Upper limits 2S down {";
    for (size_t i = 0; i <sigma_2_down.size(); ++i) {
      std::cout << sigma_2_down[i];
      if (i !=sigma_2_down.size() - 1) {
	std::cout << ", "; 
        }
    }
    std::cout << "}" << std::endl;}
   //===================================================//
      // O B S E R V E D  L I M I T  C A L C U L A T I O N //
      //===================================================//
      /*  if(obs){
       TRandom3 *rand = new TRandom3();
       double Nobs = rand->Poisson(Nbkg);
       RooDataSet *data_obs =model_0.generate(output_BDT,Nobs);
       // Generate a toyMC sample from a poissonian distribution
        TRandom3 *rand = new TRandom3(); 
        double Nobs = rand->Poisson(Nbkg);
        RooDataHist *dataHist_obs = total_model_0.generateBinned(output_BDT,Nobs);
 
      RooFitResult *fit_model_obs = total_model_1.fitTo( *dataHist_obs , RooFit::Save(), Extended(kTRUE), Strategy(2));

      Nexp_sig.setVal(Nexp_sig_th.getVal());
      Nexp_tt1.setVal(Nexp_tt1_th.getVal());
      Nexp_tt4.setVal(Nexp_tt4_th.getVal());
      Nexp_tt5.setVal(Nexp_tt5_th.getVal());
      Nexp_znn.setVal(Nexp_znn_th.getVal());
      Nexp_w.setVal(Nexp_w_th.getVal());
      Nexp_qcd.setVal(Nexp_qcd_th.getVal());
      Nexp_oth.setVal(Nexp_oth_th.getVal());
      
      
      ProposalHelper ph_obs;
      ph_obs.SetVariables((RooArgSet&)fit_model_obs->floatParsFinal());
      ph_obs.SetCovMatrix(fit_model_obs->covarianceMatrix());
      ph_obs.SetUpdateProposalParameters(kTRUE); 
      ph_obs.SetCacheSize(100);
      ProposalFunction* pf_obs = ph_obs.GetProposalFunction();

      cout << "\n>>>>>> START with observed limit calculation for mass " << mass << endl;
      cout << "\n>>>>>> Nominal value : NSignal = " << Nexp_sig.getVal() << "\n" << endl;
      cout << "\n>>>>>> denominator is : = " << denominator << "\n" << endl;
      MCMCCalculator mcmc_calculator(*data_obs, mc);
      mcmc_calculator.SetProposalFunction(*pf_obs);
      mcmc_calculator.SetConfidenceLevel(optMCMC.CL_95); 
      mcmc_calculator.SetNumIters(optMCMC.numIters);
      mcmc_calculator.SetNumBurnInSteps(optMCMC.numBurnInSteps);
      mcmc_calculator.SetLeftSideTailFraction(0.0);

      MCMCInterval* interval       = mcmc_calculator.GetInterval();
      double Nexp_Signal_upper_obs = interval->UpperLimit(Nexp_sig);
      double BR_obs                = Nexp_Signal_upper_obs / denominator;
      
      observed_limits.push_back(Nexp_Signal_upper_obs);
      observed_br_limits.push_back(BR_obs);}*/
  // Print the vector
  /*  if(obs){
    std::cout << "obs Upper limits are  observed_br_limits {";
    for (size_t i = 0; i <observed_br_limits.size(); ++i) {
      std::cout << observed_br_limits[i];
      if (i !=observed_br_limits.size() - 1) {
	std::cout << ", "; 
        }
    }
    
    std::cout << "}" << std::endl;}*/
}


