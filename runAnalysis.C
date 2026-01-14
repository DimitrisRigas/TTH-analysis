// -------------------------------------------------------------------
// Sample selector
// -------------------------------------------------------------------
enum SampleType {
  kSignal12,
  kSignal25,
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
  // NEW single-top / tt semi-leptonic / tW-like
  kTBbarQtoLNu,
  kTBbarQto2Q,
  kTTtoLNu2Q,
  kTbarWplustoNu2Q,
  kTbarWplusto4Q
};

void runAnalysis(SampleType sample = kSignal12)
{
  TString infile;

  switch (sample) {
    case kSignal12:
      infile = "TTH12Gev.root";
      std::cout << "Running on SIGNAL m=12: " << infile << std::endl;
      break;

    case kSignal25:
      infile = "TTH25Gev.root";
      std::cout << "Running on SIGNAL m=25: " << infile << std::endl;
      break;

    case kSignal60:
      infile = "TTH60Gev.root";
      std::cout << "Running on SIGNAL m=60: " << infile << std::endl;
      break;

    // ----- backgrounds -----
    case kTTbar:
      infile = "TTto2L2Nu.root";
      std::cout << "Running on TTbar background: " << infile << std::endl;
      break;

    case kDYee:
      infile = "DYto2E4JETS.root";
      std::cout << "Running on DY → ee background: " << infile << std::endl;
      break;

    case kDYmumu:
      infile = "DYto2MU4JETS.root";
      std::cout << "Running on DY → μμ background: " << infile << std::endl;
      break;

    case kVBF_Hbb:
      infile = "VBFHto2b.root";
      std::cout << "Running on VBF H→bb sample: " << infile << std::endl;
      break;

    case kGGH_Hbb:
      infile = "glugluHtobb.root";
      std::cout << "Running on ggH H→bb sample: " << infile << std::endl;
      break;

    case kTTH_Hbb:
      infile = "TTHHtobb.root";
      std::cout << "Running on ttH H→bb sample: " << infile << std::endl;
      break;

    // ----- diboson backgrounds -----
    case kWW:
      infile = "WW.root";
      std::cout << "Running on WW background: " << infile << std::endl;
      break;

    case kWZ:
      infile = "WZ.root";
      std::cout << "Running on WZ background: " << infile << std::endl;
      break;

    case kZZ:
      infile = "ZZ.root";
      std::cout << "Running on ZZ background: " << infile << std::endl;
      break;

    // ----- NEW single-top / tt semi-leptonic / tW-like backgrounds -----
    case kTBbarQtoLNu:
      infile = "TBbarQtoLNu.root";
      std::cout << "Running on TBbarQtoLNu background: " << infile << std::endl;
      break;

    case kTBbarQto2Q:
      infile = "TBbarQto2Q.root";
      std::cout << "Running on TBbarQto2Q background: " << infile << std::endl;
      break;

    case kTTtoLNu2Q:
      infile = "TTtoLNu2Q.root";
      std::cout << "Running on TTtoLNu2Q background: " << infile << std::endl;
      break;

    case kTbarWplustoNu2Q:
      infile = "TbarWplustoNu2Q.root";
      std::cout << "Running on TbarWplustoNu2Q background: " << infile << std::endl;
      break;

    case kTbarWplusto4Q:
      infile = "TbarWplusto4Q.root";
      std::cout << "Running on TbarWplusto4Q background: " << infile << std::endl;
      break;

    default:
      std::cerr << "ERROR: Unknown sample type!" << std::endl;
      return;
  }

  TFile *f = TFile::Open(infile);
  if (!f || f->IsZombie()) {
    std::cerr << "ERROR: Cannot open file " << infile << std::endl;
    return;
  }

  TTree *tree = (TTree*)f->Get("Events");
  if (!tree) {
    std::cerr << "ERROR: Tree 'Events' not found in file!" << std::endl;
    f->Close();
    return;
  }

  MyClass analysis(tree);
  analysis.Loop();

  std::cout << "Analysis finished for file: " << infile << std::endl;

  f->Close();
}
