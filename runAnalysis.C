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
  kTTH_Hbb
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

    // ----- keep backgrounds exactly as you already have -----
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
    return;
  }

  MyClass analysis(tree);
  analysis.Loop();

  std::cout << "Analysis finished for file: " << infile << std::endl;
}
