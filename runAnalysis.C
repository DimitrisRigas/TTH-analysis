// -------------------------------------------------------------------
// Sample selector
// -------------------------------------------------------------------
enum SampleType {
    kSignal,
    kTTbar,
    kDYee,
    kDYmumu
};

void runAnalysis(SampleType sample = kSignal)
{
    // -------------------------------------------------------------------
    // Choose input file based on sample
    // -------------------------------------------------------------------
    TString infile;

    switch (sample) {
        case kSignal:
            infile = "TTH12Gev.root";
            std::cout << "Running on SIGNAL sample: " << infile << std::endl;
            break;

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

        default:
            std::cerr << "ERROR: Unknown sample type!" << std::endl;
            return;
    }

    // -------------------------------------------------------------------
    // DO NOT compile MyClass here!
    // MyClass must already be compiled by:
    //   .L MyClass.C++
    // -------------------------------------------------------------------

    // -------------------------------------------------------------------
    // Open input file
    // -------------------------------------------------------------------
    TFile *f = TFile::Open(infile);
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: Cannot open file " << infile << std::endl;
        return;
    }

    // -------------------------------------------------------------------
    // Grab the TTree (NanoAOD uses "Events")
    // -------------------------------------------------------------------
    TTree *tree = (TTree*)f->Get("Events");
    if (!tree) {
        std::cerr << "ERROR: Tree 'Events' not found in file!" << std::endl;
        return;
    }

    // -------------------------------------------------------------------
    // Run the analysis
    // -------------------------------------------------------------------
    MyClass analysis(tree);
    analysis.Loop();

    std::cout << "Analysis finished for file: " << infile << std::endl;
}
