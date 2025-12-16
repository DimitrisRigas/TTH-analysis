// DO NOT include MyClass.h here — MyClass is already compiled by ROOT

void runAnalysis(bool runSignal = true)
{
    // -------------------------------------------------------------------
    // Choose input file based on flag
    // -------------------------------------------------------------------
    TString infile;

    if (runSignal) {
        infile = "TTH12Gev.root";   // <-- SIGNAL FILE
        std::cout << "Running on SIGNAL sample: " << infile << std::endl;
    } else {
        infile = "TTto2L2Nu.root";   // <-- BACKGROUND FILE
        std::cout << "Running on BACKGROUND sample: " << infile << std::endl;
    }

    // -------------------------------------------------------------------
    // DO NOT compile MyClass here!
    // MyClass must already be compiled by .L MyClass.C++
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
