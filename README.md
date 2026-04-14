# Higgs Boson (ttH) Event Analysis with Machine Learning

## Overview

This project focuses on Higgs boson production in association with a top–antitop quark pair (ttH) using CMS data from CERN. It combines event analysis, machine learning classification, visualization, and statistical fitting workflows to study signal-background separation and explore potential signatures of new physics.

In addition to ML-based classification, the repository includes tools for template fitting, RooFit-based studies, and limit-setting workflows, reflecting an end-to-end high-energy physics analysis pipeline.
The project emphasizes transferable data science skills such as high-dimensional data analysis, feature engineering,
model evaluation, and neural-based classification.

## Dataset
- Source: CMS Open Data (CERN)
- Type: High-energy physics collision events
- Characteristics: Large-scale, noisy, high-dimensional scientific data

## Methods

- Data preprocessing and event selection
- Feature engineering based on physics-motivated observables
- Exploratory data analysis and visualization
- Machine learning classification for signal vs background separation
- Template construction and BDT-based fitting
- RooFit-based statistical modeling
- Fit diagnostics and limit-setting studies
- Model validation and performance evaluation

## Repository Contents

- `runAnalysis.C` – main analysis workflow entry point
- `MyClass.C`, `MyClass.h`, `MyClass2.C` – helper classes and analysis logic
- `TMVAClassification.C`, `TMVClassification.C` – TMVA / ML classification workflows
- `myplot.C`, `reconplot.C`, `overlay_signal_allbkgs.C`, `overlay_signal_keybkgs.C` – plotting and visualization scripts
- `binopt.C` – binning optimization studies
- `roofit_example.C`, `roofit_bdt_templates_fit.C` – RooFit-based modeling and template fitting
- `fit_and_limits.cc`, `fit_and_limits1.cc`, `fit_and_limits2.cc`, `fit_and_limits3.cc` – fitting and statistical limit studies
- `logbook.txt` – project notes and development log
- `requirements.txt` – Python dependencies
- `Instructions.txt` – project-specific usage notes. Its in progress. Not final
  
## Tools & Technologies
- C++
- root
- 

## Results
- Successful classification of signal vs background events
- Identification of non-linear correlations using neural-based models
- Framework extensible to anomaly detection and new physics searches
-Pseudodata Production
-Successfull model fitting and limit extractions

## Key Takeaways
- Applied machine learning to real-world scientific data
- Worked with complex, high-dimensional datasets
- Built end-to-end data analysis pipelines
- Bridged physics domain knowledge with data-driven modeling




