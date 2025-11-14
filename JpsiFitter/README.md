# JpsiFitter

Collection of methods for the fit to the dilepton invariant mass distributions. The class is based on RooFit. Code based on `alice/O2Physics/PWGDQ/Macros/`

- DQFitter.py: class collecting all the methods. It allows to fit single invariant mass spectra or perform a multi-trial fit
- fit_library: directory containing a collection of PDFs that are used in dilepton analyses. To add another PDF the user has to follow the same template as those already included
- tutorial.py: simple script which generates the tutorial sample and fits it
- configFit.json: configuration of the fit. Contains all the parameters and PDFs which should be used in the fit

## Tutorial

For a deeper tutorial, see `Jpsi-Jets-Analysis/README.md`