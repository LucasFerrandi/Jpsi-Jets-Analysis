# Tutorial on running J/Psi jets analysis

(tutorial created on April 11, 2025)

1. **Install O2Physics (see O2 Documentation)**  
2.  **Enter O2Physics environment:**
  - Using alias:`alice` & `source bashrc_alienv`
  - or `alienv enter O2Physics/latest`
3. **The input AODs must be downloaded and prepared**
4. **Compile my task:**
  - Can be done with AliBuild or with Ninja
  - To compile with Ninja:
  	- alias: `ninjajpsi`
  	- or `cd ~/alice/sw/BUILD/O2Physics-latest/O2Physics/` & `ninja PWGJE/Tasks/install -j 3`
5. **Run my task:**
  - To run, a json file must be provided containing:
  	- the AOD file as input or a txt file containing a list of AODs paths (in this case the string must start with "@")
  	- Parameters for the jet finder, such as:
    		- Vertex z cut
    		- event selection
    		- p_T, eta, phi
    		- Jet algorithm
  	- Parameters for JPsi Task, such as:
    		- p_T bounds (ex.: 5, 7, 15 and 35 GeV)
  	- If a json file doesn't exist yet, run without it and it will be automatically generated
  - A script was created to run this, which can be executed as following:
  	- cd into ~/alice/O2Physics/PWGJE/Tasks/JPsiWorkDir:
    		- `cd $JPsiDir`
  	- Run the script providing the json file:
    		- `./RunJPsiFragmentation.sh dpl-config.json`
    		- (Or `o2-analysis-je-jet-jpsi-fragmentation -b --configuration json:<JSON> | o2-analysis-je-jet-finder-dielectron-data-charged -b --configuration json://<JSON>`)
  - A `AnalysisResults.root` should've been created with jet spectra and z-vs-mass histograms for each pT range
6. **Run the fitter:**
  - `cd ~/alice/O2Physics/PWGDQ/Macros/`
  - a json file must me passed with
  	- Input root file
  	- Histogram to be fitted
  	- Signal function
  	- Bkg function
  	- etc.
  - `tutorial.py` runs `DQFitter` for every x-projection of the histo:
  	- `python tutorial.py WithoutPsi2sconfigFit_MuensterWorkshop.json --run_fit_projections`
  	- A root file should've been created in output/ containing:
    	  - A copy of the input histogram and inclusive mass distribution
    	  - Results of the fits for every x projection
    	  - For each fit range, a histogram compiling these results