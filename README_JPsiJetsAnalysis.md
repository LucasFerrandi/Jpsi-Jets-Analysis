# Tutorial on running J/Psi jets analysis

(created on April 11, 2025. Aliases aimed to my pc)

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

          - Example: `@CompleteAO2DFilesSubJobs_24Gb_LHC24am_pass1.txt` 
     - Parameters for the jet finder, such as:

          - Vertex z cut

          - event selection
          - p_T, eta, phi
          - Jetalgorithm
     - Parameters for JPsi Task, such as:
          p_T bounds (ex.: 5, 7, 15 and 35 GeV)
- If a json file doesn't exist yet, run without it and it will be automatically generated
- A script was created to run this, which can be executed as following:
- cd into ~/alice/O2Physics/PWGJE/Tasks/JPsiWorkDir:
     `cd $JPsiDir`
- Run the script providing the json file:
     `./RunJPsiFragmentation.sh dpl-config.json`

     -  (Or `o2-analysis-je-jet-jpsi-fragmentation -b --configuration json:<JSON> | o2-analysis-je-jet-finder-dielectron-data-charged -b --configuration json://<JSON>`)
- A `AnalysisResults.root` should've been created with jet spectra and z-vs-mass histograms for each pT range
6. **Run the fitter:**
- `cd ~/alice/O2Physics/PWGDQ/Macros/`
- a json file must me passed with
     - Input root file
     - Histogram to be fitted
     - Signal function
     - Bkg function
     - etc.
- Run `tutorial.py`, which runs `DQFitter` for every x-projection of the histo:
     - `python tutorial.py configFit_z_Xi.json --run_fit_projections`
- A root file should've been created in output/ containing:
     - A copy of the input histogram and inclusive mass distribution
     -  Results of the fits for every x projection
     - For each fit range, a histogram compiling these results

## Monte Carlo
Analysis done in `/PWGJE/Tasks/JPsiWorkDir/JPsiMC/`
- Given J/Psi Monte Carlo datasets, Skim them using TableMakerMC_withAssoc
     - Provide reco and gen level MCSignals to it. Ex.: `eFromJpsi,eFromPromptJpsi,eFromNonpromptJpsi,Jpsi,nonPromptJpsi,promptJpsi,allBeautyHadrons,Bplus,protonPrimary, everythingFromEverythingFromBeauty`
     - Output:
          - AO2D.root with reduced tables for reco-level tracks (mainly electrons) and gen-level (particles matching provided MCSignals or matching reco tracks)
          - AnalysisResults.root with general statistics
- run dqEfficiency_withAssoc on reduced AO2D
     - analysis on pairs of electrons
     - For `same-event-pairing`, provide gen and reco MCsignals, such as`"cfgBarrelMCGenSignals": "promptJpsi,nonPromptJpsi"` and`"cfgBarrelMCRecSignals": "eePrimaryFromPromptJPsi,eePrimaryFromNonPromptJPsi"`
     - Output:
          - AnalysisResults.root with histograms for dielectron which passes each MCSignal
          - AnalysisResults_Trees.root with dielectron trees (tables). Can be `dielectronsAll`, for example

## Machine Learning and Efficiency
1. Create a Python virtual environment inside ALIEnv (`alice/Hipe4MLenv`, for example)
2. Install Hipe4ML in it
3. Enter ALIEnv: `alice` & `source bashrc_alienv`
4. Enter Hipe4ML venv: `source Hipe4MLev/bin/activate`
5. Efficiency calculated in `alice/EfficiencyAndML/Hipe4MCYuanjing`
     - Open with `code` and select `Hipe4MLenv` Python kernel
6. BDT trained in the same folder
7. With the model trained and converted to ONNX, do the inference on data using MLResponse