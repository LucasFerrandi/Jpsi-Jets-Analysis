# J/psi-jets Analysis

(created on April 11, 2025. Aliases aimed to my pc)

Repository for J/psi-jets analysis within ALICE's Run 3. In general terms, this framework is meant for:

- Prepare ALICE's data for analysis via data skimming
- Reconstruct J/psi mesons from dielectron decays.
- Find J/psi jets using anti-kt clustering algorithm
- Calculate the momentum-fraction distribution of J/psi in jets and other related observables
- Separate prompt from non-prompt J/psis
- Apply similar methods to Monte Carlo (MC) datasets. With this,
     - Calculate reconstruction efficiencies
     - Calculate bayesian-unfolding corrections
- Apply those corrections on measured distributions
- Compare the results with theoretical models

Below, a detailed tutorial on running things locally.

## Contents

- [Preparing the Framework](#preparing-the-framework)
- [Calculation of Momentum Fraction](#calculation-of-momentum-fraction)
- [Monte Carlo](#monte-carlo)
- [Efficiency Calculation](#efficiency-calculation)
- [Separation of Prompt and Non-Prompt J/psi](#separation-of-prompt-and-non-prompt-jpsi)

## Preparing the Framework

The usage of this repository generally demands ALICE's O2 analysis framework

1. Install O2Physics (see [O2 Documentation](https://aliceo2group.github.io/analysis-framework/docs/))
2. Recomendations:

     - add this to your .bashrc in order to enter ALICE environment with `alice`:

     ```sh

     O2PhysicsBranch='O2Physics/latest'
     alice() {
     cd ~/alice/
     alienv enter $O2PhysicsBranch
     }

     ```

     - Add `bashrc_alienv` to `alice/` folder
     - Add `Jpsi-Jets-Analysis/` repository to `alice/`

3. Enter O2Physics environment:
     - Using alias: `alice` & `source bashrc_alienv`
     - Or `alienv enter O2Physics/latest`

## Calculation of Momentum Fraction

This is mainly done by `jpsiFragmentationFunction.cxx` task.

1. The input AODs must be downloaded and prepared
2. Compile my task:
     - Can be done with AliBuild or with Ninja
     - To compile with Ninja (maybe `direnv allow` is needed!):
          - alias: `ninjajpsi`
          - or `cd ~/alice/sw/BUILD/O2Physics-latest/O2Physics/` & `ninja PWGJE/Tasks/install -j 3`
3. Run my task:
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
     - cd into `~/alice/Jpsi-Jets-Analysis/JpsiWorkDir` or
          `cd $JpsiDir` (alias)
     - Run the script providing the json file:
          - Make it executable: `chmod +x runJPsiFragmentation.sh`
          - Execute it: `./runJPsiFragmentation.sh dpl-config.json`

          - (Or `o2-analysis-je-jet-jpsi-fragmentation -b --configuration json:<JSON> | o2-analysis-je-jet-finder-dielectron-data-charged -b --configuration json://<JSON>`)
     - A `AnalysisResults.root` should've been created with jet spectra and z-vs-mass histograms for each pT range
4. Run the fitter:
     - `cd ~/alice/Jpsi-Jets-Analysis/JpsiFitter`
     - a json file must me passed with
          - Input root file
          - Histogram to be fitted
          - Signal function
          - Bkg function
          - etc.
     - Run `tutorial.py` (or similar ones), which runs `DQFitter` for every x-projection of the histo:
          - `python tutorial.py configFit_z_Xi.json --run_fit_projections`
     - A root file should've been created in `output/`, containing:
          - A copy of the input histogram and inclusive mass distribution
          - Results of the fits for every x projection
          - For each fit range, a histogram compiling these results

## Monte Carlo

Analysis done in `~/alice/Jpsi-Jets-Analysis/JpsiWorkDir/MC`

- Given J/psi Monte Carlo datasets, Skim them using TableMakerMC_withAssoc. Usually using Hyperloop.
     - Provide reco and gen level MCSignals to it. Ex.: `eFromJpsi,eFromPromptJpsi,eFromNonpromptJpsi,Jpsi,nonPromptJpsi,promptJpsi,allBeautyHadrons,Bplus,protonPrimary, everythingFromEverythingFromBeauty`
     - Output:
          - AO2D.root with reduced tables for reco-level tracks (mainly electrons) and gen-level (particles matching provided MCSignals or matching reco tracks)
          - AnalysisResults.root with general statistics
- run dqEfficiency_withAssoc on reduced AO2D (done in `~/alice/Jpsi-Jets-Analysis/JpsiWorkDir/MC/DQEfficiency`)
     - analysis on pairs of electrons
     - For `same-event-pairing`, provide gen and reco MCsignals, such as`"cfgBarrelMCGenSignals": "promptJpsi,nonPromptJpsi"` and`"cfgBarrelMCRecSignals": "eePrimaryFromPromptJPsi,eePrimaryFromNonPromptJPsi"`
     - Important: in order to produce `dielectronAll` table, one must:
          - Produce the reduced dataset containing `ReducedTracksBarrelInfo`
          - Enable `processBarrelOnlyWithCollSkimmed`, `cfgFlatTables` and `fgUseKFVertexing`
          - Add `kalman-filter` histograms to `cfgAddSEPHistogram`
     - Output:
          - AnalysisResults.root with histograms for dielectron which passes each MCSignal
          - AnalysisResults_Trees.root with dielectron trees (tables). Can be `dielectronsAll`, for example

## Efficiency Calculation

Done by the matching between MC-truth-level and MC-reconstructed-level J/psis
1. Create a Python virtual environment inside ALIEnv (`alice/Hipe4MLenv`, for example)
2. Install Hipe4ML in it
3. Enter ALIEnv: `alice` & `source bashrc_alienv`
4. Enter Hipe4ML venv: `source Hipe4MLev/bin/activate`
5. Efficiency calculated in `Jpsi-Jets-Analysis/efficiencyAndML/efficienciesJpsi.ipynb`
     - Open with `code` and select `Hipe4MLenv` Python kernel

## Separation of Prompt and Non-Prompt J/psi

This is done using Boosted Decision Trees

1. Produce a MC J/psi tree (such as one containing `dielectronall` produced using `dqEfficiency`)
2. Enter Hipe4ML virtual environment such as described in [Efficiency section](#efficiency-calculation)
3. run `promptSeparation.ipynb`
4. TODO (25-11-14): With the model trained and converted to ONNX, do the inference on data using MLResponse framework