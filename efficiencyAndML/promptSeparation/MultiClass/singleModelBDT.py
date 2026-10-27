import json
import os
from importlib import reload

import matplotlib.pyplot as plt
import numpy as np
import onnx
import pandas as pd
import xgboost as xgb
from hipe4ml import analysis_utils, plot_utils
from hipe4ml.model_handler import ModelHandler
from hipe4ml.tree_handler import TreeHandler
from hipe4ml_converter.h4ml_converter import H4MLConverter

# from hipe4ml.analysis_utils import train_test_generator, bdt_efficiency_array

# from parquet_to_root import parquet_to_root # Apparently parquet_to_root is missing from the repo
# Parameters for training

PT_BINS = [(0, 100)]  # pT bins
# CENT_BINS = [(0, 100)]  # centrality bins #(-1, 100) for pp? Not used?
CENT_Type = "kCentFT0C" # kCentFT0C, kCentFT0A, kCentFT0M

DATA_FILE = "/home/ferrandi/alice/Jpsi-Jets-Analysis/efficiencyAndML/promptSeparation/MultiClass/input/bkg1AOD.parquet.gzip"
PROMPTMC_FILE = "/home/ferrandi/alice/Jpsi-Jets-Analysis/efficiencyAndML/promptSeparation/MultiClass/input/prompt1AOD.parquet.gzip"
NONPROMPTMC_FILE = "/home/ferrandi/alice/Jpsi-Jets-Analysis/efficiencyAndML/promptSeparation/MultiClass/input/nonPrompt1AOD.parquet.gzip"

# Output directory
OUT_DIR = "output/MultiClassMoreStat/myVarsNoPIDOptimized/TEMP1AOD"
DEFAULT_OUT_DIR = "output/singleModelDefaultDir"

def runSingleModel(
    features=None,
    pt_bin=None,
    model_name="default_model",
    out_dir=None,
    data_file=DATA_FILE,
    promptmc_file=PROMPTMC_FILE,
    nonpromptmc_file=NONPROMPTMC_FILE,
    n_estimators=(20, 1000),
    max_depth=(2, 6),
    learning_rate=(0.01, 0.1),
    n_jobs=10,
    timeout=60,
    n_trials=100,
    random_state=42,
    use_previous_model=False,
    parquet_as_input=True,
    save_root_as_parquet=False,
    plot_fewer_shap_feats=False,
    plotFigures=True,
    ):  
    
    """
    Train and evaluate one BDT feature set in one pT interval.

    ---
    ``use_previous_model``: If True, the programs verifies whether a already-trained .pkl model exists in out_dir. If so, it will be used and training will be skipped. If not, training will continue normally for that case.
    ``parquet_as_input``: If True, pass parquet files instead of root files
    ``save_root_as_parquet``: If true, receives root files as input and converts and saves them to parquet
    ``plot_fewer_shap_feats``: If True, SHAP will plot only a few variables. Useful when too many variables (SHAP does not plot too many variables)
    ``plotFigures`` controls interactive display. All figures are still generated
    and saved when it is false. Pass a unique ``out_dir`` for each model/bin.
    """

    if out_dir is None:
        if features is None and pt_bin is None and model_name == "default":
            out_dir = DEFAULT_OUT_DIR
        else:
            selected_bin = tuple(pt_bin) if pt_bin is not None else (0, 100)
            out_dir = os.path.join(
                DEFAULT_OUT_DIR, str(model_name),
                f"pt{selected_bin[0]}_{selected_bin[1]}",
            )
    out_dir = str(out_dir)
    fig_dir = os.path.join(out_dir, "fig")

    PT_BINS = [tuple(pt_bin)] if pt_bin is not None else [(0, 100)]
    CENT_BINS = [(0, 100)]
    pre_cut = '1.8 < fMass < 5.0 and -0.4 < fDcaXY1 < 0.4 and -0.4 < fDcaZ1 < 0.4 and -0.4 < fDcaXY2 < 0.4 and -0.4 < fDcaZ2 < 0.4'
    tree_name = None if parquet_as_input else 'O2rtdielectronall'
    if len(PT_BINS) != 1:
        raise ValueError("runSingleModel processes exactly one pT bin per call")
    if PT_BINS[0][0] >= PT_BINS[0][1]:
        raise ValueError(f"Invalid pT bin: {PT_BINS[0]}")
    if features is not None and not features:
        raise ValueError(f"Feature set '{model_name}' cannot be empty")
    os.makedirs(fig_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True) 

    pre_cut = '1.8 < fMass < 5.0 and -0.4 < fDcaXY1 < 0.4 and -0.4 < fDcaZ1 < 0.4 and -0.4 < fDcaXY2 < 0.4 and -0.4 < fDcaZ2 < 0.4'

    if parquet_as_input:
        tree_name = None
    promptH_origin = TreeHandler(promptmc_file, tree_name, folder_name="DF*")
    nonpromptH_origin = TreeHandler(nonpromptmc_file, tree_name, folder_name="DF*")
    dataH_origin = TreeHandler(data_file, tree_name, folder_name="DF*")

    # for i, (ptmin, ptmax) in enumerate(PT_BINS):
    ptmin, ptmax = PT_BINS[0]

    print(f"Processing pT bin {ptmin} - {ptmax} GeV/c")
    tag = f'pt{ptmin}_{ptmax}'
    cut_pt = f'fPt >= {ptmin} and fPt < {ptmax}'
    cut = f'{pre_cut} and {cut_pt}'


    # Prepare datasets
    promptH = promptH_origin.get_subset(frac=1.0)
    promptH.apply_preselections(cut)
    promptH.apply_preselections('fMcDecision == 1 and 2.4 < fMass < 3.2')
    nonpromptH = nonpromptH_origin.get_subset(frac=1.0)
    nonpromptH.apply_preselections(cut)
    nonpromptH.apply_preselections('fMcDecision == 2 and 2.4 < fMass < 3.2')
    dataH = dataH_origin.get_subset(frac=1.0)
    dataH = dataH.get_subset(cut)
    bkgH = dataH.get_subset('1.8 < fMass < 2.4 or 3.2 < fMass < 5.0', frac=1.0) # size=promptH.get_n_cand() * 2) # Neet to optimize size ratio

    class_counts = {
        "background": len(bkgH),
        "prompt": len(promptH),
        "nonprompt": len(nonpromptH),
    }
    empty_classes = [name for name, count in class_counts.items() if count == 0]
    if empty_classes:
        raise ValueError(
            f"Cannot train model '{model_name}' for pT bin {ptmin}-{ptmax}: "
            f"empty class(es) after selections: {', '.join(empty_classes)}. "
            f"Candidate counts are {class_counts}; check PT_BINS, input files, "
            "and the existing selection cuts."
        )

    if not parquet_as_input and save_root_as_parquet:
        promptH_origin.write_df_to_parquet_files("prompt1AOD", "../input/")
        nonpromptH_origin.write_df_to_parquet_files("nonPrompt1AOD", "../input/")
        dataH_origin.write_df_to_parquet_files("bkg1AOD", "../input/")

    # Remove variables not used for training and plot distributions and correlations

    # "Default":
    # branches_to_remove = ["fSign", "fFilterMap", "fMcDecision", 
    #     "fDCAxyzTrk0KF", "fDCAxyzTrk1KF", "fDCAxyzBetweenTrksKF",
    #     "fDCAxyTrk0KF", "fDCAxyTrk1KF", "fDCAxyBetweenTrksKF",
    #     "fDeviationTrk0KF", "fDeviationTrk1KF", "fDeviationxyTrk0KF", "fDeviationxyTrk1KF",
    #     "fMassKFGeo", "fChi2OverNDFKFGeo", "fDecayLengthKFGeo", "fDecayLengthOverErrKFGeo",
    #     "fDecayLengthXYKFGeo", "fDecayLengthXYOverErrKFGeo", "fPseudoproperDecayTimeKFGeo", "fPseudoproperDecayTimeErrKFGeo",
    #     "fCosPAKFGeo", "fPairDCAxyz", "fPairDCAxy",
    #     "fDeviationPairKF", "fDeviationxyPairKF",
    #     "fMassKFGeoTop", "fChi2OverNDFKFGeoTop"]
    # plus
    # branches_to_remove_for_training = ["fEta", "fPhi",
    #     "fPt1", "fEta1", "fPhi1", "fTPCNClsFound1",
    #     "fDcaXY1", "fDcaZ1",
    #     "fTPCNSigmaPr1", "fTPCSignal1", "fTPCNSigmaPi1",
    #     "fTOFBeta1", "fTOFNSigmaPi1", "fTOFNSigmaPr1", "fTOFNSigmaEl1",
    #     "fPt2", "fEta2", "fPhi2", "fTPCNClsFound2",
    #     "fDcaXY2", "fDcaZ2",
    #     "fTPCNSigmaPr2", "fTPCSignal2", "fTPCNSigmaPi2",
    #     "fTOFBeta2", "fTOFNSigmaPi2", "fTOFNSigmaPr2", "fTOFNSigmaEl2",
    #     "fLz", "fLxy", "fTauz", "fTauxy"]     

    # ALL: All Features:  []'fMass', 'fPt', 'fEta', 'fPhi', 'fSign', 'fFilterMap', 'fMcDecision', 'fPt1', 'fEta1', 'fPhi1', 'fITSClusterMap1', 'fITSChi2NCl1', 'fTPCNClsCR1', 'fTPCNClsFound1', 'fTPCChi2NCl1', 'fDcaXY1', 'fDcaZ1', 'fTPCSignal1', 'fTPCNSigmaEl1', 'fTPCNSigmaPi1', 'fTPCNSigmaPr1', 'fTOFBeta1', 'fTOFNSigmaEl1', 'fTOFNSigmaPi1', 'fTOFNSigmaPr1', 'fPt2', 'fEta2', 'fPhi2', 'fITSClusterMap2', 'fITSChi2NCl2', 'fTPCNClsCR2', 'fTPCNClsFound2', 'fTPCChi2NCl2', 'fDcaXY2', 'fDcaZ2', 'fTPCSignal2', 'fTPCNSigmaEl2', 'fTPCNSigmaPi2', 'fTPCNSigmaPr2', 'fTOFBeta2', 'fTOFNSigmaEl2', 'fTOFNSigmaPi2', 'fTOFNSigmaPr2', 'fDCAxyzTrk0KF', 'fDCAxyzTrk1KF', 'fDCAxyzBetweenTrksKF', 'fDCAxyTrk0KF', 'fDCAxyTrk1KF', 'fDCAxyBetweenTrksKF', 'fDeviationTrk0KF', 'fDeviationTrk1KF', 'fDeviationxyTrk0KF', 'fDeviationxyTrk1KF', 'fMassKFGeo', 'fChi2OverNDFKFGeo', 'fDecayLengthOverErrKFGeo', 'fDecayLengthXYOverErrKFGeo', 'fPseudoproperDecayTimeKFGeo', 'fPseudoproperDecayTimeErrKFGeo', 'fCosPAKFGeo', 'fPairDCAxyz', 'fPairDCAxy', 'fDeviationPairKF', 'fDeviationxyPairKF', 'fMassKFGeoTop', 'fChi2OverNDFKFGeoTop', 'fTauz', 'fTauxy', 'fLz', 'fLxy']

    available_features = promptH.get_var_names()
    features_for_train = list(available_features)

    first_branches_to_remove = ["fPt", "fMass", "fSign", "fFilterMap", "fMcDecision", "fMassKFGeo", "fMassKFGeoTop", "fEta", "fPhi", "fPt1", "fEta1", "fPhi1", "fPt2", "fEta2", "fPhi2", 'fDecayLengthKFGeo', 'fDecayLengthXYKFGeo']

    PID_vars = ["fTPCNSigmaEl1", "fTPCNSigmaPi1", "fTPCNSigmaPr1", "fTOFNSigmaEl1", "fTOFNSigmaPi1", "fTOFNSigmaPr1", "fTPCNSigmaEl2", "fTPCNSigmaPi2", "fTPCNSigmaPr2", "fTOFNSigmaEl2", "fTOFNSigmaPi2", "fTOFNSigmaPr2", "fTPCSignal1", "fTOFBeta1", "fTPCSignal2", "fTOFBeta2" ]


    #Current: Current vars:  ['fITSClusterMap1', 'fITSChi2NCl1', 'fTPCNClsCR1', 'fTPCNClsFound1', 'fTPCChi2NCl1', 'fDcaXY1', 'fDcaZ1', 'fITSClusterMap2', 'fITSChi2NCl2', 'fTPCNClsCR2', 'fTPCNClsFound2', 'fTPCChi2NCl2', 'fDcaXY2', 'fDcaZ2', 'fDCAxyzTrk0KF', 'fDCAxyzTrk1KF', 'fDCAxyzBetweenTrksKF', 'fDCAxyTrk0KF', ', 'fDecayLengthOverErrKFGeo', 'fPseudoproperDecayTimeKFGeo', 'fPseudoproperDecayTimeErrKFGeo', 'fCosPAKFGeo', 'fDeviationPairKF', 'fDeviationxyPairKF', 'fTauz'', 'fLz', 'fLxy']

    further_branches_to_remove = ['fTPCNClsCR1', 'fTPCNClsCR2', 'fTPCNClsFound1', 'fTPCNClsFound2', 'fDCAxyTrk0KF', 'fDCAxyTrk1KF', 'fDeviationTrk0KF', 'fDeviationTrk1KF', 'fDecayLengthXYKFGeo', 'fDecayLengthXYOverErrKFGeo', 'fDecayLengthKFGeo', 'fChi2OverNDFKFGeoTop', 'fChi2OverNDFKFGeo', 'fPairDCAxy', 'fPairDCAxyz', 'fDeviationxyTrk0KF', 'fDeviationxyTrk1KF', 'fDCAxyBetweenTrksKF', 'fTauxy']


    branches_to_remove_for_training = first_branches_to_remove + PID_vars + further_branches_to_remove
    total_branches_to_remove = branches_to_remove_for_training

    for b in total_branches_to_remove:
        if b in features_for_train:
            features_for_train.remove(b)

    if features is not None:
        missing_features = sorted(set(features) - set(available_features))
        if missing_features:
            raise ValueError(f"Feature set '{model_name}' contains unavailable variables: {missing_features}")
        if len(features) != len(set(features)):
            raise ValueError(f"Feature set '{model_name}' contains duplicate variables")
        features_for_train = list(features)


    # vars_to_draw = promptH.get_var_names()
    Nvars = len(features_for_train)
    NFigs = max(1, (Nvars + 11) // 12)
    NVarPerPlot = Nvars // NFigs + (Nvars % NFigs > 0)
    leg_labels = ['background', 'prompt', 'nonprompt']


    for i in range(NFigs):
        start_idx = i * NVarPerPlot
        end_idx = min(start_idx + NVarPerPlot, Nvars)
        vars_subset = features_for_train[start_idx:end_idx]

        plot_utils.plot_distr([bkgH, promptH, nonpromptH], vars_subset, bins=50, labels=leg_labels, log=True, density=True, figsize=(12, 7), alpha=0.3, grid=True)
        plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
        plt.savefig(f'{fig_dir}/varsDistr_{tag}_part{i+1}.pdf')
        if plotFigures:
            plt.show()
        plt.close()

    Nvars = len(features_for_train)
    corr_leg_labels=[]
    for leg in leg_labels:
        corr_leg_labels.append(leg+"-vars Correlations")
    corr_plots1 = plot_utils.plot_corr([bkgH, promptH, nonpromptH], features_for_train, corr_leg_labels)
    corr_plots1[0].savefig(f'{fig_dir}/corr_bkg_{tag}.pdf', bbox_inches='tight')
    corr_plots1[1].savefig(f'{fig_dir}/corr_prompt_{tag}.pdf', bbox_inches='tight')
    corr_plots1[2].savefig(f'{fig_dir}/corr_nonprompt_{tag}.pdf', bbox_inches='tight')

    # Train model
    # model_clf = xgb.XGBClassifier()

    # model_hdl = ModelHandler(model_clf, features_for_train)
    # train_test_data = analysis_utils.train_test_generator([bkgH, promptH, nonpromptH], [0, 1, 2], test_size=0.2, random_state=random_state)
    # hyper_pars_ranges = {'n_estimators': n_estimators, 'max_depth': max_depth, 'learning_rate': learning_rate}
    # model_hdl.optimize_params_optuna(train_test_data, hyper_pars_ranges, cross_val_scoring='roc_auc_ovo', timeout=timeout, n_jobs=n_jobs, n_trials=n_trials, direction='maximize')

    train_test_data = analysis_utils.train_test_generator([bkgH, promptH, nonpromptH], [0, 1, 2], test_size=0.2, random_state=random_state)

    modelPath = f"{out_dir}/modelBDT_{tag}.pkl"

    if use_previous_model and os.path.isfile(modelPath):
        model_hdl = ModelHandler()
        model_hdl.load_model_handler(modelPath)
        print(f"Loaded previous model: {modelPath}")
    else:
        if use_previous_model:
            print(f"Model {modelPath} not found! Continuing to training.")
        model_clf = xgb.XGBClassifier()
        model_hdl = ModelHandler(model_clf, features_for_train)
        hyper_pars_ranges = {'n_estimators': n_estimators, 'max_depth': max_depth, 'learning_rate': learning_rate}
        model_hdl.optimize_params_optuna(train_test_data, hyper_pars_ranges, cross_val_scoring='roc_auc_ovo', timeout=timeout, n_jobs=n_jobs, n_trials=n_trials, direction='maximize')
        model_hdl.train_test_model(train_test_data, multi_class_opt="ovo")


    # Plot model predictions for both training and test samples
    leg_labels = ['bkg', 'prompt', 'nonprompt']
    ml_out_fig = plot_utils.plot_output_train_test(model_hdl, train_test_data, bins=100, 
                                                output_margin=False, labels=leg_labels,
                                                logscale=True, density=True)

    figLabels = ['bkg', 'prompt', 'nonprompt']
    for i, fig in enumerate(ml_out_fig):
        fig.savefig(f'{fig_dir}/prediction_{figLabels[i]}_{tag}.pdf', bbox_inches='tight')
    if plotFigures:
        plt.show()
    for fig in ml_out_fig:
        plt.close(fig)


    reload(analysis_utils)

    y_pred_test = model_hdl.predict(train_test_data[2], output_margin=False)
    y_pred_train = model_hdl.predict(train_test_data[0], output_margin=False)

    # Evaluate the efficiency over this model's full configured pT interval.
    pT_dist = train_test_data[2]['fPt']
    pTmax = min(ptmax, float(np.max(pT_dist)))
    var_name = r"$p_{\mathrm{T}}^{\mathrm{J}/\psi}$ (GeV/$c$)" # pT
    bkgCut = 0.8
    NpTpoints = 20
    NnonpCuts = 6
    NfigColumns = 2

    fig, axes = plt.subplots(NnonpCuts//NfigColumns + (NnonpCuts % NfigColumns > 0), NfigColumns, figsize=(20, 20))
    fig.suptitle(
        rf"BDT efficiency with background score $<$ {bkgCut}. Model for ${ptmin} \leq p_\mathrm{{T}}^{{\mathrm{{J}}/\psi}} < {ptmax}$",
        fontsize=24,
        y=1
    )
    axes = axes.flatten()
    efficiency_cuts = np.round(np.linspace(0, 0.95, NnonpCuts), 2)
    efficiency_arrays = {}

    for i, nonpCut in enumerate(efficiency_cuts):
        pt_centers, prompt_eff, prompt_err, nonprompt_eff, nonprompt_err =  analysis_utils.bdt_efficiency_vs_var(
                train_test_data[3],
                y_pred_test,
                pT_dist,
                nonpCut,
                bkgCut,
                NpTpoints,
                var_max=pTmax
            )
        efficiency_arrays[f"pt_centers_{i}"] = pt_centers
        efficiency_arrays[f"prompt_eff_{i}"] = prompt_eff
        efficiency_arrays[f"prompt_err_{i}"] = prompt_err
        efficiency_arrays[f"nonprompt_eff_{i}"] = nonprompt_eff
        efficiency_arrays[f"nonprompt_err_{i}"] = nonprompt_err

        reload(plot_utils)
        plot_utils.plot_bdt_efficiency_vs_var(pt_centers, prompt_eff, prompt_err, nonprompt_eff, nonprompt_err, nonpCut, bkgCut, var_name, fig, axes[i])
    plt.tight_layout()
    fig.savefig(f'{fig_dir}/eff_vs_Pt_nonPcuts_{tag}.pdf', bbox_inches='tight')
    if plotFigures:
        plt.show()
    plt.close(fig)


    # This significance may be not representative! Since: 1. signal is extracted from MC and 2. bkg from side bands of mass distribution
    y_pred_test = model_hdl.predict(train_test_data[2], output_margin=False) # For each candidate, [bkg_score, prompt_score, nonp_score]
    colors = ['darkorange', 'darkgreen']

    # bdtEff_testSample, thresholds_test, bdtCont_testSample, bdtSignif_testSample = analysis_utils.bdt_efficiency_array(train_test_data[3], y_pred_test, n_points=20, calculate_contamination=True)

    # bdt_eff_figs = plot_utils.plot_bdt_eff(thresholds_test, bdtEff_testSample, bdtCont_testSample, bdtSignif_testSample, figLabels, contaminColors=colors)

    # bdtEff_testSample, thresholds_test, bdtCont_testSample, bdtSignif_testSample = (
    #     analysis_utils.bdt_efficiency_array(
    #         train_test_data[3],
    #         y_pred_test,
    #         n_points=81,
    #         calculate_contamination=True
    #     )
    # )

    # bdt_eff_figs = plot_utils.plot_bdt_eff(
    #     thresholds_test,
    #     bdtEff_testSample,
    #     bdtCont_testSample,
    #     bdtSignif_testSample,
    #     figLabels,
    #     contaminColors=colors,
    #     y_truth=train_test_data[3],
    #     y_score=y_pred_test,
    #     bkg_keep_lower=True,
    #     chosen_bkg_rejec = 0.40
    # )

    # bdtEff_testSample, thresholds_test, *extra = analysis_utils.bdt_efficiency_array(
    #     train_test_data[3], y_pred_test, n_points=81, calculate_contamination=True, calculate_bkg_efficiency=True
    # )
    # bdt_eff_figs = plot_utils.plot_bdt_eff(thresholds_test, bdtEff_testSample, bkg_efficiencies=extra[-3:])
    # plt.show()

    bdtEff_testSample, thresholds_test, bdtCont_testSample, bdtSignif_testSample, \
        prompt_eff, nonprompt_eff, bkg_rejection = analysis_utils.bdt_efficiency_array(
            train_test_data[3],
            y_pred_test,
            n_points=81,
            calculate_contamination=True,
            calculate_bkg_efficiency=True,
        )

    bdt_eff_figs = plot_utils.plot_bdt_eff(
        thresholds_test,
        bdtEff_testSample,
        bdtCont_testSample,
        bdtSignif_testSample,
        figLabels,
        contaminColors=colors,
        chosen_bkg_rejec=0.40,
        bkg_efficiencies=(prompt_eff, nonprompt_eff, bkg_rejection),
    )
    for i, fig in enumerate(bdt_eff_figs[:-1]):
        fig.savefig(f'{fig_dir}/BDT_Efficiencies_{figLabels[i]}_{tag}.pdf', bbox_inches='tight')
    if bdt_eff_figs:
        bdt_eff_figs[-1].savefig(f'{fig_dir}/BDT_BkgRejection_{tag}.pdf', bbox_inches='tight')

    comparison_path = os.path.join(out_dir, f"comparison_{tag}.npz")
    np.savez_compressed(
        comparison_path,
        y_test=np.asarray(train_test_data[3]),
        y_score_test=np.asarray(y_pred_test),
        y_train=np.asarray(train_test_data[1]),
        y_score_train=np.asarray(y_pred_train),
        pt_test=np.asarray(pT_dist),
        efficiency_cuts=efficiency_cuts,
        **efficiency_arrays,
    )
    with open(os.path.join(out_dir, f"comparison_{tag}.json"), "w") as comparison_file:
        json.dump({"model": model_name, "features": features_for_train, "pt_bin": [ptmin, ptmax]}, comparison_file, indent=2)

    bkgRejecDic = {'inclusive':0.4} # todo: decide this value using expected significance from data. One value for each (z, pTjet)?

    y_pred_train = model_hdl.predict(train_test_data[0], output_margin=False)
    y_pred_test = model_hdl.predict(train_test_data[2], output_margin=False)
    roc_train_test_fig = plot_utils.plot_roc_train_test(train_test_data[3], y_pred_test,
                                                        train_test_data[1], y_pred_train, None, leg_labels, 'macro', multi_class_opt="ovr")
    roc_train_test_fig.savefig(f'{fig_dir}/roc_train_test_{tag}.pdf')
    plt.close(roc_train_test_fig)

    # learn_curve can take long!
    #  learn_curve_fig = plot_utils.plot_learning_curves(model_hdl, train_test_data, 4, figLabels)
    # learn_curve_fig.savefig(f'{fig_dir}/learning_curves_{tag}.pdf')

    chosenFeatures = None # If None, all training vars will be plotted

    if plot_fewer_shap_feats:
        dcaVars = [var for var in features_for_train if ('dca' in var.lower() or 'length' in var.lower() or 'deviation' in var.lower()or 'lx' in var.lower() or 'KF' in var or 'tau' in var.lower() or 'cos' in var.lower() or var=="fChi2OverNDFKFGeoTop") and ('mass' not in var.lower() and 'xyz' not in var.lower())] 
        chosenFeatures=dcaVars

    shap_figs = plot_utils.plot_feature_imp(train_test_data[2], train_test_data[3], model_hdl, labels=leg_labels, features=chosenFeatures)
    # only show the SHAP summary plot
    # plt.close(shap_figs[0])
    # plt.close(shap_figs[1])
    # plt.close(shap_figs[2])
    shap_figs[3].savefig(f'{fig_dir}/shapSummary_{tag}.pdf', bbox_inches='tight')
    if plotFigures:
        plt.show()
    plt.close('all')

    # Save parquet, root and pkl
    parquet_path = f'{out_dir}/bdt_{tag}.parquet.gzip'
    df = dataH.get_data_frame()
    # parquet_to_root(parquet_path, f'{out_dir}/bdt_{tag}.root', treename='parquettree', verbose=False) # Apparently parquet_to_root is missing from the repo
    model_hdl.dump_model_handler(modelPath)

    df = dataH.get_data_frame()
    y_pred = model_hdl.model.predict_proba(df[features_for_train]) # Sums to 1 for each candidate!
    df['background_score'] = y_pred[:, 0] # prompt
    df['prompt_score'] = y_pred[:, 1] # non-prompt
    df['nonprompt_score'] = y_pred[:, 2] # background
    df.to_parquet(parquet_path, index=False, compression='gzip')

    score_fig, score_ax = plt.subplots(figsize=(12, 7))
    score_ax.hist(df['prompt_score'], bins=100, alpha=0.6, log=True, label='prompt')
    score_ax.hist(df['nonprompt_score'], bins=100, alpha=0.6, log=True, label='nonprompt')
    score_ax.hist(df['background_score'], bins=100, alpha=0.6, log=True, label='background')
    score_ax.set_xlabel('BDT score')
    score_ax.set_ylabel('Counts')
    score_ax.legend()
    score_fig.savefig(f'{fig_dir}/score_distributions_{tag}.pdf', bbox_inches='tight')
    if plotFigures:
        plt.show()
    plt.close(score_fig)

    def infer_n_features(mh: ModelHandler, input_table=None):
        try:
            return mh.model.n_features_in_
        except AttributeError:
            if input_table is not None:
                return pd.read_parquet(input_table).shape[1]
            raise ValueError("Cannot determine n_features. Provide input table.")


    def export_to_onnx(model_handler, prefix, feature_names, input_table=None):
        n_features = infer_n_features(model_handler, input_table)
        print(f"Inferred number of input features: {n_features}")

        converter = H4MLConverter(model_handler)
        onnx_model = converter.convert_model_onnx(1)

        for i, name in enumerate(feature_names):
            meta = onnx_model.metadata_props.add()
            meta.key = f"feature_{i}"
            meta.value = name

        out_path = f"{prefix}.onnx"
        converter.dump_model_onnx(out_path)
        print(f"Exported ONNX: {out_path}")

        model = onnx.load(out_path)
        for input_tensor in model.graph.input:
            dims = [d.dim_value if d.HasField("dim_value") else "None" for d in input_tensor.type.tensor_type.shape.dim]
            print(f"  - {input_tensor.name}: shape={dims}")


    # Export to ONNX
    features_onnx = [f.replace("f", "k", 1) if f.startswith("f") else f for f in features_for_train]
    export_to_onnx(model_hdl, prefix=f"{out_dir}/modelBDT_{tag}", feature_names=features_onnx, input_table=parquet_path)

    print("Training and evaluation completed. Start JSON export.")

    cut_vals = {"background": 0.3, "prompt": 0.5, "nonprompt": 0.4} # TODO: Update with optimized cut values

    json_dict = {
        "TestCut": {
            "type": "MultiClass",
            "title": f"JpsiBDTModel_{tag}",
            "inputFeatures": features_onnx,
            "modelFiles": [],
            "cent": CENT_Type
        }
    }

    CENT_BINS = [(0, 100)]
    for centmin, centmax in CENT_BINS:
        cent_key = f"AddCentCut-Cent{centmin:02d}{centmax:02d}"
        json_dict["TestCut"][cent_key] = {
            "centMin": centmin,
            "centMax": centmax
        }
        for i, (ptmin, ptmax) in enumerate(PT_BINS):
            pt_key = f"AddPtCut-pTBin{i+1}"
            model_name_onnx = f"modelBDT_{tag}.onnx"

            json_dict["TestCut"]["modelFiles"].append(model_name_onnx)
            json_dict["TestCut"][cent_key][pt_key] = {
                "pTMin": ptmin,
                "pTMax": ptmax,
                "AddMLCut-background": {"var": "kBdtBackground", "cut": cut_vals["background"], "exclude": True},
                "AddMLCut-prompt": {"var": "kBdtPrompt", "cut": cut_vals["prompt"], "exclude": True},
                "AddMLCut-nonprompt": {"var": "kBdtNonprompt", "cut": cut_vals["nonprompt"], "exclude": False}
            }

    json_output_path = os.path.join(out_dir, "bdt_config_multiclass.json")
    with open(json_output_path, "w") as f_json:
        json.dump(json_dict, f_json, indent=2)

    with open(json_output_path, "r") as f_json_read:
        data = json.load(f_json_read)

    hyperloop_str = json.dumps(data, separators=(',', ':'))
    print("Hyperloop JSON string:\n")
    print(hyperloop_str)

    escaped_str = json.dumps(hyperloop_str)
    print("\nLocal JSON string:\n")
    print(escaped_str)

    print(f"\nJSON saved to {json_output_path}")


    return {
        "model": model_name,
        "features": list(features_for_train),
        "pt_bin": [ptmin, ptmax],
        "output_dir": out_dir,
        "comparison_file": comparison_path,
        "model_handler": model_hdl,
    }


if __name__ == "__main__":
    runSingleModel()
