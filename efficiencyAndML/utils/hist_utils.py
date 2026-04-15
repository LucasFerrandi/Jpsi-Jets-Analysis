import ROOT
import numpy as np
import pandas as pd
import ROOT
import os
import glob

pTstr = "p_{T}(GeV/c)"

def test_utils():
 print("hist_utils.py imported successfully :)")
 print("In order to update it more dinamically, use 'reload' method from 'importlib' library.")

def clean_name(name):
    """
    Clean axis names by removing spaces, "MC", and standardizing common variable names.
    """
    
    if "pt" in name.lower() or "p_{t}" in name.lower():
        return pTstr
    if "eta" in name.lower():
        return "eta"
    if "phi" in name.lower():
        return "phi"
    name = name.replace("MC", "")
    name = name.replace(" ", "")
    name = name.replace("(GeV/c)", "")
    # name = name.replace("#", "")
    return name

def hist_to_df(hist, clean_names=True, isLabeledHist=False):  # For pyROOT histograms
    """
    Convert a TH1 or TH2 histogram to a pandas DataFrame.
    If clean_names=True, removes spaces and 'MC' from axis names.
    """

    if hist.ClassName().startswith("TH1"):
        # print("Converting TH1 histogram to DataFrame...")
        nbins = hist.GetNbinsX()
        edges = np.array([hist.GetBinLowEdge(i) for i in range(1, nbins+2)])
        values = np.array([hist.GetBinContent(i) for i in range(1, nbins+1)])

        x_name = hist.GetXaxis().GetTitle()
        if clean_names:
            x_name = clean_name(x_name)
        
        if isLabeledHist: # For histograms with string labels on x-axis
            labels = []
            for i in range(1, nbins+1):
                label = hist.GetXaxis().GetBinLabel(i)
                labels.append(label)
            return pd.DataFrame({
                f"{x_name}Label": labels,
                "counts": values
            })
        else:
            return pd.DataFrame({
                f"{x_name}_left": edges[:-1],
                f"{x_name}_right": edges[1:],
                "counts": values
            })

    elif hist.ClassName().startswith("TH2"):
        # print("Converting TH2 histogram to DataFrame...")
        nx = hist.GetNbinsX()
        ny = hist.GetNbinsY()
        x_edges = np.array([hist.GetXaxis().GetBinLowEdge(i) for i in range(1, nx+2)])
        y_edges = np.array([hist.GetYaxis().GetBinLowEdge(i) for i in range(1, ny+2)])

        data = []
        for i in range(1, nx+1):
            for j in range(1, ny+1):
                x_left, x_right = x_edges[i-1], x_edges[i]
                y_left, y_right = y_edges[j-1], y_edges[j]
                count = hist.GetBinContent(i, j)
                data.append([x_left, x_right, y_left, y_right, count])

        x_name = hist.GetXaxis().GetTitle()
        y_name = hist.GetYaxis().GetTitle()
        if clean_names:
            x_name = clean_name(x_name)
            y_name = clean_name(y_name)
        
        if isLabeledHist: # For histograms with string labels on x-axis
            labels_x = []
            labels_y = []
            for i in range(1, nx+1):
                for j in range(1, ny+1):
                    label_x = hist.GetXaxis().GetBinLabel(i)
                    label_y = hist.GetYaxis().GetBinLabel(j)
                    labels_x.append(label_x)
                    labels_y.append(label_y)

            # for i in range(1, nx+1):
            #     label_x = hist.GetXaxis().GetBinLabel(i)
            #     for j in range(1, ny+1):
            #         labels_x.append(label_x)
            # for j in range(1, ny+1):
            #     label_y = hist.GetYaxis().GetBinLabel(j)
            #     for i in range(1, nx+1):
            #         labels_y.append(label_y)
            columnsNames = [f"{x_name}Label", f"{y_name}Label", "counts"]

            df = pd.DataFrame(data={
                f"{x_name}XLabel": labels_x,
                f"{y_name}YLabel": labels_y,
                "counts": [row[4] for row in data]
            })

            return df

        columnsNames = [f"{x_name}_left", f"{x_name}_right", f"{y_name}_left", f"{y_name}_right", "counts"]

        return pd.DataFrame(data, columns=columnsNames) # if isLabeledHist, this line is irrelevant

    else:
        raise NotImplementedError(f"Only TH1 and TH2 supported, but got {hist.ClassName()}")
    
def rebin_df(df, new_edges, x_left="x_left", x_right="x_right", counts="counts"):
    """
    Rebins a DataFrame with columns [x_left, x_right, counts] into new bins defined by new_edges.
    """

    new_counts = []
    
    for i in range(len(new_edges)-1):
        bin_left = new_edges[i]
        bin_right = new_edges[i+1]
        # Select all original bins that overlap the new bin
        mask = (df[x_right] > bin_left) & (df[x_left] < bin_right)
        
        # Sum the counts of overlapping bins (simple sum; more precise weighting can be applied if needed)
        new_counts.append(df.loc[mask, counts].sum())
        
    return pd.DataFrame({
        "bin_left": new_edges[:-1],
        "bin_right": new_edges[1:],
        "counts": new_counts
    })

def df_to_root_graph(df, graph_name="graph", graph_title="Histogram with errors", yerr=None):
    """
    Convert a pandas DataFrame with 'bin_left', 'bin_right', 'counts' to a ROOT TGraphErrors
    with horizontal error bars corresponding to bin widths.
    """
    
    df = df.dropna(subset=["counts"]) # Exlude bins with NaN counts
    # print(df)
    bin_left = df["bin_left"].values
    bin_right = df["bin_right"].values
    counts = df["counts"].values

    x = (bin_left + bin_right) / 2.0           # bin centers
    y = counts
    ex = (bin_right - bin_left) / 2.0          # horizontal errors = half bin width
    # ey = np.sqrt(counts)                        # vertical errors (Poisson), can set to 0 if not needed
    ey=np.zeros_like(y)  # No vertical errors, can be set to zero
    if yerr is not None:
        ey = yerr              # vertical errors from DataFrame if provided
    graph = ROOT.TGraphErrors(len(x),
                              x.astype(float), y.astype(float),
                              ex.astype(float), ey.astype(float))
    graph.SetName(graph_name)
    graph.SetTitle(graph_title)
    graph.SetMarkerStyle(20)
    # graph.SetMarkerSize(0.8)

    return graph

def EfficiencyPlot(a, b, path, graph_title, legText, ptBins = [0, 1, 2, 3, 4, 5, 6, 7, 8, 12]):
    """
    Apllies cut on |eta| < 0.9 and make plot of a/b vs pt
    Expected that a = reconstructed and b = generated
    """
    root_file = ROOT.TFile.Open(path)    
    output_dir_name = 'analysis-same-event-pairing/output;1'
    output_dir = root_file.Get(output_dir_name)
    if not output_dir:
        print(f"{output_dir_name} not found in the ROOT file.")

    aList = None
    bList = None
    for JPsiList in output_dir:
        # print("JPsiList: ", JPsiList.GetName())
        if JPsiList.GetName() == a:
            aList = JPsiList
        elif JPsiList.GetName() == b:
            bList = JPsiList
    if not aList:
        print(f"{a} not found in the output directory!")
    if not bList:
        print(f"{b} not found in the output directory!")

    df_a = None
    df_b = None

    for hist in aList:
        if "Truth" in a:
            if hist.GetName() == "EtaMC_PtMC":
                df_a = hist_to_df(hist, clean_names=True)
        else:
            if hist.GetName() == "Eta_Pt":
                df_a = hist_to_df(hist, clean_names=True)
    for hist in bList:
        if "Truth" in b:
            if hist.GetName() == "EtaMC_PtMC":
                df_b = hist_to_df(hist, clean_names=True)
        else:
            if hist.GetName() == "Eta_Pt":
                df_b = hist_to_df(hist, clean_names=True)
    
    # filter out |eta| > 0.9
    etaCut = 0.9
    etaCutPlus = etaCut + 0.0001# To account for floating-point precision issues
    df_aFilt = df_a[~((df_a["eta_left"] < -etaCutPlus) | (df_a["eta_right"] > etaCutPlus))]
    df_bFilt = df_b[~((df_b["eta_left"] < -etaCutPlus) | (df_b["eta_right"] > etaCutPlus))]

    df_aFiltReb = rebin_df(df_aFilt, ptBins, x_left=f"{pTstr}_left", x_right=f"{pTstr}_right", counts="counts")
    df_bFiltReb = rebin_df(df_bFilt, ptBins, x_left=f"{pTstr}_left", x_right=f"{pTstr}_right", counts="counts")

    df_eff = df_aFiltReb.copy()
    df_eff['counts'] = df_aFiltReb['counts'] / df_bFiltReb['counts']

    # Statistical Uncertainty
    eff_statUnc = np.sqrt(df_eff['counts'] * (1 - df_eff['counts']) / df_bFiltReb['counts']) # Binomial uncertainty
    graph = df_to_root_graph(df_eff, "efficiency_graph", graph_title, yerr=np.array(eff_statUnc))
    ROOT.SetOwnership(graph, True)
    graph.GetXaxis().SetTitle(f"{pTstr} (GeV/c)")
    graph.GetYaxis().SetTitle("A x #epsilon")
    graph.SetMarkerColor(ROOT.kRed-2)
    graph.SetLineColor(ROOT.kRed-2)
    graph.SetLineWidth(2)

    legend = ROOT.TLegend(0.55, 0.7, 0.9, 0.9)
    legend.AddEntry(legend,f"pp @ 13.6 TeV, |#eta|<{etaCut}", "")
    legend.AddEntry(graph, legText, "p")
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.03)

    output_dir.Delete()
    root_file.Close()
    root_file.Delete()

    return graph, legend

def PIDEfficiencyPlot(a, b, path, graph_title, legText, ptBins = [1, 2, 3, 4, 5, 6, 7, 8, 20], output_dir_name='analysis-p-i-d-efficiency;1'):
    """
    Apllies cut on |eta| < 0.9 and make plot of a/b vs pt (e.g. hWeightedJpsiEffPtEta / hMatchedJpsiPtEta)
    Expected that a = Weighted efficiency and b = simple spectrum (e.g. pT)
    """
    root_file = ROOT.TFile.Open(path)
    output_dir = root_file.Get(output_dir_name)
    if not output_dir:
        print(f"{output_dir_name} not found in the ROOT file.")

    df_a = None
    df_b = None

    hWeighted = output_dir.Get(a)
    df_a = hist_to_df(hWeighted, clean_names=True)
    hDist = output_dir.Get(b)
    df_b = hist_to_df(hDist, clean_names=True)

    if not hWeighted or not hDist:
        raise RuntimeError("Histograms not found")
    # filter out |eta| > 0.9
    etaCut = 0.9
    etaCutPlus = etaCut + 0.0001# To account for floating-point precision issues
    df_aFilt = df_a[~((df_a["eta_left"] < -etaCutPlus) | (df_a["eta_right"] > etaCutPlus))]
    df_bFilt = df_b[~((df_b["eta_left"] < -etaCutPlus) | (df_b["eta_right"] > etaCutPlus))]

    # Finding the pT column name dynamically
    for col in df_aFilt.columns:
        if "pt" in col.lower() or "p_{t}" in col.lower() or "{p}_{t}" in col.lower():
            pTstr = col.replace("_left", "").replace("_right", "")
            break

    # Filter ou pT < 1 Gev. Because of the cut in electron pt
    ptCut = 1.0
    ptCutPlus = ptCut + 0.0001
    df_aFilt = df_aFilt[(df_aFilt[f"{pTstr}_right"] > ptCutPlus)]
    df_bFilt = df_bFilt[(df_bFilt[f"{pTstr}_right"] > ptCutPlus)]

    # print("df_aFilt: \n", df_aFilt)
    # print("df_bFilt: \n", df_bFilt)

    df_aFiltReb = rebin_df(df_aFilt, ptBins, x_left=f"{pTstr}_left", x_right=f"{pTstr}_right", counts="counts")
    df_bFiltReb = rebin_df(df_bFilt, ptBins, x_left=f"{pTstr}_left", x_right=f"{pTstr}_right", counts="counts")

    # print("df_aFiltReb: \n", df_aFiltReb)
    # print("df_bFiltReb: \n", df_bFiltReb)

    df_eff = df_aFiltReb.copy()
    df_eff['counts'] = df_aFiltReb['counts'] / df_bFiltReb['counts']

    # print(df_eff)


    # Statistical Uncertainty
    eff_statUnc = np.sqrt(df_eff['counts'] * (1 - df_eff['counts']) / df_bFiltReb['counts']) # Binomial uncertainty
    graph = df_to_root_graph(df_eff, "efficiency_graph", graph_title, yerr=np.array(eff_statUnc))
    ROOT.SetOwnership(graph, True)
    graph.GetXaxis().SetTitle(pTstr)
    graph.GetYaxis().SetTitle("#epsilon_{PID}")
    graph.SetMarkerColor(ROOT.kRed-2)
    graph.SetLineColor(ROOT.kRed-2)
    graph.SetLineWidth(2)

    legend = ROOT.TLegend(0.55, 0.7, 0.9, 0.9)
    legend.AddEntry(legend,f"pp @ 13.6 TeV, |#eta|<{etaCut}", "")
    legend.AddEntry(graph, legText, "p")
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.03)

    output_dir.Delete()
    root_file.Close()
    root_file.Delete()

    return graph, legend

def PurityPlot(a, b, path, graph_title, legText, ptBins = [0, 1, 2, 3, 4, 5, 6, 7, 8, 12], legPos=[0.55, 0.7, 0.9, 0.9]):
    """
    Apllies cut on |eta| < 0.9 and make plot of a/(a+b) vs pt
    Expected that a = correct association and b = incorrect association
    """
    root_file = ROOT.TFile.Open(path)    
    output_dir_name = 'analysis-same-event-pairing/output;1'
    output_dir = root_file.Get(output_dir_name)
    if not output_dir:
        print(f"{output_dir_name} not found in the ROOT file.")

    aList = None
    bList = None
    for JPsiList in output_dir:
        # print("JPsiList: ", JPsiList.GetName())
        if JPsiList.GetName() == a:
            aList = JPsiList
        elif JPsiList.GetName() == b:
            bList = JPsiList
    if not aList:
        print(f"{a} not found in the output directory!")
    if not bList:
        print(f"{b} not found in the output directory!")

    df_a = None
    df_b = None

    for hist in aList:
        if "Truth" in a: # For purity, never expected "Truth"
            if hist.GetName() == "EtaMC_PtMC":
                df_a = hist_to_df(hist, clean_names=True)
        else:
            if hist.GetName() == "Eta_Pt":
                df_a = hist_to_df(hist, clean_names=True)
    for hist in bList:
        if "Truth" in b:
            if hist.GetName() == "EtaMC_PtMC":
                df_b = hist_to_df(hist, clean_names=True)
        else:
            if hist.GetName() == "Eta_Pt":
                df_b = hist_to_df(hist, clean_names=True)
    
    # filter out |eta| > 0.9
    etaCut = 0.9
    etaCutPlus = etaCut + 0.0001# To account for floating-point precision issues
    df_aFilt = df_a[~((df_a["eta_left"] < -etaCutPlus) | (df_a["eta_right"] > etaCutPlus))]
    df_bFilt = df_b[~((df_b["eta_left"] < -etaCutPlus) | (df_b["eta_right"] > etaCutPlus))]

    df_aFiltReb = rebin_df(df_aFilt, ptBins, x_left=f"{pTstr}_left", x_right=f"{pTstr}_right", counts="counts")
    df_bFiltReb = rebin_df(df_bFilt, ptBins, x_left=f"{pTstr}_left", x_right=f"{pTstr}_right", counts="counts")

    df_purity = df_aFiltReb.copy()
    nTotAssocs = df_aFiltReb['counts'] + df_bFiltReb['counts']
    df_purity['counts'] = df_aFiltReb['counts'] / nTotAssocs

    # Statistical Uncertainty
    purity_statUnc = np.sqrt(df_purity['counts'] * (1 - df_purity['counts']) / nTotAssocs) # Binomial uncertainty
    graphPurity = df_to_root_graph(df_purity, "purity_graph", graph_title, yerr=np.array(purity_statUnc))
    ROOT.SetOwnership(graphPurity, True)
    graphPurity.GetXaxis().SetTitle(pTstr)
    graphPurity.GetYaxis().SetTitle("P")
    graphPurity.SetMarkerColor(ROOT.kRed-2)
    graphPurity.SetLineColor(ROOT.kRed-2)
    graphPurity.SetLineWidth(2)

    legPurity = ROOT.TLegend(*legPos)
    legPurity.AddEntry(legPurity,f"pp @ 13.6 TeV, |#eta|<{etaCut}", "")
    legPurity.AddEntry(graphPurity, legText, "p")
    legPurity.SetFillStyle(0)
    legPurity.SetBorderSize(0)
    legPurity.SetTextSize(0.03)

    output_dir.Delete()
    root_file.Close()
    root_file.Delete()

    return graphPurity, legPurity

def EfficiencyxMultPlot(a, b, path, graph_title, legText, ptBins = [0, 1, 2, 3, 4, 5, 6, 7, 8, 12]):
    """
    Apllies cut on |eta| < 0.9 and make plot of a/b vs pt
    Expected that a = reconstructed and b = generated
    """
    root_file = ROOT.TFile.Open(path)    
    output_dir_name = 'analysis-same-event-pairing/output;1'
    output_dir = root_file.Get(output_dir_name)
    if not output_dir:
        print(f"{output_dir_name} not found in the ROOT file.")

    aList = None
    bList = None
    for JPsiList in output_dir:
        # print("JPsiList: ", JPsiList.GetName())
        if JPsiList.GetName() == a:
            aList = JPsiList
        elif JPsiList.GetName() == b:
            bList = JPsiList
    if not aList:
        print(f"{a} not found in the output directory!")
    if not bList:
        print(f"{b} not found in the output directory!")

    df_a = None
    df_b = None

    for hist in aList:
        if "Truth" in a:
            if hist.GetName() == "EtaMC_PtMC":
                df_a = hist_to_df(hist, clean_names=True)
        else:
            if hist.GetName() == "Eta_Pt":
                df_a = hist_to_df(hist, clean_names=True)
    for hist in bList:
        if "Truth" in b:
            if hist.GetName() == "EtaMC_PtMC":
                df_b = hist_to_df(hist, clean_names=True)
        else:
            if hist.GetName() == "Eta_Pt":
                df_b = hist_to_df(hist, clean_names=True)
    
    # filter out |eta| > 0.9
    etaCut = 0.9
    etaCutPlus = etaCut + 0.0001# To account for floating-point precision issues
    df_aFilt = df_a[~((df_a["eta_left"] < -etaCutPlus) | (df_a["eta_right"] > etaCutPlus))]
    df_bFilt = df_b[~((df_b["eta_left"] < -etaCutPlus) | (df_b["eta_right"] > etaCutPlus))]

    df_aFiltReb = rebin_df(df_aFilt, ptBins, x_left=f"{pTstr}_left", x_right=f"{pTstr}_right", counts="counts")
    df_bFiltReb = rebin_df(df_bFilt, ptBins, x_left=f"{pTstr}_left", x_right=f"{pTstr}_right", counts="counts")

    df_eff = df_aFiltReb.copy()
    df_eff['counts'] = df_aFiltReb['counts'] / df_bFiltReb['counts']

    # Statistical Uncertainty
    eff_statUnc = np.sqrt(df_eff['counts'] * (1 - df_eff['counts']) / df_bFiltReb['counts']) # Binomial uncertainty
    graph = df_to_root_graph(df_eff, "efficiency_graph", graph_title, yerr=np.array(eff_statUnc))
    ROOT.SetOwnership(graph, True)
    graph.GetXaxis().SetTitle(pTstr)
    graph.GetYaxis().SetTitle("A x #epsilon")
    graph.SetMarkerColor(ROOT.kRed-2)
    graph.SetLineColor(ROOT.kRed-2)
    graph.SetLineWidth(2)

    legend = ROOT.TLegend(0.55, 0.7, 0.9, 0.9)
    legend.AddEntry(legend,f"pp @ 13.6 TeV, |#eta|<{etaCut}", "")
    legend.AddEntry(graph, legText, "p")
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.03)

    output_dir.Delete()
    root_file.Close()
    root_file.Delete()

    return graph, legend

# LEGEND_SIZE_DEFAULT = 0.02
# LEGEND_COORDS_DEFAULT = [0.2, 0.6, 0.4, 0.9]
# def createLegend(legLines, legCoords=LEGEND_COORDS_DEFAULT, textSize=LEGEND_SIZE_DEFAULT, isMC=False):
#     legend = TLegend(*legCoords)
#     legend.SetFillStyle(0)
#     legend.SetBorderSize(0)
#     legend.SetTextSize(textSize)
#     for legLine in legLines:
#         legend.AddEntry(legend, legLine, "")
#     legend.AddEntry(legend, "Data, |y| < 0.9", "")
#     if isMC:
#         legend.AddEntry(legend, "MC, |y| < 0.9", "")
#     return legend


# get_hist and join_dfs are initially from efficiencies_skimming.ipynb

def get_hist(filePath, rootDir_name, histName):

    """
    Get histogram from root file
    """
    root_file = ROOT.TFile.Open(filePath)
    rootDir = root_file.Get(rootDir_name) # can be e.g "dir/subdir/subdir"
    if not rootDir:
        print(f"{rootDir_name} not found in {filePath}.")
    mcSigHist = rootDir.FindObject(histName)
    if not mcSigHist:
        print(f"{histName} not found in {rootDir_name}")
    return mcSigHist

def join_dfs(dir_path, fileNames, rootDir_name, histName, datasetsNames=None):

    """
    Create Pandas dataframes from histograms from many root files and joins them into a single dataframe.
    """
    dfs = pd.DataFrame()
    for i,fileName in enumerate(fileNames):
        filePath = f"{dir_path}/{fileName}"
        hist = get_hist(filePath, rootDir_name, histName)
        df = hist_to_df(hist, isLabeledHist=True)
        if i==0:
            # dfs["XLabel"] = df["XLabel"]
            # dfs["YLabel"] = df["YLabel"]

            dfs["XLabel"] = df.iloc[:,0]
            if hist.ClassName().startswith("TH2"):
                dfs["YLabel"] = df.iloc[:,1]
            # print(dfs.to_string())
        if datasetsNames is not None:
             dfs[datasetsNames[i]] = df["counts"]
        else:
             dfs[fileName] = df["counts"]
    return dfs

def merge_TDirs(root_file, mergeDFs = True, target_dir_name="Merged"):
    """
    NOT FUNCTIONAL YET! Merge multiple TDirectories into a single TDirectory within the same ROOT file.
    """
    root_file = ROOT.TFile.Open(root_file, "UPDATE")
    target_dir = root_file.Get(target_dir_name)
    print("target dir type: ", type(target_dir))
    if not target_dir:
        target_dir = root_file.mkdir(target_dir_name)

    if mergeDFs:
        treesNames = []
        DFFound = False
        for keyDir in root_file.GetListOfKeys():
            dir =root_file.Get(keyDir.GetName())
            print("dir: ", dir)
            dirID = 0
            if keyDir.GetName().startswith("DF"):
                if dirID == 0:
                    for treeKey in dir.GetListOfKeys():
                        print("create TList?")
                print("DF found: ", dir.GetName())
                DFFound = True
                # target_dir.Merge(dir)
                for treeKey in dir.GetListOfKeys():
                    print("treeKey Name: ", treeKey.GetName())
                    tree = dir.Get(treeKey.GetName())
                    tree_clone = tree.Clone()
                    target_dir.cd()
                    tree_clone.Write()
        if not DFFound:
            print('No "DF..." directory found')
    #     target_dir.Write("", ROOT.TObject.kOverwrite)
    # else:
    #     print('Merging of non-"DF" objects is not implemented yet.')

def PIDEfficiencyManyMaps(mapsDirPath = "~/alice/Jpsi-Jets-Analysis/workDir/PIDEfficiency/PIDEfficiencyConverter/output/JpsiEffFromIdasMapsLHC25b14",
                          MCDataset="LHC25b14",
                          test=False,
                          weightedHist="hWeightedJpsiEffPtEta",
                          matchedHist="hMatchedJpsiPtEta"
                          ):
    base_dir = os.path.expanduser(mapsDirPath)
    root_files = sorted(glob.glob(f"{base_dir}/*.root"))
    if test:
        n = 3
        root_files = root_files[:n]

    ptBins = np.concatenate((np.arange(0, 10, 1),
                            np.arange(10, 22, 2)), axis=0)
    graphs = []

    for i, path in enumerate(root_files):
        # print(os.path.basename(path)[-8:-5])
        stdLabel = "-4 < n#sigma_{El} < 4, n#sigma_{#pi} > 2.5, n#sigma_{p} > 2.5"
        label = os.path.basename(path).replace("AnalysisResults_", "").replace(".root", "").replace("TrackBarrel_Conversions_withPID_", "")
        if "nSigmaEl-" in label:
            stdLabel = stdLabel.replace("-4 <", "-"+label[-3]+" <")
            stdLabel = stdLabel.replace("< 4", "< "+label[-1])
        if "nSigmaPi" in label:
            stdLabel = stdLabel.replace("{#pi} > 2.5", "{#pi} > "+label[-3:])
        if "nSigmaPr" in label:
            stdLabel = stdLabel.replace("{p} > 2.5", "{p} > "+label[-3:])
        label = stdLabel
        graph, _ = PIDEfficiencyPlot(
            weightedHist,
            matchedHist,
            path=path,
            graph_title="J/#psi PID Efficiency",
            legText=label,
            ptBins=ptBins
        )
        graph.SetTitle(f"J/#psi PID Efficiency - {label}")
        graph.SetMarkerStyle(20)        
        graphs.append((graph, label))
    return graphs

def averagePIDEfficiency(graphs): #graphs should be a tuple (graph, label)
    """
    Calculate graphs which is the average of all electron PID-efficiency maps (for comparision between datasets)
    """
    # Calculate graphs which is the average of all maps (for comparision between datasets)
    graphs_only = [g for g, _ in graphs]
    n_points = graphs_only[0].GetN()
    n_graphs = len(graphs_only)

    x_avg = []
    y_avg = []
    ex_avg = []
    ey_avg = []

    for i in range(n_points):
        xs = []
        ys = []
        eys = []
        exs = []

        for g in graphs_only:
            x = g.GetPointX(i)
            y = g.GetPointY(i)
            xs.append(float(x))
            ys.append(float(y))
            eys.append(g.GetErrorY(i))
            exs.append(g.GetErrorX(i))

        x_avg.append(xs[0])  # same binning
        # y_avg.append(np.mean(ys)) # simple average
        # simple_avg = np.mean(ys)
        weighted_avg = np.average(ys, weights=1/np.array(eys)**2)  # average weighted by error of each point
        # print(f"Simple average: {simple_avg}, Weighted average: {weighted_avg}")
        y_avg.append(weighted_avg)
        ex_avg.append(exs[0])

        # propagate independent errors
        # ey_avg.append(np.sqrt(np.sum(np.array(eys)**2)) / n_graphs) # simple error propagation for independent errors and simple average
        sumWeights = np.sum(1/np.array(eys)**2)
        stDevWeightedAvg = np.sqrt(1/sumWeights)  # error of the weighted average
        chi2 = np.sum(((ys - weighted_avg) / eys)**2)
        # stdDev = np.sqrt(np.sum((ys - weighted_avg)**2) / (n_graphs - 1)) # standard deviation (for simple avg, but also for weighted avg since each value is not a measurement of the same quantity)
        # stdErrorOfMean = stdDev / np.sqrt(n_graphs)
        birgeRatio = np.sqrt(chi2 / (n_graphs - 1)) # how much the dispersion differs from the expected. To account for big dispersion, since each value is not a measurement of the same quantity
        # print(f"Standard deviation of weighted average: {stDevWeightedAvg} \nBirge ratio: {birgeRatio}, \nFinal error: {stDevWeightedAvg * birgeRatio}, \nStandard error of mean: {stdErrorOfMean}")
        ey_avg.append(stDevWeightedAvg * birgeRatio)
        # ey_avg.append(stdErrorOfMean)


    avgGraph = ROOT.TGraphErrors(
        n_points,
        np.array(x_avg),
        np.array(y_avg),
        np.array(ex_avg),
        np.array(ey_avg)
    )
    avgGraph.GetXaxis().SetTitle(pTstr)
    avgGraph.GetYaxis().SetTitle("#epsilon_{PID}")
    return avgGraph

def plot_th2_with_python(a, path, output_dir_name='analysis-p-i-d-efficiency;1'):
    import matplotlib.pyplot as plt
    """
    Convert a ROOT TH2 histogram to a Python plot.
    """
    root_file = ROOT.TFile.Open(path)
    output_dir = root_file.Get(output_dir_name)
    th2 = output_dir.Get(a)
    nx = th2.GetNbinsX()
    ny = th2.GetNbinsY()

    xaxis = th2.GetXaxis()
    yaxis = th2.GetYaxis()

    # bin edges
    x_edges = np.array([xaxis.GetBinLowEdge(i+1) for i in range(nx)] +
                       [xaxis.GetBinUpEdge(nx)])
    y_edges = np.array([yaxis.GetBinLowEdge(i+1) for i in range(ny)] +
                       [yaxis.GetBinUpEdge(ny)])

    # bin contents
    z = np.zeros((nx, ny))
    for ix in range(1, nx+1):
        for iy in range(1, ny+1):
            z[ix-1, iy-1] = th2.GetBinContent(ix, iy)
            # if ix < 50 and z[ix-1, iy-1] > 0:
            if z[ix-1, iy-1] > 0:
                print(f"Bin ({ix}, {iy}): {z[ix-1, iy-1]}")                
    print(f"TH2 histogram has {nx} bins in x and {ny} bins in y.")
    z_min = np.min(z[z > 0])  # minimum non-zero value for better color scaling
    print(f"Minimum non-zero bin content: {z_min}")

    # matplotlib expects (ny, nx)
    z = z.T
    # z=z[:50,:50] # to focus on low pT region

    plt.figure()
    plt.pcolormesh(x_edges, y_edges, z, shading="auto")
    plt.xlabel(xaxis.GetTitle())
    plt.ylabel(yaxis.GetTitle())
    plt.colorbar(label="Counts")
    plt.show()
