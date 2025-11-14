import ROOT
import numpy as np
import pandas as pd
import ROOT
def test_utils():
 print("hist_utils.py imported successfully :)")
 print("In order to update it more dinamically, use 'reload' method from 'importlib' library.")

def clean_name(name):
    name = name.replace("MC", "")
    name = name.replace(" ", "")
    name = name.replace("(GeV/c)", "")
    name = name.replace("#", "")
    return name

def hist_to_df(hist, clean_names=True, isLabeledHist=False):  # For pyROOT histograms
    """
    Convert a TH1 or TH2 histogram to a pandas DataFrame.
    If clean_names=True, removes spaces and 'MC' from axis names.
    """

    if hist.ClassName().startswith("TH1"):
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
        
        if isLabeledHist:
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

            return pd.DataFrame(data={
                f"{x_name}XLabel": labels_x,
                f"{y_name}YLabel": labels_y,
                "counts": [row[4] for row in data]
            })

        columnsNames = [f"{x_name}_left", f"{x_name}_right", f"{y_name}_left", f"{y_name}_right", "counts"]
        return pd.DataFrame(data, columns=columnsNames)

    else:
        raise NotImplementedError(f"Only TH1 and TH2 supported, but got {hist.ClassName()}")
    
def rebin_df(df, new_edges, x_left="x_left", x_right="x_right", counts="counts"):
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
    # Apllies cut on |eta| < 0.9 and make plot of a/b vs pt
    # Expected that a = reconstructed and b = generated
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
    df_aFilt = df_a[~((df_a["eta_left"] < -0.9001) | (df_a["eta_right"] > 0.9001))] # 0.9001 due to floating-point error
    df_bFilt = df_b[~((df_b["eta_left"] < -0.9001) | (df_b["eta_right"] > 0.9001))]

    df_aFiltReb = rebin_df(df_aFilt, ptBins, x_left="p_{T}_left", x_right="p_{T}_right", counts="counts")
    df_bFiltReb = rebin_df(df_bFilt, ptBins, x_left="p_{T}_left", x_right="p_{T}_right", counts="counts")

    df_eff = df_aFiltReb.copy()
    df_eff['counts'] = df_aFiltReb['counts'] / df_bFiltReb['counts']

    # Statistical Uncertainty
    eff_statUnc = np.sqrt(df_eff['counts'] * (1 - df_eff['counts']) / df_bFiltReb['counts']) # Binomial uncertainty
    graph = df_to_root_graph(df_eff, "efficiency_graph", graph_title, yerr=np.array(eff_statUnc))
    ROOT.SetOwnership(graph, True)
    graph.GetXaxis().SetTitle("p_{T} (GeV/c)")
    graph.GetYaxis().SetTitle("A x #epsilon")
    graph.SetMarkerColor(ROOT.kRed-2)
    graph.SetLineColor(ROOT.kRed-2)
    graph.SetLineWidth(2)

    legend = ROOT.TLegend(0.55, 0.7, 0.9, 0.9)
    legend.AddEntry(legend,"pp @ 13.6 TeV", "")
    legend.AddEntry(graph, legText, "p")
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.03)

    output_dir.Delete()
    root_file.Close()
    root_file.Delete()

    return graph, legend

def PurityPlot(a, b, path, graph_title, legText, ptBins = [0, 1, 2, 3, 4, 5, 6, 7, 8, 12], legPos=[0.55, 0.7, 0.9, 0.9]):
    # Apllies cut on |eta| < 0.9 and make plot of a/(a+b) vs pt
    # Expected that a = correct association and b = incorrect association
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
    df_aFilt = df_a[~((df_a["eta_left"] < -0.9001) | (df_a["eta_right"] > 0.9001))] # 0.9001 due to floating-point error
    df_bFilt = df_b[~((df_b["eta_left"] < -0.9001) | (df_b["eta_right"] > 0.9001))]

    df_aFiltReb = rebin_df(df_aFilt, ptBins, x_left="p_{T}_left", x_right="p_{T}_right", counts="counts")
    df_bFiltReb = rebin_df(df_bFilt, ptBins, x_left="p_{T}_left", x_right="p_{T}_right", counts="counts")

    df_purity = df_aFiltReb.copy()
    nTotAssocs = df_aFiltReb['counts'] + df_bFiltReb['counts']
    df_purity['counts'] = df_aFiltReb['counts'] / nTotAssocs

    # Statistical Uncertainty
    purity_statUnc = np.sqrt(df_purity['counts'] * (1 - df_purity['counts']) / nTotAssocs) # Binomial uncertainty
    graphPurity = df_to_root_graph(df_purity, "purity_graph", graph_title, yerr=np.array(purity_statUnc))
    ROOT.SetOwnership(graphPurity, True)
    graphPurity.GetXaxis().SetTitle("p_{T} (GeV/c)")
    graphPurity.GetYaxis().SetTitle("P")
    graphPurity.SetMarkerColor(ROOT.kRed-2)
    graphPurity.SetLineColor(ROOT.kRed-2)
    graphPurity.SetLineWidth(2)

    legPurity = ROOT.TLegend(*legPos)
    legPurity.AddEntry(legPurity,"pp @ 13.6 TeV", "")
    legPurity.AddEntry(graphPurity, legText, "p")
    legPurity.SetFillStyle(0)
    legPurity.SetBorderSize(0)
    legPurity.SetTextSize(0.03)

    output_dir.Delete()
    root_file.Close()
    root_file.Delete()

    return graphPurity, legPurity

def EfficiencyxMultPlot(a, b, path, graph_title, legText, ptBins = [0, 1, 2, 3, 4, 5, 6, 7, 8, 12]):
    # Apllies cut on |eta| < 0.9 and make plot of a/b vs pt
    # Expected that a = reconstructed and b = generated
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
    df_aFilt = df_a[~((df_a["eta_left"] < -0.9001) | (df_a["eta_right"] > 0.9001))] # 0.9001 due to floating-point error
    df_bFilt = df_b[~((df_b["eta_left"] < -0.9001) | (df_b["eta_right"] > 0.9001))]

    df_aFiltReb = rebin_df(df_aFilt, ptBins, x_left="p_{T}_left", x_right="p_{T}_right", counts="counts")
    df_bFiltReb = rebin_df(df_bFilt, ptBins, x_left="p_{T}_left", x_right="p_{T}_right", counts="counts")

    df_eff = df_aFiltReb.copy()
    df_eff['counts'] = df_aFiltReb['counts'] / df_bFiltReb['counts']

    # Statistical Uncertainty
    eff_statUnc = np.sqrt(df_eff['counts'] * (1 - df_eff['counts']) / df_bFiltReb['counts']) # Binomial uncertainty
    graph = df_to_root_graph(df_eff, "efficiency_graph", graph_title, yerr=np.array(eff_statUnc))
    ROOT.SetOwnership(graph, True)
    graph.GetXaxis().SetTitle("p_{T} (GeV/c)")
    graph.GetYaxis().SetTitle("A x #epsilon")
    graph.SetMarkerColor(ROOT.kRed-2)
    graph.SetLineColor(ROOT.kRed-2)
    graph.SetLineWidth(2)

    legend = ROOT.TLegend(0.55, 0.7, 0.9, 0.9)
    legend.AddEntry(legend,"pp @ 13.6 TeV", "")
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