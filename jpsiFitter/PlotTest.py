import ROOT
import numpy as np

# 1. Define the function (same as before)
def crystal_ball(x, params):
    A, B, C, D, E, F = params
    t = (x - A) / B
    if C < 0:
        t = -t
    
    absAlpha = abs(C)
    absAlpha2 = abs(E)
    
    if -absAlpha <= t < absAlpha2:
        return np.exp(-0.5 * t * t)
    
    if t < -absAlpha:
        a = (D/absAlpha)**D * np.exp(-0.5 * absAlpha**2)
        b = D/absAlpha - absAlpha
        return a / (b - t)**D
    
    if t >= absAlpha2:
        c = (F/absAlpha2)**F * np.exp(-0.5 * absAlpha2**2)
        d = F/absAlpha2 - absAlpha2
        return c / (d + t)**F
    
    return 0.0

# 2. ROOT function wrapper (same as before)
class CB2Pdf:
    def __call__(self, x, par):
        return crystal_ball(x[0], [par[0], par[1], par[2], par[3], par[4], par[5]])

# 3. Plot and save to ROOT file
def plot_and_save_to_root():
    # Create ROOT file for saving
    root_file = ROOT.TFile("CrystalBall.root", "RECREATE")
    
    # Create canvas and function
    c = ROOT.TCanvas("c", "Double Crystal Ball", 800, 600)
    f = ROOT.TF1("cb2", CB2Pdf(), -5, 5, 6)
    f.SetParameters(0, 1, 1, 2, 1, 2)
    f.SetParNames("A", "B", "C", "D", "E", "F")
    f.SetLineColor(ROOT.kBlue)
    f.SetLineWidth(2)
    f.SetTitle("Double Crystal Ball Function;x;Probability Density")
    
    # Draw and save to ROOT file
    f.Draw()
    c.Write()  # Save canvas to ROOT file
    f.Write()  # Save function to ROOT file
    
    # Close the file
    root_file.Close()
    print("Plot saved to CrystalBall.root")

# 4. Execute
plot_and_save_to_root()