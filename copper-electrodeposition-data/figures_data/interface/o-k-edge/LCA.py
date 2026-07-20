
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
import scipy.signal
import textwrap
from sklearn.linear_model import LinearRegression
import matplotlib
matplotlib.use("TkAgg")  # or "QtAgg" if installed
from matplotlib import font_manager
from sklearn.metrics import r2_score


arial = font_manager.FontProperties(fname="/mnt/c/Windows/Fonts/arial.ttf", size=9)
arialbd = font_manager.FontProperties(fname="/mnt/c/Windows/Fonts/arialbd.ttf", size=9)
# Constants

yoffset=1/6
E_i = -15 #eV
E_f = 17 #eV
eV_Ry = 13.60570398 

def getShift(filename):
    raw_data = np.loadtxt(filename, delimiter=",", dtype=str)
    # Extract headers (second column onward)
    headers = raw_data[0, 1:]  # ['energy', 'fermi']
    # Create the nested dictionary
    d = {
        row[0]: {headers[i]: float(row[i + 1]) for i in range(len(headers))}
        for row in raw_data[1:]
    }
    
    D_E = (d["bulkXCH/scf.out"]["energy"]-d["atomXCH/scf.out"]["energy"])-(d["bulkGS/scf.out"]["energy"]-d["atomGS/scf.out"]["energy"])
    # Fermi is given in QE in eV while the rest of energies in Ry. We convert accordingly
    shift = D_E*eV_Ry
    text=textwrap.dedent(f'''
    ##########################################################################
    {filename}
    --------------------------------------------------------------------------
    Fermi_bulkXCH = {d["bulkXCH/scf.out"]["fermi"]} eV
    total_shift = D_E + d_exp
    D_E = (E_bulkXCH - E_atomXCH) - (E_bulkGS - E_atomGS)
    D_E = {D_E*eV_Ry} eV
    E_bulkXCH = {d["bulkXCH/scf.out"]["energy"]*eV_Ry} eV
    E_atomXCH = {d["atomXCH/scf.out"]["energy"]*eV_Ry} eV
    E_bulkGS = {d["bulkGS/scf.out"]["energy"]*eV_Ry} eV
    E_atomGS = {d["atomGS/scf.out"]["energy"]*eV_Ry} eV
    ''')
    print(text)
    return shift
def exp(filename):
    data = np.loadtxt(filename)
    energy = data[:, 0]
    intensity = data[:, 1]
    return energy, intensity

def exp2(filename):
    data = np.loadtxt(filename)
    energy = data[:, 0]
    intensity = data[:, 1]
    # Check if energy is in decreasing order
    if energy[0] > energy[-1]:
        energy = energy[::-1]
        intensity = intensity[::-1]
    mask = (energy >= 525) & (energy <= 540)
    energy_range = energy[mask]
    intensity_range = intensity[mask]
    auc = simpson(intensity_range, energy_range)
    #auc = np.trapz(intensity_range, energy_range)
    intensity = intensity / auc
    return energy, intensity

def renorm(energy, intensity, Ei, Ef):
    mask = (energy >= Ei) & (energy <= Ef)
    energy_range = energy[mask]
    intensity_range = intensity[mask]
    auc = simpson(intensity_range, energy_range)
    intensity = intensity / auc
    return intensity

def raw(filename):
    data = np.loadtxt(filename)
    energy = data[:, 0]
    frames = len(data[0,:])-1
    intensity = np.sum(data[:, 1:], axis=1) / frames
    return energy, intensity

def norm_area(filename):
    data = np.loadtxt(filename)
    energy = data[:, 0]
    frames = len(data[0,:])-1
    intensity = np.sum(data[:, 1:], axis=1) / frames
    mask = (energy >= E_i) & (energy <= E_f)
    energy_range = energy[mask]
    intensity_range = intensity[mask]
    auc = simpson(intensity_range, energy_range)
    #auc = np.trapz(intensity_range, energy_range)
    intensity = intensity / auc
    return energy, intensity

def givePeak(energy, intensity, name):
    # Find peaks
    peaks, _ = scipy.signal.find_peaks(intensity, prominence=0.1)

    # Get the first peak energy
    if peaks.size > 0:
        first_peak_energy = energy[peaks[0]]
        first_peak_intensity = intensity[peaks[0]]
        #print(f"{name} peak: {first_peak_energy:.4f} eV")
    else:
        peaks, _ = scipy.signal.find_peaks(intensity)
        if peaks.size > 0:
            first_peak_energy = energy[peaks[0]]
            first_peak_intensity = intensity[peaks[0]]
            print("WARNING: Very small peaks.")
        else:
            print("No peaks where found")
    return first_peak_energy, first_peak_intensity




# Experimental reference data
energy_watbox_exp, intensity_watbox_exp = exp("exp_watbox.dat")

peaks_watbox_exp, _ = scipy.signal.find_peaks(intensity_watbox_exp, prominence=0.01)
x0_watbox_exp = energy_watbox_exp[peaks_watbox_exp[0]]
y0_watbox_exp = intensity_watbox_exp[peaks_watbox_exp[0]]
x1_watbox_exp = energy_watbox_exp[peaks_watbox_exp[1]]
y1_watbox_exp = intensity_watbox_exp[peaks_watbox_exp[1]]
# Computational reference
energy_watbox, intensity_watbox = norm_area("polAvg_watbox.dat")
x_watbox,y_watbox  = givePeak(energy=energy_watbox, intensity=intensity_watbox, name='watbox')
shift_watbox = getShift('david_watbox.csv')
# ------------------------------------------------------------------------------------------------------------
# Experimental shift
d_exp = x1_watbox_exp - (x_watbox+shift_watbox)
print(f'd_exp = {d_exp}')
# ------------------------------------------------------------------------------------------------------------
print('watbox total shift: ' + str(shift_watbox + d_exp))
energy_watbox += shift_watbox + d_exp
x_watbox,y_watbox  = givePeak(energy=energy_watbox, intensity=intensity_watbox, name='watbox')




# Process data

energy_0V_exp, intensity_0V_exp = exp2("exp_0V_trimmed.dat")
energy_minus01V_exp, intensity_minus01V_exp = exp2("exp_-0.1V_trimmed.dat")
#energy_minus02V_exp, intensity_minus02V_exp = exp2("exp_-0.2V_trimmed.dat")
energy_minus03V_exp, intensity_minus03V_exp = exp2("exp_-0.3V_trimmed.dat")


energy_Au_111_cuso4_sulfate_allFrames, intensity_Au_111_cuso4_sulfate_allFrames = norm_area("polAvg_Au_111_cuso4_sulfate_allFrames.dat")
shift_Au_111_cuso4_sulfate_allFrames = getShift('david_Au_111_cuso4_sulfate.csv')
print('Au_111_cuso4_sulfate_allFrames total_shift: ' + str(shift_Au_111_cuso4_sulfate_allFrames + d_exp))
energy_Au_111_cuso4_sulfate_allFrames += shift_Au_111_cuso4_sulfate_allFrames+ d_exp
peaks_Au_111_cuso4_sulfate_allFrames, _ = scipy.signal.find_peaks(intensity_Au_111_cuso4_sulfate_allFrames, prominence=0.01)
x0_Au_111_cuso4_sulfate_allFrames = energy_Au_111_cuso4_sulfate_allFrames[peaks_Au_111_cuso4_sulfate_allFrames[0]]
y0_Au_111_cuso4_sulfate_allFrames = intensity_Au_111_cuso4_sulfate_allFrames[peaks_Au_111_cuso4_sulfate_allFrames[0]]
x1_Au_111_cuso4_sulfate_allFrames = energy_Au_111_cuso4_sulfate_allFrames[peaks_Au_111_cuso4_sulfate_allFrames[1]]
y1_Au_111_cuso4_sulfate_allFrames = intensity_Au_111_cuso4_sulfate_allFrames[peaks_Au_111_cuso4_sulfate_allFrames[1]]


energy_cuso4_solv_sulfate_allFrames, intensity_cuso4_solv_sulfate_allFrames = norm_area("polAvg_cuso4_solv_sulfate_allFrames.dat")
shift_cuso4_solv_sulfate_allFrames = getShift('david_cuso4_solv_sulfate.csv')
print('cuso4_solv_sulfate_allFrames total_shift: ' + str(shift_cuso4_solv_sulfate_allFrames + d_exp))
energy_cuso4_solv_sulfate_allFrames += shift_cuso4_solv_sulfate_allFrames+ d_exp
peaks_cuso4_solv_sulfate_allFrames, _ = scipy.signal.find_peaks(intensity_cuso4_solv_sulfate_allFrames, prominence=0.01)
x0_cuso4_solv_sulfate_allFrames = energy_cuso4_solv_sulfate_allFrames[peaks_cuso4_solv_sulfate_allFrames[0]]
y0_cuso4_solv_sulfate_allFrames = intensity_cuso4_solv_sulfate_allFrames[peaks_cuso4_solv_sulfate_allFrames[0]]
#x1_cuso4_solv_sulfate_allFrames = energy_cuso4_solv_sulfate_allFrames[peaks_cuso4_solv_sulfate_allFrames[1]]
#y1_cuso4_solv_sulfate_allFrames = intensity_cuso4_solv_sulfate_allFrames[peaks_cuso4_solv_sulfate_allFrames[1]]


energy_2h3o_so4_sulfate, intensity_2h3o_so4_sulfate = norm_area("polAvg_2h3o_so4_sulfate.dat")
shift_2h3o_so4_sulfate = getShift('david_2h3o_so4_sulfate.csv')
print('2h3o_so4_sulfate total_shift: ' + str(shift_2h3o_so4_sulfate + d_exp))
energy_2h3o_so4_sulfate += shift_2h3o_so4_sulfate+ d_exp
peaks_2h3o_so4_sulfate, _ = scipy.signal.find_peaks(intensity_2h3o_so4_sulfate, prominence=0.01)
x0_2h3o_so4_sulfate = energy_2h3o_so4_sulfate[peaks_2h3o_so4_sulfate[0]]
y0_2h3o_so4_sulfate = intensity_2h3o_so4_sulfate[peaks_2h3o_so4_sulfate[0]]
x1_2h3o_so4_sulfate = energy_2h3o_so4_sulfate[peaks_2h3o_so4_sulfate[1]]
y1_2h3o_so4_sulfate = intensity_2h3o_so4_sulfate[peaks_2h3o_so4_sulfate[1]]

energy_2h3o_so4_sulfate_allFrames, intensity_2h3o_so4_sulfate_allFrames = norm_area("polAvg_2h3o_so4_sulfate_allFrames.dat")
shift_2h3o_so4_sulfate_allFrames = getShift('david_2h3o_so4_sulfate.csv')
print('2h3o_so4_sulfate_allFrames total_shift: ' + str(shift_2h3o_so4_sulfate_allFrames + d_exp))
energy_2h3o_so4_sulfate_allFrames += shift_2h3o_so4_sulfate_allFrames+ d_exp
peaks_2h3o_so4_sulfate_allFrames, _ = scipy.signal.find_peaks(intensity_2h3o_so4_sulfate_allFrames, prominence=0.01)
x0_2h3o_so4_sulfate_allFrames = energy_2h3o_so4_sulfate_allFrames[peaks_2h3o_so4_sulfate_allFrames[0]]
y0_2h3o_so4_sulfate_allFrames = intensity_2h3o_so4_sulfate_allFrames[peaks_2h3o_so4_sulfate_allFrames[0]]
#x1_2h3o_so4_sulfate_allFrames = energy_2h3o_so4_sulfate_allFrames[peaks_2h3o_so4_sulfate_allFrames[1]]
#y1_2h3o_so4_sulfate_allFrames = intensity_2h3o_so4_sulfate_allFrames[peaks_2h3o_so4_sulfate_allFrames[1]]


energy_Au_111_2h3o_so4_sulfate_allFrames, intensity_Au_111_2h3o_so4_sulfate_allFrames = norm_area("polAvg_Au_111_2h3o_so4_sulfate_allFrames.dat")
shift_Au_111_2h3o_so4_sulfate_allFrames = getShift('david_Au_111_2h3o_so4_sulfate.csv')
print('Au_111_2h3o_so4_sulfate_allFrames total_shift: ' + str(shift_Au_111_2h3o_so4_sulfate_allFrames + d_exp))
energy_Au_111_2h3o_so4_sulfate_allFrames += shift_Au_111_2h3o_so4_sulfate_allFrames+ d_exp
peaks_Au_111_2h3o_so4_sulfate_allFrames, _ = scipy.signal.find_peaks(intensity_Au_111_2h3o_so4_sulfate_allFrames, prominence=0.01)
x0_Au_111_2h3o_so4_sulfate_allFrames = energy_Au_111_2h3o_so4_sulfate_allFrames[peaks_Au_111_2h3o_so4_sulfate_allFrames[0]]
y0_Au_111_2h3o_so4_sulfate_allFrames = intensity_Au_111_2h3o_so4_sulfate_allFrames[peaks_Au_111_2h3o_so4_sulfate_allFrames[0]]
#x1_Au_111_2h3o_so4_sulfate_allFrames = energy_Au_111_2h3o_so4_sulfate_allFrames[peaks_Au_111_2h3o_so4_sulfate_allFrames[1]]
#y1_Au_111_2h3o_so4_sulfate_allFrames = intensity_Au_111_2h3o_so4_sulfate_allFrames[peaks_Au_111_2h3o_so4_sulfate_allFrames[1]]


energy_Au_111_2h3o_so4_228, intensity_Au_111_2h3o_so4_228 = norm_area("polAvg_Au_111_2h3o_so4_228.dat")
shift_Au_111_2h3o_so4_228 = getShift('david_Au_111_2h3o_so4_228_allFrames.csv') 
print('Au_111_2h3o_so4_228 water_shift: ' + str(shift_Au_111_2h3o_so4_228 + d_exp))
energy_Au_111_2h3o_so4_228 += shift_Au_111_2h3o_so4_228+ d_exp
peaks_Au_111_2h3o_so4_228, _ = scipy.signal.find_peaks(intensity_Au_111_2h3o_so4_228, prominence=0.01)
x0_Au_111_2h3o_so4_228 = energy_Au_111_2h3o_so4_228[peaks_Au_111_2h3o_so4_228[0]]
y0_Au_111_2h3o_so4_228 = intensity_Au_111_2h3o_so4_228[peaks_Au_111_2h3o_so4_228[0]]
#x_Au_111_2h3o_so4_228 = energy_Au_111_2h3o_so4_228[peaks_Au_111_2h3o_so4_228[1]]
#y_Au_111_2h3o_so4_228 = intensity_Au_111_2h3o_so4_228[peaks_Au_111_2h3o_so4_228[1]]

energy_2h3o_so4_hydronium29Zundel, intensity_2h3o_so4_hydronium29Zundel = norm_area("polAvg_2h3o_so4_hydronium29Zundel.dat")
shift_2h3o_so4_hydronium29Zundel = getShift('david_2h3o_so4_hydronium29Zundel.csv') 
print('2h3o_so4_hydronium29Zundel water_shift: ' + str(shift_2h3o_so4_hydronium29Zundel + d_exp))
energy_2h3o_so4_hydronium29Zundel += shift_2h3o_so4_hydronium29Zundel+ d_exp
peaks_2h3o_so4_hydronium29Zundel, _ = scipy.signal.find_peaks(intensity_2h3o_so4_hydronium29Zundel, prominence=0.01)
x0_2h3o_so4_hydronium29Zundel = energy_2h3o_so4_hydronium29Zundel[peaks_2h3o_so4_hydronium29Zundel[0]]
y0_2h3o_so4_hydronium29Zundel = intensity_2h3o_so4_hydronium29Zundel[peaks_2h3o_so4_hydronium29Zundel[0]]
x1_2h3o_so4_hydronium29Zundel = energy_2h3o_so4_hydronium29Zundel[peaks_2h3o_so4_hydronium29Zundel[1]]
y1_2h3o_so4_hydronium29Zundel = intensity_2h3o_so4_hydronium29Zundel[peaks_2h3o_so4_hydronium29Zundel[1]]

energy_target = energy_0V_exp
target = intensity_0V_exp

#energy_target = energy_minus01V_exp
#target = intensity_minus01V_exp

#energy_target = energy_minus03V_exp
#target = intensity_minus03V_exp

Ei = 530
Ef = 538

# Renormalize the requested reference spectra to the same area between Ei and Ef
intensity_2h3o_so4_sulfate_allFrames = renorm(energy_2h3o_so4_sulfate_allFrames, intensity_2h3o_so4_sulfate_allFrames, Ei, Ef)
intensity_2h3o_so4_hydronium29Zundel = renorm(energy_2h3o_so4_hydronium29Zundel, intensity_2h3o_so4_hydronium29Zundel, Ei, Ef)
intensity_Au_111_2h3o_so4_sulfate_allFrames = renorm(energy_Au_111_2h3o_so4_sulfate_allFrames, intensity_Au_111_2h3o_so4_sulfate_allFrames, Ei, Ef)
intensity_Au_111_2h3o_so4_228 = renorm(energy_Au_111_2h3o_so4_228, intensity_Au_111_2h3o_so4_228, Ei, Ef)
intensity_cuso4_solv_sulfate_allFrames = renorm(energy_cuso4_solv_sulfate_allFrames, intensity_cuso4_solv_sulfate_allFrames, Ei, Ef)
intensity_Au_111_cuso4_sulfate_allFrames = renorm(energy_Au_111_cuso4_sulfate_allFrames, intensity_Au_111_cuso4_sulfate_allFrames, Ei, Ef)

target = renorm(energy_target, target, Ei, Ef)

# helper: interpolate a reference spectrum onto the experimental energy grid
def interp(energy_from, spectrum, energy_to):
    return np.interp(energy_to, energy_from, spectrum)

s_2h3o_so4_sulfate_allFrames = interp(energy_2h3o_so4_sulfate_allFrames, intensity_2h3o_so4_sulfate_allFrames, energy_target)
s_2h3o_so4_hydronium29Zundel = interp(energy_2h3o_so4_hydronium29Zundel, intensity_2h3o_so4_hydronium29Zundel, energy_target)
s_Au_111_2h3o_so4_sulfate_allFrames = interp(energy_Au_111_2h3o_so4_sulfate_allFrames, intensity_Au_111_2h3o_so4_sulfate_allFrames, energy_target)
s_Au_111_2h3o_so4_228 = interp(energy_Au_111_2h3o_so4_228, intensity_Au_111_2h3o_so4_228, energy_target)
s_cuso4_solv_sulfate_allFrames = interp(energy_cuso4_solv_sulfate_allFrames, intensity_cuso4_solv_sulfate_allFrames, energy_target)
s_Au_111_cuso4_sulfate_allFrames = interp(energy_Au_111_cuso4_sulfate_allFrames, intensity_Au_111_cuso4_sulfate_allFrames, energy_target)

# Also include the experimental watbox spectrum as a reference
intensity_watbox_exp = renorm(energy_watbox_exp, intensity_watbox_exp, Ei, Ef)
s_watbox_exp = interp(energy_watbox_exp, intensity_watbox_exp, energy_target)

X = np.vstack([
    s_watbox_exp,
    s_2h3o_so4_sulfate_allFrames,
    #s_2h3o_so4_hydronium29Zundel,
    #s_Au_111_2h3o_so4_sulfate_allFrames,
    #s_Au_111_2h3o_so4_228,
    #s_cuso4_solv_sulfate_allFrames,
    #s_Au_111_cuso4_sulfate_allFrames,
]).T

# === Fit using non-negative least squares ===
model = LinearRegression(positive=True, fit_intercept=False)
model.fit(X, target)
raw_coeffs = model.coef_
coeffs = raw_coeffs / raw_coeffs.sum()  # normalize to sum = 1
print("Model raw coefficients:")
print(raw_coeffs)

# === Output ===
labels = [
    "exp. water",
    "(H$_3$O)$_2$SO$_4$ sulfate signal",
    #"2h3o_so4_hydronium29Zundel",
    #r"Au/(H$_3$O)$_2$SO$_4$ sulfate signal",
    #"Au_111_2h3o_so4_228",
    #"cuso4_solv_sulfate_allFrames",
    #"Au_111_cuso4_sulfate_allFrames",
]

print("Normalized LCA Coefficients:")
for lab, c in zip(labels, coeffs):
    print(f"{lab}: {c:.3f}")

fit = X @ coeffs
r2 = r2_score(target, fit)
print(f"R² score: {r2:.4f}")

# === Plot ===
plt.figure(figsize=(3.5, 3.5))
plt.plot(energy_target, target, label="Target", color='k')



plt.plot(energy_target, X @ model.coef_, label="Fit", linestyle='--', linewidth=2)
#custom_coeffs = [0.5, 0.1, 0.4]
#plt.plot(energy_target, X @ custom_coeffs, label="Fit", linestyle='--', linewidth=2)
for lab, s in zip(labels, [s_watbox_exp, 
                           s_2h3o_so4_sulfate_allFrames, 
                           #s_2h3o_so4_hydronium29Zundel, 
                           #s_Au_111_2h3o_so4_sulfate_allFrames, 
                           #s_Au_111_2h3o_so4_228, 
                           #s_cuso4_solv_sulfate_allFrames, 
                           #s_Au_111_cuso4_sulfate_allFrames
                           ]):
    plt.plot(energy_target, s, label=lab, linewidth=1)
plt.xlabel("Energy (eV)", fontproperties=arialbd)
plt.ylabel("Intensity (a.u.)", fontproperties=arialbd)
plt.xticks(fontproperties=arialbd)
plt.yticks(fontproperties=arialbd)
plt.legend(fontsize=7)
plt.title("0.0 V - LCA", fontproperties=arialbd)
plt.tight_layout()
plt.savefig('minus0V_LCA.svg', format='svg', bbox_inches='tight')
plt.show()


"""
from scipy.optimize import minimize

# Objective function: least-squares residual
def objective(coeffs):
    return np.sum((X @ coeffs - target) ** 2)

# Initial guess: 0.6 h3o, 0.1 so4, 0.3 h2o
x0 = np.array([0.6, 0.1, 0.3])

# Bounds: all coefficients ≥ 0
bounds = [(0, None)] * 3

# Optional: Enforce coefficients sum to 1
constraints = [{'type': 'eq', 'fun': lambda c: np.sum(c) - 1}]

res = minimize(objective, x0, bounds=bounds, constraints=constraints)


coeffs = res.x


coeffs = res.x / res.x.sum()               # normalize to sum = 1
print("lsq_linear raw coefficients:", res.x)
"""

"""
# TEST
# === Define energy grid ===
energy = np.linspace(525, 550, 101)  # 101 points from 525 to 550 eV

# === Synthetic reference spectra using Gaussian peaks ===
def gaussian(x, mu, sigma, amp):
    return amp * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

# Hydronium: sharp peak at 533 eV
hydronium = gaussian(energy, mu=533, sigma=0.8, amp=1.0)

# Sulfate: broader peak at 534 eV
sulfate = gaussian(energy, mu=534, sigma=1.5, amp=0.8)

# Water: broader and weaker peak at 535 eV
water = gaussian(energy, mu=535, sigma=2.0, amp=0.6)

# === Create target spectrum as a known mixture ===
# 50% hydronium, 30% sulfate, 20% water
target = 0.5 * hydronium + 0.3 * sulfate + 0.2 * water

def norm(y, x):
    return y / np.trapz(y, x)

h3o = norm(hydronium, energy)
so4 = norm(sulfate, energy)
h2o = norm(water, energy)
target = norm(target, energy)



X = np.vstack([h3o, so4, h2o]).T  # shape: (N, 3)

# === Fit using non-negative least squares ===
model = LinearRegression(positive=True)
model.fit(X, target)
coeffs = model.coef_ / model.coef_.sum()  # normalize to sum = 1

# === Output ===
print("Normalized LCA Coefficients:")
print(f"Hydronium: {coeffs[0]:.3f}")
print(f"Sulfate  : {coeffs[1]:.3f}")
print(f"Water    : {coeffs[2]:.3f}")

# === Plot ===
plt.plot(energy, target, label="Target", color='k')
plt.plot(energy, X @ coeffs, label="Fit", linestyle='--', linewidth=5)
plt.plot(energy, h3o, label="Hydronium")
plt.plot(energy, so4, label="Sulfate")
plt.plot(energy, h2o, label="Water")
plt.xlabel("Energy (eV)")
plt.ylabel("Intensity (a.u.)")
plt.legend()
plt.title("LCA using scikit-learn")
plt.tight_layout()
plt.show()

"""