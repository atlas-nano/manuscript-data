
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
import scipy.signal
import textwrap
from scipy.interpolate import interp1d
from matplotlib import font_manager
import matplotlib
matplotlib.use("TkAgg")  # or "QtAgg" if installed

arial = font_manager.FontProperties(fname="/mnt/c/Windows/Fonts/arial.ttf", size=9)
arialbd = font_manager.FontProperties(fname="/mnt/c/Windows/Fonts/arialbd.ttf", size=9)
# Constants

yoffset=1/4
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
    intensity = intensity / auc
    return energy, intensity

def norm_area_allframes(filename):
    """Load a file with multiple frames (columns), normalize each frame by its area
    in the energy window [E_i, E_f], then return energy, mean intensity and std dev
    across frames. This preserves per-frame variation for shading in plots.
    """
    data = np.loadtxt(filename)
    energy = data[:, 0]
    #if data.ndim == 1:
    #    # single-column file (unexpected), fall back to norm_area behavior
    #    frames = 1
    #    intensities = data[:, 1:]
    #else:
    #    intensities = data[:, 1:]
    #    frames = intensities.shape[1]
    #
    mean_intensity = np.mean(data[:, 1:], axis=1)
    std_intensity = np.std(data[:, 1:], axis=1)
    # Mask for normalization window
    mask = (energy >= E_i) & (energy <= E_f)
    energy_range = energy[mask]
    mean_intensity_range = mean_intensity[mask]
    
    auc = simpson(mean_intensity_range, energy_range)
    print(f"auc: {auc}")
    print(f"std : {np.max(std_intensity)}")

    mean_intensity = mean_intensity / auc
    std_intensity = std_intensity / auc

    return energy, mean_intensity, std_intensity

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
"""
# Experimental reference data
energy_bulkCuO_exp, intensity_bulkCuO_exp = exp("exp_bulkCuO.dat")
# Reverse the order of this data
energy_bulkCuO_exp = energy_bulkCuO_exp[::-1]
intensity_bulkCuO_exp = intensity_bulkCuO_exp[::-1]
x_bulkCuO_exp,y_bulkCuO_exp  = givePeak(energy=energy_bulkCuO_exp, intensity=intensity_bulkCuO_exp, name='bulkCuO')
# COmputational reference
energy_bulkCuO, intensity_bulkCuO = norm_area("polAvg_bulkCuO.dat")
x_bulkCuO,y_bulkCuO  = givePeak(energy=energy_bulkCuO, intensity=intensity_bulkCuO, name='bulkCuO')
shift_bulkCuO = getShift('david_bulkCuO.csv')
# ------------------------------------------------------------------------------------------------------------
# Experimental shift
d_exp = x_bulkCuO_exp - (x_bulkCuO+shift_bulkCuO)
print(f'd_exp = {d_exp}')
# ------------------------------------------------------------------------------------------------------------
print('bulkCuO total shift: ' + str(shift_bulkCuO + d_exp))
energy_bulkCuO += shift_bulkCuO + d_exp
x_bulkCuO,y_bulkCuO  = givePeak(energy=energy_bulkCuO, intensity=intensity_bulkCuO, name='bulkCuO')
"""

# Experimental reference data
energy_watbox_exp, intensity_watbox_exp = exp("exp_watbox.dat")

peaks_watbox_exp, _ = scipy.signal.find_peaks(intensity_watbox_exp, prominence=0.01)
x0_watbox_exp = energy_watbox_exp[peaks_watbox_exp[0]]
y0_watbox_exp = intensity_watbox_exp[peaks_watbox_exp[0]]
x1_watbox_exp = energy_watbox_exp[peaks_watbox_exp[1]]
y1_watbox_exp = intensity_watbox_exp[peaks_watbox_exp[1]]
# Computational reference
# Use per-frame normalization for watbox so we can show mean +/- std shading
energy_watbox, intensity_watbox_mean, intensity_watbox_std = norm_area_allframes("polAvg_watbox.dat")

x_watbox,y_watbox  = givePeak(energy=energy_watbox, intensity=intensity_watbox_mean, name='watbox')
shift_watbox = getShift('david_watbox.csv')
# ------------------------------------------------------------------------------------------------------------
# Experimental shift
d_exp = x1_watbox_exp - (x_watbox+shift_watbox)
print(f'd_exp = {d_exp}')
# ------------------------------------------------------------------------------------------------------------
print('watbox total shift: ' + str(shift_watbox + d_exp))
energy_watbox += shift_watbox + d_exp
x_watbox += shift_watbox + d_exp
x0_watbox = 534.4
# Use mean when sampling annotation point
y0_watbox = intensity_watbox_mean[np.argmin(np.abs(energy_watbox - 534.4))]
x1_watbox = x_watbox
y1_watbox = y_watbox

# Experimental data

energy_0V_exp, intensity_0V_exp = exp("exp_0V.dat")
x_0V_exp, y_0V_exp = givePeak(energy=energy_0V_exp, intensity=intensity_0V_exp, name='0V')

energy_minus01V_exp, intensity_minus01V_exp = exp("exp_-0.1V.dat")
peaks_0V_exp, _ = scipy.signal.find_peaks(intensity_0V_exp, prominence=0.05)
x0_0V_exp = energy_0V_exp[peaks_0V_exp[0]]
y0_0V_exp = intensity_0V_exp[peaks_0V_exp[0]]
x1_0V_exp = energy_0V_exp[peaks_0V_exp[1]]
y1_0V_exp = intensity_0V_exp[peaks_0V_exp[1]]

energy_minus02V_exp, intensity_minus02V_exp = exp("exp_-0.2V.dat")
#peaks_minus02V_exp, _ = scipy.signal.find_peaks(intensity_minus02V_exp, prominence=0.15)
#x_minus02V_exp = energy_minus02V_exp[peaks_minus02V_exp[0]]
#y_minus02V_exp = intensity_minus02V_exp[peaks_minus02V_exp[0]]
x_minus02V_exp, y_minus02V_exp = givePeak(energy=energy_minus02V_exp, intensity=intensity_minus02V_exp, name='minus02V')

energy_minus03V_exp, intensity_minus03V_exp = exp("exp_-0.3V.dat")
x_minus03V_exp, y_minus03V_exp = givePeak(energy=energy_minus03V_exp, intensity=intensity_minus03V_exp, name='minus03V')


# Process data


energy_Au_111_cuso4_total, intensity_Au_111_cuso4_total = norm_area("polAvg_Au_111_cuso4_total.dat")
shift_Au_111_cuso4_total = getShift('david_Au_111_cuso4_total.csv')
print('Au_111_cuso4_total total_shift: ' + str(shift_Au_111_cuso4_total + d_exp))
energy_Au_111_cuso4_total += shift_Au_111_cuso4_total+ d_exp
x_Au_111_cuso4_total,y_Au_111_cuso4_total  = givePeak(energy=energy_Au_111_cuso4_total, intensity=intensity_Au_111_cuso4_total, name='Au_111_cuso4_total')


energy_Au_111_cuso4_sulfate, intensity_Au_111_cuso4_sulfate = norm_area("polAvg_Au_111_cuso4_sulfate.dat")
shift_Au_111_cuso4_sulfate = getShift('david_Au_111_cuso4_sulfate.csv')
print('Au_111_cuso4_sulfate total_shift: ' + str(shift_Au_111_cuso4_sulfate + d_exp))
energy_Au_111_cuso4_sulfate += shift_Au_111_cuso4_sulfate+ d_exp
x_Au_111_cuso4_sulfate,y_Au_111_cuso4_sulfate  = givePeak(energy=energy_Au_111_cuso4_sulfate, intensity=intensity_Au_111_cuso4_sulfate, name='Au_111_cuso4_sulfate')

energy_Au_111_cuso4_sulfate_allFrames, intensity_Au_111_cuso4_sulfate_allFrames_mean, intensity_Au_111_cuso4_sulfate_allFrames_std = norm_area_allframes("polAvg_Au_111_cuso4_sulfate_allFrames.dat")
shift_Au_111_cuso4_sulfate_allFrames = getShift('david_Au_111_cuso4_sulfate.csv')
print('Au_111_cuso4_sulfate_allFrames total_shift: ' + str(shift_Au_111_cuso4_sulfate_allFrames + d_exp))
energy_Au_111_cuso4_sulfate_allFrames += shift_Au_111_cuso4_sulfate_allFrames+ d_exp
peaks_Au_111_cuso4_sulfate_allFrames, _ = scipy.signal.find_peaks(intensity_Au_111_cuso4_sulfate_allFrames_mean, prominence=0.01)
x0_Au_111_cuso4_sulfate_allFrames = energy_Au_111_cuso4_sulfate_allFrames[peaks_Au_111_cuso4_sulfate_allFrames[0]]
y0_Au_111_cuso4_sulfate_allFrames = intensity_Au_111_cuso4_sulfate_allFrames_mean[peaks_Au_111_cuso4_sulfate_allFrames[0]]
#x1_Au_111_cuso4_sulfate_allFrames = energy_Au_111_cuso4_sulfate_allFrames[peaks_Au_111_cuso4_sulfate_allFrames[1]]
#y1_Au_111_cuso4_sulfate_allFrames = intensity_Au_111_cuso4_sulfate_allFrames_mean[peaks_Au_111_cuso4_sulfate_allFrames[1]]

energy_cuso4_solv_sulfate, intensity_cuso4_solv_sulfate = norm_area("polAvg_cuso4_solv_sulfate.dat")
shift_cuso4_solv_sulfate = getShift('david_cuso4_solv_sulfate.csv')
print('cuso4_solv_sulfate total_shift: ' + str(shift_cuso4_solv_sulfate + d_exp))
energy_cuso4_solv_sulfate += shift_cuso4_solv_sulfate+ d_exp
x_cuso4_solv_sulfate,y_cuso4_solv_sulfate  = givePeak(energy=energy_cuso4_solv_sulfate, intensity=intensity_cuso4_solv_sulfate, name='cuso4_solv_sulfate')

energy_cuso4_solv_sulfate_allFrames, intensity_cuso4_solv_sulfate_allFrames_mean, intensity_cuso4_solv_sulfate_allFrames_std = norm_area_allframes("polAvg_cuso4_solv_sulfate_allFrames.dat")
shift_cuso4_solv_sulfate_allFrames = getShift('david_cuso4_solv_sulfate.csv')
print('cuso4_solv_sulfate_allFrames total_shift: ' + str(shift_cuso4_solv_sulfate_allFrames + d_exp))
energy_cuso4_solv_sulfate_allFrames += shift_cuso4_solv_sulfate_allFrames+ d_exp
peaks_cuso4_solv_sulfate_allFrames, _ = scipy.signal.find_peaks(intensity_cuso4_solv_sulfate_allFrames_mean, prominence=0.01)
x0_cuso4_solv_sulfate_allFrames = energy_cuso4_solv_sulfate_allFrames[peaks_cuso4_solv_sulfate_allFrames[0]]
y0_cuso4_solv_sulfate_allFrames = intensity_cuso4_solv_sulfate_allFrames_mean[peaks_cuso4_solv_sulfate_allFrames[0]]
#x1_cuso4_solv_sulfate_allFrames = energy_cuso4_solv_sulfate_allFrames[peaks_cuso4_solv_sulfate_allFrames[1]]
#y1_cuso4_solv_sulfate_allFrames = intensity_cuso4_solv_sulfate_allFrames_mean[peaks_cuso4_solv_sulfate_allFrames[1]]

energy_cuso4_vac, intensity_cuso4_vac = norm_area("polAvg_cuso4_vac.dat")
shift_cuso4_vac = getShift('david_cuso4_vac.csv')
print('cuso4_vac total_shift: ' + str(shift_cuso4_vac + d_exp))
energy_cuso4_vac += shift_cuso4_vac+ d_exp
x_cuso4_vac,y_cuso4_vac  = givePeak(energy=energy_cuso4_vac, intensity=intensity_cuso4_vac, name='cuso4_vac')

energy_so4_vac, intensity_so4_vac = norm_area("polAvg_so4_vac.dat")
shift_so4_vac = getShift('david_so4_vac.csv')
print('so4_vac total_shift: ' + str(shift_so4_vac + d_exp))
energy_so4_vac += shift_so4_vac+ d_exp
x_so4_vac,y_so4_vac  = givePeak(energy=energy_so4_vac, intensity=intensity_so4_vac, name='so4_vac')


energy_2h3o_so4_sulfate, intensity_2h3o_so4_sulfate = norm_area("polAvg_2h3o_so4_sulfate.dat")
shift_2h3o_so4_sulfate = getShift('david_2h3o_so4_sulfate.csv')
print('2h3o_so4_sulfate total_shift: ' + str(shift_2h3o_so4_sulfate + d_exp))
energy_2h3o_so4_sulfate += shift_2h3o_so4_sulfate+ d_exp
peaks_2h3o_so4_sulfate, _ = scipy.signal.find_peaks(intensity_2h3o_so4_sulfate, prominence=0.01)
x0_2h3o_so4_sulfate = energy_2h3o_so4_sulfate[peaks_2h3o_so4_sulfate[0]]
y0_2h3o_so4_sulfate = intensity_2h3o_so4_sulfate[peaks_2h3o_so4_sulfate[0]]
x1_2h3o_so4_sulfate = energy_2h3o_so4_sulfate[peaks_2h3o_so4_sulfate[1]]
y1_2h3o_so4_sulfate = intensity_2h3o_so4_sulfate[peaks_2h3o_so4_sulfate[1]]

energy_2h3o_so4_sulfate_allFrames, intensity_2h3o_so4_sulfate_allFrames_mean, intensity_2h3o_so4_sulfate_allFrames_std = norm_area_allframes("polAvg_2h3o_so4_sulfate_allFrames.dat")
shift_2h3o_so4_sulfate_allFrames = getShift('david_2h3o_so4_sulfate.csv')
print('2h3o_so4_sulfate_allFrames total_shift: ' + str(shift_2h3o_so4_sulfate_allFrames + d_exp))
energy_2h3o_so4_sulfate_allFrames += shift_2h3o_so4_sulfate_allFrames+ d_exp
peaks_2h3o_so4_sulfate_allFrames, _ = scipy.signal.find_peaks(intensity_2h3o_so4_sulfate_allFrames_mean, prominence=0.01)
x0_2h3o_so4_sulfate_allFrames = energy_2h3o_so4_sulfate_allFrames[peaks_2h3o_so4_sulfate_allFrames[0]]
y0_2h3o_so4_sulfate_allFrames = intensity_2h3o_so4_sulfate_allFrames_mean[peaks_2h3o_so4_sulfate_allFrames[0]]
#x1_2h3o_so4_sulfate_allFrames = energy_2h3o_so4_sulfate_allFrames[peaks_2h3o_so4_sulfate_allFrames[1]]
#y1_2h3o_so4_sulfate_allFrames = intensity_2h3o_so4_sulfate_allFrames_mean[peaks_2h3o_so4_sulfate_allFrames[1]]


energy_Au_111_2h3o_so4_allFrames_waterAndHydronium, intensity_Au_111_2h3o_so4_allFrames_waterAndHydronium = norm_area("polAvg_Au_111_2h3o_so4_allFrames_waterAndHydronium.dat")
shift_Au_111_2h3o_so4_allFrames_waterAndHydronium = getShift('david_Au_111_2h3o_so4_hydronium256.csv')
print('Au_111_2h3o_so4_allFrames_waterAndHydronium total_shift: ' + str(shift_Au_111_2h3o_so4_allFrames_waterAndHydronium + d_exp))
energy_Au_111_2h3o_so4_allFrames_waterAndHydronium += shift_Au_111_2h3o_so4_allFrames_waterAndHydronium+ d_exp
peaks_Au_111_2h3o_so4_allFrames_waterAndHydronium, _ = scipy.signal.find_peaks(intensity_Au_111_2h3o_so4_allFrames_waterAndHydronium, prominence=0.01)
x0_Au_111_2h3o_so4_allFrames_waterAndHydronium = energy_Au_111_2h3o_so4_allFrames_waterAndHydronium[peaks_Au_111_2h3o_so4_allFrames_waterAndHydronium[0]]
y0_Au_111_2h3o_so4_allFrames_waterAndHydronium = intensity_Au_111_2h3o_so4_allFrames_waterAndHydronium[peaks_Au_111_2h3o_so4_allFrames_waterAndHydronium[0]]
#x1_Au_111_2h3o_so4_allFrames_waterAndHydronium = energy_Au_111_2h3o_so4_allFrames_waterAndHydronium[peaks_Au_111_2h3o_so4_allFrames_waterAndHydronium[1]]
#y1_Au_111_2h3o_so4_allFrames_waterAndHydronium = intensity_Au_111_2h3o_so4_allFrames_waterAndHydronium[peaks_Au_111_2h3o_so4_allFrames_waterAndHydronium[1]]

energy_Au_111_2h3o_so4_sulfate, intensity_Au_111_2h3o_so4_sulfate = norm_area("polAvg_Au_111_2h3o_so4_sulfate.dat")
shift_Au_111_2h3o_so4_sulfate = getShift('david_Au_111_2h3o_so4_sulfate.csv')
print('Au_111_2h3o_so4_sulfate total_shift: ' + str(shift_Au_111_2h3o_so4_sulfate + d_exp))
energy_Au_111_2h3o_so4_sulfate += shift_Au_111_2h3o_so4_sulfate+ d_exp
peaks_Au_111_2h3o_so4_sulfate, _ = scipy.signal.find_peaks(intensity_Au_111_2h3o_so4_sulfate, prominence=0.01)
x0_Au_111_2h3o_so4_sulfate = energy_Au_111_2h3o_so4_sulfate[peaks_Au_111_2h3o_so4_sulfate[0]]
y0_Au_111_2h3o_so4_sulfate = intensity_Au_111_2h3o_so4_sulfate[peaks_Au_111_2h3o_so4_sulfate[0]]
x1_Au_111_2h3o_so4_sulfate = energy_Au_111_2h3o_so4_sulfate[peaks_Au_111_2h3o_so4_sulfate[1]]
y1_Au_111_2h3o_so4_sulfate = intensity_Au_111_2h3o_so4_sulfate[peaks_Au_111_2h3o_so4_sulfate[1]]

energy_Au_111_2h3o_so4_sulfate_allFrames, intensity_Au_111_2h3o_so4_sulfate_allFrames_mean, intensity_Au_111_2h3o_so4_sulfate_allFrames_std = norm_area_allframes("polAvg_Au_111_2h3o_so4_sulfate_allFrames.dat")
shift_Au_111_2h3o_so4_sulfate_allFrames = getShift('david_Au_111_2h3o_so4_sulfate.csv')
print('Au_111_2h3o_so4_sulfate_allFrames total_shift: ' + str(shift_Au_111_2h3o_so4_sulfate_allFrames + d_exp))
energy_Au_111_2h3o_so4_sulfate_allFrames += shift_Au_111_2h3o_so4_sulfate_allFrames+ d_exp
peaks_Au_111_2h3o_so4_sulfate_allFrames, _ = scipy.signal.find_peaks(intensity_Au_111_2h3o_so4_sulfate_allFrames_mean, prominence=0.01)
x0_Au_111_2h3o_so4_sulfate_allFrames = energy_Au_111_2h3o_so4_sulfate_allFrames[peaks_Au_111_2h3o_so4_sulfate_allFrames[0]]
y0_Au_111_2h3o_so4_sulfate_allFrames = intensity_Au_111_2h3o_so4_sulfate_allFrames_mean[peaks_Au_111_2h3o_so4_sulfate_allFrames[0]]
#x1_Au_111_2h3o_so4_sulfate_allFrames = energy_Au_111_2h3o_so4_sulfate_allFrames[peaks_Au_111_2h3o_so4_sulfate_allFrames[1]]
#y1_Au_111_2h3o_so4_sulfate_allFrames = intensity_Au_111_2h3o_so4_sulfate_allFrames_mean[peaks_Au_111_2h3o_so4_sulfate_allFrames[1]]

energy_Au_111_cuso4_water, intensity_Au_111_cuso4_water = norm_area("polAvg_Au_111_cuso4_water.dat")
shift_Au_111_cuso4_water = getShift('david_Au_111_cuso4_water.csv')
print('Au_111_cuso4_water water_shift: ' + str(shift_Au_111_cuso4_water + d_exp))
energy_Au_111_cuso4_water += shift_Au_111_cuso4_water+ d_exp
x_Au_111_cuso4_water,y_Au_111_cuso4_water  = givePeak(energy=energy_Au_111_cuso4_water, intensity=intensity_Au_111_cuso4_water, name='Au_111_cuso4_water')

energy_Au_111_2h3o_so4_hydronium256, intensity_Au_111_2h3o_so4_hydronium256 = norm_area("polAvg_Au_111_2h3o_so4_hydronium256.dat")
shift_Au_111_2h3o_so4_hydronium256 = getShift('david_Au_111_2h3o_so4_hydronium256.csv')
print('Au_111_2h3o_so4_hydronium256 water_shift: ' + str(shift_Au_111_2h3o_so4_hydronium256 + d_exp))
energy_Au_111_2h3o_so4_hydronium256 += shift_Au_111_2h3o_so4_hydronium256+ d_exp

peaks_Au_111_2h3o_so4_hydronium256, _ = scipy.signal.find_peaks(intensity_Au_111_2h3o_so4_hydronium256, prominence=0.01)
x0_Au_111_2h3o_so4_hydronium256 = energy_Au_111_2h3o_so4_hydronium256[peaks_Au_111_2h3o_so4_hydronium256[0]]
y0_Au_111_2h3o_so4_hydronium256 = intensity_Au_111_2h3o_so4_hydronium256[peaks_Au_111_2h3o_so4_hydronium256[0]]
x_Au_111_2h3o_so4_hydronium256 = energy_Au_111_2h3o_so4_hydronium256[peaks_Au_111_2h3o_so4_hydronium256[1]]
y_Au_111_2h3o_so4_hydronium256 = intensity_Au_111_2h3o_so4_hydronium256[peaks_Au_111_2h3o_so4_hydronium256[1]]

energy_Au_111_2h3o_so4_hydronium_allFrames, intensity_Au_111_2h3o_so4_hydronium_allFrames_mean, intensity_Au_111_2h3o_so4_hydronium_allFrames_std = norm_area_allframes("polAvg_Au_111_2h3o_so4_hydronium_allFrames.dat")
shift_Au_111_2h3o_so4_hydronium_allFrames = getShift('david_Au_111_2h3o_so4_hydronium256.csv')
print('Au_111_2h3o_so4_hydronium_allFrames water_shift: ' + str(shift_Au_111_2h3o_so4_hydronium_allFrames + d_exp))
energy_Au_111_2h3o_so4_hydronium_allFrames += shift_Au_111_2h3o_so4_hydronium_allFrames+ d_exp

peaks_Au_111_2h3o_so4_hydronium_allFrames, _ = scipy.signal.find_peaks(intensity_Au_111_2h3o_so4_hydronium_allFrames_mean, prominence=0.01)
x0_Au_111_2h3o_so4_hydronium_allFrames = energy_Au_111_2h3o_so4_hydronium_allFrames[peaks_Au_111_2h3o_so4_hydronium_allFrames[0]]
y0_Au_111_2h3o_so4_hydronium_allFrames = intensity_Au_111_2h3o_so4_hydronium_allFrames_mean[peaks_Au_111_2h3o_so4_hydronium_allFrames[0]]
#x_Au_111_2h3o_so4_hydronium_allFrames = energy_Au_111_2h3o_so4_hydronium_allFrames[peaks_Au_111_2h3o_so4_hydronium_allFrames[1]]
#y_Au_111_2h3o_so4_hydronium_allFrames = intensity_Au_111_2h3o_so4_hydronium_allFrames_mean[peaks_Au_111_2h3o_so4_hydronium_allFrames[1]]

# *********** THIS IS THE MOST UPDATED HYDRONIUM ALL FRAMES ONLY LOOKING AT O 228 ***************
energy_Au_111_2h3o_so4_228, intensity_Au_111_2h3o_so4_228_mean, intensity_Au_111_2h3o_so4_228_std = norm_area_allframes("polAvg_Au_111_2h3o_so4_228.dat")
shift_Au_111_2h3o_so4_228 = getShift('david_Au_111_2h3o_so4_228_allFrames.csv') 
print('Au_111_2h3o_so4_228 water_shift: ' + str(shift_Au_111_2h3o_so4_228 + d_exp))
energy_Au_111_2h3o_so4_228 += shift_Au_111_2h3o_so4_228+ d_exp
peaks_Au_111_2h3o_so4_228, _ = scipy.signal.find_peaks(intensity_Au_111_2h3o_so4_228_mean, prominence=0.01)
x0_Au_111_2h3o_so4_228 = energy_Au_111_2h3o_so4_228[peaks_Au_111_2h3o_so4_228[0]]
y0_Au_111_2h3o_so4_228 = intensity_Au_111_2h3o_so4_228_mean[peaks_Au_111_2h3o_so4_228[0]]
#x_Au_111_2h3o_so4_228 = energy_Au_111_2h3o_so4_228[peaks_Au_111_2h3o_so4_228[1]]
#y_Au_111_2h3o_so4_228 = intensity_Au_111_2h3o_so4_228_mean[peaks_Au_111_2h3o_so4_228[1]]

energy_2h3o_so4_hydronium29Zundel, intensity_2h3o_so4_hydronium29Zundel = norm_area("polAvg_2h3o_so4_hydronium29Zundel.dat")
shift_2h3o_so4_hydronium29Zundel = getShift('david_2h3o_so4_hydronium29Zundel.csv') 
print('2h3o_so4_hydronium29Zundel water_shift: ' + str(shift_2h3o_so4_hydronium29Zundel + d_exp))
energy_2h3o_so4_hydronium29Zundel += shift_2h3o_so4_hydronium29Zundel+ d_exp
peaks_2h3o_so4_hydronium29Zundel, _ = scipy.signal.find_peaks(intensity_2h3o_so4_hydronium29Zundel, prominence=0.01)
x0_2h3o_so4_hydronium29Zundel = energy_2h3o_so4_hydronium29Zundel[peaks_2h3o_so4_hydronium29Zundel[0]]
y0_2h3o_so4_hydronium29Zundel = intensity_2h3o_so4_hydronium29Zundel[peaks_2h3o_so4_hydronium29Zundel[0]]
x1_2h3o_so4_hydronium29Zundel = energy_2h3o_so4_hydronium29Zundel[peaks_2h3o_so4_hydronium29Zundel[1]]
y1_2h3o_so4_hydronium29Zundel = intensity_2h3o_so4_hydronium29Zundel[peaks_2h3o_so4_hydronium29Zundel[1]]

# Use per-frame normalization so we can plot mean +/- std shading for this dataset
energy_2h3o_so4_hydronium_allFrames, intensity_2h3o_so4_hydronium_allFrames_mean, intensity_2h3o_so4_hydronium_allFrames_std = norm_area_allframes("polAvg_2h3o_so4_hydronium_allFrames.dat")
shift_2h3o_so4_hydronium_allFrames = getShift('david_2h3o_so4_hydronium29Zundel.csv') 
print('2h3o_so4_hydronium_allFrames water_shift: ' + str(shift_2h3o_so4_hydronium_allFrames + d_exp))
energy_2h3o_so4_hydronium_allFrames += shift_2h3o_so4_hydronium_allFrames+ d_exp
peaks_2h3o_so4_hydronium_allFrames, _ = scipy.signal.find_peaks(intensity_2h3o_so4_hydronium_allFrames_mean, prominence=0.01)
x0_2h3o_so4_hydronium_allFrames = energy_2h3o_so4_hydronium_allFrames[peaks_2h3o_so4_hydronium_allFrames[0]]
y0_2h3o_so4_hydronium_allFrames = intensity_2h3o_so4_hydronium_allFrames_mean[peaks_2h3o_so4_hydronium_allFrames[0]]
#x1_2h3o_so4_hydronium_allFrames = energy_2h3o_so4_hydronium_allFrames[peaks_2h3o_so4_hydronium_allFrames[1]]
#y1_2h3o_so4_hydronium_allFrames = intensity_2h3o_so4_hydronium_allFrames[peaks_2h3o_so4_hydronium_allFrames[1]]


# Create a common energy grid (you can adjust resolution as needed)
energy_common = np.linspace(min(energy_Au_111_cuso4_sulfate.min(), energy_Au_111_cuso4_water.min(), energy_Au_111_2h3o_so4_hydronium256.min()), 
                            max(energy_Au_111_cuso4_sulfate.max(), energy_Au_111_cuso4_water.max(), energy_Au_111_2h3o_so4_hydronium_allFrames.max()), 
                            num=1000)

# Interpolate both spectra onto the common grid
interp_sulfate_func = interp1d(energy_Au_111_cuso4_sulfate, intensity_Au_111_cuso4_sulfate, kind='linear', bounds_error=False, fill_value=0)
interp_water_func = interp1d(energy_Au_111_cuso4_water, intensity_Au_111_cuso4_water, kind='linear', bounds_error=False, fill_value=0)
interp_hydronium_func = interp1d(energy_Au_111_2h3o_so4_hydronium_allFrames, intensity_Au_111_2h3o_so4_hydronium_allFrames_mean, kind='linear', bounds_error=False, fill_value=0)

intensity_sulfate_interp = interp_sulfate_func(energy_common)
intensity_water_interp = interp_water_func(energy_common)
intensity_hydronium_interp = interp_hydronium_func(energy_common)

#weight_sulfate= 0.10
#weight_water= 0.60
#weight_hydronium=  0.30

weight_sulfate= 1
weight_water= 1
weight_hydronium=  1
# Sum the interpolated intensities
intensity_combined = intensity_sulfate_interp*weight_sulfate + intensity_water_interp*weight_water + intensity_hydronium_interp*weight_hydronium



# Optional: renormalize the total spectrum
#auc_total = simpson(intensity_total, energy_common)
#intensity_total /= auc_total

# Final result
energy_Au_111_cuso4_combined = energy_common
intensity_Au_111_cuso4_combined = intensity_combined
x_Au_111_cuso4_combined,y_Au_111_cuso4_combined  = givePeak(energy=energy_Au_111_cuso4_combined, intensity=intensity_Au_111_cuso4_combined, name='Au_111_cuso4_total')


plt.figure(figsize=(5.0, 3.5))

plt.plot(energy_0V_exp, intensity_0V_exp/5 + yoffset*-2, label='exp 0 V', linewidth=2, color='#515151')
if x0_0V_exp is not None and y0_0V_exp is not None:
    plt.annotate(f'{x0_0V_exp:.1f}', (x0_0V_exp-1, y0_0V_exp/5+yoffset*-2),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black')
if x1_0V_exp is not None and y_0V_exp is not None:
    plt.annotate(f'{x1_0V_exp:.1f}', (x1_0V_exp, y_0V_exp/5+yoffset*-2),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black')

#plt.plot(energy_Au_111_2h3o_so4_hydronium_allFrames, intensity_Au_111_2h3o_so4_hydronium_allFrames*weight_hydronium+yoffset*0, label='Au_111_2h3o_so4_hydronium_allFrames', linewidth=2, color='#515151')
#if x0_Au_111_2h3o_so4_hydronium_allFrames is not None and y0_Au_111_2h3o_so4_hydronium_allFrames is not None:
#    plt.annotate(f'{x0_Au_111_2h3o_so4_hydronium_allFrames:.1f}', (x0_Au_111_2h3o_so4_hydronium_allFrames, y0_Au_111_2h3o_so4_hydronium_allFrames*weight_hydronium+yoffset*0),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black')
#if x1_Au_111_2h3o_so4_hydronium_allFrames is not None and y1_Au_111_2h3o_so4_hydronium_allFrames is not None:
#    plt.annotate(f'{x1_Au_111_2h3o_so4_hydronium_allFrames:.1f}', (x1_Au_111_2h3o_so4_hydronium_allFrames, y1_Au_111_2h3o_so4_hydronium_allFrames*weight_hydronium+yoffset*0),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black')

#plt.plot(energy_Au_111_2h3o_so4_hydronium256, intensity_Au_111_2h3o_so4_hydronium256*weight_hydronium+yoffset*0, label='Au_111_2h3o_so4_hydronium256', linewidth=2, color='purple')
#if x0_Au_111_2h3o_so4_hydronium256 is not None and y0_Au_111_2h3o_so4_hydronium256 is not None:
#    plt.annotate(f'{x0_Au_111_2h3o_so4_hydronium256:.1f}', (x0_Au_111_2h3o_so4_hydronium256, y0_Au_111_2h3o_so4_hydronium256*weight_hydronium+yoffset*0),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black')
#if x_Au_111_2h3o_so4_hydronium256 is not None and y_Au_111_2h3o_so4_hydronium256 is not None:
#    plt.annotate(f'{x_Au_111_2h3o_so4_hydronium256:.1f}', (x_Au_111_2h3o_so4_hydronium256, y_Au_111_2h3o_so4_hydronium256*weight_hydronium+yoffset*0),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black')
mean_Au_111_2h3o_so4_228 = intensity_Au_111_2h3o_so4_228_mean * weight_hydronium * 2 + yoffset*0
std_Au_111_2h3o_so4_228 = intensity_Au_111_2h3o_so4_228_std * weight_hydronium * 2 + yoffset*0
plt.fill_between(energy_Au_111_2h3o_so4_228, mean_Au_111_2h3o_so4_228 - std_Au_111_2h3o_so4_228, mean_Au_111_2h3o_so4_228 + std_Au_111_2h3o_so4_228, color='#515151', alpha=0.15)
plt.plot(energy_Au_111_2h3o_so4_228, mean_Au_111_2h3o_so4_228, label='Au_111_2h3o_so4_hydronium228Zundel', linewidth=2, color='#515151')
if x0_Au_111_2h3o_so4_228 is not None and y0_Au_111_2h3o_so4_228 is not None:
    plt.annotate(f'{x0_Au_111_2h3o_so4_228:.1f}', (x0_Au_111_2h3o_so4_228, y0_Au_111_2h3o_so4_228*weight_hydronium*2+yoffset*0),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black')
#plt.plot(energy_2h3o_so4_hydronium29Zundel, intensity_2h3o_so4_hydronium29Zundel*weight_hydronium*2+yoffset*1.5, label='2h3o_so4_hydronium29Zundel', linewidth=2, color='black')
#if x0_2h3o_so4_hydronium29Zundel is not None and y0_2h3o_so4_hydronium29Zundel is not None:
#    plt.annotate(f'{x0_2h3o_so4_hydronium29Zundel:.1f}', (x0_2h3o_so4_hydronium29Zundel, y0_2h3o_so4_hydronium29Zundel*weight_hydronium*2+yoffset*1.5),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black')
#if x1_2h3o_so4_hydronium29Zundel is not None and y1_2h3o_so4_hydronium29Zundel is not None:
#    plt.annotate(f'{x1_2h3o_so4_hydronium29Zundel:.1f}', (x1_2h3o_so4_hydronium29Zundel, y1_2h3o_so4_hydronium29Zundel*weight_hydronium*2+yoffset*1.5),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black')
# Plot mean +/- std shading for hydronium all-frames
mean_2h3o_so4_hydronium_all = intensity_2h3o_so4_hydronium_allFrames_mean * weight_hydronium * 2 + yoffset*1.5
std_2h3o_so4_hydronium_all = intensity_2h3o_so4_hydronium_allFrames_std * weight_hydronium * 2
plt.fill_between(energy_2h3o_so4_hydronium_allFrames, mean_2h3o_so4_hydronium_all - std_2h3o_so4_hydronium_all, mean_2h3o_so4_hydronium_all + std_2h3o_so4_hydronium_all, color='#515151', alpha=0.15)
plt.plot(energy_2h3o_so4_hydronium_allFrames, mean_2h3o_so4_hydronium_all, label='2h3o_so4_hydronium_allFrames', linewidth=2, color='#515151')
if x0_2h3o_so4_hydronium_allFrames is not None and y0_2h3o_so4_hydronium_allFrames is not None:
    plt.annotate(f'{x0_2h3o_so4_hydronium_allFrames:.1f}', (x0_2h3o_so4_hydronium_allFrames, y0_2h3o_so4_hydronium_allFrames*weight_hydronium*2+yoffset*1.5),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black')
#if x1_2h3o_so4_hydronium_allFrames is not None and y1_2h3o_so4_hydronium_allFrames is not None:
#    plt.annotate(f'{x1_2h3o_so4_hydronium_allFrames:.1f}', (x1_2h3o_so4_hydronium_allFrames, y1_2h3o_so4_hydronium_allFrames*weight_hydronium*2+yoffset*1.5),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black')
#plt.plot(energy_Au_111_cuso4_sulfate, intensity_Au_111_cuso4_sulfate*weight_sulfate+yoffset*1.5, label='Au_111_cuso4_sulfate', linewidth=2, color='#F14040', linestyle='--')
#if x_Au_111_cuso4_sulfate is not None and y_Au_111_cuso4_sulfate is not None:
#    plt.annotate(f'{x_Au_111_cuso4_sulfate:.1f}', (x_Au_111_cuso4_sulfate, y_Au_111_cuso4_sulfate*weight_sulfate+yoffset*1.5),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black')
mean_Au_111_cuso4_sulfate_all = intensity_Au_111_cuso4_sulfate_allFrames_mean * weight_sulfate + yoffset*3.3
std_Au_111_cuso4_sulfate_all = intensity_Au_111_cuso4_sulfate_allFrames_std * weight_sulfate
plt.fill_between(energy_Au_111_cuso4_sulfate_allFrames, mean_Au_111_cuso4_sulfate_all - std_Au_111_cuso4_sulfate_all, mean_Au_111_cuso4_sulfate_all + std_Au_111_cuso4_sulfate_all, color='#515151', alpha=0.15)
plt.plot(energy_Au_111_cuso4_sulfate_allFrames, mean_Au_111_cuso4_sulfate_all, label='Au_111_cuso4_sulfate_allFrames', linewidth=2, color='#515151')
if x0_Au_111_cuso4_sulfate_allFrames is not None and y0_Au_111_cuso4_sulfate_allFrames is not None:
    plt.annotate(f'{x0_Au_111_cuso4_sulfate_allFrames:.1f}', (x0_Au_111_cuso4_sulfate_allFrames, y0_Au_111_cuso4_sulfate_allFrames*weight_sulfate+yoffset*3.3),
    textcoords="offset points", xytext=(0,2),
    ha='center', fontproperties=arial, color='black')
if len(peaks_Au_111_cuso4_sulfate_allFrames) > 1:
    x1 = energy_Au_111_cuso4_sulfate_allFrames[peaks_Au_111_cuso4_sulfate_allFrames[1]]
    y1 = intensity_Au_111_cuso4_sulfate_allFrames_mean[peaks_Au_111_cuso4_sulfate_allFrames[1]]
    plt.annotate(f'{x1:.1f}', (x1, y1*weight_sulfate+yoffset*3.3), textcoords="offset points", xytext=(0,2), ha='center', fontproperties=arial, color='black')

#plt.plot(energy_cuso4_solv_sulfate, intensity_cuso4_solv_sulfate+yoffset*4, label='cuso4_solv_sulfate', linewidth=2, color='black', )
#if x_cuso4_solv_sulfate is not None and y_cuso4_solv_sulfate is not None:
#    plt.annotate(f'{x_cuso4_solv_sulfate:.1f}', (x_cuso4_solv_sulfate, y_cuso4_solv_sulfate+yoffset*4),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontsize=7, color='black') 

mean_cuso4_solv_sulfate_all = intensity_cuso4_solv_sulfate_allFrames_mean * weight_sulfate + yoffset*5
std_cuso4_solv_sulfate_all = intensity_cuso4_solv_sulfate_allFrames_std * weight_sulfate
plt.fill_between(energy_cuso4_solv_sulfate_allFrames, mean_cuso4_solv_sulfate_all - std_cuso4_solv_sulfate_all, mean_cuso4_solv_sulfate_all + std_cuso4_solv_sulfate_all, color='#515151', alpha=0.15)
plt.plot(energy_cuso4_solv_sulfate_allFrames, mean_cuso4_solv_sulfate_all, label='cuso4_solv_sulfate_allFrames', linewidth=2, color='#515151')
if x0_cuso4_solv_sulfate_allFrames is not None and y0_cuso4_solv_sulfate_allFrames is not None:
    plt.annotate(f'{x0_cuso4_solv_sulfate_allFrames:.1f}', (x0_cuso4_solv_sulfate_allFrames, y0_cuso4_solv_sulfate_allFrames*weight_sulfate+yoffset*5),
    textcoords="offset points", xytext=(0,2),
    ha='center', fontproperties=arial, color='black')
#if x1_cuso4_solv_sulfate_allFrames is not None and y1_cuso4_solv_sulfate_allFrames is not None:
#    plt.annotate(f'{x1_cuso4_solv_sulfate_allFrames:.1f}', (x1_cuso4_solv_sulfate_allFrames, y1_cuso4_solv_sulfate_allFrames*weight_sulfate+yoffset*4),
#    textcoords="offset points", xytext=(0,2),
#    ha='center', fontproperties=arial, color='black')

#plt.plot(energy_cuso4_vac, intensity_cuso4_vac+yoffset*5, label='cuso4_vac', linewidth=2, color='black', )
#if x_cuso4_vac is not None and y_cuso4_vac is not None:
#    plt.annotate(f'{x_cuso4_vac:.1f}', (x_cuso4_vac, y_cuso4_vac+yoffset*5),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontsize=7, color='black') 

#plt.plot(energy_so4_vac, intensity_so4_vac+yoffset*6, label='so4_vac', linewidth=2, color='black', )
#if x_so4_vac is not None and y_so4_vac is not None:
#    plt.annotate(f'{x_so4_vac:.1f}', (x_so4_vac, y_so4_vac+yoffset*6),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontsize=7, color='black') 

#/2
#plt.plot(energy_Au_111_cuso4_water, intensity_Au_111_cuso4_water*weight_water+yoffset*8.2, label='Au_111_cuso4_water', linewidth=2, color='black')
#if x_Au_111_cuso4_water is not None and y_Au_111_cuso4_water is not None:
#    plt.annotate(f'{x_Au_111_cuso4_water:.1f}', (x_Au_111_cuso4_water, y_Au_111_cuso4_water*weight_water+yoffset*8.2),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black')

#plt.plot(energy_2h3o_so4_sulfate, intensity_2h3o_so4_sulfate*weight_sulfate+yoffset*6.5, label='2h3o_so4_sulfate', linewidth=2, color='black')
#if x0_2h3o_so4_sulfate is not None and y0_2h3o_so4_sulfate is not None:
#    plt.annotate(f'{x0_2h3o_so4_sulfate:.1f}', (x0_2h3o_so4_sulfate, y0_2h3o_so4_sulfate*weight_sulfate+yoffset*6.5),
#    textcoords="offset points", xytext=(0,2),
#    ha='center', fontproperties=arial, color='black')
#if x1_2h3o_so4_sulfate is not None and y1_2h3o_so4_sulfate is not None:
#    plt.annotate(f'{x1_2h3o_so4_sulfate:.1f}', (x1_2h3o_so4_sulfate, y1_2h3o_so4_sulfate*weight_sulfate+yoffset*6.5),
#    textcoords="offset points", xytext=(0,2),
#    ha='center', fontproperties=arial, color='black')
# Plot mean +/- std as shaded band for the all-frames 2h3o sulfate dataset
mean_2h3o_so4_sulfate_all = intensity_2h3o_so4_sulfate_allFrames_mean * weight_sulfate + yoffset*6.5
std_2h3o_so4_sulfate_all = intensity_2h3o_so4_sulfate_allFrames_std * weight_sulfate
# shaded region (mean +/- 1 sigma)
plt.fill_between(energy_2h3o_so4_sulfate_allFrames, mean_2h3o_so4_sulfate_all - std_2h3o_so4_sulfate_all, mean_2h3o_so4_sulfate_all + std_2h3o_so4_sulfate_all, color='#515151', alpha=0.15)
plt.plot(energy_2h3o_so4_sulfate_allFrames, mean_2h3o_so4_sulfate_all, label='2h3o_so4_sulfate_allFrames', linewidth=2, color='#515151')
if x0_2h3o_so4_sulfate_allFrames is not None and y0_2h3o_so4_sulfate_allFrames is not None:
    plt.annotate(f'{x0_2h3o_so4_sulfate_allFrames:.1f}', (x0_2h3o_so4_sulfate_allFrames, y0_2h3o_so4_sulfate_allFrames*weight_sulfate+yoffset*6.5),
    textcoords="offset points", xytext=(0,2),
    ha='center', fontproperties=arial, color='black')
#if x1_2h3o_so4_sulfate_allFrames is not None and y1_2h3o_so4_sulfate_allFrames is not None:
#    plt.annotate(f'{x1_2h3o_so4_sulfate_allFrames:.1f}', (x1_2h3o_so4_sulfate_allFrames, y1_2h3o_so4_sulfate_allFrames*weight_sulfate+yoffset*6.5),
#    textcoords="offset points", xytext=(0,2),
#    ha='center', fontproperties=arial, color='black')

#plt.plot(energy_Au_111_2h3o_so4_sulfate, intensity_Au_111_2h3o_so4_sulfate*weight_sulfate+yoffset*10.5, label='Au_111_2h3o_so4_sulfate', linewidth=2, color='black')
#if x0_Au_111_2h3o_so4_sulfate is not None and y0_Au_111_2h3o_so4_sulfate is not None:
#    plt.annotate(f'{x0_Au_111_2h3o_so4_sulfate:.1f}', (x0_Au_111_2h3o_so4_sulfate, y0_Au_111_2h3o_so4_sulfate*weight_sulfate+yoffset*10.5),
#    textcoords="offset points", xytext=(0,2),
#    ha='center', fontproperties=arial, color='black')
#if x1_Au_111_2h3o_so4_sulfate is not None and y1_Au_111_2h3o_so4_sulfate is not None:
#    plt.annotate(f'{x1_Au_111_2h3o_so4_sulfate:.1f}', (x1_Au_111_2h3o_so4_sulfate, y1_Au_111_2h3o_so4_sulfate*weight_sulfate+yoffset*10.5),
#    textcoords="offset points", xytext=(0,2),
#    ha='center', fontproperties=arial, color='black')

mean_Au_111_2h3o_so4_sulfate_all = intensity_Au_111_2h3o_so4_sulfate_allFrames_mean * weight_sulfate + yoffset*7.5
std_Au_111_2h3o_so4_sulfate_all = intensity_Au_111_2h3o_so4_sulfate_allFrames_std * weight_sulfate
plt.fill_between(energy_Au_111_2h3o_so4_sulfate_allFrames, mean_Au_111_2h3o_so4_sulfate_all - std_Au_111_2h3o_so4_sulfate_all, mean_Au_111_2h3o_so4_sulfate_all + std_Au_111_2h3o_so4_sulfate_all, color='#515151', alpha=0.15)
plt.plot(energy_Au_111_2h3o_so4_sulfate_allFrames, mean_Au_111_2h3o_so4_sulfate_all, label='Au_111_2h3o_so4_sulfate_allFrames', linewidth=2, color='#515151')
if x0_Au_111_2h3o_so4_sulfate_allFrames is not None and y0_Au_111_2h3o_so4_sulfate_allFrames is not None:
    plt.annotate(f'{x0_Au_111_2h3o_so4_sulfate_allFrames:.1f}', (x0_Au_111_2h3o_so4_sulfate_allFrames, y0_Au_111_2h3o_so4_sulfate_allFrames*weight_sulfate+yoffset*7.5),
    textcoords="offset points", xytext=(0,2),
    ha='center', fontproperties=arial, color='black')
#if x1_Au_111_2h3o_so4_sulfate_allFrames is not None and y1_Au_111_2h3o_so4_sulfate_allFrames is not None:
#    plt.annotate(f'{x1_Au_111_2h3o_so4_sulfate_allFrames:.1f}', (x1_Au_111_2h3o_so4_sulfate_allFrames, y1_Au_111_2h3o_so4_sulfate_allFrames*weight_sulfate+yoffset*10.5),
#    textcoords="offset points", xytext=(0,2),
#    ha='center', fontproperties=arial, color='black')

#/2
#plt.plot(energy_Au_111_cuso4_water, intensity_Au_111_cuso4_water*weight_water+yoffset*12.2, label='Au_111_cuso4_water', linewidth=2, color='black')
#if x_Au_111_cuso4_water is not None and y_Au_111_cuso4_water is not None:
#    plt.annotate(f'{x_Au_111_cuso4_water:.1f}', (x_Au_111_cuso4_water, y_Au_111_cuso4_water*weight_water+yoffset*12.2),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black')

# Plot watbox mean +/- std as shaded band
mean_watbox = intensity_watbox_mean * weight_water + yoffset*9.0
std_watbox = intensity_watbox_std * weight_water
plt.fill_between(energy_watbox, mean_watbox - std_watbox, mean_watbox + std_watbox, color='#515151', alpha=0.15)
plt.plot(energy_watbox, mean_watbox, label='water', linewidth=2, color='#515151')
if x0_watbox is not None and y0_watbox is not None:
    plt.annotate(f'{x0_watbox:.1f}', (x0_watbox, y0_watbox*weight_water+yoffset*9.0),
    textcoords="offset points", xytext=(0,2),
    ha='center', fontproperties=arial, color='black')
if x1_watbox is not None and y1_watbox is not None:
    plt.annotate(f'{x1_watbox:.1f}', (x1_watbox, y1_watbox*weight_water+yoffset*9.0),
    textcoords="offset points", xytext=(0,2),
    ha='center', fontproperties=arial, color='black')
exp_factor_watbox = 0.03
plt.plot(energy_watbox_exp, intensity_watbox_exp*exp_factor_watbox+yoffset*10.0, label='watbox exp', linewidth=2, color='#515151')
if x0_watbox_exp is not None and y0_watbox_exp is not None:
    plt.annotate(f'{x0_watbox_exp:.1f}', (x0_watbox_exp, y0_watbox_exp*exp_factor_watbox+yoffset*10.0),
    textcoords="offset points", xytext=(0,2),
    ha='center', fontproperties=arial, color='black')
if x1_watbox_exp is not None and y1_watbox_exp is not None:
    plt.annotate(f'{x1_watbox_exp:.1f}', (x1_watbox_exp, y1_watbox_exp*exp_factor_watbox+yoffset*10.0),
    textcoords="offset points", xytext=(0,2),
    ha='center', fontproperties=arial, color='black')
#plt.plot(energy_Au_111_cuso4_combined, intensity_Au_111_cuso4_combined+yoffset*3, label='Au_111_cuso4_combined', linewidth=2, color='#F14040')
#if x_Au_111_cuso4_combined is not None and y_Au_111_cuso4_combined is not None:
#    plt.annotate(f'{x_Au_111_cuso4_combined:.1f}', (x_Au_111_cuso4_combined, y_Au_111_cuso4_combined+yoffset*3),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black')


#plt.plot(energy_Au_111_cuso4_total, intensity_Au_111_cuso4_total+yoffset*3, label='Au_111_cuso4_total', linewidth=2, color='#F14040')
#if x_Au_111_cuso4_total is not None and y_Au_111_cuso4_total is not None:
#    plt.annotate(f'{x_Au_111_cuso4_total:.1f}', (x_Au_111_cuso4_total, y_Au_111_cuso4_total+yoffset*3),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arialbd, color='black') 

#plt.plot(energy_Au_111_2h3o_so4_allFrames_waterAndHydronium, intensity_Au_111_2h3o_so4_allFrames_waterAndHydronium+yoffset*15.0, label='Au_111_2h3o_so4_allFrames_waterAndHydronium', linewidth=2, color='black')
#if x0_Au_111_2h3o_so4_allFrames_waterAndHydronium is not None and y0_Au_111_2h3o_so4_allFrames_waterAndHydronium is not None:
#    plt.annotate(f'{x0_Au_111_2h3o_so4_allFrames_waterAndHydronium:.1f}', (x0_Au_111_2h3o_so4_allFrames_waterAndHydronium, y0_Au_111_2h3o_so4_allFrames_waterAndHydronium+yoffset*15.0),
#    textcoords="offset points", xytext=(0,2),
#    ha='center', fontproperties=arial, color='black')
#if x1_Au_111_2h3o_so4_allFrames_waterAndHydronium is not None and y1_Au_111_2h3o_so4_allFrames_waterAndHydronium is not None:
#    plt.annotate(f'{x1_Au_111_2h3o_so4_allFrames_waterAndHydronium:.1f}', (x1_Au_111_2h3o_so4_allFrames_waterAndHydronium, y1_Au_111_2h3o_so4_allFrames_waterAndHydronium+yoffset*5.0),
#    textcoords="offset points", xytext=(0,2),
#    ha='center', fontproperties=arial, color='black')


# Plot reference 

#plt.plot(energy_bulkCuO, intensity_bulkCuO+yoffset*14, label='bulkCuO', linewidth=2, color='midnightblue', )
#if x_bulkCuO is not None and y_bulkCuO is not None:
#    plt.annotate(f'{x_bulkCuO:.1f}', (x_bulkCuO, y_bulkCuO+yoffset*14),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arialbd, color='black') 

# Plot settings
#plt.figure(figsize=(4, 3.5))
plt.xlabel('Energy (eV)', fontproperties=arialbd)
#tick_positions = np.arange(928, 945, 1)
#tick_labels = [str(tick) if tick % 2 == 0 else '' for tick in tick_positions]
#plt.xticks(tick_positions, tick_labels,fontproperties=arial)
plt.ylabel('Intensity (a.u.)', fontproperties=arialbd)
plt.gca().set_yticklabels([])  # Removes y-axis labels, keeps tick marks
plt.xticks(fontproperties=arialbd)
plt.xlim(525, 550)
# Add vertical dotted reference lines at 534.9 and 536.4 eV
plt.axvline(534.9, color='black', linestyle=':', linewidth=1)
plt.axvline(536.4, color='black', linestyle=':', linewidth=1)
#plt.grid(True)

# Place each line's label to the right of the curve instead of using a legend
ax = plt.gca()
xmin, xmax = ax.get_xlim()
# small horizontal offset from the right border (in data units)
dx = 0.25

#for line in ax.get_lines():
#    label = line.get_label()
#    # Skip automatically generated labels like '_line0'
#    if label.startswith("_"):
#        continue
#    xdata = np.array(line.get_xdata())
#    ydata = np.array(line.get_ydata())
#    if xdata.size == 0:
#        continue
#    # Choose an x position near the right edge but inside the axes
#    x_text = min(xmax - dx, xdata.max())
#    # Interpolate y at x_text (fall back to last y if outside)
#    try:
#        y_text = np.interp(x_text, xdata, ydata)
#    except Exception:
#        y_text = ydata[-1]
#    # Slight vertical offset to avoid overlapping the line
#    voff = 0.05
#    #ax.text(x_text + 0.05, y_text + voff, label, va='center', ha='left', fontproperties=arial, color=line.get_color())
#    ax.text(x_text - 5 , y_text + voff, label, va='center', ha='left', fontproperties=arial, color=line.get_color())

#plt.tight_layout()

# Saving or showing plot
plt.savefig('compare0V.svg', format='svg', bbox_inches='tight')
#plt.savefig('xasM2.png', format='png', dpi=300, bbox_inches='tight')
plt.show()