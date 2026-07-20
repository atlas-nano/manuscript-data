
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
import scipy.signal
import textwrap
from matplotlib import font_manager

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
    mean_intensity = np.mean(data[:, 1:], axis=1)
    std_intensity = np.std(data[:, 1:], axis=1)
    # Mask for normalization window
    mask = (energy >= E_i) & (energy <= E_f)
    energy_range = energy[mask]
    mean_intensity_range = mean_intensity[mask]
    
    auc = simpson(mean_intensity_range, energy_range)
    print(f"auc: {auc}")
    print(f"std : {np.max(std_intensity)}")
    
    # Normalize both mean and std by the area under the curve
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

# Experimental reference data
energy_bulkCu_exp, intensity_bulkCu_exp = exp("exp_bulkCu.dat")
x_bulkCu_exp,y_bulkCu_exp  = givePeak(energy=energy_bulkCu_exp, intensity=intensity_bulkCu_exp, name='bulkCu')
# COmputational reference
energy_bulkCu, intensity_bulkCu = norm_area("polAvg_bulkCu.dat")
x_bulkCu,y_bulkCu  = givePeak(energy=energy_bulkCu, intensity=intensity_bulkCu, name='bulkCu')
shift_bulkCu = getShift('david_bulkCu.csv')
# ------------------------------------------------------------------------------------------------------------
# Experimental shift
d_exp = x_bulkCu_exp - (x_bulkCu+shift_bulkCu)
print(f'd_exp = {d_exp}')
# ------------------------------------------------------------------------------------------------------------
print('bulkCu total shift: ' + str(shift_bulkCu + d_exp))
energy_bulkCu += shift_bulkCu + d_exp
x_bulkCu,y_bulkCu  = givePeak(energy=energy_bulkCu, intensity=intensity_bulkCu, name='bulkCu')

# Experimental data
energy_minus01V_exp, intensity_minus01V_exp = exp("exp_-0.1V.dat")
x_minus01V_exp,y_minus01V_exp  = givePeak(energy=energy_minus01V_exp, intensity=intensity_minus01V_exp, name='minus01V')

energy_minus02V_exp, intensity_minus02V_exp = exp("exp_-0.2V.dat")
peaks_minus02V_exp, _ = scipy.signal.find_peaks(intensity_minus02V_exp, prominence=0.15)
x_minus02V_exp_0 = energy_minus02V_exp[peaks_minus02V_exp[0]]
y_minus02V_exp_0 = intensity_minus02V_exp[peaks_minus02V_exp[0]]
x_minus02V_exp_1 = energy_minus02V_exp[peaks_minus02V_exp[1]]
y_minus02V_exp_1 = intensity_minus02V_exp[peaks_minus02V_exp[1]]
#x_minus02V_exp, y_minus02V_exp = givePeak(energy=energy_minus02V_exp, intensity=intensity_minus02V_exp, name='minus02V')

energy_minus03V_exp, intensity_minus03V_exp = exp("exp_-0.3V.dat")
x_minus03V_exp, y_minus03V_exp = givePeak(energy=energy_minus03V_exp, intensity=intensity_minus03V_exp, name='minus03V')

energy_0V_exp, intensity_0V_exp = exp("exp_0V.dat")
x_0V_exp, y_0V_exp = givePeak(energy=energy_0V_exp, intensity=intensity_0V_exp, name='0V')


# Process data

energy_cuso4_solv, intensity_cuso4_solv = norm_area("polAvg_cuso4_solv.dat")
shift_cuso4_solv = getShift('david_cuso4_solv.csv')
print('cuso4_solv total_shift: ' + str(shift_cuso4_solv + d_exp))
energy_cuso4_solv += shift_cuso4_solv+ d_exp
x_cuso4_solv,y_cuso4_solv  = givePeak(energy=energy_cuso4_solv, intensity=intensity_cuso4_solv, name='cuso4_solv')

energy_Au_111_Cu_atom, intensity_Au_111_Cu_atom = norm_area("polAvg_Au_111_Cu_atom.dat")
shift_Au_111_Cu_atom = getShift('david_Au_111_Cu_atom.csv')
print('Au_111_Cu_atom total_shift: ' + str(shift_Au_111_Cu_atom + d_exp))
energy_Au_111_Cu_atom += shift_Au_111_Cu_atom+ d_exp
x_Au_111_Cu_atom,y_Au_111_Cu_atom  = givePeak(energy=energy_Au_111_Cu_atom, intensity=intensity_Au_111_Cu_atom, name='Au_111_Cu_atom')

energy_bulkCu2O_exp, intensity_bulkCu2O_exp = exp("exp_bulkCu2O.dat")
# Ensure energy is increasing (some experimental files are in descending energy order)
if energy_bulkCu2O_exp[0] > energy_bulkCu2O_exp[-1]:
    energy_bulkCu2O_exp = energy_bulkCu2O_exp[::-1]
    intensity_bulkCu2O_exp = intensity_bulkCu2O_exp[::-1]
peaks_bulkCu2O_exp, _ = scipy.signal.find_peaks(intensity_bulkCu2O_exp, prominence=0.1)
# store peak positions (may be multiple)
x_bulkCu2O_exp_0 = energy_bulkCu2O_exp[peaks_bulkCu2O_exp[0]]
y_bulkCu2O_exp_0 = intensity_bulkCu2O_exp[peaks_bulkCu2O_exp[0]]
x_bulkCu2O_exp_1 = energy_bulkCu2O_exp[peaks_bulkCu2O_exp[1]]
y_bulkCu2O_exp_1 = intensity_bulkCu2O_exp[peaks_bulkCu2O_exp[1]]

energy_bulkCuO_exp, intensity_bulkCuO_exp = exp("exp_bulkCuO.dat")
# Ensure energy is increasing (some experimental files are in descending energy order)
if energy_bulkCuO_exp[0] > energy_bulkCuO_exp[-1]:
    energy_bulkCuO_exp = energy_bulkCuO_exp[::-1]
    intensity_bulkCuO_exp = intensity_bulkCuO_exp[::-1]
peaks_bulkCuO_exp, _ = scipy.signal.find_peaks(intensity_bulkCuO_exp, prominence=0.1)
# store peak positions (may be multiple)
x_bulkCuO_exp_0 = energy_bulkCuO_exp[peaks_bulkCuO_exp[0]]
y_bulkCuO_exp_0 = intensity_bulkCuO_exp[peaks_bulkCuO_exp[0]]
x_bulkCuO_exp_1 = energy_bulkCuO_exp[peaks_bulkCuO_exp[1]]
y_bulkCuO_exp_1 = intensity_bulkCuO_exp[peaks_bulkCuO_exp[1]]


# Load per-frame data so we can compute mean and std and shade the
# standard deviation like other multi-frame datasets
energy_bulkCu2O, intensity_bulkCu2O_mean, intensity_bulkCu2O_std = norm_area_allframes("polAvg_bulkCu2O.dat")
shift_bulkCu2O = getShift('david_bulkCu2O.csv')
print('bulkCu2O total_shift: ' + str(shift_bulkCu2O + d_exp))
energy_bulkCu2O += shift_bulkCu2O+ d_exp
peaks_bulkCu2O, _ = scipy.signal.find_peaks(intensity_bulkCu2O_mean, prominence=0.01)
if peaks_bulkCu2O.size > 0:
    x_bulkCu2O_0 = energy_bulkCu2O[peaks_bulkCu2O[0]]
    y_bulkCu2O_0 = intensity_bulkCu2O_mean[peaks_bulkCu2O[0]]
    # guard against index error if fewer peaks
    if peaks_bulkCu2O.size > 2:
        x_bulkCu2O_1 = energy_bulkCu2O[peaks_bulkCu2O[2]]
        y_bulkCu2O_1 = intensity_bulkCu2O_mean[peaks_bulkCu2O[2]]
    else:
        x_bulkCu2O_1 = None
        y_bulkCu2O_1 = None
else:
    x_bulkCu2O_0 = None
    y_bulkCu2O_0 = None
    x_bulkCu2O_1 = None
    y_bulkCu2O_1 = None

# Load per-frame mean and std for bulkCuO_mag so we can shade the std
energy_bulkCuO_mag, intensity_bulkCuO_mag_mean, intensity_bulkCuO_mag_std = norm_area_allframes("polAvg_bulkCuO_mag.dat")
shift_bulkCuO_mag = getShift('david_bulkCuO_mag.csv')
print('bulkCuO_mag total_shift: ' + str(shift_bulkCuO_mag + d_exp))
energy_bulkCuO_mag += shift_bulkCuO_mag+ d_exp
peaks_bulkCuO_mag, _ = scipy.signal.find_peaks(intensity_bulkCuO_mag_mean, prominence=0.05)
if peaks_bulkCuO_mag.size > 0:
    x_bulkCuO_mag_0 = energy_bulkCuO_mag[peaks_bulkCuO_mag[0]]
    y_bulkCuO_mag_0 = intensity_bulkCuO_mag_mean[peaks_bulkCuO_mag[0]]
    if peaks_bulkCuO_mag.size > 1:
        x_bulkCuO_mag_1 = energy_bulkCuO_mag[peaks_bulkCuO_mag[1]]
        y_bulkCuO_mag_1 = intensity_bulkCuO_mag_mean[peaks_bulkCuO_mag[1]]
    else:
        x_bulkCuO_mag_1 = None
        y_bulkCuO_mag_1 = None
else:
    x_bulkCuO_mag_0 = None
    y_bulkCuO_mag_0 = None
    x_bulkCuO_mag_1 = None
    y_bulkCuO_mag_1 = None

energy_Au_111_cuso4, intensity_Au_111_cuso4 = norm_area("polAvg_Au_111_cuso4.dat")
shift_Au_111_cuso4 = getShift('david_Au_111_cuso4.csv')
print('Au_111_cuso4 total_shift: ' + str(shift_Au_111_cuso4 + d_exp))
energy_Au_111_cuso4 += shift_Au_111_cuso4+ d_exp
peaks_Au_111_cuso4, _ = scipy.signal.find_peaks(intensity_Au_111_cuso4)
x0_Au_111_cuso4 = energy_Au_111_cuso4[peaks_Au_111_cuso4[0]]
y0_Au_111_cuso4 = intensity_Au_111_cuso4[peaks_Au_111_cuso4[0]]
x_Au_111_cuso4,y_Au_111_cuso4  = givePeak(energy=energy_Au_111_cuso4, intensity=intensity_Au_111_cuso4, name='Au_111_cuso4')


energy_Au_111_cu_ads_so4_diemac_1_80, intensity_Au_111_cu_ads_so4_diemac_1_80 = norm_area("polAvg_Au_111_cu_ads_so4_diemac_1_80.dat")
shift_Au_111_cu_ads_so4_diemac_1_80 = getShift('david_Au_111_cu_ads_so4.csv')
print('Au_111_cu_ads_so4_diemac_1_80 total_shift: ' + str(shift_Au_111_cu_ads_so4_diemac_1_80 + d_exp))
energy_Au_111_cu_ads_so4_diemac_1_80 += shift_Au_111_cu_ads_so4_diemac_1_80+ d_exp
peaks_Au_111_cu_ads_so4_diemac_1_80, _ = scipy.signal.find_peaks(intensity_Au_111_cu_ads_so4_diemac_1_80)
x0_Au_111_cu_ads_so4_diemac_1_80 = energy_Au_111_cu_ads_so4_diemac_1_80[peaks_Au_111_cu_ads_so4_diemac_1_80[0]]
y0_Au_111_cu_ads_so4_diemac_1_80 = intensity_Au_111_cu_ads_so4_diemac_1_80[peaks_Au_111_cu_ads_so4_diemac_1_80[0]]
x_Au_111_cu_ads_so4_diemac_1_80,y_Au_111_cu_ads_so4_diemac_1_80  = givePeak(energy=energy_Au_111_cu_ads_so4_diemac_1_80, intensity=intensity_Au_111_cu_ads_so4_diemac_1_80, name='Au_111_cu_ads_so4_diemac_1_80')

energy_Au_111_Cu_layer, intensity_Au_111_Cu_layer = norm_area("polAvg_Au_111_Cu_layer.dat")
shift_Au_111_Cu_layer = getShift('david_Au_111_Cu_layer.csv')
print('Au_111_Cu_layer total_shift: ' + str(shift_Au_111_Cu_layer + d_exp))
energy_Au_111_Cu_layer += shift_Au_111_Cu_layer+ d_exp
x_Au_111_Cu_layer,y_Au_111_Cu_layer  = givePeak(energy=energy_Au_111_Cu_layer, intensity=intensity_Au_111_Cu_layer, name='Au_111_Cu_layer')


energy_Au_111_CuO_Layer_0, intensity_Au_111_CuO_Layer_0 = norm_area("polAvg_Au_111_CuO_Layer_0.dat")
shift_Au_111_CuO_Layer_0 = getShift('david_Au_111_CuO_Layer_0.csv')
print('Au_111_CuO_Layer_0 total_shift: ' + str(shift_Au_111_CuO_Layer_0 + d_exp))
energy_Au_111_CuO_Layer_0 += shift_Au_111_CuO_Layer_0+ d_exp
x_Au_111_CuO_Layer_0,y_Au_111_CuO_Layer_0  = givePeak(energy=energy_Au_111_CuO_Layer_0, intensity=intensity_Au_111_CuO_Layer_0, name='Au_111_CuO_Layer_0')

energy_Au_111_CuO_Layer_0_solv, intensity_Au_111_CuO_Layer_0_solv = norm_area("polAvg_Au_111_CuO_Layer_0_solv.dat")
shift_Au_111_CuO_Layer_0_solv = getShift('david_Au_111_CuO_Layer_0_solv.csv')
print('Au_111_CuO_Layer_0_solv total_shift: ' + str(shift_Au_111_CuO_Layer_0_solv + d_exp))
energy_Au_111_CuO_Layer_0_solv += shift_Au_111_CuO_Layer_0_solv+ d_exp
x_Au_111_CuO_Layer_0_solv,y_Au_111_CuO_Layer_0_solv  = givePeak(energy=energy_Au_111_CuO_Layer_0_solv, intensity=intensity_Au_111_CuO_Layer_0_solv, name='Au_111_CuO_Layer_0_solv')

energy_Au_111_CuO_Layer_0_solv_mag, intensity_Au_111_CuO_Layer_0_solv_mag_mean, intensity_Au_111_CuO_Layer_0_solv_mag_std = norm_area_allframes("polAvg_Au_111_CuO_Layer_0_solv_mag.dat")
shift_Au_111_CuO_Layer_0_solv_mag = getShift('david_Au_111_CuO_Layer_0_solv_mag_cu1.csv')
print('Au_111_CuO_Layer_0_solv_mag total_shift: ' + str(shift_Au_111_CuO_Layer_0_solv_mag + d_exp))
energy_Au_111_CuO_Layer_0_solv_mag += shift_Au_111_CuO_Layer_0_solv_mag + d_exp
x_Au_111_CuO_Layer_0_solv_mag, y_Au_111_CuO_Layer_0_solv_mag = givePeak(energy=energy_Au_111_CuO_Layer_0_solv_mag, intensity=intensity_Au_111_CuO_Layer_0_solv_mag_mean, name='Au_111_CuO_Layer_0_solv_mag')

energy_Au_111_CuO_Layer_1_solv_mag, intensity_Au_111_CuO_Layer_1_solv_mag_mean, intensity_Au_111_CuO_Layer_1_solv_mag_std = norm_area_allframes("polAvg_Au_111_CuO_Layer_1_solv_mag.dat")
shift_Au_111_CuO_Layer_1_solv_mag = getShift('david_Au_111_CuO_Layer_1_solv_mag.csv')
print('Au_111_CuO_Layer_1_solv_mag total_shift: ' + str(shift_Au_111_CuO_Layer_1_solv_mag + d_exp))
energy_Au_111_CuO_Layer_1_solv_mag += shift_Au_111_CuO_Layer_1_solv_mag + d_exp
x_Au_111_CuO_Layer_1_solv_mag, y_Au_111_CuO_Layer_1_solv_mag = givePeak(energy=energy_Au_111_CuO_Layer_1_solv_mag, intensity=intensity_Au_111_CuO_Layer_1_solv_mag_mean, name='Au_111_CuO_Layer_1_solv_mag')


energy_Au_111_CuO_Layer_0_solv_mag_u8, intensity_Au_111_CuO_Layer_0_solv_mag_u8 = norm_area("polAvg_Au_111_CuO_Layer_0_solv_mag_u8.dat")
shift_Au_111_CuO_Layer_0_solv_mag_u8 = getShift('david_Au_111_CuO_Layer_0_solv_mag_u8.csv')
print('Au_111_CuO_Layer_0_solv_mag_u8 total_shift: ' + str(shift_Au_111_CuO_Layer_0_solv_mag_u8 + d_exp))
energy_Au_111_CuO_Layer_0_solv_mag_u8 += shift_Au_111_CuO_Layer_0_solv_mag_u8+ d_exp
peaks_Au_111_CuO_Layer_0_solv_mag_u8, _ = scipy.signal.find_peaks(intensity_Au_111_CuO_Layer_0_solv_mag_u8, prominence=0.02)
x_Au_111_CuO_Layer_0_solv_mag_u8_0 = energy_Au_111_CuO_Layer_0_solv_mag_u8[peaks_Au_111_CuO_Layer_0_solv_mag_u8[0]]
y_Au_111_CuO_Layer_0_solv_mag_u8_0 = intensity_Au_111_CuO_Layer_0_solv_mag_u8[peaks_Au_111_CuO_Layer_0_solv_mag_u8[0]]
x_Au_111_CuO_Layer_0_solv_mag_u8_1 = energy_Au_111_CuO_Layer_0_solv_mag_u8[peaks_Au_111_CuO_Layer_0_solv_mag_u8[1]]
y_Au_111_CuO_Layer_0_solv_mag_u8_1 = intensity_Au_111_CuO_Layer_0_solv_mag_u8[peaks_Au_111_CuO_Layer_0_solv_mag_u8[1]]
#x_Au_111_CuO_Layer_0_solv_mag_u8,y_Au_111_CuO_Layer_0_solv_mag_u8  = givePeak(energy=energy_Au_111_CuO_Layer_0_solv_mag_u8, intensity=intensity_Au_111_CuO_Layer_0_solv_mag_u8, name='Au_111_CuO_Layer_0_solv_mag_u8')

energy_Au_111_CuO_Layer_0_solv_allFrames, intensity_Au_111_CuO_Layer_0_solv_allFrames = norm_area("polAvg_Au_111_CuO_Layer_0_solv_allFrames.dat")
shift_Au_111_CuO_Layer_0_solv_allFrames = getShift('david_Au_111_CuO_Layer_0_solv.csv')
print('Au_111_CuO_Layer_0_solv_allFrames total_shift: ' + str(shift_Au_111_CuO_Layer_0_solv_allFrames + d_exp))
energy_Au_111_CuO_Layer_0_solv_allFrames += shift_Au_111_CuO_Layer_0_solv_allFrames+ d_exp
peaks_Au_111_CuO_Layer_0_solv_allFrames, _ = scipy.signal.find_peaks(intensity_Au_111_CuO_Layer_0_solv_allFrames, prominence=0.05)
x_Au_111_CuO_Layer_0_solv_allFrames_0 = energy_Au_111_CuO_Layer_0_solv_allFrames[peaks_Au_111_CuO_Layer_0_solv_allFrames[0]]
y_Au_111_CuO_Layer_0_solv_allFrames_0 = intensity_Au_111_CuO_Layer_0_solv_allFrames[peaks_Au_111_CuO_Layer_0_solv_allFrames[0]]
x_Au_111_CuO_Layer_0_solv_allFrames_1 = energy_Au_111_CuO_Layer_0_solv_allFrames[peaks_Au_111_CuO_Layer_0_solv_allFrames[1]]
y_Au_111_CuO_Layer_0_solv_allFrames_1 = intensity_Au_111_CuO_Layer_0_solv_allFrames[peaks_Au_111_CuO_Layer_0_solv_allFrames[1]]
#x_Au_111_CuO_Layer_0_solv_allFrames,y_Au_111_CuO_Layer_0_solv_allFrames  = givePeak(energy=energy_Au_111_CuO_Layer_0_solv_allFrames, intensity=intensity_Au_111_CuO_Layer_0_solv_allFrames, name='Au_111_CuO_Layer_0_solv_allFrames')

energy_Au_111_Cu2O_Layer_0_solv, intensity_Au_111_Cu2O_Layer_0_solv_mean, intensity_Au_111_Cu2O_Layer_0_solv_std = norm_area_allframes("polAvg_Au_111_Cu2O_Layer_0_solv.dat")
shift_Au_111_Cu2O_Layer_0_solv = getShift('david_Au_111_Cu2O_Layer_0_solv.csv')
print('Au_111_Cu2O_Layer_0_solv total_shift: ' + str(shift_Au_111_Cu2O_Layer_0_solv + d_exp))
energy_Au_111_Cu2O_Layer_0_solv += shift_Au_111_Cu2O_Layer_0_solv + d_exp
x_Au_111_Cu2O_Layer_0_solv, y_Au_111_Cu2O_Layer_0_solv = givePeak(energy=energy_Au_111_Cu2O_Layer_0_solv, intensity=intensity_Au_111_Cu2O_Layer_0_solv_mean, name='Au_111_Cu2O_Layer_0_solv')

energy_Au_111_Cu2O_Layer_1_solv_mag, intensity_Au_111_Cu2O_Layer_1_solv_mag_mean, intensity_Au_111_Cu2O_Layer_1_solv_mag_std = norm_area_allframes("polAvg_Au_111_Cu2O_Layer_1_solv_mag.dat")
shift_Au_111_Cu2O_Layer_1_solv_mag = getShift('david_Au_111_Cu2O_Layer_1_solv_mag.csv')
print('Au_111_Cu2O_Layer_1_solv_mag total_shift: ' + str(shift_Au_111_Cu2O_Layer_1_solv_mag + d_exp))
energy_Au_111_Cu2O_Layer_1_solv_mag += shift_Au_111_Cu2O_Layer_1_solv_mag + d_exp
x_Au_111_Cu2O_Layer_1_solv_mag, y_Au_111_Cu2O_Layer_1_solv_mag = givePeak(energy=energy_Au_111_Cu2O_Layer_1_solv_mag, intensity=intensity_Au_111_Cu2O_Layer_1_solv_mag_mean, name='Au_111_Cu2O_Layer_1_solv_mag')

energy_Au_111_Cu2O_Layer_1_solv_newNomag, intensity_Au_111_Cu2O_Layer_1_solv_newNomag = norm_area("polAvg_Au_111_Cu2O_Layer_1_solv_newNomag.dat")
shift_Au_111_Cu2O_Layer_1_solv_newNomag = getShift('david_Au_111_Cu2O_Layer_1_solv_newNomag.csv')
print('Au_111_Cu2O_Layer_1_solv_newNomag total_shift: ' + str(shift_Au_111_Cu2O_Layer_1_solv_newNomag + d_exp))
energy_Au_111_Cu2O_Layer_1_solv_newNomag += shift_Au_111_Cu2O_Layer_1_solv_newNomag+ d_exp
peaks_Au_111_Cu2O_Layer_1_solv_newNomag, _ = scipy.signal.find_peaks(intensity_Au_111_Cu2O_Layer_1_solv_newNomag, prominence=0.03)
x_Au_111_Cu2O_Layer_1_solv_newNomag_0 = energy_Au_111_Cu2O_Layer_1_solv_newNomag[peaks_Au_111_Cu2O_Layer_1_solv_newNomag[0]]
y_Au_111_Cu2O_Layer_1_solv_newNomag_0 = intensity_Au_111_Cu2O_Layer_1_solv_newNomag[peaks_Au_111_Cu2O_Layer_1_solv_newNomag[0]]
x_Au_111_Cu2O_Layer_1_solv_newNomag_1 = energy_Au_111_Cu2O_Layer_1_solv_newNomag[peaks_Au_111_Cu2O_Layer_1_solv_newNomag[1]]
y_Au_111_Cu2O_Layer_1_solv_newNomag_1 = intensity_Au_111_Cu2O_Layer_1_solv_newNomag[peaks_Au_111_Cu2O_Layer_1_solv_newNomag[1]]

# Add 1 eV as estimate of incorrect screening of elecrtrons according to DOS and based on dry results
#energy_Au_111_CuO_Layer_0_solv += 1.0
#x_Au_111_CuO_Layer_0_solv += 1.0

#energy_Au_111_CuO_Layer_0_solv_allFrames += 1.0
#x_Au_111_CuO_Layer_0_solv_allFrames += 1.0

#energy_Au_111_CuO_Layer_1, intensity_Au_111_CuO_Layer_1 = norm_area("polAvg_Au_111_CuO_Layer_1.dat")
#shift_Au_111_CuO_Layer_1 = getShift('david_Au_111_CuO_Layer_1.csv')
#print('Au_111_CuO_Layer_1 total_shift: ' + str(shift_Au_111_CuO_Layer_1 + d_exp))
#energy_Au_111_CuO_Layer_1 += shift_Au_111_CuO_Layer_1+ d_exp
#x_Au_111_CuO_Layer_1,y_Au_111_CuO_Layer_1  = givePeak(energy=energy_Au_111_CuO_Layer_1, intensity=intensity_Au_111_CuO_Layer_1, name='Au_111_CuO_Layer_1')


#plt.plot(energy_cuso4_solv, intensity_cuso4_solv+yoffset*0, label='cuso4_solv', linewidth=2, color='#515151', )
#if x_cuso4_solv is not None and y_cuso4_solv is not None:
#    plt.annotate(f'{x_cuso4_solv:.1f}', (x_cuso4_solv, y_cuso4_solv+yoffset*0),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 


#plt.plot(energy_Au_111_Cu_atom, intensity_Au_111_Cu_atom+yoffset*0, label='Au_111_Cu_atom', linewidth=2, color='#515151')
#if x_Au_111_Cu_atom is not None and y_Au_111_Cu_atom is not None:
#    plt.annotate(f'{x_Au_111_Cu_atom:.1f}', (x_Au_111_Cu_atom, y_Au_111_Cu_atom+yoffset*0),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 

# Plot settings
plt.figure(figsize=(3.0, 3.5))

#plt.plot(energy_Au_111_cuso4, intensity_Au_111_cuso4+yoffset*0, label='Au_111_cuso4', linewidth=2, color='#515151', )
#if x0_Au_111_cuso4 is not None and y0_Au_111_cuso4 is not None:
#    plt.annotate(f'{x0_Au_111_cuso4:.1f}', (x0_Au_111_cuso4, y0_Au_111_cuso4+yoffset*0),
#    textcoords="offset points", xytext=(0,2),
#    ha='center', fontproperties=arial, color='black')
#if x_Au_111_cuso4 is not None and y_Au_111_cuso4 is not None:
#    plt.annotate(f'{x_Au_111_cuso4:.1f}', (x_Au_111_cuso4, y_Au_111_cuso4+yoffset*0),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 

exp_scaling = 7

#plt.plot(energy_0V_exp, intensity_0V_exp/exp_scaling + yoffset*-1, label='0 V', linewidth=2, color='#515151')
##if x_0V_exp is not None and y_0V_exp is not None:
##    plt.annotate(f'{x_0V_exp:.1f}', (x_0V_exp, y_0V_exp/exp_scaling + yoffset*0),
##                 textcoords="offset points", xytext=(0,2),
##                 ha='center', fontproperties=arial, color='black')
#plt.plot(energy_Au_111_cu_ads_so4_diemac_1_80, intensity_Au_111_cu_ads_so4_diemac_1_80+yoffset*-1, label='Au_111_cu_ads_so4_diemac_1_80 z = 0 A', linewidth=2, color='#515151', linestyle=":")
#if x0_Au_111_cu_ads_so4_diemac_1_80 is not None and y0_Au_111_cu_ads_so4_diemac_1_80 is not None:
#    plt.annotate(f'{x0_Au_111_cu_ads_so4_diemac_1_80:.1f}', (x0_Au_111_cu_ads_so4_diemac_1_80, y0_Au_111_cu_ads_so4_diemac_1_80+yoffset*-1),
#    textcoords="offset points", xytext=(0,2),
#    ha='center', fontproperties=arial, color='black')
#if x_Au_111_cu_ads_so4_diemac_1_80 is not None and y_Au_111_cu_ads_so4_diemac_1_80 is not None:
#    plt.annotate(f'{x_Au_111_cu_ads_so4_diemac_1_80:.1f}', (x_Au_111_cu_ads_so4_diemac_1_80, y_Au_111_cu_ads_so4_diemac_1_80+yoffset*-1),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 

#plt.plot(energy_minus01V_exp, intensity_minus01V_exp/exp_scaling+yoffset*0, label='- 0.1 V', linewidth=2, color='#F14040')
##if x_minus01V_exp is not None and y_minus01V_exp is not None:
##    plt.annotate(f'{x_minus01V_exp:.1f}', (x_minus01V_exp, y_minus01V_exp/exp_scaling+yoffset*1),
##    textcoords="offset points", xytext=(0,2), 
##    ha='center', fontproperties=arial, color='black') 
#plt.plot(energy_Au_111_cu_ads_so4_diemac_1_80, intensity_Au_111_cu_ads_so4_diemac_1_80+yoffset*0, label='Au_111_cu_ads_so4_diemac_1_80 z = 0 A', linewidth=2, color='#F14040',  linestyle=":")
#if x0_Au_111_cu_ads_so4_diemac_1_80 is not None and y0_Au_111_cu_ads_so4_diemac_1_80 is not None:
#    plt.annotate(f'{x0_Au_111_cu_ads_so4_diemac_1_80:.1f}', (x0_Au_111_cu_ads_so4_diemac_1_80, y0_Au_111_cu_ads_so4_diemac_1_80+yoffset*0),
#    textcoords="offset points", xytext=(0,2),
#    ha='center', fontproperties=arial, color='black')
#if x_Au_111_cu_ads_so4_diemac_1_80 is not None and y_Au_111_cu_ads_so4_diemac_1_80 is not None:
#    plt.annotate(f'{x_Au_111_cu_ads_so4_diemac_1_80:.1f}', (x_Au_111_cu_ads_so4_diemac_1_80, y_Au_111_cu_ads_so4_diemac_1_80+yoffset*0),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 
    
plt.plot(energy_minus02V_exp, intensity_minus02V_exp/exp_scaling + yoffset*0, label='-0.2 V', color='#1A6FDF', linewidth=2)
if x_minus02V_exp_0 is not None and y_minus02V_exp_0 is not None:
    plt.annotate(f'{x_minus02V_exp_0:.1f}', (x_minus02V_exp_0, y_minus02V_exp_0/exp_scaling + yoffset*0),
                 textcoords="offset points", xytext=(0,2),
                 ha='center', fontproperties=arial, color='black')
if x_minus02V_exp_1 is not None and y_minus02V_exp_1 is not None:
    plt.annotate(f'{x_minus02V_exp_1:.1f}', (x_minus02V_exp_1, y_minus02V_exp_1/exp_scaling + yoffset*0),
                 textcoords="offset points", xytext=(0,2),
                 ha='center', fontproperties=arial, color='black')
#plt.plot(energy_Au_111_CuO_Layer_0, intensity_Au_111_CuO_Layer_0+yoffset*2, label='Au_111_CuO_Layer_0', linewidth=2, color='#1A6FDF',  linestyle=":")
#if x_Au_111_CuO_Layer_0 is not None and y_Au_111_CuO_Layer_0 is not None:
#    plt.annotate(f'{x_Au_111_CuO_Layer_0:.1f}', (x_Au_111_CuO_Layer_0, y_Au_111_CuO_Layer_0+yoffset*2),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 
#plt.plot(energy_Au_111_CuO_Layer_0_solv, intensity_Au_111_CuO_Layer_0_solv+yoffset*1.2, label='Au_111_CuO_Layer_0_solv', linewidth=2, color='#1A6FDF',  linestyle=":")
#if x_Au_111_CuO_Layer_0_solv is not None and y_Au_111_CuO_Layer_0_solv is not None:
#    plt.annotate(f'{x_Au_111_CuO_Layer_0_solv:.1f}', (x_Au_111_CuO_Layer_0_solv, y_Au_111_CuO_Layer_0_solv+yoffset*1.2),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black')
#plt.plot(energy_Au_111_CuO_Layer_0_solv_allFrames, intensity_Au_111_CuO_Layer_0_solv_allFrames+yoffset*1, label='Au_111_CuO_Layer_0_solv_allFrames', linewidth=2, color='#1A6FDF',)
#if x_Au_111_CuO_Layer_0_solv_allFrames_0 is not None and y_Au_111_CuO_Layer_0_solv_allFrames_0 is not None:
#    plt.annotate(f'{x_Au_111_CuO_Layer_0_solv_allFrames_0:.1f}', (x_Au_111_CuO_Layer_0_solv_allFrames_0, y_Au_111_CuO_Layer_0_solv_allFrames_0 + yoffset*1),
#                 textcoords="offset points", xytext=(0,2),
#                 ha='center', fontproperties=arial, color='black')
#if x_Au_111_CuO_Layer_0_solv_allFrames_1 is not None and y_Au_111_CuO_Layer_0_solv_allFrames_1 is not None:
#    plt.annotate(f'{x_Au_111_CuO_Layer_0_solv_allFrames_1:.1f}', (x_Au_111_CuO_Layer_0_solv_allFrames_1, y_Au_111_CuO_Layer_0_solv_allFrames_1 + yoffset*1),
#                 textcoords="offset points", xytext=(0,2),
#                 ha='center', fontproperties=arial, color='black')

plt.fill_between(energy_Au_111_CuO_Layer_0_solv_mag,
    intensity_Au_111_CuO_Layer_0_solv_mag_mean - intensity_Au_111_CuO_Layer_0_solv_mag_std + yoffset*1,
    intensity_Au_111_CuO_Layer_0_solv_mag_mean + intensity_Au_111_CuO_Layer_0_solv_mag_std + yoffset*1,
    color='#1A6FDF', alpha=0.15)
plt.plot(energy_Au_111_CuO_Layer_0_solv_mag, intensity_Au_111_CuO_Layer_0_solv_mag_mean+yoffset*1, label='Au_111_CuO_Layer_0_solv_mag', linewidth=2, color='#1A6FDF')
if x_Au_111_CuO_Layer_0_solv_mag is not None and y_Au_111_CuO_Layer_0_solv_mag is not None:
    plt.annotate(f'{x_Au_111_CuO_Layer_0_solv_mag:.1f}', (x_Au_111_CuO_Layer_0_solv_mag, y_Au_111_CuO_Layer_0_solv_mag + yoffset*1),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black')

plt.fill_between(energy_Au_111_CuO_Layer_1_solv_mag,
    intensity_Au_111_CuO_Layer_1_solv_mag_mean - intensity_Au_111_CuO_Layer_1_solv_mag_std + yoffset*2,
    intensity_Au_111_CuO_Layer_1_solv_mag_mean + intensity_Au_111_CuO_Layer_1_solv_mag_std + yoffset*2,
    color='#1A6FDF', alpha=0.15)
plt.plot(energy_Au_111_CuO_Layer_1_solv_mag, intensity_Au_111_CuO_Layer_1_solv_mag_mean+yoffset*2, label='Au_111_CuO_Layer_1_solv_mag', linewidth=2, color='#1A6FDF')
if x_Au_111_CuO_Layer_1_solv_mag is not None and y_Au_111_CuO_Layer_1_solv_mag is not None:
    plt.annotate(f'{x_Au_111_CuO_Layer_1_solv_mag:.1f}', (x_Au_111_CuO_Layer_1_solv_mag, y_Au_111_CuO_Layer_1_solv_mag + yoffset*2),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black')

#plt.plot(energy_Au_111_CuO_Layer_0_solv_mag_u8, intensity_Au_111_CuO_Layer_0_solv_mag_u8+yoffset*4, label='Au_111_CuO_Layer_0_solv_mag_u8', linewidth=2, color='#1A6FDF',)
#if x_Au_111_CuO_Layer_0_solv_mag_u8_0 is not None and y_Au_111_CuO_Layer_0_solv_mag_u8_0 is not None:
#    plt.annotate(f'{x_Au_111_CuO_Layer_0_solv_mag_u8_0:.1f}', (x_Au_111_CuO_Layer_0_solv_mag_u8_0, y_Au_111_CuO_Layer_0_solv_mag_u8_0 + yoffset*4),
#                 textcoords="offset points", xytext=(0,2),
#                 ha='center', fontproperties=arial, color='black')
#if x_Au_111_CuO_Layer_0_solv_mag_u8_1 is not None and y_Au_111_CuO_Layer_0_solv_mag_u8_1 is not None:
#    plt.annotate(f'{x_Au_111_CuO_Layer_0_solv_mag_u8_1:.1f}', (x_Au_111_CuO_Layer_0_solv_mag_u8_1, y_Au_111_CuO_Layer_0_solv_mag_u8_1 + yoffset*4),
#                 textcoords="offset points", xytext=(0,2),
#                 ha='center', fontproperties=arial, color='black')


plt.fill_between(energy_Au_111_Cu2O_Layer_0_solv,
    intensity_Au_111_Cu2O_Layer_0_solv_mean - intensity_Au_111_Cu2O_Layer_0_solv_std + yoffset*3,
    intensity_Au_111_Cu2O_Layer_0_solv_mean + intensity_Au_111_Cu2O_Layer_0_solv_std + yoffset*3,
    color='#1A6FDF', alpha=0.15)
plt.plot(energy_Au_111_Cu2O_Layer_0_solv, intensity_Au_111_Cu2O_Layer_0_solv_mean+yoffset*3, label='Au_111_Cu2O_Layer_0_solv', linewidth=2, color='#1A6FDF')
if x_Au_111_Cu2O_Layer_0_solv is not None and y_Au_111_Cu2O_Layer_0_solv is not None:
    plt.annotate(f'{x_Au_111_Cu2O_Layer_0_solv:.1f}', (x_Au_111_Cu2O_Layer_0_solv, y_Au_111_Cu2O_Layer_0_solv+yoffset*3),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black') 

plt.fill_between(energy_Au_111_Cu2O_Layer_1_solv_mag,
    intensity_Au_111_Cu2O_Layer_1_solv_mag_mean - intensity_Au_111_Cu2O_Layer_1_solv_mag_std + yoffset*4,
    intensity_Au_111_Cu2O_Layer_1_solv_mag_mean + intensity_Au_111_Cu2O_Layer_1_solv_mag_std + yoffset*4,
    color='#1A6FDF', alpha=0.15)
plt.plot(energy_Au_111_Cu2O_Layer_1_solv_mag, intensity_Au_111_Cu2O_Layer_1_solv_mag_mean+yoffset*4, label='Au_111_Cu2O_Layer_1_solv_mag', linewidth=2, color='#1A6FDF')
if x_Au_111_Cu2O_Layer_1_solv_mag is not None and y_Au_111_Cu2O_Layer_1_solv_mag is not None:
    plt.annotate(f'{x_Au_111_Cu2O_Layer_1_solv_mag:.1f}', (x_Au_111_Cu2O_Layer_1_solv_mag, y_Au_111_Cu2O_Layer_1_solv_mag+yoffset*4),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black')
#plt.plot(energy_Au_111_Cu2O_Layer_1_solv_newNomag, intensity_Au_111_Cu2O_Layer_1_solv_newNomag+yoffset*4, label='Au_111_Cu2O_Layer_1_solv_newNomag', linewidth=2, color='red',)
#if x_Au_111_Cu2O_Layer_1_solv_newNomag_0 is not None and y_Au_111_Cu2O_Layer_1_solv_newNomag_0 is not None:
#    plt.annotate(f'{x_Au_111_Cu2O_Layer_1_solv_newNomag_0:.1f}', (x_Au_111_Cu2O_Layer_1_solv_newNomag_0, y_Au_111_Cu2O_Layer_1_solv_newNomag_0 + yoffset*4),
#                 textcoords="offset points", xytext=(0,2),
#                 ha='center', fontproperties=arial, color='black')
#if x_Au_111_Cu2O_Layer_1_solv_newNomag_1 is not None and y_Au_111_Cu2O_Layer_1_solv_newNomag_1 is not None:
#    plt.annotate(f'{x_Au_111_Cu2O_Layer_1_solv_newNomag_1:.1f}', (x_Au_111_Cu2O_Layer_1_solv_newNomag_1, y_Au_111_Cu2O_Layer_1_solv_newNomag_1 + yoffset*4),
#                 textcoords="offset points", xytext=(0,2),
#                 ha='center', fontproperties=arial, color='black')

#plt.plot(energy_minus03V_exp, intensity_minus03V_exp/exp_scaling + yoffset*4, label='-0.3 V', linewidth=2, color='#37AD6B')
##if x_minus03V_exp is not None and y_minus03V_exp is not None:
##    plt.annotate(f'{x_minus03V_exp:.1f}', (x_minus03V_exp, y_minus03V_exp + yoffset*3),
##                 textcoords="offset points", xytext=(0,2),
##                 ha='center', fontproperties=arial, color='black')
#plt.plot(energy_Au_111_Cu_layer, intensity_Au_111_Cu_layer+yoffset*4, label='Au_111_Cu_layer', linewidth=2, color='#37AD6B',  linestyle=":")
#if x_Au_111_Cu_layer is not None and y_Au_111_Cu_layer is not None:
#    plt.annotate(f'{x_Au_111_Cu_layer:.1f}', (x_Au_111_Cu_layer, y_Au_111_Cu_layer+yoffset*4),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 

plt.plot(energy_bulkCu2O_exp, intensity_bulkCu2O_exp/ 7.0+yoffset*5, label='bulkCu2O exp', linewidth=2, color='#1A6FDF')
if x_bulkCu2O_exp_0 is not None and y_bulkCu2O_exp_0 is not None:
    plt.annotate(f'{x_bulkCu2O_exp_0:.1f}', (x_bulkCu2O_exp_0, y_bulkCu2O_exp_0/7.0 + yoffset*5),
                 textcoords="offset points", xytext=(0,2),
                 ha='center', fontproperties=arial, color='black')
if x_bulkCu2O_exp_1 is not None and y_bulkCu2O_exp_1 is not None:
    plt.annotate(f'{x_bulkCu2O_exp_1:.1f}', (x_bulkCu2O_exp_1, y_bulkCu2O_exp_1/7.0 + yoffset*5),
                 textcoords="offset points", xytext=(0,2),
                 ha='center', fontproperties=arial, color='black')
    
plt.fill_between(energy_bulkCu2O,
    intensity_bulkCu2O_mean - intensity_bulkCu2O_std + yoffset*6.5,
    intensity_bulkCu2O_mean + intensity_bulkCu2O_std + yoffset*6.5,
    color='#1A6FDF', alpha=0.15)
plt.plot(energy_bulkCu2O, intensity_bulkCu2O_mean+yoffset*6.5, label='bulkCu2O', linewidth=2, color='#1A6FDF')
if x_bulkCu2O_0 is not None and y_bulkCu2O_0 is not None:
    plt.annotate(f'{x_bulkCu2O_0:.1f}', (x_bulkCu2O_0, y_bulkCu2O_0 + yoffset*6.5),
                 textcoords="offset points", xytext=(0,2),
                 ha='center', fontproperties=arial, color='black')
if x_bulkCu2O_1 is not None and y_bulkCu2O_1 is not None:
    plt.annotate(f'{x_bulkCu2O_1:.1f}', (x_bulkCu2O_1, y_bulkCu2O_1 + yoffset*6.5),
                 textcoords="offset points", xytext=(0,2),
                 ha='center', fontproperties=arial, color='black')
    
plt.plot(energy_bulkCuO_exp, intensity_bulkCuO_exp/ 7.0+yoffset*7.5, label='bulkCuO exp', linewidth=2, color='#1A6FDF')
if x_bulkCuO_exp_0 is not None and y_bulkCuO_exp_0 is not None:
    plt.annotate(f'{x_bulkCuO_exp_0:.1f}', (x_bulkCuO_exp_0, y_bulkCuO_exp_0/7.0 + yoffset*7.5),
                 textcoords="offset points", xytext=(0,2),
                 ha='center', fontproperties=arial, color='black')
if x_bulkCuO_exp_1 is not None and y_bulkCuO_exp_1 is not None:
    plt.annotate(f'{x_bulkCuO_exp_1:.1f}', (x_bulkCuO_exp_1, y_bulkCuO_exp_1/7.0 + yoffset*7.5),
                 textcoords="offset points", xytext=(0,2),
                 ha='center', fontproperties=arial, color='black')

plt.fill_between(energy_bulkCuO_mag,
    intensity_bulkCuO_mag_mean - intensity_bulkCuO_mag_std + yoffset*9,
    intensity_bulkCuO_mag_mean + intensity_bulkCuO_mag_std + yoffset*9,
    color='#1A6FDF', alpha=0.15)
plt.plot(energy_bulkCuO_mag, intensity_bulkCuO_mag_mean+yoffset*9, label='bulkCuO_mag', linewidth=2, color='#1A6FDF')
if x_bulkCuO_mag_0 is not None and y_bulkCuO_mag_0 is not None:
    plt.annotate(f'{x_bulkCuO_mag_0:.1f}', (x_bulkCuO_mag_0, y_bulkCuO_mag_0 + yoffset*9),
                 textcoords="offset points", xytext=(0,2),
                 ha='center', fontproperties=arial, color='black')
if x_bulkCuO_mag_1 is not None and y_bulkCuO_mag_1 is not None:
    plt.annotate(f'{x_bulkCuO_mag_1:.1f}', (x_bulkCuO_mag_1, y_bulkCuO_mag_1 + yoffset*9),
                 textcoords="offset points", xytext=(0,2),
                 ha='center', fontproperties=arial, color='black')

# Plot reference 

#plt.plot(energy_bulkCu, intensity_bulkCu+yoffset*14, label='bulkCu', linewidth=2, color='dimgray', )
#if x_bulkCu is not None and y_bulkCu is not None:
#    plt.annotate(f'{x_bulkCu:.1f}', (x_bulkCu, y_bulkCu+yoffset*14),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 

# Plot settings
plt.xlabel('Energy (eV)', fontproperties=arialbd)
#tick_positions = np.arange(928, 945, 1)
#tick_labels = [str(tick) if tick % 2 == 0 else '' for tick in tick_positions]
#plt.xticks(tick_positions, tick_labels,fontproperties=arial)
plt.ylabel('Intensity (a.u.)', fontproperties=arialbd)
plt.gca().set_yticks([])  # Remove y-axis ticks and labels
#plt.yticks(fontproperties=arial)
plt.xticks(fontproperties=arialbd)
plt.xlim(926,960)
#plt.grid(True)

ax = plt.gca()
# Annotate each labeled line near the right edge instead of using a legend
# Compute a suitable x position a bit further inside the right axis limit so
# the text doesn't run outside the axes. Place the annotation to the left
# of the reference point (negative x offset) and right-align the text.
x_min, x_max = ax.get_xlim()
# Move labels 6% of the axis width inside the right edge (matches o-k-edge behavior)
x_label_pos = x_max - 0.06 * (x_max - x_min)
#for line in ax.get_lines():
#    label = line.get_label()
#    # skip lines without an explicit label (matplotlib labels starting with '_')
#    if not label or label.startswith('_'):
#        continue
#    xd = np.asarray(line.get_xdata())
#    yd = np.asarray(line.get_ydata())
#    # interpolate y at x_label_pos (safe if xd is monotonic)
#    try:
#        y_at = np.interp(x_label_pos, xd, yd)
#    except Exception:
#        # fallback to last point if interpolation fails or data empty
#        y_at = yd[-1] if yd.size > 0 else 0
#    # place text to the left of the anchor point so it remains inside axes
#    ax.annotate(label,
#                xy=(x_label_pos, y_at),
#                xytext=(-6, -5),
#                textcoords='offset points',
#                fontproperties=arialbd,
#                va='center',
#                ha='right')

# Make layout tight so annotations aren't clipped by figure margins
plt.tight_layout()
plt.axvline(933.0, color='black', linestyle=':', linewidth=1)
plt.axvline(952.9, color='black', linestyle=':', linewidth=1)
# Saving or showing plot
plt.savefig('compare_cu.svg', format='svg', bbox_inches='tight')
plt.show()
