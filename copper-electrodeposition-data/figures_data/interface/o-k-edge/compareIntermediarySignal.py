
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
import scipy.signal
import textwrap
from scipy.interpolate import interp1d
from matplotlib import font_manager

arial = font_manager.FontProperties(fname="/mnt/c/Windows/Fonts/arial.ttf", size=9)
arialbd = font_manager.FontProperties(fname="/mnt/c/Windows/Fonts/arialbd.ttf", size=9)
# Constants


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
energy_bulkCuO_exp, intensity_bulkCuO_exp = exp("exp_bulkCuO.dat")
# Reverse the order of this data
energy_bulkCuO_exp = energy_bulkCuO_exp[::-1]
intensity_bulkCuO_exp = intensity_bulkCuO_exp[::-1]
x_bulkCuO_exp,y_bulkCuO_exp  = givePeak(energy=energy_bulkCuO_exp, intensity=intensity_bulkCuO_exp, name='bulkCuO')
## Computational reference
#energy_bulkCuO, intensity_bulkCuO = norm_area("polAvg_bulkCuO.dat")
#x_bulkCuO,y_bulkCuO  = givePeak(energy=energy_bulkCuO, intensity=intensity_bulkCuO, name='bulkCuO')
#shift_bulkCuO = getShift('david_bulkCuO.csv')
## ------------------------------------------------------------------------------------------------------------
## Experimental shift
#d_exp = x_bulkCuO_exp - (x_bulkCuO+shift_bulkCuO)
#print(f'd_exp = {d_exp}')
## ------------------------------------------------------------------------------------------------------------
#print('bulkCuO total shift: ' + str(shift_bulkCuO + d_exp))
#energy_bulkCuO += shift_bulkCuO + d_exp
#x_bulkCuO,y_bulkCuO  = givePeak(energy=energy_bulkCuO, intensity=intensity_bulkCuO, name='bulkCuO')



# Computational reference
# Use norm_area_allframes to obtain mean and std across frames for proper shading
energy_bulkCuO_mag, intensity_bulkCuO_mag_mean, intensity_bulkCuO_mag_std = norm_area_allframes("polAvg_bulkCuO_mag.dat")
x_bulkCuO_mag,y_bulkCuO_mag  = givePeak(energy=energy_bulkCuO_mag, intensity=intensity_bulkCuO_mag_mean, name='bulkCuO_mag')
print(f'comp bulkCuO_mag: {x_bulkCuO_mag}')
shift_bulkCuO_mag = getShift('david_bulkCuO_mag.csv')
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Experimental shift
d_exp = x_bulkCuO_exp - (x_bulkCuO_mag+shift_bulkCuO_mag)
print(f'd_exp = {d_exp}')
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::
print('bulkCuO_mag total shift: ' + str(shift_bulkCuO_mag + d_exp))
energy_bulkCuO_mag += shift_bulkCuO_mag + d_exp
# peak detection on the mean intensity returned by norm_area_allframes
x_bulkCuO_mag,y_bulkCuO_mag  = givePeak(energy=energy_bulkCuO_mag, intensity=intensity_bulkCuO_mag_mean, name='bulkCuO_mag')

# Use norm_area_allframes for bulkCu2O so we can show real per-frame std as shading
energy_bulkCu2O, intensity_bulkCu2O_mean, intensity_bulkCu2O_std = norm_area_allframes("polAvg_bulkCu2O.dat")
shift_bulkCu2O = getShift('david_bulkCu2O.csv')
print('bulkCu2O total_shift: ' + str(shift_bulkCu2O + d_exp))
energy_bulkCu2O += shift_bulkCu2O+ d_exp
x_bulkCu2O,y_bulkCu2O  = givePeak(energy=energy_bulkCu2O, intensity=intensity_bulkCu2O_mean, name='bulkCu2O')


# Experimental data

energy_0V_exp, intensity_0V_exp = exp("exp_0V.dat")
x_0V_exp, y_0V_exp = givePeak(energy=energy_0V_exp, intensity=intensity_0V_exp, name='0V')

energy_minus01V_exp, intensity_minus01V_exp = exp("exp_-0.1V.dat")
x_minus01V_exp,y_minus01V_exp  = givePeak(energy=energy_minus01V_exp, intensity=intensity_minus01V_exp, name='minus01V')

energy_minus02V_exp, intensity_minus02V_exp = exp("exp_-0.2V.dat")
peaks_minus02V_exp, _ = scipy.signal.find_peaks(intensity_minus02V_exp, prominence=0.2)
x_minus02V_exp = energy_minus02V_exp[peaks_minus02V_exp[0]]
y_minus02V_exp = intensity_minus02V_exp[peaks_minus02V_exp[0]]
#x_minus02V_exp, y_minus02V_exp = givePeak(energy=energy_minus02V_exp, intensity=intensity_minus02V_exp, name='minus02V')

energy_minus03V_exp, intensity_minus03V_exp = exp("exp_-0.3V.dat")
x_minus03V_exp, y_minus03V_exp = givePeak(energy=energy_minus03V_exp, intensity=intensity_minus03V_exp, name='minus03V')


# Process data

#energy_cuso4_solv_sulfate, intensity_cuso4_solv_sulfate = norm_area("polAvg_cuso4_solv_sulfate.dat")
#shift_cuso4_solv_sulfate = getShift('david_cuso4_solv_sulfate.csv')
#print('cuso4_solv_sulfate total_shift: ' + str(shift_cuso4_solv_sulfate + d_exp))
#energy_cuso4_solv_sulfate += shift_cuso4_solv_sulfate+ d_exp
#x_cuso4_solv_sulfate,y_cuso4_solv_sulfate  = givePeak(energy=energy_cuso4_solv_sulfate, intensity=intensity_cuso4_solv_sulfate, name='cuso4_solv_sulfate')


#energy_Au_111_watbox, intensity_Au_111_watbox = norm_area("polAvg_Au_111_watbox.dat")
#shift_Au_111_watbox = getShift('david_Au_111_watbox.csv')
#print('Au_111_watbox total_shift: ' + str(shift_Au_111_watbox + d_exp))
#energy_Au_111_watbox += shift_Au_111_watbox+ d_exp
#x_Au_111_watbox,y_Au_111_watbox  = givePeak(energy=energy_Au_111_watbox, intensity=intensity_Au_111_watbox, name='Au_111_watbox')

energy_Au_111_cuso4_sulfate, intensity_Au_111_cuso4_sulfate = norm_area("polAvg_Au_111_cuso4_sulfate.dat")
shift_Au_111_cuso4_sulfate = getShift('david_Au_111_cuso4_sulfate.csv')
print('Au_111_cuso4_sulfate total_shift: ' + str(shift_Au_111_cuso4_sulfate + d_exp))
energy_Au_111_cuso4_sulfate += shift_Au_111_cuso4_sulfate+ d_exp
x_Au_111_cuso4_sulfate,y_Au_111_cuso4_sulfate  = givePeak(energy=energy_Au_111_cuso4_sulfate, intensity=intensity_Au_111_cuso4_sulfate, name='Au_111_cuso4_sulfate')

energy_Au_111_cuso4_water, intensity_Au_111_cuso4_water = norm_area("polAvg_Au_111_cuso4_water.dat")
shift_Au_111_cuso4_water = getShift('david_Au_111_cuso4_water.csv')
print('Au_111_cuso4_water water_shift: ' + str(shift_Au_111_cuso4_water + d_exp))
energy_Au_111_cuso4_water += shift_Au_111_cuso4_water+ d_exp
x_Au_111_cuso4_water,y_Au_111_cuso4_water  = givePeak(energy=energy_Au_111_cuso4_water, intensity=intensity_Au_111_cuso4_water, name='Au_111_cuso4_water')

energy_Au_111_2h3o_so4_hydronium256, intensity_Au_111_2h3o_so4_hydronium256 = norm_area("polAvg_Au_111_2h3o_so4_hydronium256.dat")
shift_Au_111_2h3o_so4_hydronium256 = getShift('david_Au_111_2h3o_so4_hydronium256.csv')
print('Au_111_2h3o_so4_hydronium256 water_shift: ' + str(shift_Au_111_2h3o_so4_hydronium256 + d_exp))
energy_Au_111_2h3o_so4_hydronium256 += shift_Au_111_2h3o_so4_hydronium256+ d_exp
x_Au_111_2h3o_so4_hydronium256,y_Au_111_2h3o_so4_hydronium256  = givePeak(energy=energy_Au_111_2h3o_so4_hydronium256, intensity=intensity_Au_111_2h3o_so4_hydronium256, name='Au_111_2h3o_so4_hydronium256')

energy_Au_111_2h3o_so4_hydronium_allFrames, intensity_Au_111_2h3o_so4_hydronium_allFrames = norm_area("polAvg_Au_111_2h3o_so4_hydronium_allFrames.dat")
shift_Au_111_2h3o_so4_hydronium_allFrames = getShift('david_Au_111_2h3o_so4_hydronium256.csv')
print('Au_111_2h3o_so4_hydronium_allFrames water_shift: ' + str(shift_Au_111_2h3o_so4_hydronium_allFrames + d_exp))
energy_Au_111_2h3o_so4_hydronium_allFrames += shift_Au_111_2h3o_so4_hydronium_allFrames+ d_exp
x_Au_111_2h3o_so4_hydronium_allFrames,y_Au_111_2h3o_so4_hydronium_allFrames  = givePeak(energy=energy_Au_111_2h3o_so4_hydronium_allFrames, intensity=intensity_Au_111_2h3o_so4_hydronium_allFrames, name='Au_111_2h3o_so4_hydronium_allFrames')


# Create a common energy grid (you can adjust resolution as needed)
energy_common = np.linspace(min(energy_Au_111_cuso4_sulfate.min(), energy_Au_111_cuso4_water.min(), energy_Au_111_2h3o_so4_hydronium_allFrames.min()), 
                            max(energy_Au_111_cuso4_sulfate.max(), energy_Au_111_cuso4_water.max(), energy_Au_111_2h3o_so4_hydronium_allFrames.max()), 
                            num=1000)

# Interpolate both spectra onto the common grid
interp_sulfate_func = interp1d(energy_Au_111_cuso4_sulfate, intensity_Au_111_cuso4_sulfate, kind='linear', bounds_error=False, fill_value=0)
interp_water_func = interp1d(energy_Au_111_cuso4_water, intensity_Au_111_cuso4_water, kind='linear', bounds_error=False, fill_value=0)
interp_hydronium_func = interp1d(energy_Au_111_2h3o_so4_hydronium_allFrames, intensity_Au_111_2h3o_so4_hydronium_allFrames, kind='linear', bounds_error=False, fill_value=0)

intensity_sulfate_interp = interp_sulfate_func(energy_common)
intensity_water_interp = interp_water_func(energy_common)
intensity_hydronium_interp = interp_hydronium_func(energy_common)

#weight_sulfate= 0.15
#weight_water= 0.80
#weight_hydronium=  0.05

#weight_sulfate= 0.10
#weight_water= 0.30
#weight_hydronium=  0.60

weight_sulfate_0V= 0.00
weight_water_0V= 0.423
weight_hydronium_0V=  0.577

weight_sulfate_minus01V= 0.00
weight_water_minus01V= 0.422
weight_hydronium_minus01V=  0.578

weight_sulfate_minus03V= 0.00
weight_water_minus03V= 0.573
weight_hydronium_minus03V=  0.427


# Final result
energy_Au_111_cuso4_combined = energy_common

intensity_Au_111_cuso4_combined_0V = intensity_sulfate_interp*weight_sulfate_0V + intensity_water_interp*weight_water_0V + intensity_hydronium_interp*weight_hydronium_0V
peaks_Au_111_cuso4_combined_0V, _ = scipy.signal.find_peaks(intensity_Au_111_cuso4_combined_0V, prominence=0.0001)
x0_Au_111_cuso4_combined_0V = energy_Au_111_cuso4_combined[peaks_Au_111_cuso4_combined_0V[0]]
y0_Au_111_cuso4_combined_0V = intensity_Au_111_cuso4_combined_0V[peaks_Au_111_cuso4_combined_0V[0]]
x1_Au_111_cuso4_combined_0V = energy_Au_111_cuso4_combined[peaks_Au_111_cuso4_combined_0V[1]]
y1_Au_111_cuso4_combined_0V = intensity_Au_111_cuso4_combined_0V[peaks_Au_111_cuso4_combined_0V[1]]

#x_Au_111_cuso4_combined_0V,y_Au_111_cuso4_combined_0V  = givePeak(energy=energy_Au_111_cuso4_combined, intensity=intensity_Au_111_cuso4_combined_0V, name='0V')

intensity_Au_111_cuso4_combined_minus01V = intensity_sulfate_interp*weight_sulfate_minus01V + intensity_water_interp*weight_water_minus01V + intensity_hydronium_interp*weight_hydronium_minus01V
x_Au_111_cuso4_combined_minus01V,y_Au_111_cuso4_combined_minus01V  = givePeak(energy=energy_Au_111_cuso4_combined, intensity=intensity_Au_111_cuso4_combined_minus01V, name='minus01V')

intensity_Au_111_cuso4_combined_minus03V = intensity_sulfate_interp*weight_sulfate_minus03V + intensity_water_interp*weight_water_minus03V + intensity_hydronium_interp*weight_hydronium_minus03V
x_Au_111_cuso4_combined_minus03V,y_Au_111_cuso4_combined_minus03V  = givePeak(energy=energy_Au_111_cuso4_combined, intensity=intensity_Au_111_cuso4_combined_minus03V, name='minus03V')


energy_Au_111_Cu_layer, intensity_Au_111_Cu_layer = norm_area("polAvg_Au_111_Cu_layer.dat")
shift_Au_111_Cu_layer = getShift('david_Au_111_Cu_layer.csv')
print('Au_111_Cu_layer total_shift: ' + str(shift_Au_111_Cu_layer + d_exp))
energy_Au_111_Cu_layer += shift_Au_111_Cu_layer+ d_exp
x_Au_111_Cu_layer,y_Au_111_Cu_layer  = givePeak(energy=energy_Au_111_Cu_layer, intensity=intensity_Au_111_Cu_layer, name='Au_111_Cu_layer')

#energy_Au_111_CuO_Layer_0_oxide, intensity_Au_111_CuO_Layer_0_oxide = norm_area("polAvg_Au_111_CuO_Layer_0_oxide.dat")
#shift_Au_111_CuO_Layer_0_oxide = getShift('david_Au_111_CuO_Layer_0_oxide.csv')
#print('Au_111_CuO_Layer_0_oxide total_shift: ' + str(shift_Au_111_CuO_Layer_0_oxide + d_exp))
#energy_Au_111_CuO_Layer_0_oxide += shift_Au_111_CuO_Layer_0_oxide+ d_exp
#x_Au_111_CuO_Layer_0_oxide,y_Au_111_CuO_Layer_0_oxide  = givePeak(energy=energy_Au_111_CuO_Layer_0_oxide, intensity=intensity_Au_111_CuO_Layer_0_oxide, name='Au_111_CuO_Layer_0_oxide')


#energy_Au_111_CuO_Layer_0_solv_oxide, intensity_Au_111_CuO_Layer_0_solv_oxide = norm_area("polAvg_Au_111_CuO_Layer_0_solv_oxide.dat")
#shift_Au_111_CuO_Layer_0_solv_oxide = getShift('david_Au_111_CuO_Layer_0_solv_oxide.csv')
#print('Au_111_CuO_Layer_0_solv_oxide total_shift: ' + str(shift_Au_111_CuO_Layer_0_solv_oxide + d_exp))
#energy_Au_111_CuO_Layer_0_solv_oxide += shift_Au_111_CuO_Layer_0_solv_oxide+ d_exp
#x_Au_111_CuO_Layer_0_solv_oxide,y_Au_111_CuO_Layer_0_solv_oxide  = givePeak(energy=energy_Au_111_CuO_Layer_0_solv_oxide, intensity=intensity_Au_111_CuO_Layer_0_solv_oxide, name='Au_111_CuO_Layer_0_solv_oxide')


energy_Au_111_CuO_Layer_0_solv_oxide_allFrames, intensity_Au_111_CuO_Layer_0_solv_oxide_allFrames = norm_area("polAvg_Au_111_CuO_Layer_0_solv_oxide_allFrames.dat")
shift_Au_111_CuO_Layer_0_solv_oxide_allFrames = getShift('david_Au_111_CuO_Layer_0_solv_oxide.csv')
print('Au_111_CuO_Layer_0_solv_oxide_allFrames total_shift: ' + str(shift_Au_111_CuO_Layer_0_solv_oxide_allFrames + d_exp))
energy_Au_111_CuO_Layer_0_solv_oxide_allFrames += shift_Au_111_CuO_Layer_0_solv_oxide_allFrames+ d_exp
x_Au_111_CuO_Layer_0_solv_oxide_allFrames,y_Au_111_CuO_Layer_0_solv_oxide_allFrames  = givePeak(energy=energy_Au_111_CuO_Layer_0_solv_oxide_allFrames, intensity=intensity_Au_111_CuO_Layer_0_solv_oxide_allFrames, name='Au_111_CuO_Layer_0_solv_oxide_allFrames')

energy_Au_111_CuO_Layer_0_solv_oxide_mag, intensity_Au_111_CuO_Layer_0_solv_oxide_mag_mean, intensity_Au_111_CuO_Layer_0_solv_oxide_mag_std = norm_area_allframes("polAvg_Au_111_CuO_Layer_0_solv_oxide_mag.dat")
shift_Au_111_CuO_Layer_0_solv_oxide_mag = getShift('david_Au_111_CuO_Layer_0_solv_oxide_mag.csv')
print('Au_111_CuO_Layer_0_solv_oxide_mag total_shift: ' + str(shift_Au_111_CuO_Layer_0_solv_oxide_mag + d_exp))
energy_Au_111_CuO_Layer_0_solv_oxide_mag += shift_Au_111_CuO_Layer_0_solv_oxide_mag+ d_exp
x_Au_111_CuO_Layer_0_solv_oxide_mag,y_Au_111_CuO_Layer_0_solv_oxide_mag  = givePeak(energy=energy_Au_111_CuO_Layer_0_solv_oxide_mag, intensity=intensity_Au_111_CuO_Layer_0_solv_oxide_mag_mean, name='Au_111_CuO_Layer_0_solv_oxide_mag')


energy_Au_111_CuO_Layer_1_solv_oxide_mag, intensity_Au_111_CuO_Layer_1_solv_oxide_mag_mean, intensity_Au_111_CuO_Layer_1_solv_oxide_mag_std = norm_area_allframes("polAvg_Au_111_CuO_Layer_1_solv_oxide_mag.dat")
shift_Au_111_CuO_Layer_1_solv_oxide_mag = getShift('david_Au_111_CuO_Layer_1_solv_oxide_mag.csv')
print('Au_111_CuO_Layer_1_solv_oxide_mag total_shift: ' + str(shift_Au_111_CuO_Layer_1_solv_oxide_mag + d_exp))
energy_Au_111_CuO_Layer_1_solv_oxide_mag += shift_Au_111_CuO_Layer_1_solv_oxide_mag+ d_exp
x_Au_111_CuO_Layer_1_solv_oxide_mag,y_Au_111_CuO_Layer_1_solv_oxide_mag  = givePeak(energy=energy_Au_111_CuO_Layer_1_solv_oxide_mag, intensity=intensity_Au_111_CuO_Layer_1_solv_oxide_mag_mean, name='Au_111_CuO_Layer_1_solv_oxide_mag')


energy_Au_111_CuO_Layer_0_solv_oxide_mag_u8, intensity_Au_111_CuO_Layer_0_solv_oxide_mag_u8 = norm_area("polAvg_Au_111_CuO_Layer_0_solv_oxide_mag_u8.dat")
shift_Au_111_CuO_Layer_0_solv_oxide_mag_u8 = getShift('david_Au_111_CuO_Layer_0_solv_oxide_mag_u8.csv')
print('Au_111_CuO_Layer_0_solv_oxide_mag_u8 total_shift: ' + str(shift_Au_111_CuO_Layer_0_solv_oxide_mag_u8 + d_exp))
energy_Au_111_CuO_Layer_0_solv_oxide_mag_u8 += shift_Au_111_CuO_Layer_0_solv_oxide_mag_u8+ d_exp
x_Au_111_CuO_Layer_0_solv_oxide_mag_u8,y_Au_111_CuO_Layer_0_solv_oxide_mag_u8  = givePeak(energy=energy_Au_111_CuO_Layer_0_solv_oxide_mag_u8, intensity=intensity_Au_111_CuO_Layer_0_solv_oxide_mag_u8, name='Au_111_CuO_Layer_0_solv_oxide_mag_u8')


energy_bulkCu2O_exp, intensity_bulkCu2O_exp = exp("exp_bulkCu2O.dat")
peaks_bulkCu2O_exp, _ = scipy.signal.find_peaks(intensity_bulkCu2O_exp, prominence=0.5)
x_bulkCu2O_exp = energy_bulkCu2O_exp[peaks_bulkCu2O_exp[0]]
y_bulkCu2O_exp = intensity_bulkCu2O_exp[peaks_bulkCu2O_exp[0]]

energy_Au_111_Cu2O_Layer_0_solv_oxide, intensity_Au_111_Cu2O_Layer_0_solv_oxide_mean, intensity_Au_111_Cu2O_Layer_0_solv_oxide_std = norm_area_allframes("polAvg_Au_111_Cu2O_Layer_0_solv_oxide.dat")
shift_Au_111_Cu2O_Layer_0_solv_oxide = getShift('david_Au_111_Cu2O_Layer_0_solv_oxide.csv')
print('Au_111_Cu2O_Layer_0_solv_oxide total_shift: ' + str(shift_Au_111_Cu2O_Layer_0_solv_oxide + d_exp))
energy_Au_111_Cu2O_Layer_0_solv_oxide += shift_Au_111_Cu2O_Layer_0_solv_oxide + d_exp
x_Au_111_Cu2O_Layer_0_solv_oxide,y_Au_111_Cu2O_Layer_0_solv_oxide  = givePeak(energy=energy_Au_111_Cu2O_Layer_0_solv_oxide, intensity=intensity_Au_111_Cu2O_Layer_0_solv_oxide_mean, name='Au_111_Cu2O_Layer_0_solv_oxide')

energy_Au_111_Cu2O_Layer_1_solv_oxide_mag, intensity_Au_111_Cu2O_Layer_1_solv_oxide_mag_mean, intensity_Au_111_Cu2O_Layer_1_solv_oxide_mag_std = norm_area_allframes("polAvg_Au_111_Cu2O_Layer_1_solv_oxide_mag.dat")
shift_Au_111_Cu2O_Layer_1_solv_oxide_mag = getShift('david_Au_111_Cu2O_Layer_1_solv_oxide_mag.csv')
print('Au_111_Cu2O_Layer_1_solv_oxide_mag total_shift: ' + str(shift_Au_111_Cu2O_Layer_1_solv_oxide_mag + d_exp))
energy_Au_111_Cu2O_Layer_1_solv_oxide_mag += shift_Au_111_Cu2O_Layer_1_solv_oxide_mag + d_exp
x_Au_111_Cu2O_Layer_1_solv_oxide_mag,y_Au_111_Cu2O_Layer_1_solv_oxide_mag  = givePeak(energy=energy_Au_111_Cu2O_Layer_1_solv_oxide_mag, intensity=intensity_Au_111_Cu2O_Layer_1_solv_oxide_mag_mean, name='Au_111_Cu2O_Layer_1_solv_oxide_mag')

energy_Au_111_Cu2O_Layer_1_solv_oxide_newNomag, intensity_Au_111_Cu2O_Layer_1_solv_oxide_newNomag = norm_area("polAvg_Au_111_Cu2O_Layer_1_solv_oxide_newNomag.dat")
shift_Au_111_Cu2O_Layer_1_solv_oxide_newNomag = getShift('david_Au_111_Cu2O_Layer_1_solv_oxide_newNomag.csv')
print('Au_111_Cu2O_Layer_1_solv_oxide_newNomag total_shift: ' + str(shift_Au_111_Cu2O_Layer_1_solv_oxide_newNomag + d_exp))
energy_Au_111_Cu2O_Layer_1_solv_oxide_newNomag += shift_Au_111_Cu2O_Layer_1_solv_oxide_newNomag+ d_exp
x_Au_111_Cu2O_Layer_1_solv_oxide_newNomag,y_Au_111_Cu2O_Layer_1_solv_oxide_newNomag  = givePeak(energy=energy_Au_111_Cu2O_Layer_1_solv_oxide_newNomag, intensity=intensity_Au_111_Cu2O_Layer_1_solv_oxide_newNomag, name='Au_111_Cu2O_Layer_1_solv_oxide_newNomag')


# Add 1 eV as estimate of incorrect screening of elecrtrons according to DOS and based on dry results
#energy_Au_111_CuO_Layer_0_solv_oxide += 1.0
#x_Au_111_CuO_Layer_0_solv_oxide += 1.0

#energy_Au_111_CuO_Layer_0_solv_oxide_allFrames += 1.0
#x_Au_111_CuO_Layer_0_solv_oxide_allFrames += 1.0

#energy_watbox, intensity_watbox = norm_area("polAvg_watbox.dat")
#shift_watbox = getShift('david_watbox.csv')
#print('watbox total_shift: ' + str(shift_watbox + d_exp))
#energy_watbox += shift_watbox+ d_exp
#x_watbox,y_watbox  = givePeak(energy=energy_watbox, intensity=intensity_watbox, name='watbox')

# Plot settings
plt.figure(figsize=(3.0, 3.5))


#plt.plot(energy_cuso4_solv_sulfate, intensity_cuso4_solv_sulfate+yoffset*0, label='cuso4_solv_sulfate', linewidth=2, color='black', )
#if x_cuso4_solv_sulfate is not None and y_cuso4_solv_sulfate is not None:
#    plt.annotate(f'{x_cuso4_solv_sulfate:.1f}', (x_cuso4_solv_sulfate, y_cuso4_solv_sulfate+yoffset*0),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 

#plt.plot(energy_Au_111_watbox, intensity_Au_111_watbox+yoffset*-0.5, label='Au_111_watbox', linewidth=2, color='#515151')
#if x_Au_111_watbox is not None and y_Au_111_watbox is not None:
#    plt.annotate(f'{x_Au_111_watbox:.1f}', (x_Au_111_watbox, y_Au_111_watbox+yoffset*-0.5),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 

#plt.plot(energy_Au_111_cuso4_sulfate, intensity_Au_111_cuso4_sulfate/3+yoffset*1, label='Au_111_cuso4_sulfate', linewidth=2, color='#F14040', linestyle='--')
#if x_Au_111_cuso4_sulfate is not None and y_Au_111_cuso4_sulfate is not None:
#    plt.annotate(f'{x_Au_111_cuso4_sulfate:.1f}', (x_Au_111_cuso4_sulfate, y_Au_111_cuso4_sulfate/3+yoffset*1),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black')

#plt.plot(energy_Au_111_cuso4_water, intensity_Au_111_cuso4_water/2+yoffset*1, label='Au_111_cuso4_water', linewidth=2, color='#F14040', linestyle=':')
#if x_Au_111_cuso4_water is not None and y_Au_111_cuso4_water is not None:
#    plt.annotate(f'{x_Au_111_cuso4_water:.1f}', (x_Au_111_cuso4_water, y_Au_111_cuso4_water/2+yoffset*1),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 

exp_scaling = 7
yoffset=1/6

#plt.plot(energy_0V_exp, intensity_0V_exp/12 + yoffset*0, label='0 V', linewidth=2, color='#515151')
##if x_0V_exp is not None and y_0V_exp is not None:
##    plt.annotate(f'{x_0V_exp:.1f}', (x_0V_exp, y_0V_exp/10 + yoffset*0),
##                 textcoords="offset points", xytext=(0,2),
##                 ha='center', fontproperties=arial, color='black')
#plt.plot(energy_Au_111_cuso4_combined, intensity_Au_111_cuso4_combined_0V+yoffset*0.2, label='Au_111_cuso4_combined_minus01V', linewidth=2, color='#515151', linestyle=":")
#if x0_Au_111_cuso4_combined_0V is not None and y0_Au_111_cuso4_combined_0V is not None:
#    plt.annotate(f'{x0_Au_111_cuso4_combined_0V:.1f}', (x0_Au_111_cuso4_combined_0V, y0_Au_111_cuso4_combined_0V+yoffset*0.2),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black')
#if x1_Au_111_cuso4_combined_0V is not None and y1_Au_111_cuso4_combined_0V is not None:
#    plt.annotate(f'{x1_Au_111_cuso4_combined_0V:.1f}', (x1_Au_111_cuso4_combined_0V, y1_Au_111_cuso4_combined_0V+yoffset*0.2),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black')
 

#plt.plot(energy_minus01V_exp, intensity_minus01V_exp/9+yoffset*1.5, label='- 0.1 V', linewidth=2, color='#F14040')
##if x_minus01V_exp is not None and y_minus01V_exp is not None:
##    plt.annotate(f'{x_minus01V_exp:.1f}', (x_minus01V_exp, y_minus01V_exp/exp_scaling+yoffset*1),
##    textcoords="offset points", xytext=(0,2), 
##    ha='center', fontproperties=arial, color='black') 
#plt.plot(energy_Au_111_cuso4_combined, intensity_Au_111_cuso4_combined_minus01V+yoffset*1.7, label='Au_111_cuso4_combined_minus01V', linewidth=2, color='#F14040', linestyle=":")
#if x_Au_111_cuso4_combined_minus01V is not None and y_Au_111_cuso4_combined_minus01V is not None:
#    plt.annotate(f'{x_Au_111_cuso4_combined_minus01V:.1f}', (x_Au_111_cuso4_combined_minus01V, y_Au_111_cuso4_combined_minus01V+yoffset*1.7),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black')


plt.plot(energy_minus02V_exp, intensity_minus02V_exp/5.0 + yoffset*0, label='exp -0.2 V', color='#1A6FDF', linewidth=2)
if x_minus02V_exp is not None and y_minus02V_exp is not None:
    print("NOT NONE")
    plt.annotate(f'{x_minus02V_exp:.1f}', (x_minus02V_exp, y_minus02V_exp/5.0 + yoffset*0),
                 textcoords="offset points", xytext=(0,2),
                 ha='center', fontproperties=arial, color='black')
#plt.plot(energy_Au_111_CuO_Layer_0_oxide, intensity_Au_111_CuO_Layer_0_oxide+yoffset*2, label='Au_111_CuO_Layer_0_oxide', linewidth=2, color='#1A6FDF', alpha = 0.5)
#if x_Au_111_CuO_Layer_0_oxide is not None and y_Au_111_CuO_Layer_0_oxide is not None:
#    plt.annotate(f'{x_Au_111_CuO_Layer_0_oxide:.1f}', (x_Au_111_CuO_Layer_0_oxide, y_Au_111_CuO_Layer_0_oxide+yoffset*2),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black')
#plt.plot(energy_Au_111_CuO_Layer_0_solv_oxide, intensity_Au_111_CuO_Layer_0_solv_oxide+yoffset*3.5, label='Au_111_CuO_Layer_0_solv_oxide', linewidth=2, color='#1A6FDF', linestyle=":")
#if x_Au_111_CuO_Layer_0_solv_oxide is not None and y_Au_111_CuO_Layer_0_solv_oxide is not None:
#    plt.annotate(f'{x_Au_111_CuO_Layer_0_solv_oxide:.1f}', (x_Au_111_CuO_Layer_0_solv_oxide, y_Au_111_CuO_Layer_0_solv_oxide+yoffset*3.5),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black')

#plt.plot(energy_Au_111_CuO_Layer_0_solv_oxide_allFrames, intensity_Au_111_CuO_Layer_0_solv_oxide_allFrames+yoffset*1, label='Au_111_CuO_Layer_0_solv_oxide_allFrames', linewidth=2, color='#1A6FDF')
#if x_Au_111_CuO_Layer_0_solv_oxide_allFrames is not None and y_Au_111_CuO_Layer_0_solv_oxide_allFrames is not None:
#    plt.annotate(f'{x_Au_111_CuO_Layer_0_solv_oxide_allFrames:.1f}', (x_Au_111_CuO_Layer_0_solv_oxide_allFrames, y_Au_111_CuO_Layer_0_solv_oxide_allFrames+yoffset*1),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black')

plt.fill_between(energy_Au_111_CuO_Layer_0_solv_oxide_mag, 
    intensity_Au_111_CuO_Layer_0_solv_oxide_mag_mean - intensity_Au_111_CuO_Layer_0_solv_oxide_mag_std + yoffset*1,
    intensity_Au_111_CuO_Layer_0_solv_oxide_mag_mean + intensity_Au_111_CuO_Layer_0_solv_oxide_mag_std + yoffset*1,
    color='#1A6FDF', alpha=0.15)
plt.plot(energy_Au_111_CuO_Layer_0_solv_oxide_mag, intensity_Au_111_CuO_Layer_0_solv_oxide_mag_mean+yoffset*1, label='Au_111_CuO_Layer_0_solv_oxide_mag_u4.9', linewidth=2, color='#1A6FDF')
if x_Au_111_CuO_Layer_0_solv_oxide_mag is not None and y_Au_111_CuO_Layer_0_solv_oxide_mag is not None:
    plt.annotate(f'{x_Au_111_CuO_Layer_0_solv_oxide_mag:.1f}', (x_Au_111_CuO_Layer_0_solv_oxide_mag, y_Au_111_CuO_Layer_0_solv_oxide_mag+yoffset*1),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black')

plt.fill_between(energy_Au_111_CuO_Layer_1_solv_oxide_mag,
    intensity_Au_111_CuO_Layer_1_solv_oxide_mag_mean - intensity_Au_111_CuO_Layer_1_solv_oxide_mag_std + yoffset*2,
    intensity_Au_111_CuO_Layer_1_solv_oxide_mag_mean + intensity_Au_111_CuO_Layer_1_solv_oxide_mag_std + yoffset*2,
    color='#1A6FDF', alpha=0.15)
plt.plot(energy_Au_111_CuO_Layer_1_solv_oxide_mag, intensity_Au_111_CuO_Layer_1_solv_oxide_mag_mean+yoffset*2, label='Au_111_CuO_Layer_1_solv_oxide_mag_u4.9', linewidth=2, color='#1A6FDF')
if x_Au_111_CuO_Layer_1_solv_oxide_mag is not None and y_Au_111_CuO_Layer_1_solv_oxide_mag is not None:
    plt.annotate(f'{x_Au_111_CuO_Layer_1_solv_oxide_mag:.1f}', (x_Au_111_CuO_Layer_1_solv_oxide_mag, y_Au_111_CuO_Layer_1_solv_oxide_mag+yoffset*2),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black')

#plt.plot(energy_Au_111_CuO_Layer_0_solv_oxide_mag_u8, intensity_Au_111_CuO_Layer_0_solv_oxide_mag_u8+yoffset*4, label='Au_111_CuO_Layer_0_solv_oxide_mag_u8', linewidth=2, color='#1A6FDF')
#if x_Au_111_CuO_Layer_0_solv_oxide_mag_u8 is not None and y_Au_111_CuO_Layer_0_solv_oxide_mag_u8 is not None:
#    plt.annotate(f'{x_Au_111_CuO_Layer_0_solv_oxide_mag_u8:.1f}', (x_Au_111_CuO_Layer_0_solv_oxide_mag_u8, y_Au_111_CuO_Layer_0_solv_oxide_mag_u8+yoffset*4),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black')

plt.fill_between(energy_Au_111_Cu2O_Layer_0_solv_oxide,
    intensity_Au_111_Cu2O_Layer_0_solv_oxide_mean - intensity_Au_111_Cu2O_Layer_0_solv_oxide_std + yoffset*3,
    intensity_Au_111_Cu2O_Layer_0_solv_oxide_mean + intensity_Au_111_Cu2O_Layer_0_solv_oxide_std + yoffset*3,
    color='#1A6FDF', alpha=0.15)
plt.plot(energy_Au_111_Cu2O_Layer_0_solv_oxide, intensity_Au_111_Cu2O_Layer_0_solv_oxide_mean+yoffset*3, label='Au_111_Cu2O_Layer_0_solv_oxide', linewidth=2, color='#1A6FDF')
if x_Au_111_Cu2O_Layer_0_solv_oxide is not None and y_Au_111_Cu2O_Layer_0_solv_oxide is not None:
    plt.annotate(f'{x_Au_111_Cu2O_Layer_0_solv_oxide:.1f}', (x_Au_111_Cu2O_Layer_0_solv_oxide, y_Au_111_Cu2O_Layer_0_solv_oxide+yoffset*3),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black') 

plt.fill_between(energy_Au_111_Cu2O_Layer_1_solv_oxide_mag,
    intensity_Au_111_Cu2O_Layer_1_solv_oxide_mag_mean - intensity_Au_111_Cu2O_Layer_1_solv_oxide_mag_std + yoffset*4.5,
    intensity_Au_111_Cu2O_Layer_1_solv_oxide_mag_mean + intensity_Au_111_Cu2O_Layer_1_solv_oxide_mag_std + yoffset*4.5,
    color='#1A6FDF', alpha=0.15)
plt.plot(energy_Au_111_Cu2O_Layer_1_solv_oxide_mag, intensity_Au_111_Cu2O_Layer_1_solv_oxide_mag_mean+yoffset*4.5, label=' ', linewidth=2, color='#1A6FDF')
if x_Au_111_Cu2O_Layer_1_solv_oxide_mag is not None and y_Au_111_Cu2O_Layer_1_solv_oxide_mag is not None:
    plt.annotate(f'{x_Au_111_Cu2O_Layer_1_solv_oxide_mag:.1f}', (x_Au_111_Cu2O_Layer_1_solv_oxide_mag, y_Au_111_Cu2O_Layer_1_solv_oxide_mag+yoffset*4.5),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black') 
#plt.plot(energy_Au_111_Cu2O_Layer_1_solv_oxide_newNomag, intensity_Au_111_Cu2O_Layer_1_solv_oxide_newNomag+yoffset*4.5, label='Au_111_Cu2O_Layer_1_solv_oxide_newNomag', linewidth=2, color='red', )
#if x_Au_111_Cu2O_Layer_1_solv_oxide_newNomag is not None and y_Au_111_Cu2O_Layer_1_solv_oxide_newNomag is not None:
#    plt.annotate(f'{x_Au_111_Cu2O_Layer_1_solv_oxide_newNomag:.1f}', (x_Au_111_Cu2O_Layer_1_solv_oxide_newNomag, y_Au_111_Cu2O_Layer_1_solv_oxide_newNomag+yoffset*4.5),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 

plt.plot(energy_bulkCu2O_exp, intensity_bulkCu2O_exp / 9.0 +yoffset*5.8, label='bulkCu2O exp', linewidth=2, color='#1A6FDF')
if x_bulkCu2O_exp is not None and y_bulkCu2O_exp is not None:
    plt.annotate(f'{x_bulkCu2O_exp:.1f}', (x_bulkCu2O_exp, y_bulkCu2O_exp / 9.0 +yoffset*5.8),
    textcoords="offset points", xytext=(0,2),
    ha='center', fontproperties=arial, color='black')

# Shaded envelope for bulkCu2O using per-frame std from norm_area_allframes
plt.fill_between(energy_bulkCu2O,
    (intensity_bulkCu2O_mean - intensity_bulkCu2O_std) + yoffset*7.5,
    (intensity_bulkCu2O_mean + intensity_bulkCu2O_std) + yoffset*7.5,
    color='#1A6FDF', alpha=0.15)
plt.plot(energy_bulkCu2O, intensity_bulkCu2O_mean + yoffset*7.5, label='bulkCu2O exp', linewidth=2, color='#1A6FDF')
if x_bulkCu2O is not None and y_bulkCu2O is not None:
    plt.annotate(f'{x_bulkCu2O:.1f}', (x_bulkCu2O, y_bulkCu2O + yoffset*7.5),
    textcoords="offset points", xytext=(0,2),
    ha='center', fontproperties=arial, color='black')

plt.plot(energy_bulkCuO_exp, intensity_bulkCuO_exp / 8.0 +yoffset*8.5, label='bulkCuO exp', linewidth=2, color='#1A6FDF')
if x_bulkCuO_exp is not None and y_bulkCuO_exp is not None:
    plt.annotate(f'{x_bulkCuO_exp:.1f}', (x_bulkCuO_exp, y_bulkCuO_exp / 8.0 +yoffset*8.5),
    textcoords="offset points", xytext=(0,2),
    ha='center', fontproperties=arial, color='black')

# Shaded envelope for bulkCuO_mag using per-frame std from norm_area_allframes
plt.fill_between(energy_bulkCuO_mag,
    (intensity_bulkCuO_mag_mean - intensity_bulkCuO_mag_std) + yoffset*10.5,
    (intensity_bulkCuO_mag_mean + intensity_bulkCuO_mag_std) + yoffset*10.5,
    color='#1A6FDF', alpha=0.15)
plt.plot(energy_bulkCuO_mag, intensity_bulkCuO_mag_mean + yoffset*10.5, label='bulkCuO_mag exp', linewidth=2, color='#1A6FDF')
if x_bulkCuO_mag is not None and y_bulkCuO_mag is not None:
    plt.annotate(f'{x_bulkCuO_mag:.1f}', (x_bulkCuO_mag, y_bulkCuO_mag + yoffset*10.5),
    textcoords="offset points", xytext=(0,2),
    ha='center', fontproperties=arial, color='black')



#plt.plot(energy_minus03V_exp, intensity_minus03V_exp/6 + yoffset*7.0, label='-0.3 V', linewidth=2, color='#37AD6B')
##if x_minus03V_exp is not None and y_minus03V_exp is not None:
##    plt.annotate(f'{x_minus03V_exp:.1f}', (x_minus03V_exp, y_minus03V_exp/6 + yoffset*3),
##                 textcoords="offset points", xytext=(0,2),
##                 ha='center', fontproperties=arial, color='black')
#plt.plot(energy_Au_111_cuso4_combined, intensity_Au_111_cuso4_combined_minus03V+yoffset*7.5, label='Au_111_cuso4_combined_minus03V', linewidth=2, color='#37AD6B', linestyle=":")
#if x_Au_111_cuso4_combined_minus03V is not None and y_Au_111_cuso4_combined_minus03V is not None:
#    plt.annotate(f'{x_Au_111_cuso4_combined_minus03V:.1f}', (x_Au_111_cuso4_combined_minus03V, y_Au_111_cuso4_combined_minus03V+yoffset*7.5),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black')



#plt.plot(energy_Au_111_Cu_layer, intensity_Au_111_Cu_layer+yoffset*7.5, label='Au_111_Cu_layer', linewidth=2, color='#37AD6B', linestyle=":")
#if x_Au_111_Cu_layer is not None and y_Au_111_Cu_layer is not None:
#    plt.annotate(f'{x_Au_111_Cu_layer:.1f}', (x_Au_111_Cu_layer, y_Au_111_Cu_layer+yoffset*7.5),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 

# Plot reference 

#plt.plot(energy_bulkCuO, intensity_bulkCuO+yoffset*14, label='bulkCuO', linewidth=2, color='midnightblue', )
#if x_bulkCuO is not None and y_bulkCuO is not None:
#    plt.annotate(f'{x_bulkCuO:.1f}', (x_bulkCuO, y_bulkCuO+yoffset*14),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 

# Plot settings
plt.xlabel('Energy (eV)', fontproperties=arialbd)
#tick_positions = np.arange(928, 945, 1)
#tick_labels = [str(tick) if tick % 2 == 0 else '' for tick in tick_positions]
#plt.xticks(tick_positions, tick_labels,fontproperties=arial)
plt.ylabel('Intensity (a.u.)', fontproperties=arialbd)
#plt.yticks(fontproperties=arial)
plt.gca().set_yticklabels([])  # Removes y-axis labels, keeps tick marks
plt.xticks(fontproperties=arialbd)
plt.xlim(525, 550)
plt.axvline(532.7, color='black', linestyle=':', linewidth=1)
#plt.xlim(E_i+d_exp, E_f+d_exp)
#plt.grid(True)
ax = plt.gca()
# Annotate each labeled line near the right edge instead of using a legend
# Compute a suitable x position a bit further inside the right axis limit so
# the text doesn't run outside the axes. Place the annotation to the left
# of the reference point (negative x offset) and right-align the text.
x_min, x_max = ax.get_xlim()
# Move labels 6% of the axis width inside the right edge (was 1% before)
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
#                xytext=(5, 10),
#                textcoords='offset points',
#                fontproperties=arialbd,
#                va='center',
#                ha='right')

# Make layout tight so annotations aren't clipped by figure margins
plt.tight_layout()
# no legend box

# Saving or showing plot
plt.savefig('compare-o.svg', format='svg')
plt.show()
