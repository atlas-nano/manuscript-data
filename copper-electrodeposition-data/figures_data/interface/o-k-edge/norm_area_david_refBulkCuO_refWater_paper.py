
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
def renorm(energy, intensity, Ei, Ef):
    mask = (energy >= Ei) & (energy <= Ef)
    energy_range = energy[mask]
    intensity_range = intensity[mask]
    auc = simpson(intensity_range, energy_range)
    intensity = intensity / auc
    return intensity
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
energy_bulkCuO_exp, intensity_bulkCuO_exp = exp("exp_bulkCuO.dat")
# Reverse the order of this data
energy_bulkCuO_exp = energy_bulkCuO_exp[::-1]
intensity_bulkCuO_exp = intensity_bulkCuO_exp[::-1]
x_bulkCuO_exp,y_bulkCuO_exp  = givePeak(energy=energy_bulkCuO_exp, intensity=intensity_bulkCuO_exp, name='bulkCuO')
# Computational reference
energy_bulkCuO_mag, intensity_bulkCuO_mag = norm_area("polAvg_bulkCuO_mag.dat")
x_bulkCuO_mag,y_bulkCuO_mag  = givePeak(energy=energy_bulkCuO_mag, intensity=intensity_bulkCuO_mag, name='bulkCuO_mag')
print(f'comp bulkCuO_mag: {x_bulkCuO_mag}')
shift_bulkCuO_mag = getShift('david_bulkCuO_mag.csv')
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Experimental shift
d_exp = x_bulkCuO_exp - (x_bulkCuO_mag+shift_bulkCuO_mag)
print(f'd_exp = {d_exp}')
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::
print('bulkCuO_mag total shift: ' + str(shift_bulkCuO_mag + d_exp))
energy_bulkCuO_mag += shift_bulkCuO_mag + d_exp
x_bulkCuO_mag,y_bulkCuO_mag  = givePeak(energy=energy_bulkCuO_mag, intensity=intensity_bulkCuO_mag, name='bulkCuO_mag')






# Second reference: water ****************************
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
d_exp_watbox = x1_watbox_exp - (x_watbox+shift_watbox)
print(f'd_exp = {d_exp}')
# ------------------------------------------------------------------------------------------------------------
print('watbox total shift: ' + str(shift_watbox + d_exp))
energy_watbox += shift_watbox + d_exp
x_watbox,y_watbox  = givePeak(energy=energy_watbox, intensity=intensity_watbox, name='watbox')





# Experimental data

energy_0V_exp, intensity_0V_exp = exp("exp_0V.dat")
x_0V_exp, y_0V_exp = givePeak(energy=energy_0V_exp, intensity=intensity_0V_exp, name='0V')

energy_minus01V_exp, intensity_minus01V_exp = exp("exp_-0.1V.dat")
x_minus01V_exp,y_minus01V_exp  = givePeak(energy=energy_minus01V_exp, intensity=intensity_minus01V_exp, name='minus01V')

energy_minus02V_exp, intensity_minus02V_exp = exp("exp_-0.2V.dat")
#peaks_minus02V_exp, _ = scipy.signal.find_peaks(intensity_minus02V_exp, prominence=0.15)
#x_minus02V_exp = energy_minus02V_exp[peaks_minus02V_exp[0]]
#y_minus02V_exp = intensity_minus02V_exp[peaks_minus02V_exp[0]]
x_minus02V_exp, y_minus02V_exp = givePeak(energy=energy_minus02V_exp, intensity=intensity_minus02V_exp, name='minus02V')

energy_minus03V_exp, intensity_minus03V_exp = exp("exp_-0.3V.dat")
x_minus03V_exp, y_minus03V_exp = givePeak(energy=energy_minus03V_exp, intensity=intensity_minus03V_exp, name='minus03V')

energy_bulkCu2O_exp, intensity_bulkCu2O_exp = exp("exp_bulkCu2O.dat")
peaks_bulkCu2O_exp, _ = scipy.signal.find_peaks(intensity_bulkCu2O_exp, prominence=0.5)
x_bulkCu2O_exp = energy_bulkCu2O_exp[peaks_bulkCu2O_exp[0]]
y_bulkCu2O_exp = intensity_bulkCu2O_exp[peaks_bulkCu2O_exp[0]]

# Process data

#energy_Au_111_cuso4_sulfate, intensity_Au_111_cuso4_sulfate = norm_area("polAvg_Au_111_cuso4_sulfate.dat")
#shift_Au_111_cuso4_sulfate = getShift('david_Au_111_cuso4_sulfate.csv')
#print('Au_111_cuso4_sulfate total_shift: ' + str(shift_Au_111_cuso4_sulfate + d_exp_watbox))
#energy_Au_111_cuso4_sulfate += shift_Au_111_cuso4_sulfate+ d_exp_watbox
#x_Au_111_cuso4_sulfate,y_Au_111_cuso4_sulfate  = givePeak(energy=energy_Au_111_cuso4_sulfate, intensity=intensity_Au_111_cuso4_sulfate, name='Au_111_cuso4_sulfate')
#
#energy_Au_111_cuso4_water, intensity_Au_111_cuso4_water = norm_area("polAvg_Au_111_cuso4_water.dat")
#shift_Au_111_cuso4_water = getShift('david_Au_111_cuso4_water.csv')
#print('Au_111_cuso4_water water_shift: ' + str(shift_Au_111_cuso4_water + d_exp_watbox))
#energy_Au_111_cuso4_water += shift_Au_111_cuso4_water+ d_exp_watbox
#x_Au_111_cuso4_water,y_Au_111_cuso4_water  = givePeak(energy=energy_Au_111_cuso4_water, intensity=intensity_Au_111_cuso4_water, name='Au_111_cuso4_water')
#
#energy_Au_111_2h3o_so4_hydronium_allFrames, intensity_Au_111_2h3o_so4_hydronium_allFrames = norm_area("polAvg_Au_111_2h3o_so4_hydronium_allFrames.dat")
#shift_Au_111_2h3o_so4_hydronium_allFrames = getShift('david_Au_111_2h3o_so4_hydronium256.csv')
#print('Au_111_2h3o_so4_hydronium_allFrames water_shift: ' + str(shift_Au_111_2h3o_so4_hydronium_allFrames + d_exp_watbox))
#energy_Au_111_2h3o_so4_hydronium_allFrames += shift_Au_111_2h3o_so4_hydronium_allFrames+ d_exp_watbox
#x_Au_111_2h3o_so4_hydronium_allFrames,y_Au_111_2h3o_so4_hydronium_allFrames  = givePeak(energy=energy_Au_111_2h3o_so4_hydronium_allFrames, intensity=intensity_Au_111_2h3o_so4_hydronium_allFrames, name='Au_111_2h3o_so4_hydronium_allFrames')


# Computational Au_111_2h3o_so4_sulfate_allFrames (load, shift, peak)
#energy_Au_111_2h3o_so4_sulfate_allFrames, intensity_Au_111_2h3o_so4_sulfate_allFrames = norm_area("polAvg_Au_111_2h3o_so4_sulfate_allFrames.dat")
#shift_Au_111_2h3o_so4_sulfate_allFrames = getShift('david_Au_111_2h3o_so4_sulfate.csv')
#print('Au_111_2h3o_so4_sulfate_allFrames total_shift: ' + str(shift_Au_111_2h3o_so4_sulfate_allFrames + d_exp))
#energy_Au_111_2h3o_so4_sulfate_allFrames += shift_Au_111_2h3o_so4_sulfate_allFrames + d_exp_watbox
#x_Au_111_2h3o_so4_sulfate_allFrames,y_Au_111_2h3o_so4_sulfate_allFrames  = givePeak(energy=energy_Au_111_2h3o_so4_sulfate_allFrames, intensity=intensity_Au_111_2h3o_so4_sulfate_allFrames, name='Au_111_2h3o_so4_sulfate_allFrames')


energy_2h3o_so4_sulfate_allFrames, intensity_2h3o_so4_sulfate_allFrames = norm_area("polAvg_2h3o_so4_sulfate_allFrames.dat")
shift_2h3o_so4_sulfate_allFrames = getShift('david_2h3o_so4_sulfate.csv')
print('2h3o_so4_sulfate_allFrames total_shift: ' + str(shift_2h3o_so4_sulfate_allFrames + d_exp_watbox))
energy_2h3o_so4_sulfate_allFrames += shift_2h3o_so4_sulfate_allFrames+ d_exp_watbox
peaks_2h3o_so4_sulfate_allFrames, _ = scipy.signal.find_peaks(intensity_2h3o_so4_sulfate_allFrames, prominence=0.01)
x0_2h3o_so4_sulfate_allFrames = energy_2h3o_so4_sulfate_allFrames[peaks_2h3o_so4_sulfate_allFrames[0]]
y0_2h3o_so4_sulfate_allFrames = intensity_2h3o_so4_sulfate_allFrames[peaks_2h3o_so4_sulfate_allFrames[0]]
#x1_2h3o_so4_sulfate_allFrames = energy_2h3o_so4_sulfate_allFrames[peaks_2h3o_so4_sulfate_allFrames[1]]
#y1_2h3o_so4_sulfate_allFrames = intensity_2h3o_so4_sulfate_allFrames[peaks_2h3o_so4_sulfate_allFrames[1]]



# ******************************************************

# Create a common energy grid using experimental watbox and computed sulfate-allFrames
energy_common = np.linspace(min(energy_watbox_exp.min(), energy_2h3o_so4_sulfate_allFrames.min()),
                            max(energy_watbox_exp.max(), energy_2h3o_so4_sulfate_allFrames.max()),
                            num=1000)
Ei = 530
Ef = 538
# Renormalize the requested reference spectra to the same area between Ei and Ef
intensity_2h3o_so4_sulfate_allFrames = renorm(energy_2h3o_so4_sulfate_allFrames, intensity_2h3o_so4_sulfate_allFrames, Ei, Ef)
intensity_watbox_exp = renorm(energy_watbox_exp, intensity_watbox_exp, Ei, Ef)

# helper: interpolate a reference spectrum onto the experimental energy grid
def interp(energy_from, spectrum, energy_to):
    return np.interp(energy_to, energy_from, spectrum)

# Interpolate both spectra onto the common grid
#interp_watbox_func = interp1d(energy_watbox_exp, intensity_watbox_exp, kind='linear', bounds_error=False, fill_value=0)
#interp_sulfate_allframes_func = interp1d(energy_2h3o_so4_sulfate_allFrames, intensity_2h3o_so4_sulfate_allFrames, kind='linear', bounds_error=False, fill_value=0)

intensity_watbox_exp_interp = interp(energy_watbox_exp, intensity_watbox_exp, energy_common)
intensity_2h3o_so4_sulfate_allFrames_interp = interp(energy_2h3o_so4_sulfate_allFrames, intensity_2h3o_so4_sulfate_allFrames, energy_common)




# Weights (0V) from the LCA analysis
weight_water_exp_0V = 0.446
weight_sulfate_0V = 0.554



intensity_Au_111_cuso4_combined_0V = intensity_watbox_exp_interp*weight_water_exp_0V + intensity_2h3o_so4_sulfate_allFrames_interp*weight_sulfate_0V
peaks_Au_111_cuso4_combined_0V, _ = scipy.signal.find_peaks(intensity_Au_111_cuso4_combined_0V, prominence=0.0001)
x0_Au_111_cuso4_combined_0V = energy_common[peaks_Au_111_cuso4_combined_0V[0]]
y0_Au_111_cuso4_combined_0V = intensity_Au_111_cuso4_combined_0V[peaks_Au_111_cuso4_combined_0V[0]]
x1_Au_111_cuso4_combined_0V = energy_common[peaks_Au_111_cuso4_combined_0V[1]]
y1_Au_111_cuso4_combined_0V = intensity_Au_111_cuso4_combined_0V[peaks_Au_111_cuso4_combined_0V[1]]

intensity_Au_111_cuso4_combined_minus01V = intensity_watbox_exp_interp*weight_water_exp_0V + intensity_2h3o_so4_sulfate_allFrames_interp*weight_sulfate_0V
peaks_Au_111_cuso4_combined_minus01V, _ = scipy.signal.find_peaks(intensity_Au_111_cuso4_combined_minus01V, prominence=0.0001)
x0_Au_111_cuso4_combined_minus01V = energy_common[peaks_Au_111_cuso4_combined_minus01V[0]]
y0_Au_111_cuso4_combined_minus01V = intensity_Au_111_cuso4_combined_minus01V[peaks_Au_111_cuso4_combined_minus01V[0]]
x1_Au_111_cuso4_combined_minus01V = energy_common[peaks_Au_111_cuso4_combined_minus01V[1]]
y1_Au_111_cuso4_combined_minus01V = intensity_Au_111_cuso4_combined_minus01V[peaks_Au_111_cuso4_combined_minus01V[1]]

intensity_Au_111_cuso4_combined_minus03V = intensity_watbox_exp_interp*weight_water_exp_0V + intensity_2h3o_so4_sulfate_allFrames_interp*weight_sulfate_0V
peaks_Au_111_cuso4_combined_minus03V, _ = scipy.signal.find_peaks(intensity_Au_111_cuso4_combined_minus03V, prominence=0.0001)
x0_Au_111_cuso4_combined_minus03V = energy_common[peaks_Au_111_cuso4_combined_minus03V[0]]
y0_Au_111_cuso4_combined_minus03V = intensity_Au_111_cuso4_combined_minus03V[peaks_Au_111_cuso4_combined_minus03V[0]]
x1_Au_111_cuso4_combined_minus03V = energy_common[peaks_Au_111_cuso4_combined_minus03V[1]]
y1_Au_111_cuso4_combined_minus03V = intensity_Au_111_cuso4_combined_minus03V[peaks_Au_111_cuso4_combined_minus03V[1]]


energy_Au_111_Cu_layer, intensity_Au_111_Cu_layer = norm_area("polAvg_Au_111_Cu_layer.dat")
shift_Au_111_Cu_layer = getShift('david_Au_111_Cu_layer.csv')
print('Au_111_Cu_layer total_shift: ' + str(shift_Au_111_Cu_layer + d_exp))
energy_Au_111_Cu_layer += shift_Au_111_Cu_layer+ d_exp
x_Au_111_Cu_layer,y_Au_111_Cu_layer  = givePeak(energy=energy_Au_111_Cu_layer, intensity=intensity_Au_111_Cu_layer, name='Au_111_Cu_layer')

energy_Au_111_CuO_Layer_0_oxide, intensity_Au_111_CuO_Layer_0_oxide = norm_area("polAvg_Au_111_CuO_Layer_0_oxide.dat")
shift_Au_111_CuO_Layer_0_oxide = getShift('david_Au_111_CuO_Layer_0_oxide.csv')
print('Au_111_CuO_Layer_0_oxide total_shift: ' + str(shift_Au_111_CuO_Layer_0_oxide + d_exp))
energy_Au_111_CuO_Layer_0_oxide += shift_Au_111_CuO_Layer_0_oxide+ d_exp
x_Au_111_CuO_Layer_0_oxide,y_Au_111_CuO_Layer_0_oxide  = givePeak(energy=energy_Au_111_CuO_Layer_0_oxide, intensity=intensity_Au_111_CuO_Layer_0_oxide, name='Au_111_CuO_Layer_0_oxide')

energy_Au_111_CuO_Layer_0_solv_oxide, intensity_Au_111_CuO_Layer_0_solv_oxide = norm_area("polAvg_Au_111_CuO_Layer_0_solv_oxide.dat")
shift_Au_111_CuO_Layer_0_solv_oxide = getShift('david_Au_111_CuO_Layer_0_solv_oxide.csv')
print('Au_111_CuO_Layer_0_solv_oxide total_shift: ' + str(shift_Au_111_CuO_Layer_0_solv_oxide + d_exp))
energy_Au_111_CuO_Layer_0_solv_oxide += shift_Au_111_CuO_Layer_0_solv_oxide+ d_exp
x_Au_111_CuO_Layer_0_solv_oxide,y_Au_111_CuO_Layer_0_solv_oxide  = givePeak(energy=energy_Au_111_CuO_Layer_0_solv_oxide, intensity=intensity_Au_111_CuO_Layer_0_solv_oxide, name='Au_111_CuO_Layer_0_solv_oxide')

energy_Au_111_CuO_Layer_0_solv_oxide_allFrames, intensity_Au_111_CuO_Layer_0_solv_oxide_allFrames = norm_area("polAvg_Au_111_CuO_Layer_0_solv_oxide_allFrames.dat")
shift_Au_111_CuO_Layer_0_solv_oxide_allFrames = getShift('david_Au_111_CuO_Layer_0_solv_oxide.csv')
print('Au_111_CuO_Layer_0_solv_oxide_allFrames total_shift: ' + str(shift_Au_111_CuO_Layer_0_solv_oxide_allFrames + d_exp))
energy_Au_111_CuO_Layer_0_solv_oxide_allFrames += shift_Au_111_CuO_Layer_0_solv_oxide_allFrames+ d_exp
x_Au_111_CuO_Layer_0_solv_oxide_allFrames,y_Au_111_CuO_Layer_0_solv_oxide_allFrames  = givePeak(energy=energy_Au_111_CuO_Layer_0_solv_oxide_allFrames, intensity=intensity_Au_111_CuO_Layer_0_solv_oxide_allFrames, name='Au_111_CuO_Layer_0_solv_oxide_allFrames')

energy_Au_111_Cu2O_Layer_1_solv_oxide_mag, intensity_Au_111_Cu2O_Layer_1_solv_oxide_mag = norm_area("polAvg_Au_111_Cu2O_Layer_1_solv_oxide_mag.dat")
shift_Au_111_Cu2O_Layer_1_solv_oxide_mag = getShift('david_Au_111_Cu2O_Layer_1_solv_oxide_mag.csv')
print('Au_111_Cu2O_Layer_1_solv_oxide_mag total_shift: ' + str(shift_Au_111_Cu2O_Layer_1_solv_oxide_mag + d_exp))
energy_Au_111_Cu2O_Layer_1_solv_oxide_mag += shift_Au_111_Cu2O_Layer_1_solv_oxide_mag+ d_exp
x_Au_111_Cu2O_Layer_1_solv_oxide_mag,y_Au_111_Cu2O_Layer_1_solv_oxide_mag  = givePeak(energy=energy_Au_111_Cu2O_Layer_1_solv_oxide_mag, intensity=intensity_Au_111_Cu2O_Layer_1_solv_oxide_mag, name='Au_111_Cu2O_Layer_1_solv_oxide_mag')



# Plot settings
plt.figure(figsize=(2.5, 4.0))


exp_scaling = 7
yoffset=1/6
interp_factor = 1/3

plt.plot(energy_0V_exp, intensity_0V_exp/13 + yoffset*0, label='0 V', linewidth=2, color='#515151')
#if x_0V_exp is not None and y_0V_exp is not None:
#    plt.annotate(f'{x_0V_exp:.1f}', (x_0V_exp, y_0V_exp/10 + yoffset*0),
#                 textcoords="offset points", xytext=(0,2),
#                 ha='center', fontproperties=arial, color='black')
plt.plot(energy_common, intensity_Au_111_cuso4_combined_0V*interp_factor+yoffset*0.5, label='Au_111_cuso4_combined_minus01V', linewidth=2, color='#515151', alpha=0.5)
if x0_Au_111_cuso4_combined_0V is not None and y0_Au_111_cuso4_combined_0V is not None:
    plt.annotate(f'{x0_Au_111_cuso4_combined_0V:.1f}', (x0_Au_111_cuso4_combined_0V-2, y0_Au_111_cuso4_combined_0V*interp_factor+yoffset*0.5),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black')
if x1_Au_111_cuso4_combined_0V is not None and y1_Au_111_cuso4_combined_0V is not None:
    plt.annotate(f'{x1_Au_111_cuso4_combined_0V:.1f}', (x1_Au_111_cuso4_combined_0V+1, y1_Au_111_cuso4_combined_0V*interp_factor+yoffset*0.5),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black')
 

plt.plot(energy_minus01V_exp, intensity_minus01V_exp/11+yoffset*1.5, label='- 0.1 V', linewidth=2, color='#F14040')
#if x_minus01V_exp is not None and y_minus01V_exp is not None:
#    plt.annotate(f'{x_minus01V_exp:.1f}', (x_minus01V_exp, y_minus01V_exp/exp_scaling+yoffset*1),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 
plt.plot(energy_common, intensity_Au_111_cuso4_combined_minus01V*interp_factor+yoffset*1.9, label='Au_111_cuso4_combined_minus01V', color='#F14040', alpha=0.5, linewidth=2)
if x0_Au_111_cuso4_combined_minus01V is not None and y0_Au_111_cuso4_combined_minus01V is not None:
    plt.annotate(f'{x0_Au_111_cuso4_combined_minus01V:.1f}', (x0_Au_111_cuso4_combined_minus01V-2, y0_Au_111_cuso4_combined_minus01V*interp_factor+yoffset*1.9),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black')
if x1_Au_111_cuso4_combined_minus01V is not None and y1_Au_111_cuso4_combined_minus01V is not None:
    plt.annotate(f'{x1_Au_111_cuso4_combined_minus01V:.1f}', (x1_Au_111_cuso4_combined_minus01V+1, y1_Au_111_cuso4_combined_minus01V*interp_factor+yoffset*1.9),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black')


plt.plot(energy_minus02V_exp, intensity_minus02V_exp/5.0 + yoffset*3.1, label='-0.2 V', color='#1A6FDF', linewidth=2)
#if x_minus02V_exp is not None and y_minus02V_exp is not None:
#    plt.annotate(f'{x_minus02V_exp:.1f}', (x_minus02V_exp, y_minus02V_exp/exp_scaling + yoffset*2),
#                 textcoords="offset points", xytext=(0,2),
#                 ha='center', fontproperties=arial, color='black')

plt.plot(energy_Au_111_Cu2O_Layer_1_solv_oxide_mag, intensity_Au_111_Cu2O_Layer_1_solv_oxide_mag+yoffset*3.5, label='Au_111_Cu2O_Layer_1_solv_oxide_mag', linewidth=2, color='#1A6FDF', alpha=0.5)
if x_Au_111_Cu2O_Layer_1_solv_oxide_mag is not None and y_Au_111_Cu2O_Layer_1_solv_oxide_mag is not None:
    plt.annotate(f'{x_Au_111_Cu2O_Layer_1_solv_oxide_mag:.1f}', (x_Au_111_Cu2O_Layer_1_solv_oxide_mag-2, y_Au_111_Cu2O_Layer_1_solv_oxide_mag+yoffset*3.5),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black') 
plt.plot(energy_bulkCu2O_exp, intensity_bulkCu2O_exp / 12.0 +yoffset*3.5, label='bulkCu2O exp', linewidth=2, color='#1A6FDF', linestyle=':')
if x_bulkCu2O_exp is not None and y_bulkCu2O_exp is not None:
    plt.annotate(f'{x_bulkCu2O_exp:.1f}', (x_bulkCu2O_exp+3, y_bulkCu2O_exp / 12.0 +yoffset*(3.5-0.2)),
    textcoords="offset points", xytext=(0,2),
    ha='center', fontproperties=arial, color='black')


plt.plot(energy_minus03V_exp, intensity_minus03V_exp/8 + yoffset*7.0, label='-0.3 V', linewidth=2, color='#37AD6B')
#if x_minus03V_exp is not None and y_minus03V_exp is not None:
#    plt.annotate(f'{x_minus03V_exp:.1f}', (x_minus03V_exp, y_minus03V_exp/6 + yoffset*3),
#                 textcoords="offset points", xytext=(0,2),
#                 ha='center', fontproperties=arial, color='black')
plt.plot(energy_common, intensity_Au_111_cuso4_combined_minus03V*interp_factor+yoffset*7.5, label='Au_111_cuso4_combined_minus03V', linewidth=2, color='#37AD6B', alpha=0.5)
if x0_Au_111_cuso4_combined_minus03V is not None and y0_Au_111_cuso4_combined_minus03V is not None:
    plt.annotate(f'{x0_Au_111_cuso4_combined_minus03V:.1f}', (x0_Au_111_cuso4_combined_minus03V-2, y0_Au_111_cuso4_combined_minus03V*interp_factor+yoffset*7.5),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black')
if x1_Au_111_cuso4_combined_minus03V is not None and y1_Au_111_cuso4_combined_minus03V is not None:
    plt.annotate(f'{x1_Au_111_cuso4_combined_minus03V:.1f}', (x1_Au_111_cuso4_combined_minus03V+1, y1_Au_111_cuso4_combined_minus03V*interp_factor+yoffset*7.5),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black')



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
#plt.xlim(E_i+d_exp, E_f+d_exp)
#plt.grid(True)
#handles, labels = plt.gca().get_legend_handles_labels()
#plt.legend(handles[::-1], labels[::-1],loc='upper left', bbox_to_anchor=(1, 1))
#plt.tight_layout()

# Saving or showing plot
plt.savefig('main-o.svg', format='svg')
#plt.savefig('xasM2.png', format='png', dpi=300, bbox_inches='tight')
plt.show()
