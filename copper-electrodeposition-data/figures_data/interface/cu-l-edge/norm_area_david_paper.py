
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

def givePeak2(energy, intensity, prominence, index):
    # Find peaks
    peaks, _ = scipy.signal.find_peaks(intensity, prominence)
    first_peak_energy = energy[peaks[index]]
    first_peak_intensity = intensity[peaks[index]]
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
x_minus02V_exp = energy_minus02V_exp[peaks_minus02V_exp[0]]
y_minus02V_exp = intensity_minus02V_exp[peaks_minus02V_exp[0]]
#x_minus02V_exp, y_minus02V_exp = givePeak(energy=energy_minus02V_exp, intensity=intensity_minus02V_exp, name='minus02V')

energy_minus03V_exp, intensity_minus03V_exp = exp("exp_-0.3V.dat")
x_minus03V_exp, y_minus03V_exp = givePeak(energy=energy_minus03V_exp, intensity=intensity_minus03V_exp, name='minus03V')

energy_0V_exp, intensity_0V_exp = exp("exp_0V.dat")
x_0V_exp, y_0V_exp = givePeak(energy=energy_0V_exp, intensity=intensity_0V_exp, name='0V')


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

energy_Au_111_CuO_Layer_0_solv_mag, intensity_Au_111_CuO_Layer_0_solv_mag = norm_area("polAvg_Au_111_CuO_Layer_0_solv_mag.dat")
shift_Au_111_CuO_Layer_0_solv_mag = getShift('david_Au_111_CuO_Layer_0_solv_mag_cu1.csv')
print('Au_111_CuO_Layer_0_solv_mag total_shift: ' + str(shift_Au_111_CuO_Layer_0_solv_mag + d_exp))
energy_Au_111_CuO_Layer_0_solv_mag += shift_Au_111_CuO_Layer_0_solv_mag+ d_exp
x_Au_111_CuO_Layer_0_solv_mag,y_Au_111_CuO_Layer_0_solv_mag  = givePeak(energy=energy_Au_111_CuO_Layer_0_solv_mag, intensity=intensity_Au_111_CuO_Layer_0_solv_mag, name='Au_111_CuO_Layer_0_solv_mag')

energy_Au_111_CuO_Layer_0_solv_mag_u8, intensity_Au_111_CuO_Layer_0_solv_mag_u8 = norm_area("polAvg_Au_111_CuO_Layer_0_solv_mag_u8.dat")
shift_Au_111_CuO_Layer_0_solv_mag_u8 = getShift('david_Au_111_CuO_Layer_0_solv_mag_u8.csv')
print('Au_111_CuO_Layer_0_solv_mag_u8 total_shift: ' + str(shift_Au_111_CuO_Layer_0_solv_mag_u8 + d_exp))
energy_Au_111_CuO_Layer_0_solv_mag_u8 += shift_Au_111_CuO_Layer_0_solv_mag_u8+ d_exp
x_Au_111_CuO_Layer_0_solv_mag_u8,y_Au_111_CuO_Layer_0_solv_mag_u8  = givePeak(energy=energy_Au_111_CuO_Layer_0_solv_mag_u8, intensity=intensity_Au_111_CuO_Layer_0_solv_mag_u8, name='Au_111_CuO_Layer_0_solv_mag_u8')

energy_Au_111_CuO_Layer_0_solv_allFrames, intensity_Au_111_CuO_Layer_0_solv_allFrames = norm_area("polAvg_Au_111_CuO_Layer_0_solv_allFrames.dat")
shift_Au_111_CuO_Layer_0_solv_allFrames = getShift('david_Au_111_CuO_Layer_0_solv.csv')
print('Au_111_CuO_Layer_0_solv_allFrames total_shift: ' + str(shift_Au_111_CuO_Layer_0_solv_allFrames + d_exp))
energy_Au_111_CuO_Layer_0_solv_allFrames += shift_Au_111_CuO_Layer_0_solv_allFrames+ d_exp
x_Au_111_CuO_Layer_0_solv_allFrames,y_Au_111_CuO_Layer_0_solv_allFrames  = givePeak(energy=energy_Au_111_CuO_Layer_0_solv_allFrames, intensity=intensity_Au_111_CuO_Layer_0_solv_allFrames, name='Au_111_CuO_Layer_0_solv_allFrames')

energy_Au_111_Cu2O_Layer_1_solv_mag, intensity_Au_111_Cu2O_Layer_1_solv_mag = norm_area("polAvg_Au_111_Cu2O_Layer_1_solv_mag.dat")
shift_Au_111_Cu2O_Layer_1_solv_mag = getShift('david_Au_111_Cu2O_Layer_1_solv_mag.csv')
print('Au_111_Cu2O_Layer_1_solv_mag total_shift: ' + str(shift_Au_111_Cu2O_Layer_1_solv_mag + d_exp))
energy_Au_111_Cu2O_Layer_1_solv_mag += shift_Au_111_Cu2O_Layer_1_solv_mag+ d_exp
x_Au_111_Cu2O_Layer_1_solv_mag_0,y_Au_111_Cu2O_Layer_1_solv_mag_0 = givePeak2(energy=energy_Au_111_Cu2O_Layer_1_solv_mag, intensity=intensity_Au_111_Cu2O_Layer_1_solv_mag, prominence=0.01, index=0)
x_Au_111_Cu2O_Layer_1_solv_mag_1,y_Au_111_Cu2O_Layer_1_solv_mag_1 = givePeak2(energy=energy_Au_111_Cu2O_Layer_1_solv_mag, intensity=intensity_Au_111_Cu2O_Layer_1_solv_mag, prominence=0.03, index=4)


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
plt.figure(figsize=(2.5, 4.0))

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

plt.plot(energy_0V_exp, intensity_0V_exp/exp_scaling + yoffset*-1, label='0 V', linewidth=2, color='#515151')
#if x_0V_exp is not None and y_0V_exp is not None:
#    plt.annotate(f'{x_0V_exp:.1f}', (x_0V_exp, y_0V_exp/exp_scaling + yoffset*0),
#                 textcoords="offset points", xytext=(0,2),
#                 ha='center', fontproperties=arial, color='black')
plt.plot(energy_Au_111_cu_ads_so4_diemac_1_80, intensity_Au_111_cu_ads_so4_diemac_1_80+yoffset*-1, label='Au_111_cu_ads_so4_diemac_1_80 z = 0 A', linewidth=2, color='#515151', alpha=0.5)
if x0_Au_111_cu_ads_so4_diemac_1_80 is not None and y0_Au_111_cu_ads_so4_diemac_1_80 is not None:
    plt.annotate(f'{x0_Au_111_cu_ads_so4_diemac_1_80:.1f}', (x0_Au_111_cu_ads_so4_diemac_1_80, y0_Au_111_cu_ads_so4_diemac_1_80+yoffset*-1),
    textcoords="offset points", xytext=(0,2),
    ha='center', fontproperties=arial, color='black')
if x_Au_111_cu_ads_so4_diemac_1_80 is not None and y_Au_111_cu_ads_so4_diemac_1_80 is not None:
    plt.annotate(f'{x_Au_111_cu_ads_so4_diemac_1_80:.1f}', (x_Au_111_cu_ads_so4_diemac_1_80, y_Au_111_cu_ads_so4_diemac_1_80+yoffset*-1),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black') 

plt.plot(energy_minus01V_exp, intensity_minus01V_exp/exp_scaling+yoffset*0, label='- 0.1 V', linewidth=2, color='#F14040')
#if x_minus01V_exp is not None and y_minus01V_exp is not None:
#    plt.annotate(f'{x_minus01V_exp:.1f}', (x_minus01V_exp, y_minus01V_exp/exp_scaling+yoffset*1),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 
plt.plot(energy_Au_111_cu_ads_so4_diemac_1_80, intensity_Au_111_cu_ads_so4_diemac_1_80+yoffset*0, label='Au_111_cu_ads_so4_diemac_1_80 z = 0 A', linewidth=2, color='#F14040',  alpha=0.5)
if x0_Au_111_cu_ads_so4_diemac_1_80 is not None and y0_Au_111_cu_ads_so4_diemac_1_80 is not None:
    plt.annotate(f'{x0_Au_111_cu_ads_so4_diemac_1_80:.1f}', (x0_Au_111_cu_ads_so4_diemac_1_80, y0_Au_111_cu_ads_so4_diemac_1_80+yoffset*0),
    textcoords="offset points", xytext=(0,2),
    ha='center', fontproperties=arial, color='black')
if x_Au_111_cu_ads_so4_diemac_1_80 is not None and y_Au_111_cu_ads_so4_diemac_1_80 is not None:
    plt.annotate(f'{x_Au_111_cu_ads_so4_diemac_1_80:.1f}', (x_Au_111_cu_ads_so4_diemac_1_80, y_Au_111_cu_ads_so4_diemac_1_80+yoffset*0),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black') 
    
plt.plot(energy_minus02V_exp, intensity_minus02V_exp/exp_scaling + yoffset*1.2, label='-0.2 V', color='#1A6FDF', linewidth=2)
#if x_minus02V_exp is not None and y_minus02V_exp is not None:
#    plt.annotate(f'{x_minus02V_exp:.1f}', (x_minus02V_exp, y_minus02V_exp/exp_scaling + yoffset*2),
#                 textcoords="offset points", xytext=(0,2),
#                 ha='center', fontproperties=arial, color='black')
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
#plt.plot(energy_Au_111_CuO_Layer_0_solv_allFrames, intensity_Au_111_CuO_Layer_0_solv_allFrames+yoffset*1.2, label='Au_111_CuO_Layer_0_solv_allFrames', linewidth=2, color='#1A6FDF',  linestyle=":")
#if x_Au_111_CuO_Layer_0_solv_allFrames is not None and y_Au_111_CuO_Layer_0_solv_allFrames is not None:
#    plt.annotate(f'{x_Au_111_CuO_Layer_0_solv_allFrames:.1f}', (x_Au_111_CuO_Layer_0_solv_allFrames, y_Au_111_CuO_Layer_0_solv_allFrames+yoffset*1.2),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 
#plt.plot(energy_Au_111_CuO_Layer_0_solv_mag, intensity_Au_111_CuO_Layer_0_solv_mag+yoffset*1.2, label='Au_111_CuO_Layer_0_solv_mag', linewidth=2, color='black',  linestyle=":")
#if x_Au_111_CuO_Layer_0_solv_mag is not None and y_Au_111_CuO_Layer_0_solv_mag is not None:
#    plt.annotate(f'{x_Au_111_CuO_Layer_0_solv_mag:.1f}', (x_Au_111_CuO_Layer_0_solv_mag, y_Au_111_CuO_Layer_0_solv_mag+yoffset*1.2),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black')
#plt.plot(energy_Au_111_CuO_Layer_0_solv_mag_u8, intensity_Au_111_CuO_Layer_0_solv_mag_u8+yoffset*1.2, label='Au_111_CuO_Layer_0_solv_mag_u8', linewidth=2, color='black',  linestyle=":")
#if x_Au_111_CuO_Layer_0_solv_mag_u8 is not None and y_Au_111_CuO_Layer_0_solv_mag_u8 is not None:
#    plt.annotate(f'{x_Au_111_CuO_Layer_0_solv_mag_u8:.1f}', (x_Au_111_CuO_Layer_0_solv_mag_u8, y_Au_111_CuO_Layer_0_solv_mag_u8+yoffset*1.2),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black')
plt.plot(energy_bulkCu2O_exp, intensity_bulkCu2O_exp/ 16.0+yoffset*1.4, label='bulkCu2O exp', linewidth=2, color='#1A6FDF', linestyle=':')
#if x_bulkCu2O_exp_0 is not None and y_bulkCu2O_exp_0 is not None:
#    plt.annotate(f'{x_bulkCu2O_exp_0:.1f}', (x_bulkCu2O_exp_0, y_bulkCu2O_exp_0/7.0 + yoffset*1.0),
#                 textcoords="offset points", xytext=(0,2),
#                 ha='center', fontproperties=arial, color='black')
#if x_bulkCu2O_exp_1 is not None and y_bulkCu2O_exp_1 is not None:
#    plt.annotate(f'{x_bulkCu2O_exp_1:.1f}', (x_bulkCu2O_exp_1, y_bulkCu2O_exp_1/7.0 + yoffset*1.0),
#                 textcoords="offset points", xytext=(0,2),
#                 ha='center', fontproperties=arial, color='black')

plt.plot(energy_Au_111_Cu2O_Layer_1_solv_mag, intensity_Au_111_Cu2O_Layer_1_solv_mag+yoffset*1.6, label='Au_111_Cu2O_Layer_1_solv_mag', linewidth=2, color='#1A6FDF',  alpha=0.5)
if x_Au_111_Cu2O_Layer_1_solv_mag_0 is not None and y_Au_111_Cu2O_Layer_1_solv_mag_0 is not None:
    plt.annotate(f'{x_Au_111_Cu2O_Layer_1_solv_mag_0:.1f}', (x_Au_111_Cu2O_Layer_1_solv_mag_0, y_Au_111_Cu2O_Layer_1_solv_mag_0 +yoffset*1.6),
                 textcoords="offset points", xytext=(0,2),
                 ha='center', fontproperties=arial, color='black')
if x_Au_111_Cu2O_Layer_1_solv_mag_1 is not None and y_Au_111_Cu2O_Layer_1_solv_mag_1 is not None:
    plt.annotate(f'{x_Au_111_Cu2O_Layer_1_solv_mag_1:.1f}', (x_Au_111_Cu2O_Layer_1_solv_mag_1, y_Au_111_Cu2O_Layer_1_solv_mag_1 +yoffset*1.6),
                 textcoords="offset points", xytext=(0,2),
                 ha='center', fontproperties=arial, color='black')



plt.plot(energy_minus03V_exp, intensity_minus03V_exp/exp_scaling + yoffset*4, label='-0.3 V', linewidth=2, color='#37AD6B')
#if x_minus03V_exp is not None and y_minus03V_exp is not None:
#    plt.annotate(f'{x_minus03V_exp:.1f}', (x_minus03V_exp, y_minus03V_exp + yoffset*3),
#                 textcoords="offset points", xytext=(0,2),
#                 ha='center', fontproperties=arial, color='black')
plt.plot(energy_Au_111_Cu_layer, intensity_Au_111_Cu_layer+yoffset*4, label='Au_111_Cu_layer', linewidth=2, color='#37AD6B',  alpha=0.5)
if x_Au_111_Cu_layer is not None and y_Au_111_Cu_layer is not None:
    plt.annotate(f'{x_Au_111_Cu_layer:.1f}', (x_Au_111_Cu_layer, y_Au_111_Cu_layer+yoffset*4),
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
plt.gca().set_yticklabels([])  # Removes y-axis labels, keeps tick marks
#plt.yticks(fontproperties=arial)
plt.xticks(fontproperties=arialbd)
plt.xlim(926,960)
#plt.grid(True)
#handles, labels = plt.gca().get_legend_handles_labels()
#plt.legend(handles[::-1], labels[::-1],loc='upper left', bbox_to_anchor=(1, 1))
plt.tight_layout()

# Saving or showing plot
plt.savefig('main_cu.svg', format='svg', bbox_inches='tight')
#plt.savefig('xasM2.png', format='png', dpi=300, bbox_inches='tight')
plt.show()
