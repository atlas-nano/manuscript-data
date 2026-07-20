
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


# Process data

energy_cuso4_solv_sulfate, intensity_cuso4_solv_sulfate = norm_area("polAvg_cuso4_solv_sulfate.dat")
shift_cuso4_solv_sulfate = getShift('david_cuso4_solv_sulfate.csv')
print('cuso4_solv_sulfate total_shift: ' + str(shift_cuso4_solv_sulfate + d_exp))
energy_cuso4_solv_sulfate += shift_cuso4_solv_sulfate+ d_exp
x_cuso4_solv_sulfate,y_cuso4_solv_sulfate  = givePeak(energy=energy_cuso4_solv_sulfate, intensity=intensity_cuso4_solv_sulfate, name='cuso4_solv_sulfate')

energy_Au_111_Cu_atom, intensity_Au_111_Cu_atom = norm_area("polAvg_Au_111_Cu_atom.dat")
shift_Au_111_Cu_atom = getShift('david_Au_111_Cu_atom.csv')
print('Au_111_Cu_atom total_shift: ' + str(shift_Au_111_Cu_atom + d_exp))
energy_Au_111_Cu_atom += shift_Au_111_Cu_atom+ d_exp
x_Au_111_Cu_atom,y_Au_111_Cu_atom  = givePeak(energy=energy_Au_111_Cu_atom, intensity=intensity_Au_111_Cu_atom, name='Au_111_Cu_atom')

energy_Au_111_cuso4_sulfate, intensity_Au_111_cuso4_sulfate = norm_area("polAvg_Au_111_cuso4_sulfate.dat")
shift_Au_111_cuso4_sulfate = getShift('david_Au_111_cuso4_sulfate.csv')
print('Au_111_cuso4_sulfate total_shift: ' + str(shift_Au_111_cuso4_sulfate + d_exp))
energy_Au_111_cuso4_sulfate += shift_Au_111_cuso4_sulfate+ d_exp
x_Au_111_cuso4_sulfate,y_Au_111_cuso4_sulfate  = givePeak(energy=energy_Au_111_cuso4_sulfate, intensity=intensity_Au_111_cuso4_sulfate, name='Au_111_cuso4_sulfate')

energy_Au_111_watbox, intensity_Au_111_watbox = norm_area("polAvg_Au_111_watbox.dat")
shift_Au_111_watbox = getShift('david_Au_111_watbox.csv')
print('Au_111_watbox total_shift: ' + str(shift_Au_111_watbox + d_exp))
energy_Au_111_watbox += shift_Au_111_watbox+ d_exp
x_Au_111_watbox,y_Au_111_watbox  = givePeak(energy=energy_Au_111_watbox, intensity=intensity_Au_111_watbox, name='Au_111_watbox')

energy_Au_111_Cu_ohx2_atom_hydroxide, intensity_Au_111_Cu_ohx2_atom_hydroxide = norm_area("polAvg_Au_111_Cu_ohx2_atom_hydroxide.dat")
shift_Au_111_Cu_ohx2_atom_hydroxide = getShift('david_Au_111_Cu_ohx2_atom_hydroxide.csv')
print('Au_111_Cu_ohx2_atom_hydroxide total_shift: ' + str(shift_Au_111_Cu_ohx2_atom_hydroxide + d_exp))
energy_Au_111_Cu_ohx2_atom_hydroxide += shift_Au_111_Cu_ohx2_atom_hydroxide+ d_exp
x_Au_111_Cu_ohx2_atom_hydroxide,y_Au_111_Cu_ohx2_atom_hydroxide  = givePeak(energy=energy_Au_111_Cu_ohx2_atom_hydroxide, intensity=intensity_Au_111_Cu_ohx2_atom_hydroxide, name='Au_111_Cu_ohx2_atom_hydroxide')

energy_Au_111_Cu_oh_atom_hydroxide, intensity_Au_111_Cu_oh_atom_hydroxide = norm_area("polAvg_Au_111_Cu_oh_atom_hydroxide.dat")
shift_Au_111_Cu_oh_atom_hydroxide = getShift('david_Au_111_Cu_oh_atom_hydroxide.csv')
print('Au_111_Cu_oh_atom_hydroxide total_shift: ' + str(shift_Au_111_Cu_oh_atom_hydroxide + d_exp))
energy_Au_111_Cu_oh_atom_hydroxide += shift_Au_111_Cu_oh_atom_hydroxide+ d_exp
x_Au_111_Cu_oh_atom_hydroxide,y_Au_111_Cu_oh_atom_hydroxide  = givePeak(energy=energy_Au_111_Cu_oh_atom_hydroxide, intensity=intensity_Au_111_Cu_oh_atom_hydroxide, name='Au_111_Cu_oh_atom_hydroxide')

energy_Au_111_Cu2O_atom_oxide, intensity_Au_111_Cu2O_atom_oxide = norm_area("polAvg_Au_111_Cu2O_atom_oxide.dat")
shift_Au_111_Cu2O_atom_oxide = getShift('david_Au_111_Cu2O_atom_oxide.csv')
print('Au_111_Cu2O_atom_oxide total_shift: ' + str(shift_Au_111_Cu2O_atom_oxide + d_exp))
energy_Au_111_Cu2O_atom_oxide += shift_Au_111_Cu2O_atom_oxide+ d_exp
x_Au_111_Cu2O_atom_oxide,y_Au_111_Cu2O_atom_oxide  = givePeak(energy=energy_Au_111_Cu2O_atom_oxide, intensity=intensity_Au_111_Cu2O_atom_oxide, name='Au_111_Cu2O_atom_oxide')

energy_Au_111_CuO_atom_oxide, intensity_Au_111_CuO_atom_oxide = norm_area("polAvg_Au_111_CuO_atom_oxide.dat")
shift_Au_111_CuO_atom_oxide = getShift('david_Au_111_CuO_atom_oxide.csv')
print('Au_111_CuO_atom_oxide total_shift: ' + str(shift_Au_111_CuO_atom_oxide + d_exp))
energy_Au_111_CuO_atom_oxide += shift_Au_111_CuO_atom_oxide+ d_exp
x_Au_111_CuO_atom_oxide,y_Au_111_CuO_atom_oxide  = givePeak(energy=energy_Au_111_CuO_atom_oxide, intensity=intensity_Au_111_CuO_atom_oxide, name='Au_111_CuO_atom_oxide')

energy_Au_111_Cu_layer, intensity_Au_111_Cu_layer = norm_area("polAvg_Au_111_Cu_layer.dat")
shift_Au_111_Cu_layer = getShift('david_Au_111_Cu_layer.csv')
print('Au_111_Cu_layer total_shift: ' + str(shift_Au_111_Cu_layer + d_exp))
energy_Au_111_Cu_layer += shift_Au_111_Cu_layer+ d_exp
x_Au_111_Cu_layer,y_Au_111_Cu_layer  = givePeak(energy=energy_Au_111_Cu_layer, intensity=intensity_Au_111_Cu_layer, name='Au_111_Cu_layer')

energy_Au_111_Cu2O_Layer_1_oxide, intensity_Au_111_Cu2O_Layer_1_oxide = norm_area("polAvg_Au_111_Cu2O_Layer_1_oxide.dat")
shift_Au_111_Cu2O_Layer_1_oxide = getShift('david_Au_111_Cu2O_Layer_1_oxide.csv')
print('Au_111_Cu2O_Layer_1_oxide total_shift: ' + str(shift_Au_111_Cu2O_Layer_1_oxide + d_exp))
energy_Au_111_Cu2O_Layer_1_oxide += shift_Au_111_Cu2O_Layer_1_oxide+ d_exp
x_Au_111_Cu2O_Layer_1_oxide,y_Au_111_Cu2O_Layer_1_oxide  = givePeak(energy=energy_Au_111_Cu2O_Layer_1_oxide, intensity=intensity_Au_111_Cu2O_Layer_1_oxide, name='Au_111_Cu2O_Layer_1_oxide')

energy_Au_111_Cu2O_Layer_0_oxide, intensity_Au_111_Cu2O_Layer_0_oxide = norm_area("polAvg_Au_111_Cu2O_Layer_0_oxide.dat")
shift_Au_111_Cu2O_Layer_0_oxide = getShift('david_Au_111_Cu2O_Layer_0_oxide.csv')
print('Au_111_Cu2O_Layer_0_oxide total_shift: ' + str(shift_Au_111_Cu2O_Layer_0_oxide + d_exp))
energy_Au_111_Cu2O_Layer_0_oxide += shift_Au_111_Cu2O_Layer_0_oxide+ d_exp
x_Au_111_Cu2O_Layer_0_oxide,y_Au_111_Cu2O_Layer_0_oxide  = givePeak(energy=energy_Au_111_Cu2O_Layer_0_oxide, intensity=intensity_Au_111_Cu2O_Layer_0_oxide, name='Au_111_Cu2O_Layer_0_oxide')

energy_Au_111_CuO_Layer_0_oxide, intensity_Au_111_CuO_Layer_0_oxide = norm_area("polAvg_Au_111_CuO_Layer_0_oxide.dat")
shift_Au_111_CuO_Layer_0_oxide = getShift('david_Au_111_CuO_Layer_0_oxide.csv')
print('Au_111_CuO_Layer_0_oxide total_shift: ' + str(shift_Au_111_CuO_Layer_0_oxide + d_exp))
energy_Au_111_CuO_Layer_0_oxide += shift_Au_111_CuO_Layer_0_oxide+ d_exp
x_Au_111_CuO_Layer_0_oxide,y_Au_111_CuO_Layer_0_oxide  = givePeak(energy=energy_Au_111_CuO_Layer_0_oxide, intensity=intensity_Au_111_CuO_Layer_0_oxide, name='Au_111_CuO_Layer_0_oxide')

energy_Au_111_CuO_Layer_1_oxide, intensity_Au_111_CuO_Layer_1_oxide = norm_area("polAvg_Au_111_CuO_Layer_1_oxide.dat")
shift_Au_111_CuO_Layer_1_oxide = getShift('david_Au_111_CuO_Layer_1_oxide.csv')
print('Au_111_CuO_Layer_1_oxide total_shift: ' + str(shift_Au_111_CuO_Layer_1_oxide + d_exp))
energy_Au_111_CuO_Layer_1_oxide += shift_Au_111_CuO_Layer_1_oxide+ d_exp
x_Au_111_CuO_Layer_1_oxide,y_Au_111_CuO_Layer_1_oxide  = givePeak(energy=energy_Au_111_CuO_Layer_1_oxide, intensity=intensity_Au_111_CuO_Layer_1_oxide, name='Au_111_CuO_Layer_1_oxide')

energy_watbox, intensity_watbox = norm_area("polAvg_watbox.dat")
shift_watbox = getShift('david_watbox.csv')
print('watbox total_shift: ' + str(shift_watbox + d_exp))
energy_watbox += shift_watbox+ d_exp
x_watbox,y_watbox  = givePeak(energy=energy_watbox, intensity=intensity_watbox, name='watbox')

energy_bulkCu2O, intensity_bulkCu2O = norm_area("polAvg_bulkCu2O.dat")
shift_bulkCu2O = getShift('david_bulkCu2O.csv')
print('bulkCu2O total_shift: ' + str(shift_bulkCu2O + d_exp))
energy_bulkCu2O += shift_bulkCu2O+ d_exp
x_bulkCu2O,y_bulkCu2O  = givePeak(energy=energy_bulkCu2O, intensity=intensity_bulkCu2O, name='bulkCu2O')
# Plot

#plt.plot(energy_cuso4_solv_sulfate, intensity_cuso4_solv_sulfate+yoffset*0, label='cuso4_solv_sulfate', linewidth=2, color='black', )
#if x_cuso4_solv_sulfate is not None and y_cuso4_solv_sulfate is not None:
#    plt.annotate(f'{x_cuso4_solv_sulfate:.1f}', (x_cuso4_solv_sulfate, y_cuso4_solv_sulfate+yoffset*0),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 


#plt.plot(energy_Au_111_cuso4_sulfate, intensity_Au_111_cuso4_sulfate+yoffset*2, label='Au_111_cuso4_sulfate', linewidth=2, color='red', )
#if x_Au_111_cuso4_sulfate is not None and y_Au_111_cuso4_sulfate is not None:
#    plt.annotate(f'{x_Au_111_cuso4_sulfate:.1f}', (x_Au_111_cuso4_sulfate, y_Au_111_cuso4_sulfate+yoffset*2),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 
#
#plt.plot(energy_Au_111_watbox, intensity_Au_111_watbox+yoffset*3, label='Au_111_watbox', linewidth=2, color='deepskyblue', )
#if x_Au_111_watbox is not None and y_Au_111_watbox is not None:
#    plt.annotate(f'{x_Au_111_watbox:.1f}', (x_Au_111_watbox, y_Au_111_watbox+yoffset*3),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 



#plt.plot(energy_Au_111_Cu_layer, intensity_Au_111_Cu_layer+yoffset*7, label='Au_111_Cu_layer', linewidth=2, color='darkgreen', )
#if x_Au_111_Cu_layer is not None and y_Au_111_Cu_layer is not None:
#    plt.annotate(f'{x_Au_111_Cu_layer:.1f}', (x_Au_111_Cu_layer, y_Au_111_Cu_layer+yoffset*7),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 
#
#plt.plot(energy_Au_111_Cu2O_Layer_1_oxide, intensity_Au_111_Cu2O_Layer_1_oxide+yoffset*8, label='Au_111_Cu2O_Layer_1_oxide', linewidth=2, color='purple', )
#if x_Au_111_Cu2O_Layer_1_oxide is not None and y_Au_111_Cu2O_Layer_1_oxide is not None:
#    plt.annotate(f'{x_Au_111_Cu2O_Layer_1_oxide:.1f}', (x_Au_111_Cu2O_Layer_1_oxide, y_Au_111_Cu2O_Layer_1_oxide+yoffset*8),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 
#
#plt.plot(energy_Au_111_Cu2O_Layer_0_oxide, intensity_Au_111_Cu2O_Layer_0_oxide+yoffset*9, label='Au_111_Cu2O_Layer_0_oxide', linewidth=2, color='darkorange', )
#if x_Au_111_Cu2O_Layer_0_oxide is not None and y_Au_111_Cu2O_Layer_0_oxide is not None:
#    plt.annotate(f'{x_Au_111_Cu2O_Layer_0_oxide:.1f}', (x_Au_111_Cu2O_Layer_0_oxide, y_Au_111_Cu2O_Layer_0_oxide+yoffset*9),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 
#
#plt.plot(energy_Au_111_CuO_Layer_0_oxide, intensity_Au_111_CuO_Layer_0_oxide+yoffset*10, label='Au_111_CuO_Layer_0_oxide', linewidth=2, color='slategray', )
#if x_Au_111_CuO_Layer_0_oxide is not None and y_Au_111_CuO_Layer_0_oxide is not None:
#    plt.annotate(f'{x_Au_111_CuO_Layer_0_oxide:.1f}', (x_Au_111_CuO_Layer_0_oxide, y_Au_111_CuO_Layer_0_oxide+yoffset*10),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 
#
#plt.plot(energy_Au_111_CuO_Layer_1_oxide, intensity_Au_111_CuO_Layer_1_oxide+yoffset*11, label='Au_111_CuO_Layer_1_oxide', linewidth=2, color='darkblue', )
#if x_Au_111_CuO_Layer_1_oxide is not None and y_Au_111_CuO_Layer_1_oxide is not None:
#    plt.annotate(f'{x_Au_111_CuO_Layer_1_oxide:.1f}', (x_Au_111_CuO_Layer_1_oxide, y_Au_111_CuO_Layer_1_oxide+yoffset*11),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 

#plt.plot(energy_watbox, intensity_watbox+yoffset*5, label='watbox', linewidth=2, color='darkred', )
#if x_watbox is not None and y_watbox is not None:
#    plt.annotate(f'{x_watbox:.1f}', (x_watbox, y_watbox+yoffset*5),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 

#plt.plot(energy_bulkCu2O, intensity_bulkCu2O+yoffset*6, label='bulkCu2O', linewidth=2, color='darkviolet')
#if x_bulkCu2O is not None and y_bulkCu2O is not None:
#    plt.annotate(f'{x_bulkCu2O:.1f}', (x_bulkCu2O, y_bulkCu2O+yoffset*6),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 
## Plot reference 
#
#plt.plot(energy_bulkCuO, intensity_bulkCuO+yoffset*7, label='bulkCuO', linewidth=2, color='midnightblue', )
#if x_bulkCuO is not None and y_bulkCuO is not None:
#    plt.annotate(f'{x_bulkCuO:.1f}', (x_bulkCuO, y_bulkCuO+yoffset*7),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 

plt.figure(figsize=(2.5, 3.5))

plt.plot(energy_Au_111_Cu_ohx2_atom_hydroxide, intensity_Au_111_Cu_ohx2_atom_hydroxide+yoffset*2, label='Au_111_Cu_ohx2_atom_hydroxide', linewidth=2, color='purple')
if x_Au_111_Cu_ohx2_atom_hydroxide is not None and y_Au_111_Cu_ohx2_atom_hydroxide is not None:
    plt.annotate(f'{x_Au_111_Cu_ohx2_atom_hydroxide:.1f}', (x_Au_111_Cu_ohx2_atom_hydroxide, y_Au_111_Cu_ohx2_atom_hydroxide+yoffset*2),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black') 

plt.plot(energy_Au_111_Cu_oh_atom_hydroxide, intensity_Au_111_Cu_oh_atom_hydroxide+yoffset*3, label='Au_111_Cu_oh_atom_hydroxide', linewidth=2, color='#515151', )
if x_Au_111_Cu_oh_atom_hydroxide is not None and y_Au_111_Cu_oh_atom_hydroxide is not None:
    plt.annotate(f'{x_Au_111_Cu_oh_atom_hydroxide:.1f}', (x_Au_111_Cu_oh_atom_hydroxide, y_Au_111_Cu_oh_atom_hydroxide+yoffset*3),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black') 

plt.plot(energy_Au_111_Cu2O_atom_oxide, intensity_Au_111_Cu2O_atom_oxide+yoffset*4, label='Au_111_Cu2O_atom_oxide', linewidth=2, color='#F14040')
if x_Au_111_Cu2O_atom_oxide is not None and y_Au_111_Cu2O_atom_oxide is not None:
    plt.annotate(f'{x_Au_111_Cu2O_atom_oxide:.1f}', (x_Au_111_Cu2O_atom_oxide, y_Au_111_Cu2O_atom_oxide+yoffset*4),
    textcoords="offset points", xytext=(0,2),
    ha='center', fontproperties=arial, color='black') 

plt.plot(energy_Au_111_CuO_atom_oxide, intensity_Au_111_CuO_atom_oxide+yoffset*5, label='Au_111_CuO_atom_oxide', linewidth=2, color='#1A6FDF')
if x_Au_111_CuO_atom_oxide is not None and y_Au_111_CuO_atom_oxide is not None:
    plt.annotate(f'{x_Au_111_CuO_atom_oxide:.1f}', (x_Au_111_CuO_atom_oxide, y_Au_111_CuO_atom_oxide+yoffset*5),
    textcoords="offset points", xytext=(0,2),
    ha='center', fontproperties=arial, color='black') 

#plt.plot(energy_Au_111_Cu_atom, intensity_Au_111_Cu_atom+yoffset*6, label='Au_111_Cu_atom', linewidth=2, color='#37AD6B', )
#if x_Au_111_Cu_atom is not None and y_Au_111_Cu_atom is not None:
#    plt.annotate(f'{x_Au_111_Cu_atom:.1f}', (x_Au_111_Cu_atom, y_Au_111_Cu_atom+yoffset*6),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 

plt.plot(energy_Au_111_watbox, intensity_Au_111_watbox+yoffset*6, label='Au_111_watbox', linewidth=2, color='#37AD6B', )
if x_Au_111_watbox is not None and y_Au_111_watbox is not None:
    plt.annotate(f'{x_Au_111_watbox:.1f}', (x_Au_111_watbox, y_Au_111_watbox+yoffset*6),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black') 

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
#plt.grid(True)
handles, labels = plt.gca().get_legend_handles_labels()
#plt.legend(handles[::-1], labels[::-1], bbox_to_anchor=(1, 1))
plt.tight_layout()

# Saving or showing plot
plt.savefig('o-k-moleculesXAS.svg', format='svg', bbox_inches='tight')
#plt.savefig('xasM2.png', format='png', dpi=300, bbox_inches='tight')
plt.show()
