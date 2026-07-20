
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

energy_Au_111_Cu_ohx2_atom, intensity_Au_111_Cu_ohx2_atom = norm_area("polAvg_Au_111_Cu_ohx2_atom.dat")
shift_Au_111_Cu_ohx2_atom = getShift('david_Au_111_Cu_ohx2_atom.csv')
print('Au_111_Cu_ohx2_atom total_shift: ' + str(shift_Au_111_Cu_ohx2_atom + d_exp))
energy_Au_111_Cu_ohx2_atom += shift_Au_111_Cu_ohx2_atom+ d_exp
x_Au_111_Cu_ohx2_atom,y_Au_111_Cu_ohx2_atom  = givePeak(energy=energy_Au_111_Cu_ohx2_atom, intensity=intensity_Au_111_Cu_ohx2_atom, name='Au_111_Cu_ohx2_atom')

energy_Au_111_Cu_oh_atom, intensity_Au_111_Cu_oh_atom = norm_area("polAvg_Au_111_Cu_oh_atom.dat")
shift_Au_111_Cu_oh_atom = getShift('david_Au_111_Cu_oh_atom.csv')
print('Au_111_Cu_oh_atom total_shift: ' + str(shift_Au_111_Cu_oh_atom + d_exp))
energy_Au_111_Cu_oh_atom += shift_Au_111_Cu_oh_atom+ d_exp
x_Au_111_Cu_oh_atom,y_Au_111_Cu_oh_atom  = givePeak(energy=energy_Au_111_Cu_oh_atom, intensity=intensity_Au_111_Cu_oh_atom, name='Au_111_Cu_oh_atom')

energy_Au_111_Cu2O_atom, intensity_Au_111_Cu2O_atom = norm_area("polAvg_Au_111_Cu2O_atom.dat")
shift_Au_111_Cu2O_atom = getShift('david_Au_111_Cu2O_atom.csv')
print('Au_111_Cu2O_atom total_shift: ' + str(shift_Au_111_Cu2O_atom + d_exp))
energy_Au_111_Cu2O_atom += shift_Au_111_Cu2O_atom+ d_exp
x_Au_111_Cu2O_atom,y_Au_111_Cu2O_atom  = givePeak(energy=energy_Au_111_Cu2O_atom, intensity=intensity_Au_111_Cu2O_atom, name='Au_111_Cu2O_atom')

energy_Au_111_CuO_atom, intensity_Au_111_CuO_atom = norm_area("polAvg_Au_111_CuO_atom.dat")
shift_Au_111_CuO_atom = getShift('david_Au_111_CuO_atom.csv')
print('Au_111_CuO_atom total_shift: ' + str(shift_Au_111_CuO_atom + d_exp))
energy_Au_111_CuO_atom += shift_Au_111_CuO_atom+ d_exp
x_Au_111_CuO_atom,y_Au_111_CuO_atom  = givePeak(energy=energy_Au_111_CuO_atom, intensity=intensity_Au_111_CuO_atom, name='Au_111_CuO_atom')

energy_Au_111_Cu_layer, intensity_Au_111_Cu_layer = norm_area("polAvg_Au_111_Cu_layer.dat")
shift_Au_111_Cu_layer = getShift('david_Au_111_Cu_layer.csv')
print('Au_111_Cu_layer total_shift: ' + str(shift_Au_111_Cu_layer + d_exp))
energy_Au_111_Cu_layer += shift_Au_111_Cu_layer+ d_exp
x_Au_111_Cu_layer,y_Au_111_Cu_layer  = givePeak(energy=energy_Au_111_Cu_layer, intensity=intensity_Au_111_Cu_layer, name='Au_111_Cu_layer')

energy_Au_111_Cu2O_Layer_1, intensity_Au_111_Cu2O_Layer_1 = norm_area("polAvg_Au_111_Cu2O_Layer_1.dat")
shift_Au_111_Cu2O_Layer_1 = getShift('david_Au_111_Cu2O_Layer_1.csv')
print('Au_111_Cu2O_Layer_1 total_shift: ' + str(shift_Au_111_Cu2O_Layer_1 + d_exp))
energy_Au_111_Cu2O_Layer_1 += shift_Au_111_Cu2O_Layer_1+ d_exp
x_Au_111_Cu2O_Layer_1,y_Au_111_Cu2O_Layer_1  = givePeak(energy=energy_Au_111_Cu2O_Layer_1, intensity=intensity_Au_111_Cu2O_Layer_1, name='Au_111_Cu2O_Layer_1')

energy_Au_111_Cu2O_Layer_0, intensity_Au_111_Cu2O_Layer_0 = norm_area("polAvg_Au_111_Cu2O_Layer_0.dat")
shift_Au_111_Cu2O_Layer_0 = getShift('david_Au_111_Cu2O_Layer_0.csv')
print('Au_111_Cu2O_Layer_0 total_shift: ' + str(shift_Au_111_Cu2O_Layer_0 + d_exp))
energy_Au_111_Cu2O_Layer_0 += shift_Au_111_Cu2O_Layer_0+ d_exp
x_Au_111_Cu2O_Layer_0,y_Au_111_Cu2O_Layer_0  = givePeak(energy=energy_Au_111_Cu2O_Layer_0, intensity=intensity_Au_111_Cu2O_Layer_0, name='Au_111_Cu2O_Layer_0')

energy_Au_111_CuO_Layer_0, intensity_Au_111_CuO_Layer_0 = norm_area("polAvg_Au_111_CuO_Layer_0.dat")
shift_Au_111_CuO_Layer_0 = getShift('david_Au_111_CuO_Layer_0.csv')
print('Au_111_CuO_Layer_0 total_shift: ' + str(shift_Au_111_CuO_Layer_0 + d_exp))
energy_Au_111_CuO_Layer_0 += shift_Au_111_CuO_Layer_0+ d_exp
x_Au_111_CuO_Layer_0,y_Au_111_CuO_Layer_0  = givePeak(energy=energy_Au_111_CuO_Layer_0, intensity=intensity_Au_111_CuO_Layer_0, name='Au_111_CuO_Layer_0')

energy_Au_111_CuO_Layer_1, intensity_Au_111_CuO_Layer_1 = norm_area("polAvg_Au_111_CuO_Layer_1.dat")
shift_Au_111_CuO_Layer_1 = getShift('david_Au_111_CuO_Layer_1.csv')
print('Au_111_CuO_Layer_1 total_shift: ' + str(shift_Au_111_CuO_Layer_1 + d_exp))
energy_Au_111_CuO_Layer_1 += shift_Au_111_CuO_Layer_1+ d_exp
x_Au_111_CuO_Layer_1,y_Au_111_CuO_Layer_1  = givePeak(energy=energy_Au_111_CuO_Layer_1, intensity=intensity_Au_111_CuO_Layer_1, name='Au_111_CuO_Layer_1')

energy_Au_111_Cu_Layer_0, intensity_Au_111_Cu_Layer_0 = norm_area("polAvg_Au_111_Cu_Layer_0.dat")
shift_Au_111_Cu_Layer_0 = getShift('david_Au_111_Cu_Layer_0.csv')
print('Au_111_Cu_Layer_0 total_shift: ' + str(shift_Au_111_Cu_Layer_0 + d_exp))
energy_Au_111_Cu_Layer_0 += shift_Au_111_Cu_Layer_0+ d_exp
x_Au_111_Cu_Layer_0,y_Au_111_Cu_Layer_0  = givePeak(energy=energy_Au_111_Cu_Layer_0, intensity=intensity_Au_111_Cu_Layer_0, name='Au_111_Cu_Layer_0')

energy_bulkCu2O, intensity_bulkCu2O = norm_area("polAvg_bulkCu2O.dat")
shift_bulkCu2O = getShift('david_bulkCu2O.csv')
print('bulkCu2O total_shift: ' + str(shift_bulkCu2O + d_exp))
energy_bulkCu2O += shift_bulkCu2O+ d_exp
x_bulkCu2O,y_bulkCu2O  = givePeak(energy=energy_bulkCu2O, intensity=intensity_bulkCu2O, name='bulkCu2O')

energy_bulkCuO, intensity_bulkCuO = norm_area("polAvg_bulkCuO.dat")
shift_bulkCuO = getShift('david_bulkCuO.csv')
print('bulkCuO total_shift: ' + str(shift_bulkCuO + d_exp))
energy_bulkCuO += shift_bulkCuO+ d_exp
x_bulkCuO,y_bulkCuO  = givePeak(energy=energy_bulkCuO, intensity=intensity_bulkCuO, name='bulkCuO')
# Plot

#plt.plot(energy_cuso4_solv, intensity_cuso4_solv+yoffset*0, label='cuso4_solv', linewidth=2, color='black', )
#if x_cuso4_solv is not None and y_cuso4_solv is not None:
#    plt.annotate(f'{x_cuso4_solv:.1f}', (x_cuso4_solv, y_cuso4_solv+yoffset*0),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 
#


#plt.plot(energy_Au_111_cuso4, intensity_Au_111_cuso4+yoffset*2, label='Au_111_cuso4', linewidth=2, color='red', )
#if x0_Au_111_cuso4 is not None and y0_Au_111_cuso4 is not None:
#    plt.annotate(f'{x0_Au_111_cuso4:.1f}', (x0_Au_111_cuso4, y0_Au_111_cuso4+yoffset*2),
#    textcoords="offset points", xytext=(0,2),
#    ha='center', fontproperties=arial, color='black')
#if x_Au_111_cuso4 is not None and y_Au_111_cuso4 is not None:
#    plt.annotate(f'{x_Au_111_cuso4:.1f}', (x_Au_111_cuso4, y_Au_111_cuso4+yoffset*2),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 

plt.figure(figsize=(2.5, 3.5))

plt.plot(energy_Au_111_Cu_ohx2_atom, intensity_Au_111_Cu_ohx2_atom+yoffset*1, label='Au_111_Cu_ohx2_atom', linewidth=2, color='purple', )
if x_Au_111_Cu_ohx2_atom is not None and y_Au_111_Cu_ohx2_atom is not None:
    plt.annotate(f'{x_Au_111_Cu_ohx2_atom:.1f}', (x_Au_111_Cu_ohx2_atom, y_Au_111_Cu_ohx2_atom+yoffset*1),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black') 

plt.plot(energy_Au_111_Cu_oh_atom, intensity_Au_111_Cu_oh_atom+yoffset*2, label='Au_111_Cu_oh_atom', linewidth=2, color='#515151',)
if x_Au_111_Cu_oh_atom is not None and y_Au_111_Cu_oh_atom is not None:
    plt.annotate(f'{x_Au_111_Cu_oh_atom:.1f}', (x_Au_111_Cu_oh_atom, y_Au_111_Cu_oh_atom+yoffset*2),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black') 

plt.plot(energy_Au_111_Cu2O_atom, intensity_Au_111_Cu2O_atom+yoffset*3, label='Au_111_Cu2O_atom', linewidth=2, color='#F14040', )
if x_Au_111_Cu2O_atom is not None and y_Au_111_Cu2O_atom is not None:
    plt.annotate(f'{x_Au_111_Cu2O_atom:.1f}', (x_Au_111_Cu2O_atom, y_Au_111_Cu2O_atom+yoffset*3),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black') 

plt.plot(energy_Au_111_CuO_atom, intensity_Au_111_CuO_atom+yoffset*4, label='Au_111_CuO_atom', linewidth=2, color='#1A6FDF', )
if x_Au_111_CuO_atom is not None and y_Au_111_CuO_atom is not None:
    plt.annotate(f'{x_Au_111_CuO_atom:.1f}', (x_Au_111_CuO_atom, y_Au_111_CuO_atom+yoffset*4),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black') 

plt.plot(energy_Au_111_Cu_atom, intensity_Au_111_Cu_atom+yoffset*5, label='Au_111_Cu_atom', linewidth=2, color='#37AD6B', )
if x_Au_111_Cu_atom is not None and y_Au_111_Cu_atom is not None:
    plt.annotate(f'{x_Au_111_Cu_atom:.1f}', (x_Au_111_Cu_atom, y_Au_111_Cu_atom+yoffset*5),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black') 


#plt.plot(energy_Au_111_Cu_layer, intensity_Au_111_Cu_layer+yoffset*6, label='Au_111_Cu_layer', linewidth=2, color='blue', linestyle='--')
#if x_Au_111_Cu_layer is not None and y_Au_111_Cu_layer is not None:
#    plt.annotate(f'{x_Au_111_Cu_layer:.1f}', (x_Au_111_Cu_layer, y_Au_111_Cu_layer+yoffset*6),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 
#
#plt.plot(energy_Au_111_Cu2O_Layer_1, intensity_Au_111_Cu2O_Layer_1+yoffset*7, label='Au_111_Cu2O_Layer_1', linewidth=2, color='darkgreen', )
#if x_Au_111_Cu2O_Layer_1 is not None and y_Au_111_Cu2O_Layer_1 is not None:
#    plt.annotate(f'{x_Au_111_Cu2O_Layer_1:.1f}', (x_Au_111_Cu2O_Layer_1, y_Au_111_Cu2O_Layer_1+yoffset*7),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 
#
#plt.plot(energy_Au_111_Cu2O_Layer_0, intensity_Au_111_Cu2O_Layer_0+yoffset*8, label='Au_111_Cu2O_Layer_0', linewidth=2, color='purple', )
#if x_Au_111_Cu2O_Layer_0 is not None and y_Au_111_Cu2O_Layer_0 is not None:
#    plt.annotate(f'{x_Au_111_Cu2O_Layer_0:.1f}', (x_Au_111_Cu2O_Layer_0, y_Au_111_Cu2O_Layer_0+yoffset*8),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 
#
#plt.plot(energy_Au_111_CuO_Layer_0, intensity_Au_111_CuO_Layer_0+yoffset*9, label='Au_111_CuO_Layer_0', linewidth=2, color='chocolate', )
#if x_Au_111_CuO_Layer_0 is not None and y_Au_111_CuO_Layer_0 is not None:
#    plt.annotate(f'{x_Au_111_CuO_Layer_0:.1f}', (x_Au_111_CuO_Layer_0, y_Au_111_CuO_Layer_0+yoffset*9),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 
#
#plt.plot(energy_Au_111_CuO_Layer_1, intensity_Au_111_CuO_Layer_1+yoffset*10, label='Au_111_CuO_Layer_1', linewidth=2, color='darkorange', )
#if x_Au_111_CuO_Layer_1 is not None and y_Au_111_CuO_Layer_1 is not None:
#    plt.annotate(f'{x_Au_111_CuO_Layer_1:.1f}', (x_Au_111_CuO_Layer_1, y_Au_111_CuO_Layer_1+yoffset*10),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 
#
#plt.plot(energy_Au_111_Cu_Layer_0, intensity_Au_111_Cu_Layer_0+yoffset*11, label='Au_111_Cu_Layer_0', linewidth=2, color='darkblue', )
#if x_Au_111_Cu_Layer_0 is not None and y_Au_111_Cu_Layer_0 is not None:
#    plt.annotate(f'{x_Au_111_Cu_Layer_0:.1f}', (x_Au_111_Cu_Layer_0, y_Au_111_Cu_Layer_0+yoffset*11),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 

#plt.plot(energy_bulkCu2O, intensity_bulkCu2O+yoffset*4, label='bulkCu2O', linewidth=2, color='darkred', )
#if x_bulkCu2O is not None and y_bulkCu2O is not None:
#    plt.annotate(f'{x_bulkCu2O:.1f}', (x_bulkCu2O, y_bulkCu2O+yoffset*4),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 
#
#plt.plot(energy_bulkCuO, intensity_bulkCuO+yoffset*5, label='bulkCuO', linewidth=2, color='darkviolet', linestyle='--')
#if x_bulkCuO is not None and y_bulkCuO is not None:
#    plt.annotate(f'{x_bulkCuO:.1f}', (x_bulkCuO, y_bulkCuO+yoffset*5),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 
## Plot reference 
#
#plt.plot(energy_bulkCu, intensity_bulkCu+yoffset*6, label='bulkCu', linewidth=2, color='dimgray', )
#if x_bulkCu is not None and y_bulkCu is not None:
#    plt.annotate(f'{x_bulkCu:.1f}', (x_bulkCu, y_bulkCu+yoffset*6),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 

# Plot settings
#plt.figure(figsize=(4, 3.5))
plt.xlabel('Energy (eV)', fontproperties=arialbd)
#tick_positions = np.arange(928, 945, 1)
#tick_labels = [str(tick) if tick % 2 == 0 else '' for tick in tick_positions]
#plt.xticks(tick_positions, tick_labels,fontproperties=arial)
plt.ylabel('Intensity (a.u.)', fontproperties=arialbd)
plt.gca().set_yticklabels([])  # Removes y-axis labels, keeps tick marks
plt.xticks(fontproperties=arialbd)
plt.xlim(926,944)
#plt.grid(True)
handles, labels = plt.gca().get_legend_handles_labels()
#plt.legend(handles[::-1], labels[::-1], bbox_to_anchor=(1, 1))
plt.tight_layout()

# Saving or showing plot
plt.savefig('figS7-moleculesXAS.svg', format='svg', bbox_inches='tight')
#plt.savefig('xasM2.png', format='png', dpi=300, bbox_inches='tight')
plt.show()
