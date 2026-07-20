
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
import scipy.signal
import textwrap
from matplotlib import font_manager
import matplotlib
matplotlib.use("TkAgg")

arial = font_manager.FontProperties(fname="/mnt/c/Windows/Fonts/arial.ttf", size=9)
arialbd = font_manager.FontProperties(fname="/mnt/c/Windows/Fonts/arialbd.ttf", size=9)

#E_f = 17 #eV
E_i = -25 #eV
E_f = 25 #eV

E_i_exp = 525 #eV
E_f_exp = 545 #eV From Juan

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
    :::::::::::::::::::::::::::::::::::::
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

    # Check if energy is in decreasing order
    if energy[0] > energy[-1]:
        energy = energy[::-1]
        intensity = intensity[::-1]
    mask = (energy >= E_i_exp) & (energy <= E_f_exp)
    energy_range = energy[mask]
    intensity_range = intensity[mask]
    auc = simpson(intensity_range, energy_range)
    intensity = intensity / auc
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
    peaks, _ = scipy.signal.find_peaks(intensity, prominence=0.05)

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
x_bulkCuO_exp,y_bulkCuO_exp  = givePeak(energy=energy_bulkCuO_exp, intensity=intensity_bulkCuO_exp, name='bulkCuO')
print(f'exp bulkCuO: {x_bulkCuO_exp}')
## Computational reference


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


#energy_bulkCuO_mag_u8, intensity_bulkCuO_mag_u8 = norm_area("polAvg_bulkCuO_mag_u8.dat")
#x_bulkCuO_mag_u8,y_bulkCuO_mag_u8  = givePeak(energy=energy_bulkCuO_mag_u8, intensity=intensity_bulkCuO_mag_u8, name='bulkCuO_mag_u8')
#print(f'comp bulkCuO_mag_u8: {x_bulkCuO_mag_u8}')
#shift_bulkCuO_mag_u8 = getShift('david_bulkCuO_mag_u8.csv')
## ::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Experimental shift
#d_exp = x_bulkCuO_exp - (x_bulkCuO_mag_u8+shift_bulkCuO_mag_u8)
#print(f'd_exp = {d_exp}')
## ::::::::::::::::::::::::::::::::::::::::::::::::::::::
#print('bulkCuO_mag_u8 total shift: ' + str(shift_bulkCuO_mag_u8 + d_exp))
#energy_bulkCuO_mag_u8 += shift_bulkCuO_mag_u8 + d_exp
#x_bulkCuO_mag_u8,y_bulkCuO_mag_u8  = givePeak(energy=energy_bulkCuO_mag_u8, intensity=intensity_bulkCuO_mag_u8, name='bulkCuO_mag_u8')


energy_bulkCuO, intensity_bulkCuO = norm_area("polAvg_bulkCuO.dat")
#x_bulkCuO,y_bulkCuO  = givePeak(energy=energy_bulkCuO, intensity=intensity_bulkCuO, name='bulkCuO')
#print(f'comp bulkCuO: {x_bulkCuO}')
shift_bulkCuO = getShift('david_bulkCuO.csv')
## ::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Experimental shift
#d_exp = x_bulkCuO_exp - (x_bulkCuO+shift_bulkCuO)
#print(f'd_exp = {d_exp}')
## ::::::::::::::::::::::::::::::::::::::::::::::::::::::
#print('bulkCuO total shift: ' + str(shift_bulkCuO + d_exp))
energy_bulkCuO += shift_bulkCuO + d_exp
x_bulkCuO,y_bulkCuO  = givePeak(energy=energy_bulkCuO, intensity=intensity_bulkCuO, name='bulkCuO')


# Process data

energy_watbox, intensity_watbox = norm_area("polAvg_watbox.dat")
shift_watbox = getShift('david_watbox.csv')
print('watbox total_shift: ' + str(shift_watbox + d_exp))
energy_watbox += shift_watbox+ d_exp
x_watbox,y_watbox  = givePeak(energy=energy_watbox, intensity=intensity_watbox, name='watbox')
energy_watbox_exp, intensity_watbox_exp = exp("exp_watbox.dat")
peaks_watbox_exp, _ = scipy.signal.find_peaks(intensity_watbox_exp, prominence=0.005)
x0_watbox_exp = energy_watbox_exp[peaks_watbox_exp[0]]
y0_watbox_exp = intensity_watbox_exp[peaks_watbox_exp[0]]
x_watbox_exp = energy_watbox_exp[peaks_watbox_exp[1]]
y_watbox_exp = intensity_watbox_exp[peaks_watbox_exp[1]]
#x_watbox_exp,y_watbox_exp  = givePeak(energy=energy_watbox_exp, intensity=intensity_watbox_exp, name='watbox_exp')

energy_bulkCu2O, intensity_bulkCu2O = norm_area("polAvg_bulkCu2O.dat")
shift_bulkCu2O = getShift('david_bulkCu2O.csv')
print('bulkCu2O total_shift: ' + str(shift_bulkCu2O + d_exp))
energy_bulkCu2O += shift_bulkCu2O+ d_exp
x_bulkCu2O,y_bulkCu2O  = givePeak(energy=energy_bulkCu2O, intensity=intensity_bulkCu2O, name='bulkCu2O')
print(f"Comp Bulk Cu2O peak: {x_bulkCu2O}")
energy_bulkCu2O_exp, intensity_bulkCu2O_exp = exp("exp_bulkCu2O.dat")
x_bulkCu2O_exp,y_bulkCu2O_exp  = givePeak(energy=energy_bulkCu2O_exp, intensity=intensity_bulkCu2O_exp, name='bulkCu2O_exp')

energy_bulkCuSO4, intensity_bulkCuSO4 = norm_area("polAvg_cuso4_bulk.dat")
shift_bulkCuSO4 = getShift('david_cuso4_bulk.csv')
print('bulkCuSO4 total_shift: ' + str(shift_bulkCuSO4 + d_exp))
energy_bulkCuSO4 += shift_bulkCuSO4+ d_exp
x_bulkCuSO4,y_bulkCuSO4  = givePeak(energy=energy_bulkCuSO4, intensity=intensity_bulkCuSO4, name='bulkCuSO4')
energy_bulkCuSO4_exp, intensity_bulkCuSO4_exp = exp("exp_bulkCuSO4.dat")
peaks_bulkCuSO4_exp, _ = scipy.signal.find_peaks(intensity_bulkCuSO4_exp, prominence=0.01)
x_bulkCuSO4_exp = energy_bulkCuSO4_exp[peaks_bulkCuSO4_exp[0]]
y_bulkCuSO4_exp = intensity_bulkCuSO4_exp[peaks_bulkCuSO4_exp[0]]
#x_bulkCuSO4_exp,y_bulkCuSO4_exp  = givePeak(energy=energy_bulkCuSO4_exp, intensity=intensity_bulkCuSO4_exp, name='bulkCuSO4_exp')
# Plot
plt.figure(figsize=(2.5, 3.5))
yoffset=1/3
exp_factor = 1.2
theory_factor = 1.5

# I will show water in other plot
plt.plot(energy_watbox, intensity_watbox*theory_factor+yoffset*-0.2, label='watbox None', linewidth=1.5, color='#515151', linestyle=':', alpha=0)
#if x_watbox is not None and y_watbox is not None:
#    plt.annotate(f'{x_watbox:.1f}', (x_watbox, y_watbox+yoffset*0),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 
plt.plot(energy_watbox_exp, intensity_watbox_exp*exp_factor+yoffset*0, label='watbox exp', linewidth=2, color='#515151', alpha=0)
#if x0_watbox_exp is not None and y0_watbox_exp is not None:
#    plt.annotate(f'{x0_watbox_exp:.1f}', (x0_watbox_exp, y0_watbox_exp*exp_factor+yoffset*0),
#    textcoords="offset points", xytext=(-9,2), 
#    ha='center', fontproperties=arial, color='black') 
#if x_watbox_exp is not None and y_watbox_exp is not None:
#    plt.annotate(f'{x_watbox_exp:.1f}', (x_watbox_exp, y_watbox_exp*exp_factor+yoffset*0),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 

plt.plot(energy_bulkCu2O, intensity_bulkCu2O*theory_factor+yoffset*0.8, label='bulkCu2O None', linewidth=1.5, color='#F14040', linestyle=':')
#if x_bulkCu2O is not None and y_bulkCu2O is not None:
#    plt.annotate(f'{x_bulkCu2O:.1f}', (x_bulkCu2O, y_bulkCu2O+yoffset*1),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 
plt.plot(energy_bulkCu2O_exp, intensity_bulkCu2O_exp*exp_factor+yoffset*1, label='bulkCu2O exp', linewidth=2, color='#F14040')
if x_bulkCu2O_exp is not None and y_bulkCu2O_exp is not None:
    plt.annotate(f'{x_bulkCu2O_exp:.1f}', (x_bulkCu2O_exp, y_bulkCu2O_exp*exp_factor+yoffset*1),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black') 

# Plot reference 
#plt.plot(energy_bulkCuO, intensity_bulkCuO*theory_factor+yoffset*1.8, label='bulkCuO None', linewidth=1.5, color='#1A6FDF', linestyle=':')
#if x_bulkCuO is not None and y_bulkCuO is not None:
#    plt.annotate(f'{x_bulkCuO:.1f}', (x_bulkCuO, y_bulkCuO+yoffset*2),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black')

plt.plot(energy_bulkCuO_mag, intensity_bulkCuO_mag*theory_factor+yoffset*1.8, label='bulkCuO_mag', linewidth=1.5, color='#1A6FDF', linestyle=':')
#plt.plot(energy_bulkCuO_mag_u8, intensity_bulkCuO_mag_u8*theory_factor+yoffset*1.8, label='bulkCuO_mag_u8', linewidth=1.5, color='#1A6FDF', linestyle=':')
#if x_bulkCuO_mag is not None and y_bulkCuO_mag is not None:
#    plt.annotate(f'{x_bulkCuO_mag:.1f}', (x_bulkCuO_mag, y_bulkCuO_mag+yoffset*2),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black')
plt.plot(energy_bulkCuO_exp, intensity_bulkCuO_exp*exp_factor+yoffset*2, label='bulkCuO exp', linewidth=2, color='#1A6FDF')
if x_bulkCuO_exp is not None and y_bulkCuO_exp is not None:
    plt.annotate(f'{x_bulkCuO_exp:.1f}', (x_bulkCuO_exp, y_bulkCuO_exp*exp_factor+yoffset*2),
    textcoords="offset points", xytext=(0,2), 
    ha='center', fontproperties=arial, color='black')

plt.plot(energy_bulkCuSO4, intensity_bulkCuSO4*theory_factor+yoffset*2.8, label='bulkCuSO4', linewidth=1.5, color='#37AD6B', linestyle=':')
#if x_bulkCuSO4 is not None and y_bulkCuSO4 is not None:
#    plt.annotate(f'{x_bulkCuSO4:.1f}', (x_bulkCuSO4, y_bulkCuSO4+yoffset*3),
#    textcoords="offset points", xytext=(0,2), 
#    ha='center', fontproperties=arial, color='black') 
plt.plot(energy_bulkCuSO4_exp, intensity_bulkCuSO4_exp*exp_factor+yoffset*3, label='bulkCuSO4 exp', linewidth=2, color='#37AD6B')
if x_bulkCuSO4_exp is not None and y_bulkCuSO4_exp is not None:
    plt.annotate(f'{x_bulkCuSO4_exp:.1f}', (x_bulkCuSO4_exp, y_bulkCuSO4_exp*exp_factor+yoffset*3),
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
#plt.yticks(fontproperties=arial)
plt.xticks(fontproperties=arialbd)
plt.xlim(E_i_exp, E_f_exp)
#plt.grid(True)
#handles, labels = plt.gca().get_legend_handles_labels()
#plt.legend(handles[::-1], labels[::-1], fontsize=15,loc='upper left', bbox_to_anchor=(1, 1))
plt.tight_layout()

# Saving or showing plot
plt.savefig('figS-bulk-o.svg', format='svg', bbox_inches='tight')
#plt.savefig('xasM2.png', format='png', dpi=300, bbox_inches='tight')
plt.show()
