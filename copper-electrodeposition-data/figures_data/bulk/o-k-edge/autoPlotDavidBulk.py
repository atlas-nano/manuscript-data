import argparse
import sys
import textwrap

startOfScript ="""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
import scipy.signal
import textwrap

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

endOfScript = """
# Plot settings
#plt.figure(figsize=(4, 3.5))
plt.xlabel('Energy (eV)', fontsize=7)
#tick_positions = np.arange(928, 945, 1)
#tick_labels = [str(tick) if tick % 2 == 0 else '' for tick in tick_positions]
#plt.xticks(tick_positions, tick_labels,fontsize=7)
plt.ylabel('Intensity (a.u.)', fontsize=7)
plt.yticks(fontsize=7)
plt.xticks(fontsize=7)
plt.xlim(E_i+d_exp, E_f+d_exp)
plt.grid(True)
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1], fontsize=7,loc='upper left', bbox_to_anchor=(1, 1))
plt.tight_layout()

# Saving or showing plot
#plt.savefig('xasM2.svg', format='svg', bbox_inches='tight')
#plt.savefig('xasM2.png', format='png', dpi=300, bbox_inches='tight')
plt.show()
"""
colors = ["black", "red", "blue", "deepskyblue", "deepskyblue", "blue", "blue", "darkgreen", "purple", "darkorange", "darkorange"]
ynumber = [0,1,2,3,4,5,6,7,8,9,10]

def create_plot(basenames, reference, suffix=None, method="norm_area"):
    if suffix is None:
        new_file = f"{method}_david.py"
        
        with open(new_file, "w") as file:
            file.write(startOfScript)
            file.write(textwrap.dedent(f"""
            # Experimental reference data
            energy_{reference}_exp, intensity_{reference}_exp = exp("exp_{reference}.dat")
            x_{reference}_exp,y_{reference}_exp  = givePeak(energy=energy_{reference}_exp, intensity=intensity_{reference}_exp, name='{reference}')
            # Computational reference
            energy_{reference}, intensity_{reference} = norm_area("polAvg_{reference}.dat")
            x_{reference},y_{reference}  = givePeak(energy=energy_{reference}, intensity=intensity_{reference}, name='{reference}')
            shift_{reference} = getShift('david_{reference}.csv')
            # ------------------------------------------------------------------------------------------------------------
            # Experimental shift
            d_exp = x_{reference}_exp - (x_{reference}+shift_{reference})
            print(f'd_exp = {{d_exp}}')
            # ------------------------------------------------------------------------------------------------------------
            print('{reference} total shift: ' + str(shift_{reference} + d_exp))
            energy_{reference} += shift_{reference} + d_exp
            x_{reference},y_{reference}  = givePeak(energy=energy_{reference}, intensity=intensity_{reference}, name='{reference}')
            """)
            )
            file.write('\n# Process data\n')
            for b in basenames:
                file.write(textwrap.dedent(f"""
                    energy_{b}, intensity_{b} = {method}("polAvg_{b}.dat")
                    shift_{b} = getShift('david_{b}.csv')
                    print('{b} total_shift: ' + str(shift_{b} + d_exp))
                    energy_{b} += shift_{b}+ d_exp
                    x_{b},y_{b}  = givePeak(energy=energy_{b}, intensity=intensity_{b}, name='{b}')
                    energy_{b}_exp, intensity_{b}_exp = exp("exp_{b}.dat")
                    """)
                )
            file.write('# Plot\n')
            file.write('# Plot reference \n')
            file.write(textwrap.dedent(f"""
                plt.plot(energy_{reference}, intensity_{reference}+yoffset*{ynumber[0]}, label='{reference} {suffix}', linewidth=1, color='{colors[0]}', linestyle='--')
                if x_{reference} is not None and y_{reference} is not None:
                    plt.annotate(f'{{x_{reference}:.1f}}', (x_{reference}, y_{reference}+yoffset*{ynumber[0]}),
                    textcoords="offset points", xytext=(0,2), 
                    ha='center', fontsize=7, color='black')
                plt.plot(energy_{reference}_exp, intensity_{reference}_exp/13+yoffset*{ynumber[0]}, label='{reference} exp', linewidth=1, color='{colors[0]}')
                """)
                )
            i = 1
            for b in basenames:
                file.write(textwrap.dedent(f"""
                    plt.plot(energy_{b}, intensity_{b}+yoffset*{ynumber[i]}, label='{b} {suffix}', linewidth=1, color='{colors[i]}', linestyle='--')
                    if x_{b} is not None and y_{b} is not None:
                        plt.annotate(f'{{x_{b}:.1f}}', (x_{b}, y_{b}+yoffset*{ynumber[i]}),
                        textcoords="offset points", xytext=(0,2), 
                        ha='center', fontsize=7, color='black') 
                    plt.plot(energy_{b}_exp, intensity_{b}_exp/13+yoffset*{ynumber[i]}, label='{b} exp', linewidth=1, color='{colors[i]}')
                    """)
                    )
                i = i+1 
            file.write(endOfScript)
            print(f"{method}_david.py created successfully!")
    else:
        new_file = f"{method}_{suffix}_david.py"
        with open(new_file, "w") as file:
            file.write(startOfScript)
            file.write(textwrap.dedent(f"""
                # Experimental reference data
                energy_{reference}_exp, intensity_{reference}_exp = exp("exp_{reference}.dat")
                x_{reference}_exp,y_{reference}_exp  = givePeak(energy=energy_{reference}_exp, intensity=intensity_{reference}_exp, name='{reference}')
                # Computational reference
                energy_{reference}_{suffix}, intensity_{reference}_{suffix} = norm_area("polAvg_{reference}_{suffix}.dat")
                x_{reference}_{suffix},y_{reference}_{suffix}  = givePeak(energy=energy_{reference}_{suffix}, intensity=intensity_{reference}_{suffix}, name='{reference}')
                shift_{reference} = getShift('david_{reference}.csv')
                # ------------------------------------------------------------------------------------------------------------
                # Experimental shift
                d_exp = x_{reference}_exp - (x_{reference}_{suffix}+shift_{reference})
                print(f'd_exp = {{d_exp}}')
                # ------------------------------------------------------------------------------------------------------------
                print('{reference} total shift: ' + str(shift_{reference} + d_exp))
                energy_{reference}_{suffix} += shift_{reference} + d_exp
                x_{reference}_{suffix},y_{reference}_{suffix}  = givePeak(energy=energy_{reference}_{suffix}, intensity=intensity_{reference}_{suffix}, name='{reference}')
                """)
            )
            file.write('\n# Process data\n')
            for b in basenames:
                file.write(textwrap.dedent(f"""
                    energy_{b}_{suffix}, intensity_{b}_{suffix} = {method}("polAvg_{b}_{suffix}.dat")
                    shift_{b} = getShift('david_{b}.csv')
                    print('{b}_{suffix} total_shift: ' + str(shift_{b} + d_exp))
                    energy_{b}_{suffix} += shift_{b}+ d_exp
                    x_{b}_{suffix},y_{b}_{suffix}  = givePeak(energy=energy_{b}_{suffix}, intensity=intensity_{b}_{suffix}, name='{b}')
                    energy_{b}_exp, intensity_{b}_exp = exp("exp_{b}.dat")
                    """)
                )
            file.write('# Plot\n')
            file.write('# Plot reference \n')
            file.write(textwrap.dedent(f"""
                plt.plot(energy_{reference}_{suffix}, intensity_{reference}_{suffix}+yoffset*{ynumber[0]}, label='{reference} {suffix}', linewidth=1, color='{colors[0]}', linestyle='--')
                if x_{reference}_{suffix} is not None and y_{reference}_{suffix} is not None:
                    plt.annotate(f'{{x_{reference}_{suffix}:.1f}}', (x_{reference}_{suffix}, y_{reference}_{suffix}+yoffset*{ynumber[0]}),
                    textcoords="offset points", xytext=(0,2), 
                    ha='center', fontsize=7, color='black')
                plt.plot(energy_{reference}_exp, intensity_{reference}_exp/13+yoffset*{ynumber[0]}, label='{reference} exp', linewidth=1, color='{colors[0]}')
                """)
                )
            i = 1
            for b in basenames:
                file.write(textwrap.dedent(f"""
                    plt.plot(energy_{b}_{suffix}, intensity_{b}_{suffix}+yoffset*{ynumber[i]}, label='{b} {suffix}', linewidth=1, color='{colors[i]}', linestyle='--')
                    if x_{b}_{suffix} is not None and y_{b}_{suffix} is not None:
                        plt.annotate(f'{{x_{b}_{suffix}:.1f}}', (x_{b}_{suffix}, y_{b}_{suffix}+yoffset*{ynumber[i]}),
                        textcoords="offset points", xytext=(0,2), 
                        ha='center', fontsize=7, color='black') 
                    plt.plot(energy_{b}_exp, intensity_{b}_exp/13+yoffset*{ynumber[i]}, label='{b} exp', linewidth=1, color='{colors[i]}')
                    """)
                    )
                i = i+1 
            file.write(endOfScript)
            print(f"{method}_{suffix}_david.py created successfully!")

def main():
    parser = argparse.ArgumentParser(
        usage="""
        python autoPlot.py  -b basename1 [basename2 ...] -r reference [-s suffix -m method]
            -b: List of basenames
            -r: reference calculation: bulkCu.
            -s: (Optional) Suffix. 
            -m: (Optional) method: raw or norm_area
        """
    )
    parser.add_argument("-s", "--suffix", required=False)
    parser.add_argument("-m", "--method", choices=["raw", "norm_area"], required=False)
    parser.add_argument("-b", "--basename", nargs='+', required=True)
    parser.add_argument("-r", "--reference", required=True)

    if len(sys.argv) == 1:  # No arguments provided
        parser.print_usage()
        sys.exit(1)

    args = parser.parse_args()

    if args.method is None:
        create_plot(basenames = args.basename, suffix = args.suffix, reference=args.reference)
    else:
        create_plot(basenames = args.basename, suffix = args.suffix, method=args.method, reference=args.reference)


if __name__ == "__main__":
    main()

