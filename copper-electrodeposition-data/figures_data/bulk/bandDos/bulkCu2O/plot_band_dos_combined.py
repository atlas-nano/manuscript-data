import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import re
from matplotlib import font_manager

arial = font_manager.FontProperties(fname="/mnt/c/Windows/Fonts/arial.ttf", size=9)
arialbd = font_manager.FontProperties(fname="/mnt/c/Windows/Fonts/arialbd.ttf", size=9)


# --- Parameters ---
dos_file = "system.dos"   # QE DOS file
fermi = 10.442             # Fermi level (eV)
bands_file = "bands.dat"
nbnd = 60                  # Number of bands (adjust if needed)

# --- Load DOS ---
dos_data = np.loadtxt(dos_file, comments="#")
# Ensure 2D shape when a single row is present
if dos_data.ndim == 1:
    dos_data = dos_data.reshape(-1, dos_data.size)
energy_dos = dos_data[:, 0]
energy_dos -= fermi  # shift energy relative to EFermi

# Heuristic detection of spin-polarization:
# - Many non-spin QE DOS files have columns: E, dos(E), Int dos(E)
#   where the last column (integrated DOS) is monotonic increasing.
# - Spin-polarized DOS usually has E, dos_up, dos_dn, ... where the
#   second and third columns are the two spin channels.
def _is_monotonic_increasing(arr):
    return np.all(np.diff(arr) >= -1e-12)

spin_polarized = False
dos_total = None
dos_up = None
dos_dw = None

ncols = dos_data.shape[1]
if ncols >= 3 and _is_monotonic_increasing(dos_data[:, -1]):
    # Treat as non-spin-polarized: columns are E, dos, int_dos
    dos_total = dos_data[:, 1]
    spin_polarized = False
else:
    # Otherwise assume spin-polarized: E, dos_up, dos_dn, ...
    # Mirror spin-down for plotting (negative values)
    if ncols >= 3:
        dos_up = dos_data[:, 1]
        dos_dw = -dos_data[:, 2]
        spin_polarized = True
    elif ncols == 2:
        # Minimal case: E, dos
        dos_total = dos_data[:, 1]
        spin_polarized = False
    else:
        raise ValueError(f"Unexpected DOS file format: {dos_file} has {ncols} columns")

# --- Load Bands ---
kpoints = []
bands = []
with open(bands_file, 'r') as f:
    lines = f.readlines()
# Skip header lines
band_lines = []
for line in lines:
    if line.strip().startswith('&plot') or line.strip().endswith('/'):
        continue
    if line.strip() == '':
        continue
    band_lines.append(line)
i = 0
while i < len(band_lines):
    kpt = [float(x) for x in band_lines[i].strip().split()]
    kpoints.append(kpt)
    i += 1
    band_vals = []
    for _ in range(nbnd//10):  # 10 bands per line
        band_vals.extend([float(x) for x in band_lines[i].strip().split()])
        i += 1
    bands.append(band_vals)
kpoints = np.array(kpoints)
bands = np.array(bands)

# Calculate k-path distance
kpath = np.zeros(len(kpoints))
for j in range(1, len(kpoints)):
    dk = np.linalg.norm(kpoints[j] - kpoints[j-1])
    kpath[j] = kpath[j-1] + dk

# --- Plot Combined ---

fig = plt.figure(figsize=(4.0, 4.0))
gs = gridspec.GridSpec(1, 2, width_ratios=[3,1], wspace=0.05)

# Band structure (left)
ax0 = fig.add_subplot(gs[0])
for i in range(bands.shape[1]):
    ax0.plot(kpath, bands[:, i] - fermi, color='b', lw=1)
ax0.axhline(0, linestyle='--', label='Fermi level', color="black")
ax0.set_ylabel("Energy (eV)", fontproperties=arialbd)
ax0.set_xlabel("k-path", fontproperties=arialbd)
#ax0.set_title("Band Structure")
ax0.grid(True, linestyle="--", alpha=0.5)

# --- High-symmetry points ---
ppout_file = 'pp.bands.out'
hs_xcoords = []
label_map = [
    'X', 'M', 'Γ', 'X', 'R', 'Γ', 'M', 'R'
]

with open(ppout_file, 'r') as f:
    for line in f:
        m = re.search(r'high-symmetry point:.*x coordinate\s+([0-9.]+)', line)
        if m:
            hs_xcoords.append(float(m.group(1)))
for x in hs_xcoords:
    ax0.axvline(x, color='k', linestyle='--', lw=0.7)
ax0.set_xticks(hs_xcoords)
ax0.set_xticklabels(label_map[:len(hs_xcoords)])

# --- Band gap shading (robust detection) ---
# Use a small tolerance so bands exactly equal to EF aren't counted as both
tol = 1e-6
bands_flat = bands.flatten()
# energies strictly below EF
below = bands_flat[bands_flat < fermi - tol]
# energies strictly above EF
above = bands_flat[bands_flat > fermi + tol]
if below.size > 0:
    vbm = below.max()
else:
    # fallback: include energies equal to EF if no strictly-below values
    le = bands_flat[bands_flat <= fermi + tol]
    vbm = le.max() if le.size > 0 else bands_flat.min()
if above.size > 0:
    cbm = above.min()
else:
    # fallback: include energies equal to EF if no strictly-above values
    ge = bands_flat[bands_flat >= fermi - tol]
    cbm = ge.min() if ge.size > 0 else bands_flat.max()
ax0.axhspan(vbm - fermi, cbm - fermi, color='gray', alpha=0.3)
print(f"VBM (eV): {vbm:.6f}")
print(f"CBM (eV): {cbm:.6f}")
print(f"Band gap (eV): {cbm-vbm:.6f}")
if cbm - vbm <= tol:
    print("Note: calculated band gap is zero (or within tolerance). Check EF and band energies.")

# DOS (right)
ax1 = fig.add_subplot(gs[1], sharey=ax0)
if spin_polarized:
    ax1.plot(dos_up, energy_dos, label="Spin Up", color='tab:gray')
    ax1.plot(dos_dw, energy_dos, label="Spin Down", color='tab:blue')
    ax1.set_title("Spin-resolved DOS")
else:
    ax1.plot(dos_total, energy_dos, label="DOS", color='tab:gray')
    #ax1.set_title("DOS")

ax1.axhline(0, linestyle="--", color="black")
ax1.set_xlabel("DOS", fontproperties=arialbd)
ax1.grid(True, linestyle="--", alpha=0.5)
#ax1.legend()
plt.setp(ax1.get_yticklabels(), visible=False)
# Align DOS plot y-limits with the band structure
ax1.set_ylim(-5,5)
plt.tight_layout()
plt.savefig("bulkCu2O_band_dos_combined.svg", format='svg')
plt.show()
