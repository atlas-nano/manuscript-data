import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import re
from matplotlib import font_manager

arial = font_manager.FontProperties(fname="/mnt/c/Windows/Fonts/arial.ttf", size=9)
arialbd = font_manager.FontProperties(fname="/mnt/c/Windows/Fonts/arialbd.ttf", size=9)

# --- Parameters ---
dos_file = "system.dos"   # QE DOS file
fermi = 12.077             # Fermi level (eV)
bands_file = "CuO_bands2.dat"
nbnd = 40                  # Number of bands (adjust if needed)

# --- Load DOS ---
dos_data = np.loadtxt(dos_file, comments="#")
energy_dos = dos_data[:,0] #- fermi
energy_dos -= fermi  # shift energy relative to EFermi
dos_up = dos_data[:,1]
dos_dw = -dos_data[:,2]  # negative for plotting spin-down

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
fig = plt.figure(figsize=(4.5, 4.0))
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
ppout_file = 'bands2/pp.bands2.out'
hs_xcoords = []
label_map = [
    'Γ', 'Y', 'H', 'C', 'E', 'M₁', 'A', 'X', 'Γ', 'Z', 'D', 'Y'
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

# --- Band gap shading ---
vbm = -1e9
cbm = 1e9
for i in range(bands.shape[0]):
    for j in range(bands.shape[1]):
        E = bands[i, j]
        if E <= fermi and E > vbm:
            vbm = E
        if E >= fermi and E < cbm:
            cbm = E
ax0.axhspan(vbm - fermi, cbm - fermi, color='gray', alpha=0.3)
print(f"VBM (eV): {vbm:.6f}")
print(f"CBM (eV): {cbm:.6f}")
print(f"Band gap (eV): {cbm-vbm:.6f}")

# DOS (right)
ax1 = fig.add_subplot(gs[1], sharey=ax0)
ax1.plot(dos_up, energy_dos, label="Spin Up", color='grey')
ax1.plot(dos_dw, energy_dos, label="Spin Down", color='b')
ax1.axhline(0, linestyle="--", color="black")
ax1.set_xlabel("DOS", fontproperties=arialbd)
#ax1.set_title("DOS")
#ax1.legend()
ax1.grid(True, linestyle="--", alpha=0.5)
plt.setp(ax1.get_yticklabels(), visible=False)
plt.ylim(-5, 5)
plt.tight_layout()
plt.savefig("bulkCuO_band_dos_combined.svg", format='svg')
plt.show()
