from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from matplotlib import font_manager

arial = font_manager.FontProperties(fname="/mnt/c/Windows/Fonts/arial.ttf", size=9)
arialbd = font_manager.FontProperties(fname="/mnt/c/Windows/Fonts/arialbd.ttf", size=9)

# Constants
# Ka1 = 10**3  # First dissociation constant for H2SO4
Ka2 = 1.03e-2 # Harris
pH = 4.5
H = 10**-pH
OH = 10**(pH-14)
Cu2 = 0.1  # From CuSO4 dissociation
S_total = Cu2  # From 1:1 CuSO4

# Define function system: [0] = SO4^2-, [1] = HSO4-
def equations(vars):
    so4, hso4 = vars
    # eq1 = Ka1 - (H * hso4 / h2so4)                       # acid-base equilibrium
    eq2 = Ka2 - (H * so4 / hso4)                       # acid-base equilibrium
    eq3 = 2*Cu2 + H - 2*so4 - hso4 - OH                 # charge balance
    return [eq2, eq3]

# Initial guesses
initial_guess = [0.09, 0.01]

# Solve
so4, hso4 = fsolve(equations, initial_guess)

# Output
print("Equilibrium concentrations (mol/L):")
print(f"[Cu²⁺]     = {Cu2:.4e} M")
print(f"[SO₄²⁻]    = {so4:.4e} M")
print(f"[HSO₄⁻]    = {hso4:.4e} M")
# print(f"[H2SO4]    = {h2so4:.4e} M")
print(f"[H⁺]       = {H:.4e} M  (from pH = {pH})")
print(f"[OH-]       = {OH:.4e} M  (from pH = {pH})")

# -------------------------------
# 2. Near-Electrode Concentration using GCS model (kB/e version)
# -------------------------------

# Physical constants
eps0 = 8.854e-12      # vacuum permittivity (F/m)
epsr = 78.5           # relative permittivity of water
NA = 6.022e23         # Avogadro's number (1/mol)
e = 1.602e-19         # elementary charge (C)
kB = 1.381e-23        # Boltzmann constant (J/K)
T = 298.15            # temperature (K)
phi_total = -0.1      # total electrode potential (V)

# 1. Ionic strength (mol/L)
I = 0.5 * (Cu2*2**2 + so4*2**2 + hso4*1**2 + H*1**2 + OH*1**2)
number_m3__mol_L = 1e3*NA  # conversion factor

# 2. Debye length (m), in terms of kB and e
lambda_D = np.sqrt((epsr * eps0 * kB * T) / (2 * e**2 * I*number_m3__mol_L))  # meters

# Capacitance (F/m²), e.g., 20 µF/cm² = 0.2 F/m²
x_2 = 6.0 * 1e-10 # m
C_s = 0.2  # F/m²
epsr_stern = 6.0  # dielectric constant in Stern layer (typical range 4–10)
C_s  = epsr_stern * eps0 / x_2 # F/m² 
print(f"x_2 (OHP) = {x_2*1e10:.4f} Å")

# Valence of electrolytes
z = 2 

# 4. Solve for diffuse potential ψ_d (Gouy–Chapman)
# σ = sqrt(8 εε₀ kBT * NA * I) sinh(e ψ_d / 2kBT)
def surface_charge(phi_2):
    prefactor = np.sqrt(8 * epsr * eps0 * kB * T * Cu2*number_m3__mol_L)
    return prefactor * np.sinh((z*e * phi_2) / (2 * kB * T))  # C/m²

def total_potential_residual(phi_2):
    # phi_2 = phi_total + dphi/dx * x_2
    # phi_2 = phi_total - (surface_charge / (epsr * eps0 ))*x_2
    return phi_total - (surface_charge(phi_2) / (epsr * eps0 ))*x_2 - phi_2

sol = root_scalar(total_potential_residual, bracket=[-1.0, 1.0], method='bisect')
phi_2 = sol.root  # potential at diffuse layer / OHP
print(f"phi_2 = {phi_2*1000:.4f} mV")
C_diffuse = epsr * eps0 / lambda_D
print(f"C_s       = {C_s:.3f} F/m²")
print(f"C_diffuse = {C_diffuse:.3f} F/m²")

# 5. Compute near-electrode concentrations
def gcs_concentration(c_bulk, z, phi):
    return c_bulk * np.exp((-z * e * phi) / (kB * T))

def phi_profile(x, phi_total, phi_2, x_2):
    """Vectorized potential profile calculation."""
    arg = np.tanh((z * e * phi_2) / (4 * kB * T)) * np.exp(-(x - x_2) / lambda_D)
    #slope = (((8*kB*T*Cu2)/(eps0*epsr))**0.5) * np.sinh((z * e * phi_2) / (2 * kB * T))
    #slope = - (2 * kB * T) / (z * e * lambda_D) * np.sinh((z * e * phi_2) / (2 * kB * T))
    slope = (phi_2 - phi_total) / x_2 
    return np.where(
        x > x_2,
        ((4 * kB * T) / (z * e)) * np.arctanh(arg),
        slope*x + phi_total #slope * (x-x_2) + phi_2
    )  
phi_x = phi_2 # phi_profile(lambda_D, phi_total, phi_2, x_2)  p # -0.085  # Use the calculated diffuse potential
Cu2_near = gcs_concentration(Cu2, +2, phi_x)
so4_near = gcs_concentration(so4, -2, phi_x)
hso4_near = gcs_concentration(hso4, -1, phi_x)
H_near = gcs_concentration(H, +1, phi_x)
pH_near = -np.log10(H_near)
OH_near = gcs_concentration(OH, -1, phi_x)

# -------------------
# Output GCS Results
# -------------------
print("\n=== GCS Concentrations at OHP ===")
print(f"Ionic strength (I)     : {I:.4e} mol/L")
print(f"Debye length (λ_D)     : {lambda_D*1e9:.2f} nm")
print(f"Potential at OHP phi_2 : {phi_2:.4f} V")
print(f"[Cu²⁺]_near            : {Cu2_near:.4e} M")
print(f"[SO₄²⁻]_near           : {so4_near:.4e} M")
print(f"[HSO₄⁻]_near           : {hso4_near:.4e} M")
print(f"[H⁺]_near              : {H_near:.4e} M")
print(f"pH_near                : {pH_near:.1f}")
print(f"[OH⁻]_near             : {OH_near:.4e} M")
 
# Enrichment/Depletion Factors
print("\n=== GCS Enrichment/Depletion Factors at OHP ===")
print(f"Cu²⁺     : {Cu2_near / Cu2:.2f}×")
print(f"SO₄²⁻    : {so4_near / so4:.2f}×")
print(f"HSO₄⁻    : {hso4_near / hso4:.2f}×")
print(f"H⁺       : {H_near / H:.2f}×")
print(f"OH⁻      : {OH_near / OH:.2f}×")


# ------------------------------
# 6. Total Oxygen Atom Analysis
# ------------------------------

# Volume near electrode 
volume_m3 = 1e-18   # m³ or 1 µm³
volume_L = volume_m3 * 1000  # liters
H2O_density = 1e6 # g/m³
H2O_FM = 18.0  # g/mol (molar mass of water)

H2O_molecules = ((volume_m3 * H2O_density) / (H2O_FM))*NA
print(f"H2O_molecules (density method): {H2O_molecules:.2e}")

# Number of molecules in volume
def num_molecules(conc_M):
    return conc_M * volume_L * NA  # molecules in volume

sum_molecules = 0

for molecule in [Cu2_near, so4_near, hso4_near, H_near, OH_near]:
    sum_molecules += num_molecules(molecule)

H2O_molecules -= sum_molecules  # Adjust water molecules after accounting for other species

# Number of oxygen atoms
O_sulfate  = 4 * num_molecules(so4_near)    # SO₄²⁻ has 4 O
O_bisulfate = 4 * num_molecules(hso4_near)  # HSO₄⁻ has 4 O
O_hydronium = 1 * num_molecules(H_near)     # H₃O⁺ has 1 O
O_water     = 1 * H2O_molecules   # H₂O has 1 O

# Total oxygen atoms
O_total = O_sulfate + O_bisulfate + O_hydronium + O_water

# -----------------------
# Output Oxygen Breakdown
# -----------------------
print("\n=== Oxygen Atom Breakdown (in 1 um³) ===")
print(f"O from SO₄²⁻     : {O_sulfate:.2e} atoms ({O_sulfate/O_total:.2%})")
print(f"O from HSO₄⁻     : {O_bisulfate:.2e} atoms ({O_bisulfate/O_total:.2%})")
print(f"O from H₃O⁺      : {O_hydronium:.2e} atoms ({O_hydronium/O_total:.2%})")
print(f"O from H₂O       : {O_water:.2e} atoms ({O_water/O_total:.2%})")
print(f"Total O atoms    : {O_total:.2e}")


# ------------------------------
# 7. Potential Profile vs. Distance (Gouy–Chapman)
# ------------------------------

# Gouy–Chapman potential profile

# Distance array (from interface, in meters)
x = np.linspace(0, 5*lambda_D, 200)  # up to 5 Debye lengths




plt.figure(figsize=(3.5 , 3.0))

plt.plot(x*1e10, phi_profile(x, phi_total , phi_2, x_2), label=r"$\phi(x)$")
# Add vertical lines for x_2 and lambda_D
plt.axvline(x_2*1e10, color='r', linestyle='--', label=r"$x_2$ (OHP)")
plt.axvline(lambda_D*1e10, color='g', linestyle='--', label=r"$\lambda_D$ (Debye length)")
plt.xlabel("Distance from interface (Å)", fontproperties=arialbd)
plt.ylabel("Potential $\phi$ (V)", fontproperties=arialbd)
plt.title("Potential Profile (GCS Model)", fontproperties=arialbd)
plt.xticks(fontproperties=arialbd)
plt.yticks(fontproperties=arialbd)
plt.grid(True)
plt.legend(prop=arial)
plt.tight_layout()
plt.savefig('phi_profile.svg', format='svg', bbox_inches='tight')
plt.show()


# Plot
phi_x_profile = phi_profile(x, phi_total, phi_2, x_2)

# Plot concentration profiles for each species
species = {
    'Cu$^{2+}$': (Cu2, +2),
    'SO$_4^{2-}$': (so4, -2),
    'HSO$_4^-$': (hso4, -1),
    'H$_3$O$^+$': (H, +1),
    'OH$^-$': (OH, -1)
}

plt.figure(figsize=(4.0 , 3.0))

for label, (c_bulk, z_val) in species.items():
    c_x = gcs_concentration(c_bulk, z_val, phi_x_profile)
    plt.plot(x*1e10, c_x, label=label)
plt.axvline(x_2*1e10, color='r', linestyle='--', label=r"$x_2$")
plt.axvline(lambda_D*1e10, color='g', linestyle='--', label=r"$\lambda_D$")
plt.xlabel("Distance from interface (Å)", fontproperties=arialbd)
ax = plt.gca()
# move y-axis label and ticks to the right side
ax.yaxis.set_label_position('right')
ax.yaxis.tick_right()
ax.set_ylabel("Concentration (mol/L)", fontproperties=arialbd,
              rotation=270, rotation_mode='anchor', labelpad=12, va='center')
plt.xticks(fontproperties=arialbd)
plt.yticks(fontproperties=arialbd)
plt.title("Concentration Profiles (GCS Model)", fontproperties=arialbd)
plt.yscale('log')
plt.grid(True, which='both', ls='--')
# Ensure the left edge of the plotting box (x-axis) starts exactly at 0 Å
ax.set_xlim(0, x[-1]*1e10)
# Place textual labels near the right end of each curve (uses data coords + small offset in points)
# move species labels to the left side of the plot (near the interface)
x_label_pos_left = x[0] + (x[1]-x[0])*4  # a small offset from x=0
for label, (c_bulk, z_val) in species.items():
    c_x = gcs_concentration(c_bulk, z_val, phi_x_profile)
    # interpolate concentration at the chosen left x position
    c_val = np.interp(x_label_pos_left, x, c_x)

    # annotate slightly to the left in display points
    plt.annotate(label,
                 xy=(x_label_pos_left*1e10, c_val),
                 xytext=(-30, 0),
                 textcoords='offset points',
                 fontproperties=arial,
                 va='center')

#plt.legend(prop=arial, loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.savefig('conc_profile.svg', format='svg', bbox_inches='tight')
plt.show()
