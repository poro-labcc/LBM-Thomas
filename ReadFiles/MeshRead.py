import meshio
import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes,mark_inset

os.chdir("/home/thomas/remote/hal/ReadFiles")
Reynolds = 300
Mark = "E"

mass = np.load(f"massRE{Reynolds}_Mrk{Mark}.npy")
mass = mass / (320*2000)

plt.figure(figsize=(8, 6))
plt.title(f"Mass over time Re={Reynolds} Maximum Velocity CBDC")
plt.xlabel("Timesteps")
plt.ylabel(r"$M/M_{O}$")
plt.grid(True)

# Gráfico principal
plt.plot(mass, label="Massa normalizada")

# Criar o inset (zoom)
ax_inset = inset_axes(plt.gca(), width="70%", height="30%", loc="lower right", borderpad=2)
ax_inset.plot(np.arange(100, len(mass)), mass[100:], color='orange')
ax_inset.set_title("Zoom", fontsize=8)
ax_inset.grid(True)
ax_inset.tick_params(axis='both', labelsize=8)

plt.savefig(f"/home/thomas/Documents/TCC/ReadFiles/MassImgs/mass_RE{Reynolds}_Mrk{Mark}.png")
plt.show()


# Carregamento e normalização dos dados
FlowIN = np.load(f"massRE{Reynolds}_Mrk{Mark}_MassIN.npy") / (212 * 0.0577)
FlowOUT = np.load(f"massRE{Reynolds}_Mrk{Mark}_MassOUT.npy") / (212 * 0.0577)
timesteps = np.arange(len(FlowIN))

# ====== Primeiro gráfico: FlowIN e FlowOUT ======
fig1, ax1 = plt.subplots(figsize=(8, 6))
crop = 5
ax1.set_title(f"Flow over time Re={Reynolds} Maximum Velocity CBDC")
ax1.set_xlabel("Timesteps")
ax1.set_ylabel(r"$Flow$")
ax_inset = inset_axes(ax1, width="70%", height="30%", loc="center right", borderpad=2)
ax_inset.plot(timesteps[crop:], FlowIN[crop:], color="blue")
ax_inset.plot(timesteps[crop:], FlowOUT[crop:], color="green")
ax_inset.grid(True)
ax_inset.tick_params(axis='both', labelsize=8)
ax1.grid(True)

ax1.plot(timesteps, FlowIN, label="Inlet", color="blue")
ax1.plot(timesteps, FlowOUT, label="Outlet", color="green")
ax1.legend(loc="best")

mark_inset(ax1, ax_inset, loc1=2, loc2=4, fc="none", ec="0.5",color="red")
plt.savefig(f"/home/thomas/Documents/TCC/ReadFiles/FlowIN_OUT_Re{Reynolds}_Mrk{Mark}.png")

# ====== Segundo gráfico: |FlowIN - FlowOUT| em escala log ======
fig2, ax2 = plt.subplots(figsize=(8, 6))

ax2.set_title(f"Mass conservation error Re={Reynolds}")
ax2.set_xlabel("Timesteps")
ax2.set_ylabel(r"$|Flow_{in} - Flow_{out}|$")
ax2.set_yscale("log")
ax2.grid(True, which="both", linestyle='--')

# Evitar valores zero no log
eps = 1e-16
difference = np.abs(FlowIN - FlowOUT) + eps
ax2.plot(timesteps, difference, color="red", label="|Inlet - Outlet|", linestyle='--')
ax2.legend(loc="best")

plt.savefig(f"/home/thomas/Documents/TCC/ReadFiles/FlowDifference_Re{Reynolds}_Mrk{Mark}.png")

# Mostrar ambos
plt.show()

