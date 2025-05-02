import numpy as np
import matplotlib.pyplot as plt

my_data = np.load("./plots/Cd_St_Re.npy")
# Data
sample = np.array([1.46945313477662, 1.43068243440537, 1.40284069578276, 1.38210278229675,1.3418011930088, 1.37035414841697, 1.46010720393776, 1.58046239519459])
Breuer = np.array([1.4592093082354651, 1.400324916955583, 1.3503530522226863, 1.3206258735546643, 1.3305077054328296, 1.3605852135557535, 1.490203481512407,1.5509466155996445])
#my_data = np.array([1.4698913746563362, 1.430306946610678, 1.400605948810238, 1.3799417982017983, 1.3465266733266734, 1.3749573026973025,1.4665992707292705, 1.5874802697302697])
cd_mine = np.sort(my_data,order="Re")["Cd"]

StSample = np.array([0.126306483300589, 0.131375245579568, 0.135677799607073, 0.139155206286837, 0.147760314341847, 0.140333988212181, 0.130373280943026, 0.129666011787819])
#StMyData = np.array([0.114287, 0.119141, 0.121912, 0.124683, 0.148530, 0.140609, 0.130707, 0.128726])
St_mine = np.sort(my_data,order="Re")["St"]

x = [70,80,90,100,150,200,250,300]
x_b = [65,80,100,130,165, 200,270, 300]

# Create figure and axes
fig, ax = plt.subplots(1, 2, figsize=(10, 5))  # Adjust figure size

# Main title
fig.suptitle("Comparison Wang (2015) vs. In-House LBM-BGK Code", fontsize=14, y=0.95)

# First plot: Cd vs Re
ax[0].plot(x, sample, "--", marker='o', label="Wang (2015)")
ax[0].plot(x, cd_mine, marker='o', label="In-House")
ax[0].plot(x_b, Breuer, marker='o', label="Breuer")
ax[0].legend()
ax[0].set_xlabel(r"$Re$")
ax[0].set_ylabel(r"$\overline{C_d}$")
ax[0].grid(True)
ax[0].set_box_aspect(1)

# Second plot: St vs Re
ax[1].plot(x, StSample, "--", marker='o', label="Wang (2015)")
ax[1].plot(x, St_mine, marker='o', label="In-House")  # Fixed label
ax[1].set_xlabel(r"$Re$")
ax[1].set_ylabel(r"$S_t$")
ax[1].grid(True)
ax[1].legend()
ax[1].set_box_aspect(1)

# Adjust layout
 # Leaves space for title
plt.subplots_adjust(wspace=0.3)  # Adjust spacing

plt.show()
