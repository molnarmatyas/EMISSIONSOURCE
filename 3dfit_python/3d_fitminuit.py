from scipy.stats import chi2
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Set the backend to 'Agg' for non-GUI environments
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from iminuit import Minuit

# Define the fit function
def fit_function(k, lambd, R_o, R_s, R_l, alpha):
    kx, ky, kz = k
    factor = 1. / (197. ** 2.)
    exponent = -(factor * ((R_o ** 2.) * kx**2. + (R_s ** 2.) * ky**2. + (R_l ** 2.) * kz**2.)) ** (alpha / 2.)
    return 1. + lambd * np.exp(exponent)

# Load the data from the file
data = np.loadtxt("results/q_0275_0325_3d_100M.dat", comments='!', usecols=(0, 1, 2, 3, 4))

# Extract columns from the data
kx_all, ky_all, kz_all = data[:, 0], data[:, 1], data[:, 2]
Ck_all, Ck_err_all = data[:, 3], data[:, 4]

# Filter out invalid rows
valid = (Ck_err_all > 0) & (Ck_all > 0)
kx, ky, kz, Ck, Ck_err = kx_all[valid], ky_all[valid], kz_all[valid], Ck_all[valid], Ck_err_all[valid]

#Fitting boundaries
cube_cut = False
cube = [0,80]

sphere_cut = False
sphere = [0,np.sqrt(3)*80]

edge_cut = False
edge = 0

if cube_cut:
    incut = (kx > cube[0]) & (kx < cube[1]) & (ky > cube[0]) & (ky < cube[1]) & (kz > cube[0]) & (kz < cube[1])
    kx, ky, kz, Ck, Ck_err = kx[incut], ky[incut], kz[incut], Ck[incut], Ck_err[incut]

if sphere_cut == True:
  incut = (np.sqrt(kx**2 + ky**2 + kz**2) > sphere[0]) & (np.sqrt(kx**2 + ky**2 + kz**2) < sphere[1])
  kx, ky, kz, Ck, Ck_err = kx[incut], ky[incut], kz[incut], Ck[incut], Ck_err[incut]

if edge_cut == True:
  outcut = ((kx < edge)&(ky < edge)) | ((kx < edge)&(kz < edge)) | ((ky < edge)&(kz < edge))
  incut = ~outcut
  kx, ky, kz, Ck, Ck_err = kx[incut], ky[incut], kz[incut], Ck[incut], Ck_err[incut]

# Stack kx, ky, kz into a 2D array to pass as a single variable
k_data = np.vstack((kx, ky, kz))
k_data_all = np.vstack((kx_all, ky_all, kz_all))

# FILTER C(k)
# Compute k magnitude (assuming k = [kx, ky, kz])
k_tofit = np.array([kx,ky,kz])
print("Are the same?")
print(k_tofit[0] == k_data[0])
#Ck_tofit = 
k_magnitude = np.sqrt(k_tofit[0]**2 + k_tofit[1]**2 + k_tofit[2]**2)

# Apply the filter: select indices where k < 65 MeV/c
mask = (k_magnitude < 65.) * (k_magnitude > 7.5) # 65, 5
filtered_k_data = k_data[:, mask]
filtered_Ck = Ck[mask]
filtered_Ck_err = Ck_err[mask]
print("End of Ck: ",filtered_Ck[-1])


# Define chi-squared
def chi_squared(lambd, R_o, R_s, R_l, alpha):
    model = fit_function(filtered_k_data, lambd, R_o, R_s, R_l, alpha)
    residuals = (filtered_Ck - model) / filtered_Ck_err  # Ck_err is the error on Ck_all
    return np.sum(residuals**2)

assert filtered_k_data.shape[1] == filtered_Ck.shape[0] == filtered_Ck_err.shape[0], \
    "Mismatch in filtered data shapes!"
print("Shape of filtered_k_data", filtered_k_data.shape)
print("Shape of filtered_Ck", filtered_Ck.shape)

# Perform the fit with Minuit  ---------------------
print("Starting Minuit fit")
# Initial parameters
initial_params = [0.4, 7.0, 7.0, 7.0, 1.6]
param_names = ["lambda", "R_o", "R_s", "R_l", "alpha"]

# Set up Minuit
print("Setting up Minuit")
m = Minuit(chi_squared, lambd=initial_params[0], R_o=initial_params[1], 
           R_s=initial_params[2], R_l=initial_params[3], alpha=initial_params[4])

# Optionally set parameter bounds or fix parameters
#print("Setting up parameter bounds")
#m.limits = {
#    "lambd": (0, 1),
#    "R_o": (0, None),
#    "R_s": (0, None),
#    "R_l": (0, None),
#    "alpha": (0, 2)
#}

# Perform minimization
print("Starting minimization...")
#m.migrad(ncall=None)  # Minimization step
m.migrad()
print("Minimization step done.")
# Minos error estimation
print("Starting error estimation...")
m.minos()
print("Error estimation with minos done.")

# Chi-squared per degree of freedom
ndof = len(Ck_all) - len(m.parameters)  # Degrees of freedom
print(f"Chi2/Ndof: {m.fval / ndof}")
# ---------------------------------------------------

# Print optimized parameters and their errors
print("\nOptimized parameters and errors:")
for name, value, error in zip(param_names, m.values, m.errors):
    print(f"{name} = {value:.4f} ± {error:.4f}")
print("Asymmetric errors from Minos:")
for param in m.parameters:
    err = m.merrors[param]
    print(f"{param}: -{err.lower} +{err.upper}")

# Covariance matrix
cov_matrix = m.covariance
print("Covariance matrix:")
print(cov_matrix)

# ************ PLOTTING *************

# Calculate the fitted values for plotting
lambd, R_o, R_s, R_l, alpha = m.values[0], m.values[1], m.values[2], m.values[3], m.values[4]
Ck_fit = fit_function(k_data_all, lambd, R_o, R_s, R_l, alpha)
#diagonal=(kx_all == ky_all) & (kx_all == kz_all)
#Ck_fitted_diagonal = Ck_fit[diagonal]

#Calculate average diagonal Ck data
bins_all = round(len(kx_all)**(1/3))
k_diagonal = kx_all[::bins_all**2]
Ck_average_diagonal = np.zeros(bins_all)
Ck_err_average_diagonal = np.zeros(bins_all)
Ck_fit_average_diagonal = np.zeros(bins_all)

for i in range(1,1+bins_all-2):
  n=0
  Ck_average=0
  Ck_err_average2=0
  Ck_fit_average=0
  for a in [-1,0,1]:
    for b in [-1,0,1]:
      for c in [-1,0,1]:
        Ck_here = Ck_all[(i+a)*bins_all**2 + (i+b)*bins_all + (i+c)]
        Ck_err_here = Ck_err_all[(i+a)*bins_all**2 + (i+b)*bins_all + (i+c)]
        Ck_fit_here = Ck_fit[(i+a)*bins_all**2 + (i+b)*bins_all + (i+c)]
        Ck_fit_average += Ck_fit_here
        if (Ck_here > 0) & (Ck_err_here > 0):
          n += 1
          Ck_average += Ck_here
          Ck_err_average2 += Ck_err_here**2
  if n != 0:
    Ck_average_diagonal[i] = Ck_average/n
    Ck_err_average_diagonal[i] = np.sqrt(Ck_err_average2)/n
  Ck_fit_average_diagonal[i] = Ck_fit_average/27

if (Ck_all[0] > 0) & (Ck_err_all[0] > 0):
  Ck_average_diagonal[0] = Ck_all[0]
  Ck_err_average_diagonal[0] = Ck_err_all[0]

if (Ck_all[-1] > 0) & (Ck_err_all[-1] > 0):
  Ck_average_diagonal[-1] = Ck_all[-1]
  Ck_err_average_diagonal[-1] = Ck_err_all[-1]

Ck_fit_average_diagonal[0] = Ck_fit[0]
Ck_fit_average_diagonal[-1] = Ck_fit[-1]

# Plot the data and fit
plt.errorbar(k_diagonal, Ck_average_diagonal, yerr=Ck_err_average_diagonal, fmt='o', label="Data")
plt.plot(k_diagonal, Ck_fit_average_diagonal, 'r-', label=f"Fit: $\lambda=${m.values[0]:.3f}")
plt.ylim(0, 2)
plt.xlabel("k (MeV/c)")
plt.ylabel("C(k)")
plt.legend()
plt.title("3D Fit of C(k)")
plt.savefig("result_plots/q_0275_0325_3d_100M.png")  # Save plot to file

