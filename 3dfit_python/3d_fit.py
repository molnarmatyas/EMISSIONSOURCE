from scipy.stats import chi2
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Set the backend to 'Agg' for non-GUI environments
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# Define the fit function
def fit_function(k, lambd, R_o, R_s, R_l, alpha):
    kx, ky, kz = k
    factor = 1. / (197. ** 2.)
    exponent = -(factor * ((R_o ** 2.) * kx**2. + (R_s ** 2.) * ky**2. + (R_l ** 2.) * kz**2.)) ** (alpha / 2.)
    return 1. + lambd * np.exp(exponent)

# Load the data from the file
data = np.loadtxt("results/q_0275_0325_3d_30M.dat", comments='!', usecols=(0, 1, 2, 3, 4))

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

if cube_cut == True:
  incut = (kx > cube[0]) & (kx < cube[1]) & (kx > cube[0]) & (kx < cube[1]) & (kx > cube[0]) & (kx < cube[1])
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


# Initial parameter guess for lambda, R_o, R_s, R_l, and alpha
initial_params = [0.8, 7.0, 7.0, 7.0, 1.6]  # [lambda, R_o, R_s, R_l, alpha]

# Define the chi-squared minimization function
def chi_squared(params, k_data, Ck, Ck_err):
    model = fit_function(k_data,*params)
    chi2 = np.sum(((Ck - model) / Ck_err) ** 2.)  # chi-squared
    return chi2

# Define a callback function to print chi-squared values at each step
def callback(params):
    chi2 = chi_squared(params, k_data, Ck, Ck_err)
    print(f"Current parameters: {params}, Chi-squared: {chi2}")

# Perform the fit using minimize with bounds (you can adjust bounds as needed)
result = minimize(
    chi_squared, initial_params, args=(k_data, Ck, Ck_err),
    bounds=[(0, 2), (1, 20), (1, 20), (1, 20), (0.5, 2)], method='L-BFGS-B', callback=callback, options={'disp': True}
)

# Extract fit parameters from result
lambd, R_o, R_s, R_l, alpha = result.x

# Calculate the covariance matrix (approximation using Hessian inverse)
try:
    cov_matrix = result.hess_inv.todense()  # For L-BFGS-B, hess_inv is a scipy.sparse.linalg.LinAlgError
    param_errors = np.sqrt(np.diag(cov_matrix))  # Standard errors are sqrt of diagonal elements
except AttributeError:
    param_errors = [np.nan] * len(result.x)  # If hess_inv not available, set errors to NaN

#Calculate confidence level
#chi2 = chi_squared(result.x, k_data, Ck, Ck_err)
chi2_min = result.fun
NDF = len(Ck) - len(result.x)
confidence = 1 - chi2.cdf(chi2_min, df=NDF)

# Print optimized parameters and their errors
print("\nOptimized parameters and errors:")
param_names = ["lambda", "R_o", "R_s", "R_l", "alpha"]
for name, value, error in zip(param_names, result.x, param_errors):
    print(f"{name} = {value:.4f} ± {error:.4f}")
print(f"Chi^2/NDF = {chi2_min:.4f} / {NDF}")
print(f"Confidence level = {confidence:.4f}")

# Calculate the fitted values for plotting
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
plt.plot(k_diagonal, Ck_fit_average_diagonal, 'r-', label="Fit")
plt.ylim(0, 2)
plt.xlabel("k (MeV/c)")
plt.ylabel("C(k)")
plt.legend()
plt.title("3D Fit of C(k)")
plt.savefig("result_plots/q_0275_0325_3d_30M.png")  # Save plot to file

