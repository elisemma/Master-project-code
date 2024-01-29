import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.special import gammaln
from scipy.integrate import quad
from scipy.linalg import inv
from scipy.optimize import approx_fprime
from scipy.stats import multivariate_normal
import scipy.sparse.linalg

# Read histogram data from text file
bin_centers, counts = np.loadtxt('histogram_data.txt', unpack=True)

# Define Gaussian distribution function -- normalized to 1 over entire range

def gaussian(x, mu, sigma, amp):
    gaus = amp * np.exp(-(x - mu)**2 / (2 * sigma**2))
    return gaus

def gaussianAndPolynomial(x, p0, p1, mu, sigma, amp):
    polynomial = p0 + p1*x
    # gaussian = amp / (sigma * np.sqrt(2 * np.pi)) * np.exp(-(x - mu)**2 / (2 * sigma**2))
    gaus = gaussian(x, mu, sigma, amp)
    return polynomial + gaus

# Define negative log of the Poisson likelihood
def neg_log_poisson_likelihood(params, x, y):
    p0, p1, mu, sigma, amp = params
    # Calculate expected counts based on fitted Gaussian
    lambdas = gaussianAndPolynomial(x, p0, p1, mu, sigma, amp)
    # Add a small constant to lambdas before taking the log
    epsilon = 1e-8
    # Negative log of the Poisson pmf; we use gammaln to compute the log-factorial in a stable way
    return np.sum(lambdas - y * np.log(lambdas + epsilon) + gammaln(y + 1))


# Initial guess for the parameters
initial_params = [0.0, 0.0, 1000.0, 1.0, 100.0]

bnds = ((0.0, 1.0e+02), (0.0, 0.0), (998.0, 1002.0), (0, 2.0), (0.0, 1.0e+03))
# bnds = ((0.0, 0.0), (0.0, 0.0), (999.5, 1000.5), (0, 3.0), (0.0, 2.0e+05))

# Minimize the negative log-likelihood
# result = minimize(neg_log_poisson_likelihood, initial_params, args=(bin_centers, counts), bounds=bnds, method='BFGS')

# method='BFGS' to get the Hessian matrix (used to determine the covariance matrix) cannot be used with the bounds option
result = minimize(neg_log_poisson_likelihood, initial_params, args=(bin_centers, counts), method='BFGS')

# Extract the maximum likelihood estimates (MLE)
p0_mle, p1_mle, mu_mle, sigma_mle, amp_mle = result.x
# print(f"Maximum Likelihood Estimates: mu = {mu_mle}, sigma = {sigma_mle}")
print(f"Maximum Likelihood Estimates:")
print(f"p0_mle = {p0_mle}")
print(f"p1_mle = {p1_mle}")
print(f"mu_mle = {mu_mle}")
print(f"sigma_mle = {sigma_mle}")
print(f"amp_mle = {amp_mle}")

# Generate bin points for the fitted Gaussian pdf
x_values = np.linspace(min(bin_centers), max(bin_centers), 10000)
# Multiply by total count and bin width to convert to same scale as histogram count data
y_totalFit = gaussianAndPolynomial(x_values, p0_mle, p1_mle, mu_mle, sigma_mle, amp_mle)
y_gaussian = gaussian(x_values, mu_mle, sigma_mle, amp_mle)

binWidth = bin_centers[1] - bin_centers[0]

# Use the parameters from the Maximum Likelihood Estimates
integral, integralError_numericalIntegration = quad(gaussian,
                            900.0, # Lower bound of integration
                            1100.0,  # Upper bound of integration
                            args=(mu_mle, sigma_mle, amp_mle) # MLE parameters
                           )

integral *= (1/binWidth)
integralError_numericalIntegration *= (1/binWidth) # this is the numerical integration error



# Use this function to calculate the Hessian matrix at the MLE
# H = hessian_finite_difference(result.x, lambda x: neg_log_poisson_likelihood(x, bin_centers, counts))
# Then invert the Hessian to get the covariance matrix
# cov_matrix = inv(H)
# print("Covariance matrix:")
# print(cov_matrix)

# Check the type of hess_inv and handle accordingly:
hessian_inv_object = result.hess_inv

# print("hessian_inv_object:")
# print(hessian_inv_object)


# If the inverse Hessian is represented as a LinearOperator, its dense form can be constructed as follows:
if isinstance(hessian_inv_object, scipy.sparse.linalg.LinearOperator):
    identity_matrix = np.eye(hessian_inv_object.shape[0])
    hessian_inv_dense = hessian_inv_object.matmat(identity_matrix)

    # Inverting the inverse Hessian gives the Hessian matrix itself, which is then used to compute the covariance matrix
    covariance_matrix = hessian_inv_dense
else:
    # If hess_inv is already in the desired dense format, no further action needed:
    covariance_matrix = hessian_inv_object


# print("Covariance matrix:")
# print(covariance_matrix)


# Define the number of samples for the Monte Carlo simulation
num_samples = 10000

# Sample the parameter space using the multivariate normal distribution
samples = multivariate_normal.rvs(mean=[p0_mle, p1_mle, mu_mle, sigma_mle, amp_mle], cov=covariance_matrix, size=num_samples)

# Define a function that integrates the Gaussian over the range of interest
def integrate_gaussian(params):
    p0, p1, mu, sigma, amp = params
    integral, _ = quad(gaussian, 900, 1100, args=(mu, sigma, amp))
    integral *= (1 / binWidth)  # Adjust for bin width
    return integral

# Compute the integral for each sample
integrals = np.array([integrate_gaussian(sample) for sample in samples])

# Calculate the standard deviation of these integrals, which is the error of the integral from fit uncertainties
integralUncertainty_fitUncertainties = np.std(integrals)

totalIntegralError = np.sqrt(pow(np.sqrt(integral), 2.0) + pow(integralError_numericalIntegration, 2.0) + pow(integralUncertainty_fitUncertainties, 2.0))

print(f"Integration of gaussian function: {integral} ({totalIntegralError})")

# # Generate a set of mu and sigma values drawn from normal distributions
# # centered on the MLE values and with standard deviations equal to the errors
# n_samples = 10000
# mu_samples = np.random.normal(mu_mle, mu_error, n_samples)
# sigma_samples = np.random.normal(sigma_mle, sigma_error, n_samples)

# # Now for each mu and sigma value, compute the integral and store the results
# integral_samples = np.zeros(n_samples)
# for i in range(n_samples):
#     integral_samples[i], _ = quad(gaussian,
#                                   -np.inf,
#                                   np.inf,
#                                   args=(mu_samples[i], sigma_samples[i]))

# print(f"Integration of gaussian function: {integral} ({totalIntegralError})")

# Plot histogram data and the fitted Gaussian
plt.bar(bin_centers, counts, color='gray', alpha=0.5, label='Data', width=bin_centers[1]-bin_centers[0])
plt.plot(x_values, y_totalFit, color='red', label='Total fit')
plt.plot(x_values, y_gaussian, color='blue', label='Fitted Gaussian')
plt.xlabel('Bin Centers')
plt.ylabel('Counts')
plt.legend()
plt.show()
