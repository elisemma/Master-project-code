import numpy as np
import matplotlib.pyplot as plt

# Parameters for Gaussian peak
peak_mean = 1000
peak_std = 0.8

# Generate Gaussian data
num_gaussian_points = 5000
gaussian_data = np.random.normal(peak_mean, peak_std, num_gaussian_points)

# Generate linear background data
num_linear_points = 1000
linear_data = np.random.uniform(low=min(gaussian_data), high=max(gaussian_data), size=num_linear_points)
# linear_data = linear_data * background_slope * np.random.uniform(0.0, 2000.0, num_linear_points)

# Combine Gaussian and linear data
combined_data = np.hstack((gaussian_data, linear_data))

# Create histogram
num_bins = 100
histogram_values, bin_edges = np.histogram(combined_data, bins=num_bins, density=False)

# Plotting the histogram for visualization (optional)
plt.hist(combined_data, bins=num_bins, alpha=0.75, label='Gaussian with Linear Background')
plt.legend()
plt.title('Gaussian Peak with Linear Background')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.grid(True)
plt.show()

# Saving histogram data to a text file
out_filename = 'histogram_data.txt'
hist_data_with_edges = np.vstack((bin_edges[:-1], histogram_values)).T
np.savetxt(out_filename, hist_data_with_edges, comments='', fmt='%10.5f')
