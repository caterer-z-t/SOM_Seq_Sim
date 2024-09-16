import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde
from scipy.stats import norm

# Generate sample data (Normal distribution in this case)
data = np.random.normal(size=1000)

# Create histogram of the data
hist, bin_edges = np.histogram(data, bins=50, density=True)  # Normalize histogram
x_values = (bin_edges[:-1] + bin_edges[1:]) / 2  # Get midpoints of bins

# Calculate mean and standard deviation of the data
mean, std_dev = np.mean(data), np.std(data)

# Generate points for the normal distribution curve
x_curve = np.linspace(min(data), max(data), 1000)
normal_curve = norm.pdf(x_curve, mean, std_dev)

# KDE for the data
kde = gaussian_kde(data)
x_kde = np.linspace(min(data), max(data), 1000)
kde_vals = kde(x_kde)

# Detect peaks in the histogram
peaks_hist, _ = find_peaks(hist, height=0.02)  # Adjust height threshold as needed

# Detect peaks in the KDE
peaks_kde, _ = find_peaks(kde_vals, height=0.02)  # Adjust height threshold as needed

# Create subplots side by side
fig, axs = plt.subplots(1, 3, figsize=(15, 5))

# Plot the histogram and peaks
axs[0].plot(x_values, hist, label="Histogram")
axs[0].plot(x_values[peaks_hist], hist[peaks_hist], "x", label="Peaks", color="red")
axs[0].set_title("Histogram with Peaks")
axs[0].set_xlabel("Data Value")
axs[0].set_ylabel("Density")
axs[0].legend()

# Plot the KDE and peaks
axs[1].plot(x_kde, kde_vals, label="KDE")
axs[1].plot(x_kde[peaks_kde], kde_vals[peaks_kde], "x", label="Peaks", color="red")
axs[1].set_title("KDE with Peaks")
axs[1].set_xlabel("Data Value")
axs[1].set_ylabel("Density")
axs[1].legend()

# Plot the original data with normal distribution curve
axs[2].hist(
    data, bins=30, density=True, alpha=0.5, color="skyblue", label="Data Histogram"
)
axs[2].plot(x_curve, normal_curve, label="Normal Distribution", color="red")
axs[2].set_title("Data Histogram with Normal Distribution")
axs[2].set_xlabel("Data Value")
axs[2].set_ylabel("Density")
axs[2].legend()

# Show the plots
plt.tight_layout()
plt.show()

# Output the number of peaks detected
print(f"Number of peaks detected in Histogram: {len(peaks_hist)}")
print(f"Number of peaks detected in KDE: {len(peaks_kde)}")
