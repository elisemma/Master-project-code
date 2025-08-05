
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
import math
from foil_class_A0 import Foil_A0


# Load the data (assuming it's properly formatted)
data = np.loadtxt('./test_86Zr_data.Spe')
index = np.arange(len(data))

# Quadratic calibration parameters
a = 1.0205925340309977e-07   # Quadratic term
b = 0.4024344800692546        # Linear term
c = -0.39417117772666244      # Constant term

# Calibrated energy values based on index
calibrated_index = c + b * index + a * (index ** 2)

# Specify the calibrated energy range for linear fit
energy_min = 253
energy_max = 263

# Extract the indices where the calibrated energy falls within the specified range
mask = (calibrated_index >= energy_min) & (calibrated_index <= energy_max)
x_fit = calibrated_index[mask]
y_fit = data[mask]

# Perform linear regression
slope, intercept, r_value, p_value, std_err = linregress(x_fit, y_fit)

# Create linear fit line
fit_line = slope * x_fit + intercept

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(calibrated_index, data, label='Calibrated Data', color='blue', alpha=0.6)
plt.plot(x_fit, fit_line, label='Linear Fit', color='red', linestyle='--')

# Adding titles and labels
plt.title('Linear Fit from Calibrated Energy 253 to 263')
plt.xlabel('Energy (calibrated)')
plt.ylabel('Counts')
plt.axvline(x=energy_min, color='green', linestyle=':', label='Fit Range Start')
plt.axvline(x=energy_max, color='orange', linestyle=':', label='Fit Range End')
plt.legend()
plt.grid()
# plt.show()


mask_2kev = (calibrated_index >= 254) & (calibrated_index <= 256)
blue = np.trapezoid(data[mask_2kev], calibrated_index[mask_2kev])
red = np.trapezoid(slope * np.array([254, 256]) + intercept, calibrated_index[mask_2kev])

print(blue-red) #5 counts 

spectrum_counts = 4.8
eff = 0.005
intensity = 0.9584
true_counts = spectrum_counts*intensity/eff
activity=true_counts/1679
t = 172127
A_0 = activity*math.exp(np.log(2)/(16.5*3600)*t)
print(f'{A_0}Bq')

#Gpt: So, there are 172127 seconds between '02/12/2017 19:21:00' and '02/14/2017 19:09:47'.






def calc_xs_from_A0(foil_name, reaction_product): #fungerer ikke med denne klassen

    # A0_concat_df = pd.read_csv(f'./Calculated_A0/{foil_name}_A0_by_hand.csv')


    # A0_filtered_df = A0_concat_df[A0_concat_df['Isotope']==reaction_product]
    # A0 = A0_filtered_df['A0'].iloc[0]
    # A0_unc = A0_filtered_df['A0_unc'].iloc[0]
    A0 = 4
    A0_unc = 4

    foil = Foil_A0(foil_name, reaction_product, A0, A0_unc)
    foil.assign_areal_dens_w_unc_percent()
    foil.assign_molar_mass()
    foil.calculate_decay_constant_w_unc()
    foil.assign_beam_current_w_unc()
    foil.calculate_xs_w_unc_old()

    xs = foil.calc_xs
    xs_unc = foil.calc_xs_unc
    return xs, xs_unc #[mb]



xs = calc_xs_from_A0("Zr06", "86ZR")
print(xs)