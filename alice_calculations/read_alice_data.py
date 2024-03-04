import numpy as np
import matplotlib.pyplot as plt

# def read_data(file_path, Z, A):
#     # Initialize lists to store energy and cross-section data
#     energies = []
#     cross_sections = []

#     # Open the file and read line by line
#     with open(file_path, 'r') as file:
#         # Skip lines until reaching the cross-section data
#         for line in file:
#             if 'Ebeam' in line:
#                 break
#         # Read cross-section data for the specified isotope
#         for line in file:
#             if line.strip():  # Check if line is not empty
#                 data = line.split()
#                 if int(data[1]) == Z and int(data[2]) == A:
#                     energies.append(float(data[0]))
#                     cross_sections.append(float(data[3]))
    
#     return np.array(energies), np.array(cross_sections)

def read_data(file_path, Z, A):
    # Initialize lists to store energy and cross-section data
    energies = []
    cross_sections = []

    # Open the file and read line by line
    with open(file_path, 'r') as file:
        # Skip lines until reaching the cross-section data
        for line in file:
            if 'Ebeam' in line:
                break
        # Read cross-section data for the specified isotope
        for line in file:
            if line.strip():  # Check if line is not empty
                data = line.split()
                try:
                    if int(float(data[1])) == Z and int(float(data[2])) == A:
                        energies.append(float(data[0]))
                        cross_sections.append(float(data[3]))
                except ValueError:
                    pass  # Skip lines that can't be converted to integers
    
    return np.array(energies), np.array(cross_sections)


def plot_cross_section(energies, cross_sections, Z, A):
    plt.figure(figsize=(10, 6))
    plt.plot(energies, cross_sections, label=f'Z={Z}, A={A} Cross-section')
    plt.xlabel('Energy (MeV)')
    plt.ylabel('Cross-section (mb)')
    plt.title(f'Cross-section as a function of Energy for Z={Z}, A={A}')
    plt.legend()
    plt.grid(True)
    plt.show()




E, xs = read_data('./plot_natZr_dx', 41, 90)
plot_cross_section(E, xs, 41, 90)