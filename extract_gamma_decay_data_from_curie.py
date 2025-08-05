import curie as ci 
import pandas as pd



def make_gamma_decay_DataFrame(isotope, isotopeLatexString, unit_t_half, intensityLim=1.0):
    ip = ci.Isotope(isotope)                                                                     # Defining isotope for curie
    t_half = ip.half_life(unit_t_half)                                                           # Get half life from curie
    gammaInfo = ip.gammas(I_lim=intensityLim)                                                    # Get gamma decay info from curie
    gammaInfo.drop(columns=['unc_intensity'], inplace=True)                                      # Deleting gamma unc col
    gammaInfo['Nuclide'] = [""] * len(gammaInfo)                                                 # Init col for nuclide with empty strings
    gammaInfo['Half-life'] = f"{t_half} {unit_t_half}"                                           # Fill entire col with half-life
    gammaInfo.loc[0, 'Nuclide'] = isotopeLatexString                                             # Set the nuclide only for the first row
    gammaInfo.rename(columns={'energy': r'$E_{\gamma}$ (keV)'}, inplace=True)                     # Renaming col 
    gammaInfo.rename(columns={'intensity': r'$I{\gamma}$ ($\%$)'}, inplace=True)                  # Renaming col 
    gammaInfo = gammaInfo[['Nuclide', 'Half-life', r'$E_{\gamma}$ (keV)', r'$I{\gamma}$ ($\%$)']]  # Reorder columns to match the desired output
    return gammaInfo 



# Make latex code for Zr: __________________________________________________________________________________________________

gammaInfo_86Y = make_gamma_decay_DataFrame('86Y', r'$^{86}$Y', unit_t_half = 'h')
gammaInfo_86mY = make_gamma_decay_DataFrame('86Ym', r'$^{86\text{m}}$Y', unit_t_half = 'm')
gammaInfo87Y = make_gamma_decay_DataFrame('87Y', r'$^{87}$Y', unit_t_half = 'h')
gammaInfo_87mY = make_gamma_decay_DataFrame('87Ym', r'$^{87\text{m}}$Y', unit_t_half = 'h')
gammaInfo88Y = make_gamma_decay_DataFrame('88Y', r'$^{88}$Y', unit_t_half = 'd')
gammaInfo88Zr = make_gamma_decay_DataFrame('88ZR', r'$^{88}$Zr', unit_t_half = 'd')
gammaInfo88Nb = make_gamma_decay_DataFrame('88NB', r'$^{88}$Nb', unit_t_half = 'm')
gammaInfo89Zr = make_gamma_decay_DataFrame('89ZR', r'$^{89}$Zr', unit_t_half = 'h')
gammaInfo90Nb = make_gamma_decay_DataFrame('90NB', r'$^{90}$Nb', unit_t_half = 'h')
gammaInfo_90mY = make_gamma_decay_DataFrame('90Ym', r'$^{90\text{m}}$Y', unit_t_half = 'h')
gammaInfo_92mNb = make_gamma_decay_DataFrame('92NBm', r'$^{92\text{m}}$Nb', unit_t_half = 'd')
gammaInfo95Zr = make_gamma_decay_DataFrame('95ZR', r'$^{95}$Zr', unit_t_half = 'd')
gammaInfo95Nb = make_gamma_decay_DataFrame('95NB', r'$^{95}$Nb', unit_t_half = 'd')
gammaInfo_95mNb = make_gamma_decay_DataFrame('95NBm', r'$^{95\text{m}}$Nb', unit_t_half = 'd')
gammaInfo96Nb = make_gamma_decay_DataFrame('96NB', r'$^{96}$Nb', unit_t_half = 'h')


# Concatenate DataFrames containing gamma decay info for several nuclides
gammaInfo_allZr = pd.concat([gammaInfo_86Y, gammaInfo_86mY, gammaInfo87Y, gammaInfo_87mY, gammaInfo88Y, gammaInfo88Zr, gammaInfo88Nb, gammaInfo89Zr, gammaInfo90Nb, gammaInfo_90mY, gammaInfo_92mNb, gammaInfo95Zr, gammaInfo95Nb, gammaInfo_95mNb, gammaInfo96Nb], ignore_index=True)
print(gammaInfo_allZr)

# Convert data fram into latex code 
latexcode_allZr = gammaInfo_allZr.to_latex(index=False, escape=False, float_format="%.2f")
print(latexcode_allZr)






# Make latex code for Ni, Ti og Fe: __________________________________________________________________________________________________

gammaInfo_46Sc = make_gamma_decay_DataFrame('46SC', r'$^{46}$SC', unit_t_half = 'd')
gammaInfo_47Sc = make_gamma_decay_DataFrame('47SC', r'$^{47}$SC', unit_t_half = 'd')
gammaInfo_48Sc = make_gamma_decay_DataFrame('48SC', r'$^{48}$SC', unit_t_half = 'h')
gammaInfo_48V  = make_gamma_decay_DataFrame('48V',  r'$^{48}$V',  unit_t_half = 'd')
gammaInfo_52Mn = make_gamma_decay_DataFrame('52MN', r'$^{52}$Mn', unit_t_half = 'd')
gammaInfo_54Mn = make_gamma_decay_DataFrame('54MN', r'$^{54}$Mn', unit_t_half = 'd')
gammaInfo_55Co = make_gamma_decay_DataFrame('55CO', r'$^{55}$Co', unit_t_half = 'h')
gammaInfo_56Co = make_gamma_decay_DataFrame('56CO', r'$^{56}$Co', unit_t_half = 'd')
gammaInfo_57Co = make_gamma_decay_DataFrame('57CO', r'$^{57}$Co', unit_t_half = 'd')
gammaInfo_57Ni = make_gamma_decay_DataFrame('57NI', r'$^{57}$Ni', unit_t_half = 'h')
gammaInfo_58Co = make_gamma_decay_DataFrame('58CO', r'$^{58}$Co', unit_t_half = 'd')
gammaInfo_60Co = make_gamma_decay_DataFrame('60CO', r'$^{60}$Co', unit_t_half = 'y')
gammaInfo_61Cu = make_gamma_decay_DataFrame('61CU', r'$^{61}$Cu', unit_t_half = 'h')
gammaInfo_65Ni = make_gamma_decay_DataFrame('65NI', r'$^{65}$Ni', unit_t_half = 'h')


# Concatenate DataFrames containing gamma decay info for several nuclides
gammaInfo_allNiTiFe = pd.concat([gammaInfo_46Sc, gammaInfo_47Sc, gammaInfo_48Sc, gammaInfo_48V, gammaInfo_52Mn, gammaInfo_54Mn, gammaInfo_55Co, gammaInfo_56Co, gammaInfo_57Co, gammaInfo_57Ni, gammaInfo_58Co, gammaInfo_60Co, gammaInfo_61Cu, gammaInfo_65Ni], ignore_index=True)
print(gammaInfo_allNiTiFe)

# Convert data fram into latex code 
latexcode_allNiTiFe = gammaInfo_allNiTiFe.to_latex(index=False, escape=False, float_format="%.2f")
print(latexcode_allNiTiFe)