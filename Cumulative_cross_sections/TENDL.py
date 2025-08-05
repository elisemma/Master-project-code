# from tools import *
import numpy as np
import requests
import matplotlib.pyplot as plt

# class Tendl:

#     def __init__(self, target):
#         # target = {"Ir191": 0.373, "Ir193": 0.627}
#         self.target = target

#     def tendlDeuteronData(self, productZ, productA, isomerLevel = None):
#         targetFoil = list(self.target.keys())[0][0:2]
#         product = self.product(productZ, productA)
#         fileEnding = self.tendlFileEnding(isomerLevel)
#         E = []
#         Cs = []
#         for t in self.target.keys():
#             data = self.retrieveTendlDataFromUrl(
#                 self.tendDeuteronlUrl(targetFoil, t, product, fileEnding), t
#             )
#             E.append(data[0])
#             Cs.append(data[1])
#         CsSummed = sum(Cs)
#         E, Cs = Tools().interpolate(E[0], CsSummed)
#         return E, Cs

#     def plotTendl23(self, productZ, productA, isomerLevel = None, betaFeeding = None, branchingRatio = None, parentIsomerLevel = None):
#         E, Cs = self.tendlDeuteronData(productZ, productA, isomerLevel)
#         if betaFeeding:
#             CsParent = self.correctForBetaFeeding(productZ, productA, betaFeeding, branchingRatio, parentIsomerLevel)[0]
#             Cs = Cs + CsParent
#         plt.plot(E, Cs, label='TENDL-2023', linestyle='--', color='blue')

#     def correctForBetaFeeding(self, productZ, productA, betaFeeding, branchingRatio, parentIsomerLevel):
#         if (betaFeeding  == 'beta+'):
#             parentZ = str(int(productZ)+1); parentA = productA
#         elif (betaFeeding == 'beta-'):
#             parentZ = str(int(productZ)-1); parentA = productA
#         E, Cs = self.tendlDeuteronData(parentZ, parentA, parentIsomerLevel)
#         return E, Cs*branchingRatio

#     def product(self, productZ, productA):
#         if len(productZ) <= 2:
#             productZ = '0' + productZ
#         else:
#             productZ = productZ
#         if len(productA) <= 2:
#             productA = '0' + productA
#         else:
#             productA = productA
#         return productZ + productA

#     def tendDeuteronlUrl(self, targetFoil, target, product, fileEnding):
#         return ('https://tendl.web.psi.ch/tendl_2023/deuteron_file/'
#         + targetFoil + '/' + target
#         + '/tables/residual/rp'
#         + product + fileEnding)

#     def tendlFileEnding(self, isomerLevel=None):
#         return '.tot' if isomerLevel==None else '.L'+ isomerLevel

#     def retrieveTendlDataFromUrl(self, url, target):
#         try:
#             tendlData = requests.get(url).text.split("\n")[27:] # skipping 27 first lines in tendl file
#         except:
#             raise Exception("Unable to retrieve tendl data from url: " + url)
#         tendlData = np.genfromtxt(tendlData)
#         abundance = self.target[target]
#         E = tendlData[:,0]
#         Cs = tendlData[:,1]
#         return E, Cs*abundance

# tendl = Tendl({"Ir191": 0.373, "Ir193": 0.627})
# # tendl.plotTendl23( '78', '188', betaFeeding = 'beta+', branchingRatio=1.0)
# tendl.plotTendl23( '78', '188', 'Ir_188Ir', betaFeeding = 'beta+', branchingRatio=1.0)







def tendDeuteronlUrl(targetFoil, target, product, fileEnding):
    return ('https://tendl.web.psi.ch/tendl_2023/deuteron_file/'
        + targetFoil + '/' + target
        + '/tables/residual/rp'
        + product + fileEnding)


def tendlFileEnding(isomerLevel=None):
    return '.tot' if isomerLevel==None else '.L'+ isomerLevel


def retrieveTendlDataFromUrl(url, target):
    try:
        tendlData = requests.get(url).text.split("\n")[27:] # skipping 27 first lines in tendl file
    except:
        raise Exception("Unable to retrieve tendl data from url: " + url)
    tendlData = np.genfromtxt(tendlData)
    # abundance = target[target]
    E = tendlData[:,0]
    Cs = tendlData[:,1]
    return E, Cs#*abundance


def run_tendl(targetFoil, target, product, isomer=None):
    fileEnding = tendlFileEnding(isomer)
    url = tendDeuteronlUrl(targetFoil, target, product, fileEnding)
    E, xs = retrieveTendlDataFromUrl(url, target)
    return E, xs 


E, xs = run_tendl('Zr', 'Zr090', '039086')
print('E', E)
print('xs', xs)
