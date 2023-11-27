from pyexcel.cookbook import merge_all_to_a_book
# import pyexcel.ext.xlsx # no longer required if you use pyexcel >= 0.2.2 
import glob


merge_all_to_a_book(glob.glob("/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/Master-project-code/MyGeneratedFiles/Ni_foils/AY130217_Ni02_18cm_30MeV/*.csv"), "output.xlsx")