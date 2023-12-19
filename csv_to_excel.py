from pyexcel.cookbook import merge_all_to_a_book
# import pyexcel.ext.xlsx # no longer required if you use pyexcel >= 0.2.2 
import glob


merge_all_to_a_book(glob.glob("./MyGeneratedFiles/combined_peak_data/*.csv"), "./MyGeneratedFiles/combined_peak_data/combined_peak_data.xlsx")

