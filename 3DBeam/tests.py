# from math import log10
import source.utilities.utilities as utils
#from scipy.stats import rayleigh
#import source.utilities.fatigue_utilities as fat
# from inputs import holz
import numpy as np
import inputs.DIN_Windlasten as windDIN
import inputs.wind_modelle_IEC as windIEC
# import pandas as pd
# from os.path import join as os_join
import matplotlib.pyplot as plt

import numpy as np

import tkinter as tk
from tkinter import ttk

def get_values_page1():
    height = float(entry_height.get())
    diameter1 = float(entry_diameter1.get())
    diameter2 = float(entry_diameter2.get())
    print("Height:", height)
    print("Diameter 1:", diameter1)
    print("Diameter 2:", diameter2)

def get_values_page2():
    material = entry_material.get()
    layer_structure = entry_layer_structure.get()
    print("Material:", material)
    print("Layer Structure:", layer_structure)

root = tk.Tk()
root.title("Parameter Input")

notebook = ttk.Notebook(root)

page1 = ttk.Frame(notebook)
page2 = ttk.Frame(notebook)

notebook.add(page1, text="Page 1")
notebook.add(page2, text="Page 2")

label_height = tk.Label(page1, text="Höhe:")
label_diameter1 = tk.Label(page1, text="Durchmesser 1:")
label_diameter2 = tk.Label(page1, text="Durchmesser 2:")

entry_height = tk.Entry(page1)
entry_diameter1 = tk.Entry(page1)
entry_diameter2 = tk.Entry(page1)

label_height.grid(row=0, column=0, sticky="W")
label_diameter1.grid(row=1, column=0, sticky="W")
label_diameter2.grid(row=2, column=0, sticky="W")

entry_height.grid(row=0, column=1)
entry_diameter1.grid(row=1, column=1)
entry_diameter2.grid(row=2, column=1)

button1 = tk.Button(page1, text="Eingabe übernehmen", command=get_values_page1)
button1.grid(row=3, column=0, columnspan=2, pady=10, padx=10, ipadx=100)

label_material = tk.Label(page2, text="Material:")
label_layer_structure = tk.Label(page2, text="Schichtenaufbau:")

entry_material = tk.Entry(page2)
entry_layer_structure = tk.Entry(page2)

label_material.grid(row=0, column=0, sticky="W")
label_layer_structure.grid(row=1, column=0, sticky="W")

entry_material.grid(row=0, column=1)
entry_layer_structure.grid(row=1, column=1)

button2 = tk.Button(page2, text="Eingabe übernehmen", command=get_values_page2)
button2.grid(row=3, column=0, columnspan=2, pady=10, padx=10, ipadx=100)

root.mainloop()

def nonlinear_stiffness_matrix_3d(E, Ix, Iy, Iz, L, theta_x = 0, theta_y = 0, theta_z=0):
    # Berechne die Jacobi-Matrix
    dN_dX = np.array([
        [-1.0,  1.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, -1.0,  1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, -1.0,  1.0],
        [-L/2, L/2, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, -L/2, L/2, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, -L/2, L/2]
    ])
    dN_dX_inv = np.linalg.inv(dN_dX)
    dN_dU = np.dot(dN_dX_inv, np.array([
        [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, Ix, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, Iy, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, Iz]
    ]))
    dN_dU *= E
    
    # Berechne die Steifigkeitsmatrix
    K = np.dot(dN_dU.T, dN_dU)
    return K

#print (nonlinear_stiffness_matrix_3d(2000, 10, 200, 200,5))