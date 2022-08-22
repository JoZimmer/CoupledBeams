import numpy as np
''' 
LF1: max Biegemoment + Querkraft und Torsionsmoment 
LF2: max Querkraft + zeitgleiches Biegemoment und Torsionsmoment
LF3: max Torsionsmoment + zeitgleiches Biegemoment und querkraft 
LF4: min Torsionsmoment + zeitgleiches Biegemoment und Querkraft 
'''
lasten_design_LF = {
'LF1':{
    #maximales  Biegemoment
    'Mt': 3.25E+06,
    'Md': 1.33E+08,
    'Mimp':100,
    'Nd':-1E+07,
    'Q': 1.33E+06},
'LF2':{ 
    #maximale Querkraft
    'Mt': 3.25E+06,
    'Md': 1.33E+08,
    'Mimp':100,
    'Nd':-1E+07,
    'Q': 1.33E+06},
'LF3':{ 
    #maximales positives Torsionsmoment 
    'Mt': 3.25E+06,
    'Md': 1.33E+08,
    'Mimp':100,
    'Nd':-1E+07,
    'Q': 1.33E+06},
'LF4':{ 
    #maximales negatives Torsionsmoment
    'Mt': 3.25E+06,
    'Md': 1.33E+08,
    'Mimp':100,
    'Nd':-1E+07,
    'Q': 1.33E+06},
'LF5':{
    #Vorspannung
    'Nd'
}}


