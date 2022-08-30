'''
Vielleicht irgendwann mal wirkliche Tests
Bisher mit Excel und Papier getestet
'''

nabenhöhe = 110 #m
d_achse_unten = 12
d_achse_oben = 3.4
t = 0.48

A_unten = 18.10
V_ges = 1277.3
Iz_unten = 326.24

# mit SGR aus FE Model mit Belastung nur IEA Kopf
Fy = 1.17 # MN
Mz = 3.24 # MNm
Fx = -3.64 # MN

# Fußpunkt SGR inkl. Eigengewicht:
Mz = 187 # MNm
Nx = -9.4 #MN

sigma_druck = 4.096 #MN/m² 