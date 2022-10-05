'''
SUSPA-EX Drähte
T:\Projekte\Projekte 2022\(E) Entwicklung\22-E-001 - Studie Windraeder\Literatur\Vorspannung\
        2021_Zulassung_Spannverfahren2 ABG für SUSPA Draht für Windenergieanlagen.pdf
        ETA-07_0186_GER_Spannsystem.pdf

'''

spanndraht_parameter = {
    'St_1470/1670':{'Ep':200000, # N/mm² -> für Litzen (siehe SBT: 5.37 )
                    'Ap':38.5, #mm²
                    'fp01k':1420, # N/mm²
                    'alphaT':1E-05,# Wärmedehnzahl K^-1
                    'schlupf':1, #'mm'
                    'rho1000':2.5, # % TODO von Martin bisher in der Zulassung nichts dazu gefunden -> Suche "Relaxationsklassen"
                    'Pm0_x':{}}, # N
    'St_1570/1770':{'Ep':200000, # N/mm²
                    'Ap':38.5, #mm²
                    'fp01k':1500, # N/mm²
                    'alphaT':1E-05,# Wärmedehnzahl K^-1
                    'schlupf':1, #'mm'
                    'rho1000':2.5, # % TODO von Martin bisher in der Zulassung nichts dazu gefunden
                    'Pm0_x':{}}, # N
}


anzahl_drähte = [30,36,42,48,54,60,66,72,78,84]

for stahl in spanndraht_parameter:
    for n in anzahl_drähte:
        spanndraht_parameter[stahl]['Pm0_x'][n] = round(0.85 * spanndraht_parameter[stahl]['Ap'] * n * spanndraht_parameter[stahl]['fp01k'])



