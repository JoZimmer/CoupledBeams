'''
SUSPA-EX Drähte
T:\Projekte\Projekte 2022\(E) Entwicklung\22-E-001 - Studie Windraeder\Literatur\Vorspannung\
        2021_Zulassung_Spannverfahren2 ABG für SUSPA Draht für Windenergieanlagen.pdf
        ETA-07_0186_GER_Spannsystem.pdf

'''

suspa_draht_ex = {
    'Stahlparameter':{'Ep':200000, # N/mm² -> für Litzen (siehe SBT: 5.37 )
                    'Ap':38.5, #mm²},
                    'alphaT':1E-05,# Wärmedehnzahl K^-1
                    'schlupf':1, #'mm'
                    'rho1000':2.5, # % TODO von Martin bisher in der Zulassung nichts dazu gefunden -> Suche "Relaxationsklassen"
    },
    'St_1470/1670':{  
                    'fp01k':1420, # N/mm²
                    'Pm0_n':{}}, # N
    'St_1570/1770':{
                    'fp01k':1500, # N/mm²
                    'Pm0_n':{}}, # N
}

anzahl_drähte_suspa_ex = [30,36,42,48,54,60,66,72,78,84]

for stahl in suspa_draht_ex:
    if stahl[:3] == 'St_':
        for n in anzahl_drähte_suspa_ex:
            suspa_draht_ex[stahl]['Pm0_n'][n] = round(0.85 * suspa_draht_ex['Stahlparameter']['Ap'] * n * suspa_draht_ex[stahl]['fp01k'])

monolitzen = {
    # Übernommen aus SUSPA -> muss ggf noch angepasst werden an werte aus eigener ETA
    'Stahlparameter':{'Ep':200000, # N/mm² -> für Litzen (siehe SBT: 5.37 )
                    'Ap':None, #mm²},
                    'alphaT':1E-05,# Wärmedehnzahl K^-1
                    'schlupf':1, #'mm'
                    'rho1000':2.5, # % TODO von Martin bisher in der Zulassung nichts dazu gefunden -> Suche "Relaxationsklassen"
    },
    'St_1860':{
        'fp01k':1860, # N/mm²
        'Fp01k':246000, #N
        'Pmax_n':{}, # N
        'Pm0_n':{} #
    }
}

anzahl_monolitzen = [1,2,3,4,5]

for stahl in monolitzen:
    if stahl[:3] == 'St_':
        for n in anzahl_monolitzen:
            monolitzen[stahl]['Pm0_n'][n] = round(0.85 * monolitzen[stahl]['Fp01k'] * n)
            monolitzen[stahl]['Pmax_n'][n] = round(0.9 * monolitzen[stahl]['Fp01k'] * n)






