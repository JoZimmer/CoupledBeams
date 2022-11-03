


HOLZBAU = {
    'k_mod':{'kurz': 0.9, 'mittel':0.8, 'lang': 0.7, 'ständig':0.6},
    # kdef von Martin Tabelle 3.7
    'k_def':{'NKL1':{
                      30: 0.2, # Belastungsgrad > 30%
                      15:0.1,   # Belastungsgrad < 15%
                      'sonst':0.15
            },
            'NKL2':{
                      30: 0.4, # Belastungsgrad > 30%
                      15:0.3,   # Belastungsgrad < 15%
                      'sonst':0.2
            }
    },
    'gamma_m': 1.2, 
    'k_sys': 1.2,
    'mu':0.5,  # Reibungsbeiwert
}
#
# def calculate_k_sys(Breite_BSP):
#     '''k_sys: Systembeiwert zur erhöhung des Bauteilwiderstands wegen statistischer Effekte gegenüber Brettfestigkeit'''
#     k_sys= 1.08+0.2*Breite_BSP
#     if k_sys >1.2:
#         k_sys= 1.2
#     return k_sys



charakteristische_werte = {
    'units':{'Festigkeit':'N/mm²', 'Steifigkeit':'N/m²', 'Rohdichte':'kg/m³'},
    'C24':{
        'fmk':24,
        'ft0k':14,
        'ft90k':0.4,
        'fc0k':21,
        'fc90k':2.5,
        'fvk':4.0,
        'fvxyk':5.5, # Standartwert 
        'fvtork':2.5, # Standartwerd
        'rhok':460,
        'E0mean':12000E+06,
        'E90mean':370E+06,
        'alpha':0.02, # schwindmaß je % Holzfuechteänderung
        'alphaT':4E-06, # K^-1 Temperaturausdehnungskoeffizeint
    },
    'C30':{
        'fmk':30,
        'ft0k':18,
        'ft90k':0.4,
        'fc0k':23,
        'fc90k':2.7,
        'fvk':4.0,
        'fvxyk':5.5,  # Standartwert 
        'ftornodek':2.5, 
        'fvFek':31, # Scherfestigkeit einer 45x45° Furnierebene
        'rhok':460,
        'E0mean':12000E+06,
        'E90mean':370E+06,
        'Furnierebene':True,
        'rhok_Fe':735, # aus Lechner Baubuche Furnier
        'G0mean':500, # zwischen Lechner BSH 586 und ETA Radiusholz 690
        'G4545':4227,# Aus Dissertation ermittelte Schubsteifigkeit der Furnierebene
        'alpha':0.02, # schwindmaß je % Holzfuechteänderung
        'alphaT':4E-06, # K^-1 Temperaturausdehnungskoeffizeint
       
    }
}

# Key gibt die gesamt dicke an
# 'a': ist hier nicht aktuell 
t_querlagen = 0.04
t_furnier = 0.018 # 2 Furnierlagen
lagenaufbauten = {
    'BSP_furnier':{
        36:[{'ortho':'X','ti':0.12, 'a': 0.12},
        {'ortho':'Y','ti':t_furnier, 'a': 0.05},
        {'ortho':'X','ti':0.08, 'a': 0},
        {'ortho':'Y','ti':t_furnier, 'a': 0.05},
        {'ortho':'X','ti':0.12, 'a': 0.12}],

    },

    'BSP_normal':{
    36:[{'ortho':'X','ti':0.10, 'a': 0.12},
        {'ortho':'Y','ti':t_querlagen, 'a': 0.05},
        {'ortho':'X','ti':0.08, 'a': 0},
        {'ortho':'Y','ti':t_querlagen, 'a': 0.05},
        {'ortho':'X','ti':0.10, 'a': 0.12}],

    40:[{'ortho':'X','ti':0.12, 'a': 0.12},
        {'ortho':'Y','ti':t_querlagen, 'a': 0.05},
        {'ortho':'X','ti':0.08, 'a': 0},
        {'ortho':'Y','ti':t_querlagen, 'a': 0.05},
        {'ortho':'X','ti':0.12, 'a': 0.12}],

    44:[{'ortho':'X','ti':0.14, 'a': 0.14},
        {'ortho':'Y','ti':t_querlagen, 'a': 0.05},
        {'ortho':'X','ti':0.08, 'a': 0},
        {'ortho':'Y','ti':t_querlagen, 'a': 0.05},
        {'ortho':'X','ti':0.14, 'a': 0.14}],

    48:[{'ortho':'X','ti':0.16, 'a': 0.14},
        {'ortho':'Y','ti':t_querlagen, 'a': 0.05},
        {'ortho':'X','ti':0.08, 'a': 0},
        {'ortho':'Y','ti':t_querlagen, 'a': 0.05},
        {'ortho':'X','ti':0.16, 'a': 0.14}],

    56:[{'ortho':'X','ti':0.20, 'a': 0.14},
        {'ortho':'Y','ti':t_querlagen, 'a': 0.05},
        {'ortho':'X','ti':0.08, 'a': 0},
        {'ortho':'Y','ti':t_querlagen, 'a': 0.05},
        {'ortho':'X','ti':0.20, 'a': 0.14}],

    64:[{'ortho':'X','ti':0.20, 'a': 0.14},
        {'ortho':'Y','ti':t_querlagen, 'a': 0.05},
        {'ortho':'X','ti':0.16, 'a': 0},
        {'ortho':'Y','ti':t_querlagen, 'a': 0.05},
        {'ortho':'X','ti':0.20, 'a': 0.14}],
    }
}
