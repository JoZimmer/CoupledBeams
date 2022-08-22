HOLZBAU = {
    'k_mod':{'kurz': 0.9, 'mittel':0.8, 'lang': 0.7, 'ständig':0.6},
    'gamma_m': 1.2, 
    'k_sys': 1.2
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
        'rhok':350,
        'E0mean':11000E+06,
        'E90mean':370E+06
    },
    'BSP_RFEM':{
        'fmk':24,
        'ft0k':14,
        'ft90k':0.4,
        'fc0k':21,
        'fc90k':2.5,
        'fvk':4.0,
        'rhok':460,
        'E0mean':12000E+06,
        'E90mean':370E+06
    }
}