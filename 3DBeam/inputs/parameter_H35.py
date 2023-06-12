
params_dict = {
    "model_parameter":{
        "Inhalt":"alles was zum beam gehoert",
        'dimension': '2D', # '3D' auch möglich hat aber auf die statische analyse keinen einfluss
        'n_elements': 10, # Anzahl FE Elemente
        'nacelle_mass': 13325, # kg WTN Lastdaten #14400, # kg Kleinwind,# IEA37:267910,#  Gondel Masse in kg
        'imperfektion':0.008, # m/m Schiefstellung + Fundament schief
        'E_Modul': 12000E+06,# N/m²
        'type_of_bc':'Eingespannt',#'Feder'# TODO Feder evlt Falsch gemacht
        'spring_stiffness':[0.71E+09,14.29E+09], # Federsteifigkeit in N/m bzw. Nm/rad am Boden in u (horizontal) und gamma richtung 
        # Vergeleiche: Bögl 40 GNm/rad | TimberTower 10 GNm/rad | WTN 14.29 GNm/rad hor: 0.71 GN/m vert: 0.86 GN/m
        'damping_coeff': 0.01, # D = Prozent of the critical (siehe Bemssungskonzept Kleinwind für werte und Beschreibung)
        },

    "einheiten_input":{
        "Inhalt":"Einheiten Eingabe Werte",
        "Kraft":"N", "Moment":"Nm", "Festigkeit":"N/mm²", "Länge":"m", "Masse":"kg"
    }, 

    "material":{
        "Inhalt":"parameter die dann die holz Werte aus holz.charakteristische Werte holen",
        "holzguete":"C24"
    },

    "lagenaufbau":{
        "Inhalt":"Der Input fuer die Lagenaufbau Klasse aus holz.py (dicken ist in cm alles andere m)",
        "mit_furnier":True, # wenn "typ" in nächster Zeile angegeben ist, dann fällt das hier weg
        "typ":"BSP_furnier", # "BSP_normal"
        "dicken":[20], # Muss in inputs.holz enthalten sein
        "t_querlagen":0.04,
        "t_furnier":0.009
    },
    
    "turm_geometrie":{
        "Inhalt":"Hoehe etc. Durchmesser ueber die Hoehe. Hoehe sektionen: wenn 2 werte dann range sodass es auf die nabenhoehe kommt. Die meisten parameter beziehen sich auf definition mit knick.",
        "nabenhöhe" :35, # Entspricht Turm höhe
        "d_unten_oben" :[3.0, 1.8],
        "höhe_sektionen" :[11,12],# wenn 2 dann ist es eine range 
        "transportbreite_max":3,
        # ab hier nur sachen für knicke z.B.: knicke = {'gerade':[None,10], '1Knick':[4.5,10], 'Knick mit Übergang':[4.5,8]}
        "d_knick" :None,
        "h_knick_von_oben" :70,
        "d_unten_angepasst":None, 
        "d_knick_uebergang":"automatisch",
        "n_sektionen_uebergang":1
    },
    
    "lasten":{
        "einheiten_in":{"Kraft":"N","Moment":"Nm"},#von IEA
        "kopflasten":{
            "Inhalt":"Rotorflaechen [0] / [1] dividiert und ggf mit Faktor [2] multipliziert. Zu rechnende Lastfaelle",
            "skalierungs_rotor_flaechen":[1, 1, 1],# skalierung der IEA lasten [707,13273, 2], 0.:
            # der erste lastfall muss der zu max_druck sein damit die erforderliche spannkraft berechnet werden kann
            # die @... sind keys im Kopflast dictonary
            "lastfälle":{"max_druck": "@max_Fxy","max_all":"@max_all"}, # zu berechnende
            # Im Kosys von IEA eingeben: x = Windrichtung z = Turmachse
            "lastwerte":{ '@max_Fxy':{'Fx':1.33E+05,'Fy':300, 'Fz':-1.80E+05, 'Mx':5.2E+03, 'My':5.79E+04, 'Mz':3.3E+03},# NOTE My = Mxy 
                          '@max_Mz':{'Fx':37000,'Fy':3300, 'Fz':-152700, 'Mx':69800, 'My':-66100, 'Mz':1.288E+05}, # Fx = Fxy, My = Mxy
                          '@max_all':{'Fx':1.33E+05,'Fy':300, 'Fz':-1.80E+05, 'Mx':5.2E+03, 'My':5.79E+04, 'Mz':1.288E+05}} # bisher max Fx und max Mz 
        },
        "turmlasten":{
            "Inhalt":"paramter fuer die DIN Windlasten auf den Turm. wind_nabenhoehe: anhand dessen die vb in 10 m bestimmt wird",
            "cd_zylinder":0.8, # TODO berechnen
            "terrain_kategorie":"II",
            "windzone":2,
            "wind_nabenhoehe":35 # Windgeschindigkeit auf Nabenhöhe
        }},
        
    "vorspannung":{
        "Inhalt":"Parameter der Spannglieder und alles was dazu gehoert.",
        # n_segmente_ext: Anzahl der Segmente die durch Externe Vorspannung Überdrückt sein sollen: Inklusive der Fuge Übergnag stahl holz, Exkl. der Fußfuge des n_ext.ten segments
        "n_segmente_extern":3,
        "inkl_fuß_ext":True, # ob auch die Kraft am Fuß des n_exten segments durch die externe Vorspannung abgedeckt werden soll (als n_segmente + 1 Ebenen)
        "masse_spannkanal":[0.06,0.04], # Maße in m 
        "n_draehte_suspa_ex":20,# Nabenhöhe > 100: 84
        "stahl_suspa_ex":"St_1570/1770",
        "n_monolitzen":3, # Nabenhöhe > 100:5
        "stahl_monolitzen":"St_1860",
        "spannkraft_verluste_pauschal":20 # % oder None -> dann müssen aber weitere parameter für die verlust berechnung angegeben werden
        # TODO Spannkraftsverlust parameter hier ergänzen
    },

    "nachweise":{
        "Inhalt":"sicherheiten, dauern, holzbau aus GD",
        "einwirkungsdauer":["ständig", "kurz", "egal","spannkraft"],
        # dlc entsprechend der Kopflasten q diesem entsprechend und g TODO ansich komplexer im IEC
        "sicherheitsbeiwerte": {"dlc":1.35, "wind":1.35, "g":1.35, "vorspannung":1.0, "guenstig":1.0},
        "gamma_M":1.25,
        "NKL":2,
        "k_mod":{"kurz": 0.9, "mittel":0.8, "lang": 0.7, "ständig":0.6},
        "k_def":"z.B. fuer die Detaillierten Spannkraftverluste -> siehe hier MGR: kdef abhaengig von NKL und Belastungsgrad", #TODO auch für Steifigkeits variation über Lebendauer!!
        "k_sys":1.0
    },

    "output":{
        "Excel_Datei":["output","Berechnungs_Ergebnisse.xlsx"],
        "einheiten":{
            "lasten":{"Kraft":"kN","Moment":"kNm"},
            "schnittgrößen":{"Kraft":"MN","Moment":"MNm"},# Lasten einheiten sind in 'lasten' angegeben
            "sgr/m":{"Kraft":"kN","Moment":"kNm"},# SGR pro m
            "Spannung":"N/mm²",
            "Vorspannung":"MN"},
        "ausgabe_an_knoten" :'FE',# 'Ebenen'#
        "Nachkommastellen":{
            "default":2,
            "QS_werte":2,
            "Ausnutzung":3,
            "Maße":1,
            "Vorspannung":2
        },
        "on_off":[]}# eigentlich mal was für die writer klasse mit den dataframes
}

