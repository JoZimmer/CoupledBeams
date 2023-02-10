import numpy as np

'''
Das Skript ist dazu gedacht sämtliche (relevante) Formlen aus der Norm für Windmodelle zentral zu Sammeln.
Siehe DIN EN 61400-1 Auslegungsanforderungen Kapitel 6

Gleichungs nummer aus DIN EN 61400-1 wird angegeben
'''
# DIN EN 61400-1 Tabelle 1
WEA_KLASSEN = {
    'I':{'v_ave':10, 'v_ref':50},
    'II':{'v_ave':8.5, 'v_ref':42.5},
    'III':{'v_ave':7.5, 'v_ref':37.5},
    'A+':{'I_ref':0.18},
    'A':{'I_ref':0.16},
    'B':{'I_ref':0.14},
    'C':{'I_ref':0.12},
    'S':{'individuell festlegbar'}
}


#könnte irgendwann als Klasse angelegt werden
# z_hub, v_hub als parameter?
wea = {'typ':'Enercon E-115', 'D':100, 'z_hub':130} 

def A_1(z_hub):
    '''
    Turbulenzlängenparameter in Nabenhöhe z_hub
    Gl. (5)
    '''
    if z_hub <=60:
        return 0.7*z_hub
    else:
        return 42

def rayleigh(v_hub, v_ave):
    '''
    CDF
    mittelwert der Windgeschwindigkeitsverteilung über 10 min muss als so eine Rayleigh Verteilung angenommen werden 
    Gl. (8)
    v_ave ist das Jahres mittel der Windgeschwindigkeit auf Nabenhöhe
    v_hub sind sämtliche Windgeschwindigkeiten auf Nabenhöhe
    -> somit ist P_R die Wahrscheinlichkeit, dass v_hub gemessen wird
    '''

    P_R = 1 - np.exp(-np.pi*(v_hub/(2*v_ave))**2)

    return P_R

class normale_windbedinungen():

    def rayleigh(v_hub, wea_klasse):
        '''
        mittelwert der Windgeschwindigkeitsverteilung über 10 min muss als so eine Rayleigh Verteilung angenommen werden 
        Gl. (8)
        anhand der WEA Klasse wird v_ave, das Jahres mittel der Windgeschwindigkeit auf Nabenhöhe gewählt 
        NOTE für einen spezifischen Standort kann somit das tatsächliche v_ave genutzt werden
        v_hub sind sämtliche Windgeschwindigkeiten auf Nabenhöhe
        -> somit ist P_R die Wahrscheinlichkeit, dass v_hub gemessen wird
        '''

        P_R = 1 - np.exp(-np.pi*(v_hub/(2*WEA_KLASSEN[wea_klasse]['v_ave']))**2)

        return P_R

    def v_z_nwp(v_hub, wea, z):
        ''' 
        windprfil für WEA Standartklassen
        NWP: Normales Windprofilmodell
        Gl. (9)
        '''
        return v_hub*(z/wea['z_hub'])**0.2

    def sigma_1_ntm(v_hub, wea_klasse_turb):
        ''' 
        Standartabweichung der Turbulenz im normalen Turbulenzmodell in abhängigkeit von v_hub
        Gl. (10)
        '''
        sigma_1 = WEA_KLASSEN[wea_klasse_turb]['I_ref']*(0.75*v_hub+5.6) #Standartabweichung von der mittleren Windgeschwindigkeit

        print('\nAlternativ hierzu kann auch eine Weibull verteilung von sigma angenommen werden -> Gl. (11)')

        return sigma_1

class extreme_windbedinungen():

    def v_z_ewm_50 (v_hub, wea, z, wea_klasse, modell_art):
        '''
        extreme Windbedinungen mit wiederkehrzeitrum 50 Jahre
        modell_art: ob 'stationär' oder 'turbulent'
        Gl. (13,15)
        '''

        v_50 = WEA_KLASSEN[wea_klasse]['v_ref']*(z/wea['z_hub'])**0.11

        if modell_art == 'stationär':
            print ('\nfür eine kurzzeitige Abweichung von der mittleren Windrichtung muss ein konstanter Gierfehler von +/- 15° angesetzt werden')
            return 1.4 *v_50

        elif modell_art == 'turbulent':
            return v_50

    def v_z_ewm_1 (v_hub, wea:dict, z, wea_klasse, modell_art):
        '''
        extreme Windbedinungen mit wiederkehrzeitrum 1 Jahr
        modell_art: ob 'stationär' oder 'turbulent'
        Gl. (14,16)
        '''

        v_1 = 0.8*WEA_KLASSEN[wea_klasse]['v_ref']*(z/wea['z_hub'])**0.11 #0.8 * v50

        return v_1

    def sigma_1_ewm (v_hub):
        '''
        Für das turbulente extremwindmodell
        Standartabweichung der longitudinalen Turbulenzt unter extremwindbedingungen (ewm)
        Gl. (17)
        '''
        return 0.11*v_hub

    def v_gust_eog(v_hub, z, v_z, t, wea_klasse, wea:dict):
        ''' 
        Böenwert
        v_z: Windgeschwindigkeit in abhängigkeit der Höhe z
        t: Zeitpunkt t, oder time array
        Gl. (18,19)
        '''
        v_e1 = extreme_windbedinungen.v_z_ewm_1(v_hub, wea, z)
        sigma_1 = normale_windbedinungen.sigma_1_ntm(v_hub, wea_klasse)
        turb_len = A_1(wea['z_hub'])

        v_gust = min(
            1.35*(v_e1-v_hub),
            3.3*(sigma_1/(1+0.1*(wea['D']/turb_len)))
        )

        T = 10.5

        v_z_t = v_z - 0.37*v_gust * np.sin(3*np.pi*t/T) * (1-np.cos(2*np.pi*t/T))

        if t >= 0 and t <= T:
            return v_z_t

        else:
            return v_z

    def sigma_1_etm(v_hub, wea_klasse, wea_klasse_turb):
        '''
        extremes Turbulenzmodel - normales windprofil modell muss verwendet werden (normale_windbedinungen.v_z)
        Gl. (20)
        '''
        print ('in verbindung mit dem ETM muss das normale Windprofil modell verwendet werden Gl. (9)')

        c = 2 #m/s
        I_ref = WEA_KLASSEN[wea_klasse_turb]['I_ref'] 
        v_ave = WEA_KLASSEN[wea_klasse]['v_ave'] 

        sigma_1 = c * I_ref * (0.072* (v_ave/c+3) * (v_hub/c - 4) +10)

        return sigma_1

    def teta_e_v_hub_edc(v_hub, wea, wea_klasse_turb):
        '''
        Extreme windrichtungsänderung in abhängigkeit von v_hub positiv oder negativ! (EDC)
        Für die Windgeschwindigkeit ist das normale windgeschwindigkeitsprofil anzusetzen
        Gl. 21, 22
        '''
        print('\nFür die Windgeschwindigkeit ist das normale windgeschwindigkeitsprofil anzusetzen')

        sigma_1 = normale_windbedinungen.sigma_1_ntm(v_hub,wea_klasse_turb)
        turb_len = A_1(wea['z_hub'])
        
        teta_e = min(abs(4*np.arctan(sigma_1 / (v_hub * (1+0.1*(wea['D']/turb_len))))), 180) #ist auf den Bereich +/- 180° begrenzt

        return teta_e

    def teta_e_t_edc(v_hub, t, wea, wea_klasse_turb):
        '''
        Extreme windrichtungsänderung in abhängigkeit von v_hub und der Zeit t positiv oder negativ!
        Für die Windgeschwindigkeit ist das normale windgeschwindigkeitsprofil anzusetzen
        Gl. 21, 22
        '''
        print('\nFür die Windgeschwindigkeit ist das normale windgeschwindigkeitsprofil anzusetzen')

        sigma_1 = normale_windbedinungen.sigma_1_ntm(v_hub,wea_klasse_turb)
        turb_len = A_1(wea['z_hub'])
        
        teta_e = min(abs(4*np.arctan(sigma_1 / (v_hub * (1+0.1*(wea['D']/turb_len))))), 180) #ist auf den Bereich +/- 180° begrenzt

        T = 6 #6 sekunden ist die Dauer der windrichtungsänderung

        teta_e_t = 0.5*teta_e*(1- np.cos(np.pi*t/T))

        if t <= 0:
            return 0

        elif t <= T:
            return teta_e_t

        elif t> T:
            return teta_e

    def v_cg_z_t(z, t, v_hub, wea, wea_klasse):
        '''
        extreme kohärente Böe mit Richtungsänderung (ECD)
        Kohärenz zwischen Böe und Richtungsänderung -> Überlagerung der Phase der beiden "Signale"
        teta_cg positiv oder negativ
        Gl. 23-26
        '''
        v_cg = 15

        T = 10 #Anstiegszeit

        v_z = normale_windbedinungen.v_z_nwp(v_hub, wea, z)

        v_z_t = v_z + 0.5*v_cg * (1-np.cos(np.pi * t/T))

        # diese Geschwindigkeitsänderung tritt gleichzeitig mit einer Richtungsänderung auf

        if v_hub <= 4:
            teta_cg = 180
        elif v_hub < WEA_KLASSEN[wea_klasse]['v_ref']:
            teta_cg = 720/v_hub
        else:
            print ('teta_cg nicht definiert für v_hub > v_ref?!')

        teta_cg_t = 0.5*teta_cg * (1-np.cos(np.pi*t/T))

        if t < 0:
            return v_z, 0
        elif t <= T:
            return v_z_t, teta_cg_t
        else:
            return v_z + v_cg, teta_cg

    def v_z_t_ews(z, t,wea, wea_klasse, wea_klasse_turb, v_hub, sign=1):
        '''
        Extremer Windgradient vertikal: Situation in dem extreme unterschiede zwischen windgeschwindigkeiten auftreten
        vertikal und horizontal werden nicht gleichzeitig aufgebracht

        sign: 1 oder -1 um positiven oder negativen Gradient zu erhalten
        Gl. 27,28
        '''

        alpha = 0.2
        beta = 6.4
        T = 12
        sigma_1 = normale_windbedinungen.sigma_1_ntm(v_hub, wea_klasse_turb)

        v_z = v_hub*(z/wea['z_hub'])**alpha

        v_z_t =  v_z + sign * (((z-wea['z_hub'])/wea['D']) * (2.5+0.2*beta*sigma_1 * (wea['D']/A_1(wea['z_hub'])**0.25)) * (1- np.cos(2*np.pi*t/T)) )

        if t >= 0 and t <= T:
            return v_z_t
        else: 
            return v_z

    def v_y_z_t_ews(y, z, t,wea, wea_klasse, wea_klasse_turb, v_hub, sign=1):
        '''
        Extremer Windgradient horizontal: Situation in dem extreme unterschiede zwischen windgeschwindigkeiten auftreten
        vertikal und horizontal werden nicht gleichzeitig aufgebracht

        sign: 1 oder -1 um positiven oder negativen Gradient zu erhalten
        Gl. 27,28
        '''

        alpha = 0.2
        beta = 6.4
        T = 12
        sigma_1 = normale_windbedinungen.sigma_1_ntm(v_hub, wea_klasse_turb)

        v_z = v_hub*(z/wea['z_hub'])**alpha

        v_y_z_t = v_z + sign * (y/wea['D']) * (2.5+0.2*beta*sigma_1 * (wea['D']/A_1(wea['z_hub'])**0.25)) * (1- np.cos(2*np.pi*t/T)) 

        if t >= 0 and t <= T:
            return v_y_z_t
        else: 
            return v_z







