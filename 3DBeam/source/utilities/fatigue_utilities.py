import numpy as np

def find_mean_start(mean_min, step_size):
    '''
    fetlegen eines mean bin starts anhand der step size 
    '''
    if mean_min > 0:
        start = int(mean_min/step_size) * step_size
    else:
        start = (int(mean_min/step_size) - 1) * step_size
    return start

def combine_markovs(markov_gesamt, step_size):
    '''
    markov_gesamt: dict mit mean, range N(rfcmat) als keys, values sind listen der einzelenen zu kombinierenden markovs
    step_size: der mean und range bins NOTE bei beiden die selbe bisher
    '''
    # grenzen der gesamt markov finden
    mean_min = min([min(x) for x in markov_gesamt['mean']])
    mean_max = max([max(x) for x in markov_gesamt['mean']])
    range_min = min([min(x) for x in markov_gesamt['range']])
    range_max = max([max(x) for x in markov_gesamt['range']])

    mean_bin_gesamt = list(np.arange(mean_min, mean_max+step_size, step_size))
    range_bin_gesamt = list(np.arange(range_min, range_max+step_size, step_size))

    # NOTE rfcmat wurde mit der transponierten form von mean und range bestimmt, mean_bin udn mean_range sind aber immer noch die orginale
    # mean: col, range: row
    n_cols = len(mean_bin_gesamt)
    n_rows = len(range_bin_gesamt)
    
    rfcmat_gesamt = np.zeros((n_rows, n_cols)) # TODO erst zeilen dann spalten richtig?

    # einsortieren der einzelenen matrizen in die gesamte
    for i, rfcmat in enumerate(markov_gesamt['N']):
        start_col= mean_bin_gesamt.index(markov_gesamt['mean'][i][0])
        end_col= mean_bin_gesamt.index(markov_gesamt['mean'][i][-1])
        start_row = range_bin_gesamt.index(markov_gesamt['range'][i][0])
        end_row = range_bin_gesamt.index(markov_gesamt['range'][i][-1])

        rfcmat_gesamt[start_row:end_row, start_col:end_col] += rfcmat
    
    return {'mean':mean_bin_gesamt, 'range':range_bin_gesamt, 'N':rfcmat_gesamt}

def get_R_and_abs_max(sigma1, sigma2):
    '''
    sigma könnte jede spannung sein
    R = sigma_min/sigma_max für abs(sigma_min) < abs(sigma_max)
    '''
    signed = [sigma1, sigma2]
    sigma_max = max(signed)
    sigma_min = min(signed)
    abs_val = [abs(sigma1), abs(sigma2)]
    numerator = signed[abs_val.index(min(abs_val))]
    denominator = signed[abs_val.index(max(abs_val))]

    abs_max = denominator

    return numerator/denominator, abs_max, sigma_max, sigma_min

def lgN_R_kfat_din(kfat, R, last:str):
    '''
    lgN und N als Funktion von kfat und R (Din Gleichung mit a,b einfach umgestellt)
    last: durck, zug, schub
    '''

    a = {'druck':2.0, 'zug':9.5, 'schub':6.7,'schub_Fe':6.7, 'schub_längs':6.7}
    b = {'druck':9.0, 'zug':1.1, 'schub':1.3, 'schub_Fe':1.3, 'schub_längs':1.3}

    lgN = ((1-kfat) * a[last] * (b[last]-R))/(1-R)

    N = 10**lgN
    '''
    und das jetzt alles mal anschauen für
      kfat = x * kmod/gamma_m (also ansatz Schwingung um X geringer als maximale statische last)
      R = 0.1 - 0.9 (untere Grenze Validieren bzw. argumentieren)
    '''

    return lgN, N

def lgN_R_kfat_zü(kfat, R, last:str, lgN_case = 1):
    '''
    Formeln von Züblin
    lgN: erst guess dann wirds berechnet und ggf ne andere Formel genutzt
    N_ertragbar aus Züblin Versuchen für LVL 
    Dokument: "\Relevanteste Daten von Züblin aus 16-G-015\02-Gutachten-Versuche\01-MPA-Stuttgart_LVL_QS_Vertikalfuge\01_LENO_LVL\Gutachten 902 6504 000_3-Ai Merk Timber_12.05.2014.pdf" 
    '''
    if last == 'druck':
        def get_case(lgN):
                if lgN > 0 and lgN <= 5:
                    return 1
                elif lgN > 5 and lgN <= 6:
                    return 2
                elif lgN > 6:
                    return 3

        if lgN_case == 1:
            result_lg_init = (kfat - 1) / (0.01395*R**2 + 0.004765*R - 0.06160)
        
        elif lgN_case == 2:
            result_lg_init = (kfat - (0.05494* R**2 - 0.06043*R + 1.00549)) / (0.00296*R**2 + 0.05974*R -0.0627)

        elif lgN_case == 3:
            result_lg_init = (kfat - (-0.40735* R**2 + 0.03202*R + 1.37532)) / (0.08333*R**2 + 0.04367*R -0.127)

        case = get_case(result_lg_init)
        if case == lgN_case:
            result_N = 10**result_lg_init
            return result_lg_init, result_N

        else:
            return lgN_R_kfat_zü(kfat, R, 'druck', case)

    elif last == 'schub':
        def get_case(lgN):
                if lgN > 0 and lgN <= 6:
                    return 1
                elif lgN > 6:
                    return 2
                elif lgN == float("inf"):
                    return 2
        
        if kfat >= 1.0:
            return 1, 1 # TODO das sollte anisch gar nicht vorkommen, da es bedeutet das der statische nachweis nicht eingehalten wird

        if lgN_case == 1:
            if R == 1:
                result_lg_init = float("inf") # lgN wird negativ 
            else:
                result_lg_init = (kfat - 1) / (0.009135*R**2 + 0.05962*R - 0.06875)
        
        elif lgN_case == 2:
            if R == 1:
                result_lg_init = float("inf") # divisor wird 0
            else:
                result_lg_init = (kfat - (-0.27049* R**2 +0.05410*R + 1.21640)) / (0.05358*R**2 + 0.05073*R -0.10431)

        if 'result_lg_init' not in locals():
            print()
        else:
            case = get_case(result_lg_init)
        if case == lgN_case:
            result_N = 10**result_lg_init
            return result_lg_init, result_N

        else:
            return lgN_R_kfat_zü(kfat, R, 'schub', case)

def N_von_R(fatigue_dict, bin_size=0.1):
    '''
    soll für einen Bereich bin_size von R die jeweilige gesamte Anzahl an Schwingspielen auf diesem Spannungsniveau bestimmen
    '''
    min_R = min(fatigue_dict['R'])
    max_R = max(fatigue_dict['R'])
    n_bins = round((max_R-min_R)/bin_size + 0.5) # aufrunden
    bin_edges = np.histogram(fatigue_dict['R'], bins=n_bins)[1]
    bin_edges[0] -= 0.0001

    bin_centers = np.around((bin_edges[:-1] + bin_edges[1:]) / 2, 2)

    N_R = np.zeros(len(bin_edges)-1)

    for i, R_i in enumerate(fatigue_dict['R']):
        b = np.where(bin_edges < R_i)
        bin_ist = int(b[0][-1])
        N_R[bin_ist] += fatigue_dict['N_ist'][i]

    return N_R, bin_centers

