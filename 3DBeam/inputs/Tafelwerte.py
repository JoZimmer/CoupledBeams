import xlwings as xw
import matplotlib.pyplot as plt
import numpy as np


excel_input_file = 'Tafelwerte.xlsx'

wb = xw.Book(excel_input_file)
ws = wb.sheets['Tabelle1']

row_start = 4 # Excel Zeile -1 
row_end = 45

längenverhältnisse = ws[row_start:row_end, 1].value
tafelwerte = {
                'feld1':ws[row_start:row_end, 2].value,
                'feld2':ws[row_start:row_end, 3].value,
                'feld3':ws[row_start:row_end, 4].value
}

degrees = [2,3,4]

fitted_polynomial = {}
tafelwerte_fitted = {}
for feld in tafelwerte:
    tafelwerte_fitted[feld] = []
    fitted_polynomial[feld] = []

for degree in degrees:
    for feld in tafelwerte:
        y = tafelwerte[feld]
        coeffs = np.polyfit(längenverhältnisse, y, deg=degree)
        polynomial = np.poly1d(np.polyfit(längenverhältnisse, y, deg=degree))
        fitted_polynomial[feld].append(polynomial)
        # evaluate polynomial 
        current_tafelwert_fitted = polynomial(längenverhältnisse)
        tafelwerte_fitted[feld].append(current_tafelwert_fitted)

for feld in tafelwerte:
    plt.plot(längenverhältnisse, tafelwerte[feld], label=feld)
    for d_i, degree in enumerate(degrees):
        plt.plot(längenverhältnisse, tafelwerte_fitted[feld][d_i], label=feld + ' fit^' + str(degree))

    plt.xlabel('l1:l2')
    plt.ylabel('Tafelwert')
    plt.title(feld.upper() + ' - Tafelwert SBT & Tafelwert fit')
    plt.legend()
    plt.xlim(left=1)
    #plt.show()
    plt.close()

digits = (4)
for feld in tafelwerte:
    #plt.plot(längenverhältnisse, tafelwerte[feld], label=feld)
    for d_i, degree in enumerate(degrees):
        diff = tafelwerte[feld] - tafelwerte_fitted[feld][d_i]
        poly = fitted_polynomial[feld][d_i]
        f_x_text = ''
        for c_i, coeff in enumerate(poly.coeffs):
            f_x_text += str(round(coeff, digits))
            if c_i != len(poly.coeffs)-1:
                f_x_text += 'x^' + str(degree - c_i) + '+'


        plt.plot(längenverhältnisse, diff, label=f_x_text)

    plt.xlabel('l1:l2')
    plt.ylabel('Tafelwert_SBT - Tafelwert_fit')
    plt.legend()
    plt.xlim(left=1)
    plt.title(feld.upper() + ' - Fehler Tafelwert fit')
    plt.show()

start = {
        'feld1':[4,7],
        'feld2':[4,13],
        'feld3':[4,19]
}
for feld in tafelwerte:
    for d_i, degree in enumerate(degrees):
        poly = fitted_polynomial[feld][d_i]
        for c_i, coeff in enumerate(poly.coeffs):
            row = start[feld][0] + d_i
            col = start[feld][1] + c_i
            
            ws[row, col].value = coeff


