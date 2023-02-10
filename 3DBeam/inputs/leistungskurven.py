import numpy as np
import matplotlib.pyplot as plt

'''
dict: Anlagenname: [0 = v_werte, 1 = Leistung [kW], 2 = polyfit_coeffs, 3 = ausgewertetes polynom an v_werten]
'''
def plot_fit_leistungskurve(anlagen):
    colors = ['tab:blue','tab:green','tab:orange']
    fig, ax = plt.subplots(1,1, num='Leistungkurve')
    
    for i, anlage in enumerate(anlagen):
        ax.scatter(anlagen[anlage][0],anlagen[anlage][1], marker='o', label = anlage, color=colors[i])
        #ax.plot(anlagen[anlage][0],anlagen[anlage][3], color=colors[i])
        ax.plot(anlagen[anlage][0],anlagen[anlage][4], color=colors[i])

    ax.set(xlabel='v_hub [m/s]', ylabel='Leistung [kW]', title= 'Leistungskurve')
    ax.legend()
    ax.grid()
    plt.show()

# NOTE WTN250 Kurve für Anlage mit 50 m Nabenhöhe
anlagen = { 'WTN250':[np.arange(0,26,1), [0,0,0,0,10,25,40,58,95,125,160,190,215,230,245,248,251,255,251,248,240,232,222,213,205,195]],
            'E30':[np.linspace(2,25,24,True), [0,9,15,23,45,60,97,132,181,230,272,290,300,300,300,300,300,300,300,300,300,300,300,300]]}

anlagen_dict = {'WTN250':dict(zip(anlagen['WTN250'][0],anlagen['WTN250'][1])),
                'E30':dict(zip(anlagen['E30'][0],anlagen['E30'][1]))}

degree = 5
for i,anlage in enumerate(anlagen):
    

    polynomial = np.poly1d(np.polyfit(anlagen[anlage][0], anlagen[anlage][1], deg=degree))
    anlagen[anlage].append(polynomial)

    if anlage== 'E30':
        anlagen[anlage].append(polynomial(anlagen[anlage][0][:-11])) #evaluate
        anlagen[anlage][-1]=np.append(anlagen[anlage][-1], [300]*11) #evaluate

        coeffs = np.polyfit(anlagen[anlage][0][1:-11], np.log(anlagen[anlage][1][1:-11]), deg=1)
        y_exp = coeffs[0] * np.exp(coeffs[1] * anlagen[anlage][0][:-11])
        y_exp = np.append(y_exp, [300]*11) #evaluate
    else:
        anlagen[anlage].append(polynomial(anlagen[anlage][0])) #evaluate
        coeffs = np.polyfit(anlagen[anlage][0][4:], np.log(anlagen[anlage][1][4:]), deg=1)
        y_exp = coeffs[0] * np.exp(coeffs[1] * anlagen[anlage][0])

    anlagen[anlage].append(y_exp)


#if __name__ == '__main':
#plot_fit_leistungskurve(anlagen)
