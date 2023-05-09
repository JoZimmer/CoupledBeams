from scipy.optimize import minimize
from functools import partial
import numpy as np
import matplotlib.pyplot as plt

def f(parameter):
    '''
    parameter = list aller zu optimierenden parameter
    fix_val = werte die nicht teil der optimierung sind
    -> sinnvoll für jede gewünsche optimierung eine eigene Funktion zu erstellen in diesem stil
    Alternativ: fix_val und partial weglassen und innerhalb der funktion die fixen werte definieren
    '''

    x = parameter[0] # Zahnhöhe
    y = parameter[1] # Zahnlänge, 

    fix_val1 = 3
    fix_val2 = 0.5

    return x**2 + (y-2)**2 + fix_val1 + fix_val2

x0 = [2,3]

# parameter_fix = 3
# opt_f = partial(f, parameter_fix)

# print (opt_f(x0))

res = minimize(f, x0)

print (res.x)
print (res.fun)

# werte bereich zur auswertung der Funktion definieren
x = np.linspace(-6, 6, 30)
y = np.linspace(-6, 6, 30)

# Ab hier alles nur für den 3D plot
X, Y = np.meshgrid(x, y)
Z = f([X, Y])
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(X, Y, Z, 60, cmap='viridis')
#ax.plot_surface(X, Y, Z, cmap='spring')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

# Optimierungs ergebnis als punkt plotten
ax.scatter(res.x[0], res.x[1], res.fun, color = 'red')
plt.show()
