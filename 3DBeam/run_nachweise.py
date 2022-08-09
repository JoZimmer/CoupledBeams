import xlwings as xl
import source.utilities.holz as holz
import source.utilities.utilities as utils
import numpy as np
from source.Querschnitte import nEck, KreisRing 
import inputs.lasten as lasten
from os.path import join  as os_join


excel_file = os_join (*['..','..','RFEM','0704_parametrisches_RFEM_model_12eck_save2.xlsx'])
wb = xl.Book(excel_file)
ws = wb.sheets['Formel-Parameter']

ebenen_radien = utils.read_xl_column(ws, start_cell = 'D90', end_row= 104)
# boden knoten dazu 
heights_list = utils.read_xl_column(ws, start_cell='D75', end_row=88)
heights_list.append(0)
heights_parameter = {}
# flip -> von unten nach oben
heights_parameter['absolute_höhen'] = np.flip(np.array(heights_list))
heights_parameter['hfract'] = np.flip(np.array(utils.read_xl_column(ws, start_cell='J90', end_row=104)))

# hier soll definiert werden in welchen Einheiten inputs an Querschnitt ÜBERGEBEN wird
# 
einheiten = {'Kraft':'N', 'Moment':'Nm', 'Festigkeit':'N/mm²', 'Länge':'m'}

n_ecken = 12
ecken_possible = [8,10,12]
d_achse = np.flip(np.array(ebenen_radien) *2)
t_wand = 0.4

lagen_aufbau = [{'ortho':'X','ti':0.08, 'di':d_achse + 0.03 + 0.04 + 0.04},
                {'ortho':'Y','ti':0.04, 'di':d_achse + 0.03 +0.02},
                {'ortho':'X','ti':0.06, 'di':d_achse},
                {'ortho':'Y','ti':0.04, 'di':d_achse -  0.03 - 0.02},
                {'ortho':'X','ti':0.08, 'di':d_achse - 0.03 -  0.04 -  0.04}]

kreis_ring = KreisRing(d_achse, 0.4, lagen_aufbau=lagen_aufbau,
                    holz_parameter = holz.charakteristische_werte['BSP_RFEM'], nachweis_parameter = holz.HOLZBAU, hoehen_parameter= heights_parameter, einheiten=einheiten)
nEck12_fuß = nEck(12, d_achse, t_wand, lagen_aufbau=lagen_aufbau, 
                    holz_parameter = holz.charakteristische_werte['BSP_RFEM'], nachweis_parameter = holz.HOLZBAU, hoehen_parameter= heights_parameter, einheiten=einheiten)
kreis_ring.calculate_ausnutzung_normalspannung(lasten.lasten_design_LF['LF1'])
kreis_ring.calculate_ausnutzung_schub(lasten.lasten_design_LF['LF1'])
nEck12_fuß.calculate_ausnutzung_normalspannung(lasten.lasten_design_LF['LF1'])
nEck12_fuß.calculate_ausnutzung_schub(lasten.lasten_design_LF['LF1'])
nEck12_fuß.ausnutzung_Plattenbeanspruchung_Nebentragrichtung('kurz', 3200, 1)