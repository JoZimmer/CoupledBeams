

DOF_LABELS = {'2D': ['x','y','g'], '3D':['x', 'y', 'z', 'a', 'b', 'g']}

GRAVITY = 9.81 

DOFS_PER_NODE = {'2D':3, '3D':6}

RESPONSE_DOF_MAP = {'Nx':'x', 'Qy':'y', 'Qz':'z', 'Mx':'a', 'My':'b', 'Mz':'g'}

DOF_RESPONSE_MAP = {'x':'Nx', 'y':'Qy', 'z':'Qz', 'a':'Mx', 'b':'My', 'g':'Mz'}

DOF_LOAD_MAP = {'x':'Fx', 'y':'Fy', 'z':'Fz', 'a':'Mx', 'b':'My', 'g':'Mz'}

LOAD_DOF_MAP = {'Fx':'x', 'Fy':'y', 'Fz':'z', 'Mx':'a', 'My':'b', 'Mz':'g'}

GREEK = {'y':'y','z':'z', 'x':'x','a':r'\alpha', 'b':r'\beta', 'g':r'\gamma'}

GREEK_UNICODE = {'sigma':'\u03C3', 'mu':'\u03BC', 'delta':'\u03C3', 'eta':'\u03B7', }

UNITS_LOAD_DIRECTION = {'x':'[N/m]', 'y':'[N/m]', 'z':'[N/m]', 'a':'[Nm/m]', 'b':'[Nm/m]', 'g':'[Nm/m]'}

UNITS_POINT_LOAD_DIRECTION = {'x':'[N]', 'y':'[N]', 'z':'[N]', 'a':'[Nm]', 'b':'[Nm]', 'g':'[Nm]'}

# Basis einheiten Kraft: N, Distanz: m, Druck: N/m²
# Umrechnung von basis in "key" vom dict
UNIT_SCALE = {'N':1.0,'kN':1e-3,'MN':1e-6, 'm':1, 'cm':100, 'mm':1000, 'N/mm²':1/1E+06}

MODE_CATEGORIZATION_REVERSE = {'3D':{'x':'longitudinal',
                                      'a':'torsional',
                                      'z':'sway_y',
                                      'b':'sway_y',
                                      'y':'sway_z',
                                      'g':'sway_z'}}

CAARC_MODE_DIRECTIONS = {'0_deg':{'sway_z':0, 'sway_y':1, 'torsional':2}}
