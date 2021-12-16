

DOF_LABELS = {'2D': ['y','g'], '3D':['x', 'y', 'z', 'a', 'b', 'g']}
n_dofs_node = {'2D':2, '3D':6}
dofs_of_bc = {'2D': [0,1], '3D':[0,1,2,3,4,5] }
tip_load = {'2D': -2, '3D':-5}

RESPONSE_DIRECTION_MAP = {'Qx':'x', 'Qy':'y', 'Qz':'z', 'Mx':'a', 'My':'b', 'Mz':'g'}

DIRECTION_LOAD_MAP = {'x':'Fx', 'y':'Fy', 'z':'Fz', 'a':'Mx', 'b':'My', 'g':'Mz'}

DIRECTION_RESPONSE_MAP = {'x':'Qx', 'y':'Qy', 'z':'Qz', 'a':'Mx', 'b':'My', 'g':'Mz'}

GREEK = {'y':'y','z':'z', 'x':'x','a':r'\alpha', 'b':r'\beta', 'g':r'\gamma'}

UNITS_LOAD_DIRECTION = {'x':'[N/m]', 'y':'[N/m]', 'z':'[N/m]', 'a':'[Nm/m]', 'b':'[Nm/m]', 'g':'[Nm/m]'}

UNITS_POINT_LOAD_DIRECTION = {'x':'[N]', 'y':'[N]', 'z':'[N]', 'a':'[Nm]', 'b':'[Nm]', 'g':'[Nm]'}

UNIT_SCALE = {'N':1.0,'KN':1e-3,'MN':1e-6}

MODE_CATEGORIZATION_REVERSE = {'3D':{'x':'longitudinal',
                                      'a':'torsional',
                                      'z':'sway_y',
                                      'b':'sway_y',
                                      'y':'sway_z',
                                      'g':'sway_z'}}

CAARC_MODE_DIRECTIONS = {'0_deg':{'sway_z':0, 'sway_y':1, 'torsional':2}}
