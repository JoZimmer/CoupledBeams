

DOF_LABELS = {'2D': ['x','y','g'], '3D':['x', 'y', 'z', 'a', 'b', 'g']}

GRAVITY = 9.81 

RHO_AIR = 1.225

DOFS_PER_NODE = {'2D':3, '3D':6}

RESPONSE_DOF_MAP = {'Nx':'x', 'Qy':'y', 'Qz':'z', 'Mx':'a', 'My':'b', 'Mz':'g'}

DOF_RESPONSE_MAP = {'x':'Nx', 'y':'Qy', 'z':'Qz', 'a':'Mx', 'b':'My', 'g':'Mz'}

DOF_LOAD_MAP = {'x':'Fx', 'y':'Fy', 'z':'Fz', 'a':'Mx', 'b':'My', 'g':'Mz'}

LOAD_DOF_MAP = {'Fx':'x', 'Fy':'y', 'Fz':'z', 'Mx':'a', 'My':'b', 'Mz':'g'}

GREEK = {'y':'y','z':'z', 'x':'x','a':r'\alpha', 'b':r'\beta', 'g':r'\gamma'}

# https://pythonforundergradengineers.com/unicode-characters-in-python.html 
GREEK_UNICODE = {'sigma':'\u03C3', 'mu':'\u03BC', 'delta':'\u03C3', 'eta':'\u03B7', 'tau':'\u03C4', 'gamma':'\u03B3'}

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


class VariableNamesOpenFAST():
    '''
    die verschiedenen Variablen von FAST sortiert
    '''

    displacments = ["TTDspFA_[m]","TTDspSS_[m]"]
    forces = ['YawBrMxp_[kN-m]','YawBrMyp_[kN-m]','YawBrMzp_[kN-m]']#'YawBrFxp_[kN]','YawBrFyp_[kN]', # FROM ElastoDyn
    turbine_output = ['GenPwr_[kW]','GenSpeed_[rpm]']
    control_variables = ['BldPitch1_[deg]','RotSpeed_[rpm]', 'RotTorq_[kN-m]']#,'GenTq_[kN-m]','GenPwr_[kW]','Wind1VelX_[m/s]']
    blade_tip = ['OoPDefl1_[m]']
    tower_aerodyn = ['TwN9DynP_[Pa]']
