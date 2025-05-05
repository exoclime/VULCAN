import os

# Path to VULCAN folder
VULCAN_DIR = os.path.dirname(os.path.abspath(__file__)) + '/'

# Path to FastChem folder (custom version packaged with VULCAN)
FASTCHEM_DIR = os.path.join(VULCAN_DIR, 'fastchem_vulcan') + '/'

# Thermodynamic and chemistry data
THERMO_DIR = os.path.join(VULCAN_DIR,'thermo',)

# Composition and gas properties
COM_FILE = os.path.join(VULCAN_DIR,'thermo','all_compose.txt')

# (all the nasa9 files must be placed in the folder: thermo/NASA9/)
GIBBS_FILE = os.path.join(VULCAN_DIR,'thermo','gibbs_text.txt')

# Photochemistry cross-sections
CROSS_DIR = os.path.join(VULCAN_DIR,'thermo','photo_cross') + '/'

# Symbolic chemical functions
CHEM_FUNS_FILE = os.path.join(VULCAN_DIR,'chem_funs.py')
