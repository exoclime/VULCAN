import os

# Path to VULCAN folder
VULCAN_DIR = os.path.dirname(os.path.abspath(__file__)) + '/'

# Composition and gas properties
COM_FILE = os.path.join(VULCAN_DIR,'thermo','all_compose.txt')

# (all the nasa9 files must be placed in the folder: thermo/NASA9/)
GIBBS_FILE = os.path.join(VULCAN_DIR,'thermo','gibbs_text.txt')

# Photochemistry cross-sections
CROSS_FOLDER = os.path.join(VULCAN_DIR,'thermo','photo_cross') + '/'

# Symbolic chemical functions
CHEM_FUNS_FILE = os.path.join(VULCAN_DIR,'chem_funs.py')
