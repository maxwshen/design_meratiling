import sys

PRJ_DIR = '/cluster/mshen/prj/design_meratiling/'  
SRC_DIR = PRJ_DIR + 'src/'

# toy = True
toy = False
if toy:
  PRJ_DIR += 'toy/'
#######################################################
# Note: Directories should end in / always
#######################################################
DATA_DIR = PRJ_DIR + 'data/'
OUT_PLACE = PRJ_DIR + 'out/'
RESULTS_PLACE = PRJ_DIR + 'results/'
#######################################################
#######################################################

CLEAN = False       # Values = 'ask', True, False

# which data are we using? import that data's parameters
DATA_FOLD = '2016-11-14/'

sys.path.insert(0, DATA_DIR + DATA_FOLD)
import _dataconfig as d
print 'Using data folder:\n', DATA_DIR + DATA_FOLD
DATA_DIR += DATA_FOLD
OUT_PLACE += DATA_FOLD
RESULTS_PLACE += DATA_FOLD

#######################################################
# Project-specific parameters
#######################################################

chiapetwindow = 500   #+/-, so total window is 2x