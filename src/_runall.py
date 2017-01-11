import sys, os, datetime

import _config, _clean

from mylib import util

# Import all your steps
import a_getgrnas, b_bowtie

##############################################################
##############################################################
util.shell_cp(_config.SRC_DIR + '_config.py', _config.RESULTS_PLACE)

start = datetime.datetime.now()
print start
##############################################################
##############################################################

a_res_dir = a_getgrnas.main(_config.DATA_DIR + '/2016-11-14/', 
  _config.OUT_PLACE, 
  run = True)

b_res_dir = b_bowtie.main(a_res_dir, 
  _config.OUT_PLACE,
  run = True)


##############################################################
##############################################################


print '\n\nDone.\n'
end = datetime.datetime.now()
print end, '\nTotal Time:', end - start

_clean.main()