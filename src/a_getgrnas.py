

import _config
import sys, os, fnmatch, datetime, subprocess
import numpy as np

from mylib import util
from mylib import compbio

# Default params
DEFAULT_INP_DIR = _config.DATA_DIR
NAME = util.get_fn(__file__)

# Functions
def get_grnas(out_dir):
  for i in range(len(_config.d.NMS)):
    nm, spc, chrm, pos = _config.d.NMS[i], _config.d.SPCS[i], _config.d.CHRMS[i], _config.d.POS[i]
    print '\t', nm
    sq = compbio.get_genomic_seq_twoBitToFa(spc, chrm, str(pos[0]), str(pos[1]))

    hdrs, sqs = [], []
    print '\t\tRegion Length:', len(sq), 'bp'
    prev, too_close_filtered = 0, 0
    for j in range(len(sq) - 3):
      found = False
      if sq[j : j+2] == 'GG':
        found, start, end, orient = True, j - 22, j + 1, '+'
      if sq[j : j+2] == 'CC':
        found, start, end, orient = True, j, j + 23, '+'
      if found:
        if start - prev < _config.d.MIN_DIST:
          too_close_filtered += 1
          continue
        s = sq[start : end]
        if len(s) != 23:
          continue
        hdr = '>' + nm + '__' + chrm + '__' + str(pos[0] + start) + '__' + str(pos[0] + end) + '__' + orient
        hdrs.append(hdr)
        sqs.append(s)
        prev = start

    with open(out_dir + nm + '.fa', 'w') as f:
      for h, s in zip(hdrs, sqs):
        f.write(h + '\n' + s + '\n')
    
    print '\t\tFound', len(hdrs), 'candidate gRNAs'
    print '\t\tFiltered', too_close_filtered, 'gRNAs for being within', _config.d.MIN_DIST - 1, 'bp from another gRNA'
    if len(hdrs) < 12472:
      print '\t\tWARNING: Fewer than 12472 candidates found'



  return



@util.time_dec
def main(inp_dir, out_dir, run = True):
  print NAME  
  util.ensure_dir_exists(out_dir)
  if not run:
    print '\tskipped'
    return out_dir

  # Function calls
  get_grnas(out_dir)

  return out_dir


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1], sys.argv[2])
  else:
    main(DEFAULT_INP_DIR, _config.OUT_PLACE + NAME + '/')