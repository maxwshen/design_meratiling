

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
    nm = _config.d.NMS[i]
    spc = _config.d.SPCS[i]
    chrm = _config.d.CHRMS[i]
    pos = _config.d.POS[i]
    print '\t', nm
    sq = compbio.get_genomic_seq_twoBitToFa(spc, chrm, str(pos[0]), str(pos[1]))

    RepeatMasker = compbio.RepeatMasker(spc)
    RepeatMasker.trim(chrm, pos[0], pos[1])

    headers, sqs = [], []
    print '\tRegion Length:', len(sq), 'bp'
    prev, too_close_filtered = 0, 0
    num_repeats_found = 0
    timer = util.Timer(total = len(sq))
    for j in range(len(sq) - 3):
      found = False
      if sq[j : j+2] == 'GG':
        found = True
        start, end = j - 22, j + 1
        orient = '+'
      if sq[j : j+2] == 'CC':
        found = True
        start, end = j, j + 23 
        orient = '-'
      if found:
        # filter grnas that are too close
        if start - prev < _config.d.MIN_DIST:
          too_close_filtered += 1
          continue

        # filter incomplete grnas
        s = sq[start : end]
        if len(s) != 23:
          continue

        # filter repeats
        if RepeatMasker.search(chrm, pos[0] + j - 10, pos[0] + j + 10):
          num_repeats_found += 1
          continue

        # if all is ok
        hdr = '>' + nm + '__' + chrm + '__' + str(pos[0] + start) + '__' + str(pos[0] + end) + '__' + orient
        headers.append(hdr)
        sqs.append(s)
        prev = start

      timer.update()

    with open(out_dir + nm + '.fa', 'w') as f:
      for h, s in zip(headers, sqs):
        f.write(h + '\n' + s + '\n')
    
    print '\tFound', len(headers), 'candidate gRNAs\n'
    print '\tFiltered', too_close_filtered, 'gRNAs for being within', _config.d.MIN_DIST - 1, 'bp from another gRNA'
    print '\tFiltered', num_repeats_found, 'gRNAs in repeats'
    if len(headers) < 12472:
      print '\tWARNING: Fewer than 12472 candidates found'



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