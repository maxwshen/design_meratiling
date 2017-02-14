

import _config
import sys, os, fnmatch, datetime, subprocess
import numpy as np

from mylib import util
from mylib import compbio
from mylib import plot

# Default params
DEFAULT_INP_DIR = _config.DATA_DIR
NAME = util.get_fn(__file__)

def search_region(nm, spc, chrm, startpos, endpos, RepeatMasker):
  startpos, endpos = int(startpos), int(endpos)
  sq = compbio.get_genomic_seq_twoBitToFa(spc, chrm, str(startpos), str(endpos))

  dists = []
  headers, sqs = [], []
  prev, too_close_filtered = 0, 0
  num_repeats_found = 0
  timer = util.Timer(total = len(sq))
  for j in range(len(sq) - 3):
    found = False
    if sq[j : j+2] == 'GG':
      found = True
      start, end = j - 21, j + 2
      cut_site = j - 4
      orient = '+'
    if sq[j : j+2] == 'CC':
      found = True
      start, end = j, j + 23 
      cut_site = j + 5
      orient = '-'
    if found:
      # filter grnas that are too close
      if cut_site - prev < _config.d.MIN_DIST:
        too_close_filtered += 1
        continue

      # filter incomplete grnas
      s = sq[start : end]
      if len(s) != 23:
        continue

      # filter repeats
      if RepeatMasker.search(chrm, startpos + j - 10, startpos + j + 10):
        num_repeats_found += 1
        continue

      # if all is ok
      if orient == '+':
        pass
      if orient == '-':
        start, end = end, start
        s = compbio.reverse_complement(s)

      # G-N19-NGG, G-N18-NGG, N20-NGG
      if s[0] != 'G' and s[1] == 'G':
        s = s[1:]
      elif s[0] != 'G' and s[1] != 'G':
        s = 'G' + s

      hdr = '>' + '__'.join([nm, chrm, str(startpos + start), str(startpos + end), str(startpos + cut_site), orient])
      if hdr not in headers:
        headers.append(hdr)
        sqs.append(s)
      dists.append(cut_site - prev)
      prev = cut_site

    timer.update()
  return headers, sqs, too_close_filtered, num_repeats_found, dists


def get_grnas(out_dir):
  for i in range(len(_config.d.NMS)):
    nm = _config.d.NMS[i]
    spc = _config.d.SPCS[i]
    chrm = _config.d.CHRMS[i]
    pos = _config.d.POS[i]
    chiapetfn = _config.d.CHIAPET[i]
    print '\t', nm

    RepeatMasker = compbio.RepeatMasker(spc)
    RepeatMasker.trim(chrm, pos[0] - 500000, pos[1] + 500000)

    print '\tRegion Length:', pos[1] - pos[0], 'bp'
    headers, sqs, too_close_filtered, num_repeats_found, dists = search_region(nm, spc, chrm, pos[0], pos[1], RepeatMasker)
    print '\tFound', len(headers), 'candidate gRNAs\n'

    chiapet_loci = []
    if chiapetfn != '':
      with open(chiapetfn) as f:
        for i, line in enumerate(f):
          site = int(line.strip())
          chiapet_loci.append( (site - _config.chiapetwindow, site + _config.chiapetwindow) )
    print '\tFound', len(chiapet_loci), 'chiapet sites'

    timer = util.Timer(total = len(chiapet_loci))
    for cpl in chiapet_loci:
      new_headers, new_sqs, tcf, nrf, new_dists = search_region(nm, spc, chrm, cpl[0], cpl[1], RepeatMasker)

      for (nh, ns, nd) in zip(new_headers, new_sqs, new_dists):
        if nh not in headers and nh + '_chia' not in headers:
          headers.append(nh + '_chia')
          sqs.append(ns)
          dists.append(nd)
      too_close_filtered += tcf
      num_repeats_found += nrf
      timer.update()

    # Ensure no duplicate gRNAs exist
    if len(set(headers)) != len(headers):
      print 'ERROR: Duplicate gRNAs found'
      import pdb; pdb.set_trace()
      sys.exit(0)

    # Analyze the distribution of distances
    print '\tMean Distance:', np.mean(dists), '\n\tStd Distance:', np.std(dists)
    plot.hist(dists, out_dir + nm + '_dist_hist.png')

    with open(out_dir + nm + '.fa', 'w') as f:
      for h, s in zip(headers, sqs):
        # s = 'ATATATCTTGTGGAAAGGACGAAACACC' + s[:-3] + 'GTTTAAGAGCTATGCTGGAAACAGCATAGC'
        f.write(h + '\n' + s[:-3] + '\n')
    
    print '\tFound', len(headers), 'candidate gRNAs with ChIA-PET\n'
    print '\tFiltered', too_close_filtered, 'gRNAs for being within', _config.d.MIN_DIST - 1, 'bp from another gRNA'
    print '\tFiltered', num_repeats_found, 'gRNAs in repeats'
    if len(headers) < 12270:
      print '\tWARNING: Fewer than 12270 candidates found'

    print 'Found', sum([1 for s in sqs if len(s) == 22]), '19-len gRNAs'
    print 'Found', sum([1 for s in sqs if len(s) == 23]), '20-len gRNAs'
    print 'Found', sum([1 for s in sqs if len(s) == 24]), '21-len gRNAs'


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