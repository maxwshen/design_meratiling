# trim local candidiates and add ChIA-PET linked gRNAs

import _config
import sys, os, fnmatch, datetime, subprocess
import numpy as np

import a_getgrnas
from mylib import util
from mylib import compbio

# Default params
DEFAULT_INP_DIR = _config.OUT_PLACE + 'b_bowtie/'
NAME = util.get_fn(__file__)

def trim_and_add_chiapet(name, inp_dir, out_dir):
  for i in range(len(_config.d.NMS)):
    nm, spc, chrm, pos = _config.d.NMS[i], _config.d.SPCS[i], _config.d.CHRMS[i], _config.d.POS[i]
    print '\t', nm
    if name != 'ignore' and nm != name:
      continue

    inpfn = inp_dir + nm + '.fa'

    headers, reads = compbio.read_fasta(inpfn)
    print '\t\t', len(headers)

    if len(headers) < 12270:
      remove = []
      pass
    else:   # need to remove down to 12270
      dists = []
      for j in range(len(headers) - 1):
        w1 = headers[j].split('_')
        start1 = int(w1[5])
        w2 = headers[j+1].split('_')
        start2 = int(w2[5])
        dists.append(abs(start2 - start1))

      indices = range(len(dists))
      indices = [x for (y, x) in sorted(zip(dists, indices))]

      remove = [indices[0]]
      for ind in indices[1:]:
        skip = False
        for rm in remove:
          if abs(ind - rm) == 1:
            skip = True
        if skip:
          continue
        else:
          remove.append(ind)
        if len(remove) == len(dists) - 12269:
          break

      print '\t\tReached index', indices.index(ind), 'for', len(dists) - 12269, 'removal candidates'

    new_h, new_r = [], []
    for j in range(len(headers)):
      if j not in remove:
        new_h.append(headers[j])
        new_r.append(reads[j])
    print '\t\tFiltered', len(remove), 'reads'

    # Add positive controls
    pos_control_heads, pos_control_reads = add_positive_controls()
    new_h += pos_control_heads
    new_r += pos_control_reads

    # Add negative controls
    aret = search_region(nm, spc, chrm, pos[0] - 1e6, pos[1] + 1e6)

    negch, negcr = add_negative_controls(aret)
    new_h += negch
    new_r += negcr

    # Write out
    write_files(new_h, new_r, out_dir + nm + '.fa')
  
  # End experiment loop


  return

def add_positive_controls():
  reads, headers = [], []
  with open(_config.DATA_DIR + 'egfp_nomm2_shuf70_pluscont.txt') as f:
    for i, line in enumerate(f):
      reads.append('ATATATCTTGTGGAAAGGACGAAACACC' + line.strip() + 'GTTTAAGAGCTATGCTGGAAACAGCATAGC')
      headers.append('pos_control_' + str(i))
  return headers, reads

def add_negative_controls(aret):
  # ensure gfp is not included
  gfpgrnas = []
  with open(_config.DATA_DIR + 'egfp_NGG_NNG_seq.patman_format.txt') as f:
    for i, line in enumerate(f):
      if line[0] != '>':
        gfpgrnas.append( line.strip() )

  headers, sqs = aret
  sqs += gfpgrnas
  sqs = [s[-15:-3] for s in sqs]  # use 12-mer seed region
  sqs = set(sqs)

  new_h, reads = [], []
  timer = util.Timer(total = 125)
  for i in range(125):
    ok = False
    while True:
      skip = False
      cand = ''.join(np.random.choice(['A', 'C', 'G', 'T'], 12))
      for s in sqs:
        if mismatch(cand, s) < 3:
          skip = True
          break
      if skip:
        continue
      # if ok... proceed
      reads.append('ATATATCTTGTGGAAAGGACGAAACACC' + ''.join(np.random.choice(['A', 'C', 'G', 'T'], 8)) + cand + 'GTTTAAGAGCTATGCTGGAAACAGCATAGC')
      new_h.append('neg_control_' + str(i))
      break
    timer.update()

  return new_h, reads

def mismatch(a, b):
  return sum([1 for i in range(len(a)) if a[i] != b[i]])

def search_region(nm, spc, chrm, startpos, endpos):
  startpos, endpos = int(startpos), int(endpos)
  sq = compbio.get_genomic_seq_twoBitToFa(spc, chrm, str(startpos), str(endpos))

  headers, sqs = [], []
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
      s = sq[start : end]
      if len(s) != 23:
        continue      
      if orient == '+':
        pass
      if orient == '-':
        start, end = end, start
        s = compbio.reverse_complement(s)
      if s[0] != 'G' and s[1] == 'G':
        s = s[1:]
      elif s[0] != 'G' and s[1] != 'G':
        s = 'G' + s
      hdr = '>' + '__'.join([nm, chrm, str(startpos + start), str(startpos + end), str(startpos + cut_site), orient])
      headers.append(hdr)
      sqs.append(s)
    timer.update()
  return headers, sqs

def write_files(headers, reads, out_fn):
  with open(out_fn, 'w') as f:
    for j in range(len(headers)):
      f.write('>' + headers[j] + '\n' + reads[j] + '\n')
  return

@util.time_dec
def main(inp_dir, out_dir, run = True):
  print NAME  
  util.ensure_dir_exists(out_dir)
  if not run:
    print '\tskipped'
    return out_dir

  # Function calls
  # name = 'apobec3_mm10'
  # name = 'ignore'
  name = 'apobec3b_hg38'
  trim_and_add_chiapet(name, inp_dir, out_dir)

  return out_dir


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1], sys.argv[2])
  else:
    main(DEFAULT_INP_DIR, _config.OUT_PLACE + NAME + '/')