

import _config
import sys, os, fnmatch, datetime, subprocess
import numpy as np

from mylib import util
from mylib import compbio

# Default params
DEFAULT_INP_DIR = _config.OUT_PLACE + 'a_getgrnas/'
NAME = util.get_fn(__file__)

# Functions
def bowtie(inp_dir, out_dir):
  for i in range(len(_config.d.NMS)):
    nm, spc, chrm, pos = _config.d.NMS[i], _config.d.SPCS[i], _config.d.CHRMS[i], _config.d.POS[i]
    print '\t', nm

    inpfn = inp_dir + nm + '.fa'
    bt_tool = '/cluster/mshen/tools/bowtie2-2.2.9/bowtie2'

    if spc == 'mm10':
      bt_db = '/cluster/mshen/tools/bowtie2-2.2.9/db/' + spc
    if spc == 'hg38':
      bt_db = '/cluster/mshen/tools/bowtie2-2.2.9/db/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index'

    outfn = out_dir + nm + '.sam'
    print '\t\tRunning Bowtie2...'
    ans = subprocess.check_output(bt_tool + ' -f -k ' + str(_config.d.DUP_CUTOFF + 1) + ' -x ' + bt_db + ' -U ' + inpfn + ' -S ' + outfn,
      shell = True)
    print ans

    # Parse output SAM file from bowtie2
    cts = {}
    with open(outfn) as f:
      for i, line in enumerate(f):
        if line[0] != '@':
          w = line.split()
          if w[0] not in cts:
            cts[w[0]] = 0
          cts[w[0]] += 1

    # Write a new filtered gRNA candidate file
    headers, reads = compbio.read_fasta(inpfn)
    numskip = 0
    with open(out_dir + nm + '.fa', 'w') as f:
      for h, r in zip(headers, reads):
        skip = False
        if h in cts and cts[h] >= _config.d.DUP_CUTOFF:
          skip = True
          numskip += 1
        if not skip:
          f.write('>' + h + '\n' + 'ATATATCTTGTGGAAAGGACGAAACACC' + r + 'GTTTAAGAGCTATGCTGGAAACAGCATAGC\n')
    print '\tFiltered out', numskip, 'candidates occurring at least', _config.d.DUP_CUTOFF, 'times in', spc, 'genome'
    print '\tKept', len(headers) - numskip, 'candidate gRNAs'
    if len(headers) - numskip < 12472:
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
  bowtie(inp_dir, out_dir)

  return out_dir


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1], sys.argv[2])
  else:
    main(DEFAULT_INP_DIR, _config.OUT_PLACE + NAME + '/')
