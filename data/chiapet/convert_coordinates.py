import sys, os, fnmatch, subprocess

tool = '/cluster/mshen/tools/liftOver/liftOver'
chain_fn = '/cluster/mshen/tools/liftOver/hg19ToHg38.over.chain'

curr_dir = '/cluster/mshen/prj/design_meratiling/data/chiapet/'

for fn in os.listdir(curr_dir):
  if '.bed' in fn and 'hg38' not in fn and 'leftover' not in fn:
    print fn
    inp_fn = curr_dir + fn
    out_fn = curr_dir + fn.strip('.bed') + '_hg38.bed'
    leftover_fn = curr_dir + fn.strip('.bed') + '_leftover.bed'
    subprocess.call(tool + ' ' + inp_fn + ' ' + chain_fn + ' ' + out_fn + ' ' + leftover_fn, shell = True)