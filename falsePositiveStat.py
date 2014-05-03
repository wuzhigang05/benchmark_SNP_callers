#!/usr/bin/env python
from pandas import *
import numpy as np
import os
import pdb
import sys
import re
# to solve 'Invalid DISPLAY variable' error when run this script on biocluster
# http://bit.ly/1nkAwuk
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('PDF')

import matplotlib.pyplot as plt
from matplotlib.pyplot import setp
from matplotlib_venn import venn2, venn3
import brewer2mpl as b2m
options.display.mpl_style = 'default'

def getSNPs(file):
  oo = os.popen("""grep -v '#' %s | awk '{print $2}'""" % (file)).readlines()
  return set(map(lambda x: int(x.strip()), oo))

def simpleaxis(ax):
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  ax.spines['bottom'].set_visible(False)
  ax.get_xaxis().tick_bottom()
  ax.get_yaxis().tick_left()
   
def line_yielder(file):
  """
  generator function, which yields line by line from a file

  Parameters:
  ----------
  file: either a string of a file path or a `file`-like object

  Return:
  --------
  well, this function not actually return anything but `yield` line
  """
  if hasattr(file, "read"):
    for line in file:
      yield line
  else:
    with open(file) as OUT:
      for line in OUT:
        yield line

def getQualDp(file, positions, format):
  '''get the QUAL and Read Depth data for positins'''
  quals = []
  deps = [] 
  pat = re.compile('DP=(\d+)')
  if format in ('sam', 'gatk'):
    for line in line_yielder(file):
      if not line.startswith('#'):
        fds = line.split()
        pos = int(fds[1])
        if pos in positions:
          qual = np.round(float(fds[5]), 2)
          dep = int(pat.search(fds[7]).group(1))
          quals.extend((qual,))
          deps.extend((dep,)) 
  elif format in ('vart'): 
    for line in line_yielder(file):
      if not line.startswith('#'):
        fds = line.split()
        pos = int(fds[1])
        if pos in positions:
          qual = np.round(float(fds[-1].split(':')[7]), 2)
          dep = int(fds[-1].split(':')[1])
          quals.extend((qual,))
          deps.extend((dep,)) 
  else:
    e.write("unknown format\n")
    sys.exit(1)
  return np.array(quals), np.array(deps)   

def setBoxColors(bp):
  setp(bp['boxes'][0], color='blue')
  setp(bp['boxes'][1], color='red')
  setp(bp['caps'][0], color='blue')
  setp(bp['caps'][1], color='blue')
  setp(bp['caps'][2], color='red')
  setp(bp['caps'][3], color='red')
  setp(bp['whiskers'][0], color='blue')
  setp(bp['whiskers'][1], color='blue')
  setp(bp['whiskers'][2], color='red')
  setp(bp['whiskers'][3], color='red')
  setp(bp['fliers'][0], color='blue')
  setp(bp['fliers'][1], color='blue')
  setp(bp['fliers'][2], color='red')
  setp(bp['fliers'][3], color='red')
  setp(bp['medians'][0], color='blue')
  setp(bp['medians'][1], color='red')

if __name__ == '__main__':
  o = sys.stdout
  e = sys.stderr
  theoreticalSnps = set(np.load("mutatedPos.npy") + 1)
  allSnps = theoreticalSnps 
  files = ["sam.raw.vcf", "gatk.raw.vcf", "varianttools.raw.vcf"]
 
  # analyze the DP and QUAL distribution for false positives
  samRaw = getSNPs(files[0]) 
  gatkRaw = getSNPs(files[1])
  vartRaw = getSNPs(files[2])
  samfp = set.difference(samRaw, allSnps)
  gatkfp = set.difference(gatkRaw, allSnps)
  vartfp = set.difference(vartRaw, allSnps)
  samtp = set.intersection(samRaw, allSnps)
  gatktp = set.intersection(gatkRaw, allSnps)
  varttp = set.intersection(vartRaw, allSnps)

  samqualsFp, samdpsFp = getQualDp(files[0], samfp, 'sam')  
  gatkqualsFp, gatkdpsFp = getQualDp(files[1], gatkfp, 'gatk')  
  vartqualsFp, vartdpsFp = getQualDp(files[2], vartfp, 'vart')  
  
  samqualsTp, samdpsTp = getQualDp(files[0], samtp, 'sam')  
  gatkqualsTp, gatkdpsTp = getQualDp(files[1], gatktp, 'gatk')  
  vartqualsTp, vartdpsTp = getQualDp(files[2], varttp, 'vart')  


  fig, ax = plt.subplots(nrows=2, ncols = 3, sharex = True)
  datas = [
    [samqualsTp, samqualsFp],
    [gatkqualsTp, gatkqualsFp],
    [vartqualsTp, vartqualsFp],
    [samdpsTp, samdpsFp],
    [gatkdpsTp, gatkdpsFp],
    [vartdpsTp, vartdpsFp]]


  a = []
  a.extend(ax[0])
  a.extend(ax[1])

  for x, data in zip(a, datas):
    sa = Series(data[0])
    sb = Series(data[1])
    ss = concat((sa, sb), ignore_index=True, axis=1)
    ss.columns = ["TP", "FP"]
    bp = ss.boxplot(ax=x, sym="", grid=False, fontsize='small')
    setBoxColors(bp)
    simpleaxis(x)
# reduce yticks by half
    x.set_yticks(x.get_yticks()[::2])

  ax[0][0].set_ylabel("QUAL")
  ax[1][0].set_ylabel("Read Coverage")
  degree = 0
  ax[1][0].set_xlabel("Samtools", rotation = 0)
  ax[1][1].set_xlabel("GATK", rotation = 0)
  ax[1][2].set_xlabel("VariantTools", rotation = 0)
  fig.set_tight_layout(True)
  fig.savefig("qualDpDist.pdf")


