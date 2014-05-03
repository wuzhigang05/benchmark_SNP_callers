#!/usr/bin/env python
from pandas import *
import numpy as np
import os
import argparse
import pdb
import sys
import collections
# to solve 'Invalid DISPLAY variable' error
# http://bit.ly/1nkAwuk
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('PDF')

import matplotlib.pyplot as plt
import prettyplotlib as ppl
import brewer2mpl as b2m

def getQuals(file):
  oo = os.popen("""grep -v '#' %s | awk '{print $6}'""" % (file)).readlines()
  return map(lambda x: np.round(float(x.strip()), 2), oo)

def getQualFromVariantToolsOut(file):
  oo = os.popen("""grep -v '#' %s | awk '{print $10}' | awk -F':' '{print $8}' """ % (file)).readlines()
  return map(lambda x: np.round(float(x.strip()), 2), oo)

def getDPFromVariantToolsOut(file):
  oo = os.popen("""grep -v '#' %s | awk '{print $10}' | awk -F':' '{print $2}' """ % (file)).readlines()
  return map(lambda x: np.round(float(x.strip()), 2), oo)

def getOtherTag(file, tag):
  oo = os.popen("""grep '#' -v %s | awk -F'%s' '{print $2}' | awk -F'[=;]' '{print $2}'""" % (file, tag)).readlines()
  return map(lambda x: np.round(float(x.strip()), 2), oo)

def simpleaxis(ax):
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  ax.get_xaxis().tick_bottom()
  ax.get_yaxis().tick_left()

def boxplot(ax, data, xticklabels, ylabel):
  if isinstance(data, collections.Iterable):
    ppl.boxplot(ax, data, notch=True, xticklabels=xticklabels, sym="", vert=False)
    ax.set_ylabel(ylabel)
  else:
    sys.stderr.write("data has to be iterable\n")
    sys.exit(1)
  return ax

def hist(ax, data, bins, labels=['Samtools', 'GATK', 'VariantTools'], legend=True):
  assert isinstance(data, collections.Iterable)
  colors = b2m.get_map('RdYlGn', 'Diverging', 11 if len(data) > 11 else len(data) ).mpl_colors
  ax.hist(data, bins, color=colors, label= labels)
  simpleaxis(ax)
  if legend:
    ax.legend(loc=0, frameon=False, fontsize='small')
  return ax

if __name__ == '__main__':
  o = sys.stdout
  e = sys.stderr
  parser= argparse.ArgumentParser(description="")
  parser.add_argument("sam", help="samtools raw vcf output")
  parser.add_argument("gatk", help="gatk raw vcf output")
  parser.add_argument("vart", help="varianttools raw vcf output")
  #parser.add_argument("vatt", help="varianttools raw vcf output")
  args = parser.parse_args() 

  samquals = np.array(getQuals(args.sam))
  gatkquals = np.array(getQuals(args.gatk))
  vartquals = np.array(getQualFromVariantToolsOut(args.vart))

  samDP = np.array(getOtherTag(args.sam, 'DP'))
  gatkDP = np.array(getOtherTag(args.gatk, 'DP'))
  vartDP = np.array(getDPFromVariantToolsOut(args.vart))

  #fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True)
  #boxplot(ax[0], (samquals, gatkquals, vartquals), xticklabels=['Samtools', 'GATK', 'VariantTools'], ylabel="QUAL")
  #boxplot(ax[1], (samDP, gatkDP, vartDP), xticklabels=['Samtools', 'GATK', 'VariantTools'], ylabel="Read Coverage")

  fig, ax = plt.subplots(nrows=2, ncols=1)
  #ppl.hist(ax[0], samquals, 50, ylabel="QUAL")
  #boxplot(ax[1], (samDP, gatkDP, vartDP), xticklabels=['Samtools', 'GATK', 'VariantTools'], ylabel="Read Coverage")
  #boxplot(ax[2], (samDP, gatkDP, vartDP), xticklabels=['Samtools', 'GATK', 'VariantTools'], ylabel="Read Coverage")
  #ax.hist((samquals, gatkquals, vartquals), bins=20, color=b2m.get_map('RdYlGn', 'Diverging', 3).mpl_colors, label=['Samtools', 'GATK', 'VariantTools'])
  #simpleaxis(ax)
  ##ax.set_title("QUAL distribution")
  #ax.legend(loc=0, frameon=False)

  ax[0] = hist(ax[0], (samquals, gatkquals, vartquals), 20)
  ax[0].set_title("QUAL")
  ax[0].set_ylabel("Frequency")
  ax[1] = hist(ax[1], (samDP, gatkDP, vartDP), 20, legend=False)
  ax[1].set_title("Read Coverage")
  ax[1].set_ylabel("Frequency")
  fig.savefig("qualDistribution.pdf")

