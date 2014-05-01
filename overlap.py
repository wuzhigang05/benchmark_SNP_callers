#!/usr/bin/env python
from pandas import *
import numpy as np
import os
import pdb
import sys
# to solve 'Invalid DISPLAY variable' error
# http://bit.ly/1nkAwuk
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('PDF')

import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import brewer2mpl as b2m

def sensitivity(found, total):
    try:
        assert isinstance(found, set)
    except:
        e.write("found has to be a set\n")
        sys.exit(1)
    try:
        assert isinstance(total, set)
    except:
        e.write("total has to be a set\n")
        sys.exit(1)
    return np.round(len(found & total)/float(len(total)), 2)

def fdr(found, total):
    """
    http://en.wikipedia.org/wiki/Sensitivity_and_specificity
    Actually, this function returns False Discoverate Rate 
    """ 
    try:
        assert isinstance(found, set)
    except:
        e.write("found has to be a set\n")
        sys.exit(1)
    try:
        assert isinstance(total, set)
    except:
        e.write("total has to be a set\n")
        sys.exit(1)
    fp = len(set.difference(found, total))
    return np.round(fp/float(len(found)), 5)

def getSNPs(file):
  oo = os.popen("""grep -v '#' %s | awk '{print $2}'""" % (file)).readlines()
  return set(map(lambda x: int(x.strip()), oo))

def calculate(files, P):
  snps = map(getSNPs, files)
  sensitivities = Series(map(lambda x: sensitivity(x, P), snps), index = files, name="TPR")
  specificities = Series(map(lambda x: fdr(x, P), snps), index = files, name="FDR")
  result = concat([sensitivities, specificities], axis=1)
  print result
  return result

def venn3Diagram(A, B, C, labels, figName):
    try:
        for a in [A, B, C]:
            assert isinstance(A, set)
    except:
        e.write("The first three arguments have to be a type of 'set'\n")
        sys.exit(1)
    try: 
        assert len(labels) == 3
    except:
        e.write("labels has to have a length of 3\n")
        sys.exit(1)
    try:
        assert figName.split('.')[-1] in ['pdf', 'png', 'svg', 'jpg']
    except:
        e.write("%s has to have a suffix of ['pdf', 'png', 'svg', 'jpg']\n" % figName)
        sys.exit(1)

    fig, ax = plt.subplots()
    v = venn3(subsets=[A, B, C], set_labels=labels)
    for l in v.set_labels:
        l.set_size(18)
    fig.savefig(figName)

def TPRandFDR(file, theoreticalSnps):
    snps = getSNPs(file)
    TP = snps & theoreticalSnps
    FP = snps - TP
    print "#TP: ", len(TP)
    print "TPR: ", float(len(TP))/len(theoreticalSnps)
    print "#FP: ", len(FP)
    print "FDR: ", float(len(FP))/len(snps)
    
def simpleaxis(ax):
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  ax.get_xaxis().tick_bottom()
  ax.get_yaxis().tick_left()

if __name__ == '__main__':
  o = sys.stdout
  e = sys.stderr
  theoreticalSnps = set(np.load("mutatedPos.npy") + 1)
  allSnps = theoreticalSnps 
  files = ["sam.raw.vcf", "gatk.raw.vcf", "varianttools.raw.vcf"]
  newFiles = []
  for f in files:
      tmp = [f]
      for i in range(1, 10):
        tmp.append(".".join(f.split('.')[:2] + [str(i), 'vcf']))
      newFiles += tmp 

  dpFiles = []
  for f in files:
    tmp = []
    for i in range(3, 13):
      tmp.append(".".join(f.split('.')[:2] + ['dp' + str(i), 'vcf']))
    dpFiles += tmp

  for i in range(10):
      venn3Diagram(getSNPs(newFiles[i]), getSNPs(newFiles[i+10]), getSNPs(newFiles[i+20]), 
              #[newFiles[i], newFiles[i+10], newFiles[i+20]], "venn_%d.pdf" % i)
              ["Samtools", "GATK", "VariantTools"], "venn_%d.pdf" % i)
      
      venn3Diagram(getSNPs(dpFiles[i]), getSNPs(dpFiles[i+10]), getSNPs(dpFiles[i+20]), 
              ["Samtools", "GATK", "VariantTools"], "venn_dp%d.pdf" % i)

  sams = calculate(newFiles[:10], allSnps)
  gatks = calculate(newFiles[10:20], allSnps)
  varianttools = calculate(newFiles[20:], allSnps)
  
  fig, ax = plt.subplots()
  # ax.plot(sams['FDR'], sams['TPR'], 'r-^', label='samtools/bcftools')
  # ax.plot(gatks['FDR'], gatks['TPR'], 'b-o', label='gatkUG')
  # ax.plot(varianttools['FDR'], varianttools['TPR'], 'g-s', label='Varianttools')

  #colors = b2m.get_map('RdYlGn', 'Diverging', 3).mpl_colors
  colors = ["red", "green", "blue"]
  ax.plot(sams['TPR'], sams['FDR'], 'r-^', label='Samtools', color=colors[0])
  ax.plot(gatks['TPR'], gatks['FDR'], 'b-o', label='GATKUG', color=colors[1])
  ax.plot(varianttools['TPR'], varianttools['FDR'], 'g-s', label='Varianttools', color=colors[2])

  ax.set_ylabel("False Discovery Rate")
  ax.set_xlabel("True Positive Rate")
  ax.set_title("Filtration by QUAL")
  ax.legend(loc=0, frameon=False, fontsize='small')
  simpleaxis(ax)
  fig.savefig("plot.pdf")

  sams = calculate(dpFiles[:10], allSnps)
  gatks = calculate(dpFiles[10:20], allSnps)
  varianttools = calculate(dpFiles[20:], allSnps)
  
  fig, ax = plt.subplots()
  ax.plot(sams['TPR'], sams['FDR'], 'r-^', label='Samtools', color=colors[0])
  ax.plot(gatks['TPR'], gatks['FDR'], 'b-o', label='GATK', color=colors[1])
  ax.plot(varianttools['TPR'], varianttools['FDR'], 'g-s', label='VariantTools', color=colors[2])
  ax.set_ylabel("False Discovery Rate")
  ax.set_xlabel("True Positive Rate")
  ax.set_title("Filtration by Read Coverage")
  ax.legend(loc=0, frameon=False, fontsize='small')
  simpleaxis(ax)
  fig.savefig("tprDP.pdf")
