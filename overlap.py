#!/usr/bin/env python
from pandas import *
import numpy as np
import os
import pdb
import sys
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3

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
    # for l in labels:
    #     v.get_label_by_id(l)
    fig.savefig(figName)

if __name__ == '__main__':
  #rawSnps = getSNPs("raw.vcf")
  #rawSnps = set(np.load("mutatedPos.npy"))
  #theoreticalSnps = getSNPs("vcf_chr_4.vcf")
  theoreticalSnps = set(np.load("mutatedPos.npy") + 1)
  o = sys.stdout
  e = sys.stderr
  # total possible SNPs foundable by SNP caller is the intersection of dbSNPs
  # and the positions  reporting SNP with at least 1X read coverage
  #allSnps = rawSnps & theoreticalSnps 
  allSnps = theoreticalSnps 
  # in order to calculate the false positive rate, we need to know the 
  # true negative rate, which is for all positions with read coverage but 
  # are not overallping with know dbSNPs
  #trueNegative = set.difference(rawSnps, theoreticalSnps)
  files = ["sam.raw.vcf", "gatk.raw.vcf", "varianttools.raw.vcf"]
      #"gatk_haplotypeCaller.vcf", "gatk_unifiedCaller_realigned.vcf"]
      #"gatk_haplotypeCaller.vcf", "snp.vcf"]
  results = []
  #for q in [0, 50, 100, 150, 200, 220]:
  newFiles = []
  for f in files:
      tmp = [f]
      for i in range(1, 4):
        tmp.append(".".join(f.split('.')[:2] + [str(i), 'vcf']))
      newFiles += tmp 
  for i in range(4):
      venn3Diagram(getSNPs(newFiles[i]), getSNPs(newFiles[i+4]), getSNPs(newFiles[i+8]), 
              [newFiles[i], newFiles[i+4], newFiles[i+8]], "venn_%d.pdf" % i)
  sams = calculate(newFiles[:4], allSnps)
  gatks = calculate(newFiles[4:8], allSnps)
  varianttools = calculate(newFiles[8:], allSnps)
  
  fig, ax = plt.subplots()
  ax.plot(sams['FDR'], sams['TPR'], 'r-^', label='samtools/bcftools')
  ax.plot(gatks['FDR'], gatks['TPR'], 'b-o', label='gatkUG')
  ax.plot(varianttools['FDR'], varianttools['TPR'], 'g-s', label='Varianttools')
  ax.set_xlabel("False Discovery Rate")
  ax.set_ylabel("True Positive Rate")
  ax.legend(loc=0)
  fig.savefig("plot.pdf")
  
