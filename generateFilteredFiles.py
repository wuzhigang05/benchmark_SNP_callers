#!/usr/bin/env python2.7
import sys
import pdb
import os
import argparse
from pandas import *
import re
import tempfile
import numpy as np

def line_yielder(file):
  if hasattr(file, "read"):
    for line in file:
      yield line
    else:
      with open(file) as IN:
        for line in IN:
          yield line

def callSys(cmd):
  return os.popen(cmd).readlines()

def getQualStat(vcfFile):
  return callSys("""grep '#' -v %s | tabtk num -c6 -Q""" % vcf)

def getOtherStat(vcfFile, tag):
  return callSys("""grep '#' -v %s |  awk -F'%s' '{print $2}' | awk -F'[=;]' '{print $2}' | tabtk num -c1 -Q""" % (vcfFile, tag))

def parsingFile(file):
  quals = []
  dp = []
  mq = []
  dpR = re.compile('DP=(\d+\.?\d*)')
  mqR = re.compile('MQ=(\d+\.?\d*)')
  for line in line_yielder(file):
    if not line.startswith('#'):
      f = line.split()
      quals.append(f[5])
      dp.append(dpR.search(f[7]).group(1))
      mq.append(mqR.search(f[7]).group(1))
  return concat([Series(np.array(quals).astype(float), name="QUAL"),
    Series(np.array(dp).astype(float), name="DP"),
    Series(np.array(mq).astype(float), name = "MQ")], axis=1)

if __name__ == '__main__':
  o = sys.stdout
  e = sys.stderr
  parser= argparse.ArgumentParser(description="")

  parser.add_argument("--vcfFiles", nargs= '+', 
      default = ["sam.raw.vcf", "gatk.raw.vcf"],
      help="the vcfFiles you want to know the statistics. " +
      "default = [sam.raw.vcf, gatk.raw.vcf]")
  parser.add_argument("--tags", default = "QUAL",
      help="the specifc tags, for which you want to filter on. defualt= QUAL")
  parser.add_argument("--ref", default= 'Chr4.fasta', 
      help="reference fasta file default = [Chr4.fasta]")
  args = parser.parse_args() 
  print args.vcfFiles
  results = []
  j = 0
  if args.tags == "QUAL":
    for f in args.vcfFiles:
      for i, q in enumerate(np.linspace(.1, .9, 9)):
        qualCutoff = os.popen('''grep '#' -v %s | tabtk num -c6 -q %f''' % (f, q)).readline().strip()
        cmd1 = 'java -Xmx4g -jar GenomeAnalysisTK.jar -T VariantFiltration -R %s -V %s ' % (args.ref, f)
        cmd4 = '--filterExpression "QUAL < %s " --filterName "lowQUAL" ' %  qualCutoff
        tmp = tempfile.NamedTemporaryFile()
        try:
          cmd5 = '-o %s' % tmp.name
          outputName = ".".join(f.split('.')[:2] + [str(i+1), 'vcf'])
          j += 1
          print j
          selectcmd = "java -Xmx4g -jar GenomeAnalysisTK.jar -R %s -T SelectVariants --variant %s -o %s -select 'vc.isNotFiltered()' " % (args.ref, tmp.name, outputName)
          #print cmd1 + cmd4 + cmd5
          #print selectcmd
          callSys(cmd1 + cmd4 + cmd5)
          callSys(selectcmd)
        finally:
          tmp.close()
  elif args.tags == "DP":
    for f in args.vcfFiles:
       #for i, q in enumerate(np.linspace(.1, .9, 9)):
       for i, q in enumerate(np.linspace(3, 12, 10)): 
        #qualCutoff = os.popen('''grep '#' -v %s | awk -F'DP=' '{print $2}' | awk -F';' '{print $1}' | tabtk num -q %f''' % (f, q)).readline().strip()
        qualCutoff = str(q)
        cmd1 = 'java -Xmx4g -jar GenomeAnalysisTK.jar -T VariantFiltration -R %s -V %s ' % (args.ref, f)
        cmd4 = '--filterExpression "DP < %s " --filterName "lowDP" ' %  qualCutoff
        tmp = tempfile.NamedTemporaryFile()
        try:
          cmd5 = '-o %s' % tmp.name
          outputName = ".".join(f.split('.')[:2] + ['dp' + str(i+3), 'vcf'])
          j += 1
          print j
          selectcmd = "java -Xmx4g -jar GenomeAnalysisTK.jar -R %s -T SelectVariants --variant %s -o %s -select 'vc.isNotFiltered()' " % (args.ref, tmp.name, outputName)
          #print cmd1 + cmd4 + cmd5
          #print selectcmd
          callSys(cmd1 + cmd4 + cmd5)
          callSys(selectcmd)
        finally:
          tmp.close()
  else:
    e.write("Uknown Tag Type, only support ['QUAL', 'DP']\n")
    sys.exit(1)
