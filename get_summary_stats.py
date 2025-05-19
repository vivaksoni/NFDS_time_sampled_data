import libsequence
import sys
import pandas as pd
import math
import argparse
import vcf
import numpy as np

#parsing user given constants
parser = argparse.ArgumentParser(description='Information about length of region and sample size')
parser.add_argument('-inFile', dest = 'inFile', action='store', nargs = 1, type = str, help = 'vcf file')
parser.add_argument('-outFile', dest = 'outFile', action='store', nargs = 1, type = str, help = 'output file name')
#parser.add_argument('-samples', dest = 'samples', action='store', nargs = 1, type = int, help = 'no. of samples required')
parser.add_argument('-win_size', dest = 'win_size', action='store', nargs = 1, type = int, help = 'size of sliding windows')
parser.add_argument('-step_size', dest = 'step_size', action='store', nargs = 1, type = int, help = 'step size for sliding windows')
parser.add_argument('-regionLen', dest = 'regionLen', action='store', nargs = 1, type = int, help = 'size of simulated region')

args = parser.parse_args()
#samples = args.samples[0]
inFile = args.inFile[0]
outFile = args.outFile[0]
chr_len =  args.regionLen[0]
win_size = args.win_size[0]/float(chr_len)
step_size = args.step_size[0]/float(chr_len)


#Get chimeric site counts from vcf file
def get_genotype(vcf_file, chr_len):
    vcf_reader = vcf.Reader(filename=vcf_file)
    l_data = []
    for i,record in enumerate(vcf_reader):
        #create list of genotypes
        lst = (''.join([x['GT'] for x in record.samples]).replace('|',''))
        l_data.append([record.POS/chr_len, ''.join([x['GT'] for x in record.samples]).replace('|','')])
    return(l_data)


#Function to create dictionary of polySIM summary statistics. Returns dictionary of summary stats
def get_polySIM_stats(sd):
    #Create polysim object
    ps = libsequence.PolySIM(sd)
    #Create list of methods (ie polySIM summaryStats)
    a = [method for method in dir(ps) if callable(getattr(ps, method)) if not method.startswith('_')]
    #Loop through methods, storing names as keys, and estimates as values in dict
    ss_dict = {}
    for method in a:
        ss_dict[method] = getattr(ps, method)()
        
    return(ss_dict)


#Function to create dictionary of LD stats. Returns dictionary.
#Function to create dictionary of LD stats. Returns dictionary.
def get_LD_stats(sd):
    ld = libsequence.ld(sd)
    df = pd.DataFrame(ld)
    ss_dict = {}
    try:
        ss_dict['meanrsq'] = sum(df['rsq'])/len(df['rsq'])
        ss_dict['meanD'] = sum(df['D'])/len(df['D'])
        ss_dict['meanDprime'] = sum(df['Dprime'])/len(df['Dprime'])
    except Exception:
        ss_dict['meanrsq'] = 'NA'
        ss_dict['meanD'] = 'NA'
        ss_dict['meanDprime'] = 'NA'
    return(ss_dict)

df = pd.DataFrame()
vcf_file=inFile
l_data = get_genotype(vcf_file, chr_len)
sd = libsequence.SimData(l_data)

#define sliding windows
wins = libsequence.Windows(sd,window_size=win_size,step_len=step_size,starting_pos=0.0,ending_pos=1.0)
num_wins = len(wins)
d = {}
#Loop through windows, calculating stats for each window
for k, win in enumerate(wins):
    d['window'] = k
    d = {**d, **get_polySIM_stats(win)}
    if len(win.pos()) >= 5: #LD stats are pairwise. If only 1 site exists, it'll show an error.
        d = {**d, **get_LD_stats(win)}
    else:
        d['meanrsq'] = 'NA'
        d['meanD'] = 'NA'
        d['meanDprime'] = 'NA'

    df2 = pd.DataFrame.from_dict(d, orient='index').T
    df = pd.concat([df, df2])

df.to_csv(outFile, header=True, index=False, sep='\t')
