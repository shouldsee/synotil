#!/usr/bin/env python2
import synotil.util as sutil
import pymisca.util as pyutil
import argparse
import synotil.ptn as sptn
parser = argparse.ArgumentParser()
parser.add_argument('fnames',nargs='+')
args = parser.parse_args()
fnames = args.fnames
print fnames
import sys

vals = map(sptn.getBrachy,fnames)
bnames = map(lambda x: x.split('/')[-1].rsplit('.',2)[0].replace('_','-'), fnames )
for i,d in enumerate(vals):
	d['Alias']=bnames[i]

flats = pyutil.meta2flat([ [ (x,d.pop(x)) for x in ['RunID','SampleID','Alias']] 
	for d in vals])
nonEmptyVals = [{k:v for k,v in d.items() if v is not None} for d in vals]
conds = pyutil.meta2flat(nonEmptyVals)

colNames = pyutil.paste0([flats,conds])
print colNames
#print flats,conds
#print bnames
#print vals
#sys.exit(0)
def callback(df):
	df = pyutil.filterMatch(df,key='STRG',negate=1)
	return df
df = pyutil.Table2Mat(fnames,callback=callback, valCol='TPM',idCol='Gene ID',ext='tsv')
df.columns = colNames
print df.to_csv()
 


