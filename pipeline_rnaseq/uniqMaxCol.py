#!/usr/bin/env python2

# coding: utf-8

# In[142]:


import pandas as pd
def uniqMaxCol(df=None,fname=None,idCol=8, maxCol=13):
    ''' Keep the record among duplicates by maxCol
    '''
    if df is None:
        df = pd.read_table(fname,header=None,dtype=default_dtype)
    dup = df[idCol].duplicated(keep=False)
    ddup,duniq = df.loc[dup],df.loc[~dup]
    ddup[maxCol].astype('float16',copy=False)
#     ddup[maxCol]=ddup[maxCol].astype('float16')
    def maxBycol( val,col=maxCol):
        if len(val)>1:
            signal = val[col]
            ix = signal.idxmax()
            return val.loc[[ix]]
        else:
            return val
    ddedup = ddup.groupby(idCol,as_index=0).apply(maxBycol)
    dout = duniq.append(ddedup)
    dout.columns = default_header
    return dout

import sys
lst='''feat_chrom,str
feat_source,str
feat_type,str
feat_start,int
feat_end,int
feat_score,str
feat_strand,str
feat_phase,str
feat_attributes,str
peak_chrom,str
peak_start,int
peak_end,int
peak_id,str
peak_score,str
peak_strand,str
peak_fc,str
peak_logP,str
peak_logQ,str
peak_sumLoc,str'''.splitlines()
default_header,default_dtype = zip(*[x.split(',') for x in lst])

default_dtype={i:d for i,d in enumerate(default_dtype)}
if __name__=='__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Keep the record among duplicates by maxCol')
#     parser.add_argument('integers', metavar='N', type=int, nargs='+',
#                         help='an integer for the accumulator')
    parser.add_argument('fname', metavar='input_file', type=str, nargs=1,
                        help='Input file')
    parser.add_argument('ofname', metavar='output_file', type=str, nargs='?', default=None,
                        help='Output file')
    #     parser.add_argument('ofname', metavar='input_file', type=str, nargs='1',
#                         help='Input file')
    parser.add_argument('--header', dest='incHeader', action='store_true',
                        help='include the header',default=None)
    args = parser.parse_args()
#     print args
#     print args.ofname
#     fname = sys.argv[1]
#     ofname= sys.argv[2:]
    df = uniqMaxCol(fname=args.fname[0],idCol=8,maxCol=13)
    
#     saver = lambda df,*args: 
    df.to_csv(
        path_or_buf=args.ofname or sys.stdout,
        sep='\t',index=0,
        header=args.incHeader)  
#     saver(df,args.ofname or sys.stdout)
#         df


# In[ ]:


exit(0)


# In[ ]:


get_ipython().run_cell_magic(u'bash', u'', u'jupyter nbconvert --to python uniqMaxCol.ipynb\nsed -i "1i #!/usr/bin/env python2" *.py\nchmod +x uniqMaxCol.py\n\ncp *.py $repos/BrachyPhoton/pipeline_rnaseq')


# In[144]:


get_ipython().system(u' rm -f tmp_dedup.bed')
get_ipython().system(u' ./uniqMaxCol.py tmp.bed >tmp_dedup.bed')
get_ipython().system(u' head tmp_dedup.bed -n2')
get_ipython().system(u' echo \\\\n')
get_ipython().system(u' ./uniqMaxCol.py tmp.bed --header |head -n3')
get_ipython().system(u' echo \\\\n')
get_ipython().system(u' ./uniqMaxCol.py tmp.bed tmp_dedup.bed && head tmp_dedup.bed -n2')
# ! head tmp_dedup.bed


# In[54]:


df = uniqMaxCol(fname='tmp.bed')
df.head(5)
# print saver(df)[:300]
# df.to_csv('tmp_dedup.bed',sep='\t',index=0,header=None)
# !head tmp_dedup.bed


# In[32]:


get_ipython().run_cell_magic(u'bash', u'', u'curl http://gmod.org/wiki/GFF3 | grep -E "Column [0-9]:"')


# In[13]:


get_ipython().run_cell_magic(u'bash', u'', u'. /home/feng/envs/ref/bin/activate\n. config_default.sh\nsize2sum $GSIZE')


# In[67]:


get_ipython().run_cell_magic(u'bash', u'', u'. /home/feng/envs/ref/bin/activate\nmain(){\n    local SELF=`readlink -f ${BASH_SOURCE[0]}`        \n    \n    show_help(){\n    echo " No help for you yet :X\n    Options:\n    -d dry-run\n    "\n    \n    }\n    if [[ $# -eq 0 ]] ; then\n        show_help\n        exit 0\n    fi\n    \n    ###############\n    #### ==== Parsing argument\n    CTRL=""\n    GSUM=2.7E9\n    OPTIND=1         # Reset in case getopts has been used previously in the shell.\n    while getopts "h?c:t:g:d" opt; do\n        case "$opt" in\n        h|\\?)\n            show_help\n            exit 0\n            ;;\n        c)  CTRL=$OPTARG\n            ;;\n        t)  NCORE=$OPTARG\n            ;;\n        g)  GSUM=$OPTARG\n            ;;\n        d)  local DRY=1\n            ;;\n        esac\n    done\n    shift $((OPTIND-1))\n    [ "${1:-}" = "--" ] && shift\n    #### ---- Parsing argument\n    ###############\n    echo "[Positional]: $@"\n    \n    IN=$1 ; ALI=`bname $IN`      \n    QLIM=${3:-5E-2}\n    local GSUM=${4:-2.7E9}\n    \n    \n    \n    #### Step 1: MACS2\n    PROG="macs2 callpeak"\n    PROGALI=${PROG// /_}\n    OPT="-q $QLIM -g $GSUM --keep-dup 1"\n\n    ### Use a control if specified\n    if [[ ! -z "$CTRL" ]]; then\n        OPT="$OPT -c $CTRL"\n    fi\n    \n    CMD="$PROG $OPT -t $IN -n $ALI &>$ALI.$PROGALI.log"\n    echo $CMD\n    [[ ! -z $DRY ]] || eval $CMD    \n    \n    #### Step 2: \n}\nmain -d -c hi.bam test.bam\nmain test.bam')


# 
#     5th: integer score for display calculated as int(-10*log10qvalue). Please note that currently this value might be out of the [0-1000] range defined in UCSC Encode narrowPeak format<https://genome.ucsc.edu/FAQ/FAQformat.html#format12>
#     7th: fold-change
#     8th: -log10pvalue
#     9th: -log10qvalue
#     10th: relative summit position to peak start
# 

# In[218]:


get_ipython().run_cell_magic(u'bash', u'', u". /home/feng/envs/ref/bin/activate ; . config_default.sh \n# GREF=/home/feng/ref/Brachypodium_Bd21_v3-1/annotation/gene_only.gff\n# echo $GFF\n\nIN=NA_peaks.narrowPeak\n# echo $IN\n# type bname\nALI=`bname ${IN}`\n# exit 0\n\n#### Reference annotation\nGREF=/home/feng/ref/Brachypodium_Bd21_v3-1/annotation/Bdistachyon_314_v3.1.gene.gff3\n# bedtools intersect -a $\n# echo $GTF\n# ls /home/feng/ref/Brachypodium_Bd21_v3-1/annotation\nPAD=1000\ngrep gene $GREF | bedtools sort  >tmp.ref\n# bedtools slop -i - -g $GSIZE -b $PAD\n\nhead tmp -n10\necho\n\nbedtools slop -i tmp.ref -g $GSIZE -b $PAD >tmp.ref.s${PAD}\nbedtools closest -a tmp.ref.s${PAD} -b $IN >tmp.out  \n\n#### The output of 'bedtools closest' is not dircetly accepted by bedtools\ntabCut tmp.out -9  | bedtools slop -i - -g $GSIZE -b -${PAD} >tmp.out.-9\ntabCut tmp.out 10- | paste -d$'\\t' tmp.out.-9 - >tmp.bed\n\nhead tmp.out tmp.bed -n1\n\nln -f tmp.bed ${ALI}.bed\n\n\n#### Peaks are associated to each gene")


# In[308]:


df = uniqMaxCol(fname='tmp.bed',)


# In[ ]:


import pandas as pd


# In[216]:


# %%bash 
# head tmp.`bed
# bedtools merge -i tmp.bed


# In[195]:


get_ipython().run_cell_magic(u'bash', u'', u'head tmp.out -n1\nhead tmp.bed -n1')


# In[155]:


# %%bash
# echo =======
# head tmp.ref
# echo =======
# bedtools slop -i tmp.ref.s1000 -g $GSIZE -b -1000 | head 
# # head tmp.out


# In[108]:


get_ipython().system(u'head -n 100 /home/feng/ref/Brachypodium_Bd21_v3-1/annotation/gene_all.gff')


# In[112]:


get_ipython().system(u' head -n 100 /home/feng/ref/Brachypodium_Bd21_v3-1/annotation/Bdistachyon_314_v3.1.gene.gff3')


# In[105]:


get_ipython().run_cell_magic(u'bash', u'', u'ln -s /home/feng/ref/Brachypodium_Bd21_v3-1/transcriptome/gene_index.gff /home/feng/ref/Brachypodium_Bd21_v3-1/annotation/gene_all.gff')


# In[101]:


get_ipython().run_cell_magic(u'bash', u'', u'cat /home/feng/ref/Brachypodium_Bd21_v3-1/transcriptome/gene_index.gff | grep gene >/home/feng/ref/Brachypodium_Bd21_v3-1/annotation/gene_only.gff')


# In[73]:


df[2]-df[1]


# https://github.com/taoliu/MACS

# In[9]:


get_ipython().system(u' cat HELP')

