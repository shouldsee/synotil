#!/usr/bin/env python
# Usage: python chip-summary.py summary.txt
#
# Purpose: find target genes for many ChIP-seq samples.
#
# Note: summary.txt specifies BAM files. See QUICKSTART.txt to get started.
#
# Last modified 10 Oct 2016, slcu, hui lan
# Last modified 14 July 2017, slcu, hui
# Jun 2018, Feng: Refactored to take *.bed as input.

import os, sys, subprocess, operator, glob, shutil
import pymisca.util as pyutil


def delete_file(fname):
    if os.path.exists(fname):
        os.remove(fname)

def delete_files(s):
    for file in glob.glob(s):
        delete_file(file)


def move_files(s, dst):
    for file in glob.glob(s):
        if not os.path.exists(dst.rstrip('/') + '/' + file):
            shutil.move(file, dst)
        else:
            os.remove(file)

def calc_npeak(lst, fname, pvalue, qvalue):

    n = []

    for x in lst:
        # print('...calculate number of peaks for various thresholds')
        parameter_string = 'fc=%g p=%g q=%g' % (x, pvalue, qvalue)        
        temp_small_narrowPeak = 'temp_small.narrowPeak'
#         df = pd.read_table(fname)
#         df.query()
        cmd = 'python ' + PEAK_SELECT_SCRIPT + ' ' + fname + ' ' + parameter_string + ' > ' + temp_small_narrowPeak
        # print(cmd)
        os.system(cmd)
        with open(temp_small_narrowPeak) as f:
            n.append(len(f.readlines()))

    return n
        
def plot_npeak_vs_fc(lstx, lsty, fname):

    n = len(lstx)
    max_npeak = max(lsty)
    interval = max(1,max_npeak) / 120.0
    f = open(fname, 'w')
    f.write('Fold-change (fc) versus Number of Peaks plot\n\n')
    for i in range(n):
        f.write('fc=%3.1f |' % (lstx[i]))
        npeaks = lsty[i]
        f.write('|'*int(1.0*npeaks/interval))
        f.write(' %d\n' % npeaks)
    f.close()


def get_go(fname, type):

    d = {}
    f = open(fname)
    for line in f:
        if line.startswith('GO:'):
            lst = line.split('\t')
            go_name = lst[0]
            go_description = lst[3]
            if lst[1] == type.upper():
                d[go_name] = go_description
    f.close()
    return d


def make_go_string(lst, d):

    s = ''
    for x in lst:
        s += x + '\t' + d[x] + '\n'
    return s


def go_diff3(f1, f2, f3, type):

    d1 = get_go(f1, type)
    d2 = get_go(f2, type)
    d3 = get_go(f3, type)
    uniq1 = list(set(d1.keys()) - (set(d2.keys()) | set(d3.keys())))
    uniq2 = list(set(d2.keys()) - (set(d1.keys()) | set(d3.keys())))
    uniq3 = list(set(d3.keys()) - (set(d1.keys()) | set(d2.keys())))
    s1 = make_go_string(uniq1, d1)
    s2 = make_go_string(uniq2, d2)
    s3 = make_go_string(uniq3, d3)
    return [s1, s2, s3]

def go_diff3_bpmfcc(f1, f2, f3):

    bp = go_diff3(f1, f2, f3, 'bp')
    mf = go_diff3(f1, f2, f3, 'mf')
    cc = go_diff3(f1, f2, f3, 'cc')    

    result = ['', '', '']
    for i in range(len(bp)):
        result[i] += 'BP:\n' +  bp[i]

    for i in range(len(mf)):
        result[i] += '\n'
        result[i] += 'MF:\n' +  mf[i]        

    for i in range(len(cc)):
        result[i] += '\n'
        result[i] += 'CC:\n' +  cc[i]        

    return result
    
    
def make_pairwise_comparison_html_page(c1, c2, d, dir_name, agi2genename_dict):

    d1 = d[c1]
    d2 = d[c2]
    html_page = dir_name.rstrip('/') + '/' + c1 + '_vs_' + c2 + '.html'
    text_file1 = c1 + '_vs_' + c2 + '_1.txt'
    text_file2 = c1 + '_vs_' + c2 + '_2.txt'
    text_file3 = c1 + '_vs_' + c2 + '_3.txt'    

    fname1 = c1 + '_vs_' + c2 + '_1_GOenrichment.txt'
    fname2 = c1 + '_vs_' + c2 + '_2_GOenrichment.txt'
    fname3 = c1 + '_vs_' + c2 + '_3_GOenrichment.txt'

    ufname1 = os.path.splitext(fname1)[0] + '_uniq.txt'
    ufname2 = os.path.splitext(fname2)[0] + '_uniq.txt'
    ufname3 = os.path.splitext(fname3)[0] + '_uniq.txt'

    d1minusd2 = sorted(list(set(d1.keys()) - set(d2.keys())))
    d2minusd1 = sorted(list(set(d2.keys()) - set(d1.keys())))
    d1intd2 = sorted(list(set(d1.keys()).intersection(d2.keys())))
    sz_d1minusd2 = len(d1minusd2)
    sz_d2minusd1 = len(d2minusd1)
    sz_d1intd2 = len(d1intd2)
    f = open(html_page, 'w')
    f1 = open(text_file1,'w')
    f2 = open(text_file2,'w')
    f3 = open(text_file3,'w')
    f.write('<html>')
    f.write('<head>')
    f.write('<style> body {font-family:\"HelveticaNeue-Light\", \"Helvetica Neue Light\", \"Helvetica neue\"} </style>')
    f.write('</head>')
    f.write('<body>')
    f.write('<p>Condition A: %s</p>' % (c1))
    f.write('<p>Condition B: %s</p>' % (c2))    
    f.write('<p>First column: all genes in condition A but not in condition B; second column: all genes in both conditions; third column: all genes in condition B but not in condition A.</p>')
    f.write('<table border=1>')
    f.write('<tr>')
    f.write('<td>%s <br/><i>SUBTRACT</i><br/> %s <br/>(number=%d) <br/><a href=\"%s\">Download</a><br/><a href=\"%s\">GO enrichment</a><br/><a href=\"%s\">Unqiue GO enrichment</a></td> <td>%s <br/><i>AND</i><br/> %s <br/>(number=%d) <br/><a href=\"%s\">Download</a> <br/><a href=\"%s\">GO enrichment</a><br/><a href=\"%s\">Unique GO enrichment</a></td> <td>%s <br/><i>SUBTRACT</i></br> %s <br/>(number=%d) <br/><a href=\"%s\">Download</a><br/><a href=\"%s\">GO enrichment</a><br/><a href=\"%s\">Unique GO enrichment</a></td>' % (c1,c2,sz_d1minusd2,text_file1,fname1,ufname1,c1,c2,sz_d1intd2,text_file2,fname2,ufname2,c2,c1,sz_d2minusd1,text_file3,fname3,ufname3))
    f.write('</tr>')
    for i in range(max(sz_d1minusd2, sz_d1intd2, sz_d2minusd1)):
        f.write('<tr>')
        if i < sz_d1minusd2:
            x = d1minusd2[i]
            if (x in agi2genename_dict) and (agi2genename_dict[x] != x):
                f.write('<td>%s (%s)</td>' % (x, agi2genename_dict[x]))
            else:
                f.write('<td>%s</td>' % (x))
            f1.write('%s\n' % (x))
        else:
            f.write('<td>-</td>')
        if i < sz_d1intd2:
            x = d1intd2[i]
            if (x in agi2genename_dict) and (agi2genename_dict[x] != x):
                f.write('<td>%s (%s)</td>' % (x, agi2genename_dict[x]))
            else:
                f.write('<td>%s</td>' % (x))
            f2.write('%s\n' % (x))                
        else:
            f.write('<td>-</td>')
        if i < sz_d2minusd1:
            x = d2minusd1[i]
            if (x in agi2genename_dict) and (agi2genename_dict[x] != x):
                f.write('<td>%s (%s)</td>' % (x, agi2genename_dict[x]))
            else:
                f.write('<td>%s</td>' % (x))
            f3.write('%s\n' % (x))                
        else:
            f.write('<td>-</td>')            
        f.write('<tr>')
    f.write('</table>')    
    f.write('</body>')
    f.write('</html>')    
    f.close()
    f1.close()
    f2.close()
    f3.close()

    
    cmd = 'bash ' + GO_ENRICHMENT_SCRIPT + ' ' + text_file1 + ' > ' + fname1
    print(cmd)
    os.system(cmd)

    cmd = 'bash ' + GO_ENRICHMENT_SCRIPT + ' ' + text_file2 + ' > ' + fname2
    print(cmd)
    os.system(cmd)

    cmd = 'bash ' + GO_ENRICHMENT_SCRIPT + ' ' + text_file3 + ' > ' + fname3
    print(cmd)
    os.system(cmd)    

    diff_lst = go_diff3_bpmfcc(fname1, fname2, fname3)
    with open(ufname1, 'w') as text_file:
        text_file.write(diff_lst[0])
    with open(ufname2, 'w') as text_file:
        text_file.write(diff_lst[1])
    with open(ufname3, 'w') as text_file:
        text_file.write(diff_lst[2])
    
    cmd = ' '.join(['mv', text_file1, text_file2, text_file3, fname1, fname2, fname3, ufname1, ufname2, ufname3, dir_name])
    os.system(cmd)
    return [sz_d1minusd2, sz_d2minusd1, sz_d1intd2, html_page]
    
    
def make_comparison_table(f, gene_lists, agi2genename_dict, dir_name):
    conditions = sorted(gene_lists.keys())
    size = len(conditions)
    f.write('<table border=1>')
    f.write('<tr>')
    f.write('<td>-</td>')
    for i in range(1,size):
        f.write('<td>%s</td>' % (conditions[i]))
    f.write('</tr>')
    for i in range(size-1):
        c1 = conditions[i]
        f.write('<tr>')
        f.write('<td>%s</td>' % (c1))
        for j in range(1,size):
            if i >= j:
                f.write('<td>-</td>')
            else:
                c2 = conditions[j]
                t = make_pairwise_comparison_html_page(c1, c2, gene_lists, dir_name, agi2genename_dict)
                f.write('<td><a href=\"%s\">%d,%d,%d</a></td>' % (t[3], t[0], t[2], t[1]))
        f.write('</tr>')
    f.write('</table>')            
    
    
def make_gene_list_dict(s):
    lst = s.split('\n')
    d = {}
    for x in lst:
        t = x.split()
        if len(t) > 1 and t[0] != 'AGI_code':
            gene = t[0]
            if not gene in d:
                d[gene] = t[1]
    return d


def make_goenrichment_file(name, d):

    f = open('temp_gene_list.txt', 'w')
    gene_names = d.keys()
    f.write('\n'.join(gene_names))
    f.close()
    fname = name + '_GOenrichment.txt'
    cmd = 'bash ' + GO_ENRICHMENT_SCRIPT + ' temp_gene_list.txt > ' + fname
    os.system(cmd)
    return fname


def make_goenrichment_diff(lst, dir_name):

    cmd = 'python ' + GO_ENRICHMENT_DIFF_SCRIPT
    filename = dir_name.rstrip('/') + '/enrich_diff.txt'
    for x in lst:
        cmd += ' ' + x
    cmd += ' > ' + filename
    os.system(cmd)
    return filename


def make_text_gene_list(name, d0, d):
    '''
    name - intended file name
    d0 - dictionary where key is gene id, and value is fold-change
    d  - gene id to gene name conversion, e.g., AT3G14180 to IGN
    '''
    
    fname = name + '_gene_list.txt'
    f = open(fname, 'w')
    f.write('\t'.join(['AGI_code', 'max_fold_change_in_nearby_peaks', 'gene_name']))
    f.write('\n')
    for x in sorted(d0.keys()):
        v = d0[x]
        if (x in d) and (d[x] != x): # has a gene name, show gene name
            f.write('%s\t%s\t%s\n' % (x, v, d[x]))
        else:
            f.write('%s\t%s\n' % (x, v))
    f.close()
    return fname


def make_AGI_to_gene_name_dict(fname):
    d = {}
    f = open(fname)
    for line in f:
        line = line.strip()
        lst = line.split()
        agi = lst[0]
        name = lst[1]
        if not agi in d:
            d[agi] = name
        else:
            d[agi] += d[agi] + ';' + name
    f.close()
    return d

def get_parameter_value_float(s):
    index = s.find('=')
    result = 0.0
    try:
        result = float(s[index+1:])
    except:
        result = 1.0 * int(s[index+1:])
    return result


def get_parameter_value_string(s):
    return s[s.find('=') + 1:].strip()

def get_parameter_value_int(s):
    return int(get_parameter_value_string(s))
    
def get_global_parameters(fname):
    d = {'FC':2.0, 'QVALUE':0.01, 'PVALUE':0.05, 'PAIRWISE_COMPARE':'NO', 'TITLE':'None', 'TARGET_RANGE':DEFAULT_TARGET_RANGE}
    f = open(fname)
    for line in f:
        line = line.strip()
        if line.startswith('%%'):
            l = line[2:]
            lst = l.split()
            for x in lst:
                if x.upper().startswith('FC='):
                    d['FC'] = get_parameter_value_float(x)
                if x.upper().startswith('QVALUE='):
                    d['QVALUE'] = get_parameter_value_float(x)
                if x.upper().startswith('PVALUE='):
                    d['PVALUE'] = get_parameter_value_float(x)
                if x.upper().startswith('PATH'):
                    key = x[0:x.find('=')]
                    d[key] = get_parameter_value_string(x)
                if x.upper().startswith('PAIRWISE_COMPARE'):
                    d['PAIRWISE_COMPARE'] = get_parameter_value_string(x)
                if x.upper().startswith('TITLE'):
                    d['TITLE'] = get_parameter_value_string(x)
                if x.upper().startswith('TARGET_RANGE'):
                    d['TARGET_RANGE'] = get_parameter_value_string(x)
    f.close()
    return d


def get_local_parameter(s):

    d = {}
    t = s[1:]
    t = t.strip()
    lst = t.split()
    for x in lst:
        if x.upper().startswith('FC='):
            d['FC'] = get_parameter_value_float(x)
        if x.upper().startswith('QVALUE='):
            d['QVALUE'] = get_parameter_value_float(x)
        if x.upper().startswith('PVALUE='):
            d['PVALUE'] = get_parameter_value_float(x)

    return d
    
def get_path(s, glb_dict):
    '''
    expand PATH macro with real path
    '''
    index = s.find(':')
    t = s[index+1:]
    result = ''
    lst = t.split()
    for x in lst:
        a = x.find('<')
        b = x.find('>')
        if b > a and a >= 0:
            key = x[a+1:b]
            result += glb_dict[key]
            result += x[b+1:]
        else:
            result += x
    res = glob.glob(result)
    print res
    assert res,'[glob] empty glob:"%s"'%result
    assert len(res)==1,'[glob] ambiguous: %s'%('\n'.join(res))
    return res[0]


def get_conditions(fname, glb_dict):
    '''
    get all ChIPs/INPUT in different conditions
    '''

    d = {} # will be {cond1:{'CHIP':'...', 'INPUT':'...'}, cond2 :{}}
    f = open(fname)
    for line in f:
        line = line.strip()
        if line.startswith('@'):
            key = line[1:]
            key = key.replace(' ', '-')  # replace all spaces with a dash
            if not key in d:
                d[key] = {}
                
        if line.upper().startswith('CHIP:'):
            d[key]['CHIP'] = get_path(line, glb_dict)
        if line.upper().startswith('INPUT:'):
            d[key]['INPUT'] = get_path(line, glb_dict)
        if line.upper().startswith('%') and not line.upper().startswith('%%'):
            d[key]['PEAK_SELECT_PARAM'] = get_local_parameter(line)

    f.close()
    return d


def dumpclean(obj):
    '''
    show dictionary content, recursively
    '''
    if type(obj) == dict:
        for k, v in obj.items():
            if hasattr(v, '__iter__'):
                print k
                dumpclean(v)
            else:
                print '%s : %s' % (k, v)
    elif type(obj) == list:
        for v in obj:
            if hasattr(v, '__iter__'):
                dumpclean(v)
            else:
                print v
        else:
            print obj


def make_peak_call_script(key, dict, template_file):

    d = dict[key]
    outfilename = 'pipeline_' + key + '.sh'
    f = open(outfilename, 'w')
    f.write('# bam file paths\n')
    f.write('name=\"%s\"\n' % (key.strip()))
    for k in d:
        if k.startswith('CHIP'):
            f.write('chip=\"%s\"\n' % (d[k].strip()))
        if k.startswith('INPUT'):
            f.write('input=\"%s\"\n' % (d[k].strip()))
            
    f.write('\n')
    with open(template_file) as infile:
        f.write(infile.read())
    f.close()
    return outfilename

def make_param_string(d):
    result = ''
    for k in d:
        if k == 'FC':
            result += 'fc=' + str(d[k]) + ' '
        if k == 'QVALUE':
            result += 'q=' + str(d[k]) + ' '
        if k == 'PVALUE':
            result += 'p=' + str(d[k]) + ' '
    return result.strip()


def peak_to_target_genes(peak_file, bedmap_range,dbg=0):

    cmd = 'bedmap --echo --echo-map-id-uniq --delim \'\t\' ' \
          + '--range ' + bedmap_range \
          + ' ' + peak_file \
          + ' ' + ANNOTATION_FILE + ' > a.txt'
    if dbg: 
        print cmd
    os.system(cmd)

    cmd = 'cut -f 11 a.txt > a2.txt'
    os.system(cmd)

    cmd = 'python ' + GENELOCUS_TO_GENENAME_SCRIPT + ' a2.txt ' +  GENE_DESCRIPTION + ' > b.txt'
    os.system(cmd)

    cmd = 'paste a.txt b.txt > c.txt'
    os.system(cmd)

    # delete a line in c.txt if that lines contains transposable_element_gene or Not Found
    cmd = 'sed \'/transposable_element_gene\|Not Found/d\' c.txt > d.txt'
    os.system(cmd)
    
    cmd = 'python ' + EXTRACT_AGI_CODE_AND_FC + ' d.txt'
    print cmd
    res = subprocess.check_output(cmd,shell=1)
    return res
#     subprocess.callos.system(cmd)

def program_installed(pname):
    try:
        subprocess.call([pname, "--help"], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            print('ERROR: %s not found.' % (pname))
        else:
            print('ERROR: something is wrong with %s' % (pname))
        sys.exit()

def goatools_installed():
    f = open(GO_ENRICHMENT_SCRIPT)
    lines = f.readlines()
    f.close()
    line = lines[0].strip()
    lst = line.split('=')
    goatool_path = lst[1].strip()
    if not os.path.exists(os.path.join(goatool_path, 'scripts/find_enrichment.py')):
        print('ERROR: goatools not installed.  Install it and specify its path in the first line of chip-summary/depend/script/fe.sh')
        sys.exit()
    if not os.path.exists(os.path.join(goatool_path, 'scripts/go-basic.obo')):
        print('ERROR: go-basic.obo does not exists.  Download it (http://geneontology.org/ontology/go-basic.obo) and put it in folder goatools/scripts/.')
        sys.exit()
        
def main(f,dbg= 0,reCallPeak=0,gPar= None):
    global shellexec
    def shellexec(cmd,dbg= 0 ):
        if dbg:
            print cmd
            res = 'dbg'
        else:
            res = subprocess.check_output(cmd,shell=1)
        return res
    
#     #############################################################################################

#     DEPENDENT_FILES_PATH        = '/media/pw_synology3/Software/chip-summary/'  # [path of chip-summary.py]
#     DEFAULT_TARGET_RANGE        = '3000' # [change]  a string, not a number
#     SUMMARY_FILE_NAME           = 'summary.html'
#     SUMMARY_DIR                 = 'summary'
#     PEAK_CALL_PIPELINE_TEMPLATE = os.path.join(DEPENDENT_FILES_PATH, 'depend/script/pipeline724-t.sh')
#     PEAK_SELECT_SCRIPT          = os.path.join(DEPENDENT_FILES_PATH, 'depend/script/select_peaks.py')
#     GENELOCUS_TO_GENENAME_SCRIPT= os.path.join(DEPENDENT_FILES_PATH, 'depend/script/genelocus2genename.py')
    
#     #### Slowest part to be refactored???
#     EXTRACT_AGI_CODE_AND_FC     = os.path.join(DEPENDENT_FILES_PATH, 'depend/script/extract_AGI_code_and_fold_change.py')
#     GO_ENRICHMENT_SCRIPT        = os.path.join(DEPENDENT_FILES_PATH, 'depend/script/fe.sh')  # install goatools (GO enrichment) and edit fe.sh
#     GO_ENRICHMENT_DIFF_SCRIPT   = os.path.join(DEPENDENT_FILES_PATH, 'depend/script/goterm-matrix.py')
#     AGI_TO_GENE_NAMES           = os.path.join(DEPENDENT_FILES_PATH, 'depend/data/AGI-to-gene-names.txt')
#     ANNOTATION_FILE             = os.path.join(DEPENDENT_FILES_PATH, 'depend/data/genesTAIR10.bed') # for bedmap
#     GENE_DESCRIPTION            = os.path.join(DEPENDENT_FILES_PATH, 'depend/data/gene_description_20140101.txt')
#     MAX_FOLD_CHANGE             = 10  # for number of peaks versus fold-change plot

#     #############################################################################################
    
    gPar =  gPar or get_global_parameters(f)
    condDict = get_conditions(f, gPar)
    DIR = pyutil.dict2flat(gPar)
#     os.system('mkdir -p ' + DIR); os.chdir(DIR)
    if dbg == 1:
        d = gPar,condDict
        for dd in d:
            print pyutil.ppJson(dd)

        return d
#     try:
    if 1:

        # Collect results
        os.system('mkdir -p %s'%SUMMARY_DIR)


        # make pipeline files for peak calling
        def getPeak(k):
            sname = make_peak_call_script(k, condDict, PEAK_CALL_PIPELINE_TEMPLATE)
            print('Run %s ...' % (sname))        
            res = subprocess.call(['bash', sname])
#             return '%s_peaks.narrowPeak'%k
            return res
        
        if reCallPeak:
            # check that every ChIP file is present
            for k in condDict.keys():
                chip_file = condDict[k]['CHIP']
                input_file = condDict[k]['INPUT']
                if not os.path.exists(chip_file):
                    print('%s dose not exist. STOP' % (chip_file))
                    sys.exit()
                if not os.path.exists(input_file):
                    print('%s dose not exist. STOP' % (input_file))
                    sys.exit()
            [getPeak(k) for k in condDict.keys()]
            
        npkFS = ['%s_peaks.narrowPeak'%k for k in condDict.keys()]
        peakSummary(npkFS)
        
        gene_lists = {} # a dictionary of form d  = {'condition1': {'AT1G12345':'2.3', 'AT1G12346':'1.2'} }
        

def process(k=None,npkFile=None,gPar= None,dbg = 0):
    if k is None:
        assert npkFile,'must specify one arg'
        k = npkFile.rsplit('.',1)[0].split('/')[-1]
        

#             d = condDict[k]
    outd = {'files':{}}

    if 1:
#         parameter_string = 'fc=1.5 q=0.001 p=0.01'
        parameter_string = make_param_string(gPar)            

#             fileMacs = '%s/peaks/%s_peaks.narrowPeak'%(SUMMARY_DIR,k)
    small_narrowPeak = '{key}.snpk'.format(key=k)

    ### PeakFiltering
    cmd = 'python {SCRIPT} {INFILE} {PARAM} > {OUTF}'.format(
        SCRIPT=  PEAK_SELECT_SCRIPT,
        INFILE = npkFile,
        PARAM= parameter_string,
        OUTF = small_narrowPeak,
    )
    print(cmd);os.system(cmd)
    outd['param'] = parameter_string

#     cmd = 'cat {INFILE}  > {OUTF}'.format(
# #         SCRIPT=  PEAK_SELECT_SCRIPT,
#         INFILE = npkFile,
# #         PARAM= parameter_string,
#         OUTF = small_narrowPeak,
#     )
#     print(cmd);os.system(cmd)
#     outd['param'] = parameter_string
#     raise 0
    #### Fancy Histogram
    fc_thresholds = [x * 0.1 for x in range(10, 10*MAX_FOLD_CHANGE, 2)]
    npeak_lst = calc_npeak(fc_thresholds, k + '_peaks.narrowPeak', gPar['PVALUE'], gPar['QVALUE'])            
    plotName = SUMMARY_DIR.strip('/') + '/' + 'npeaks_vs_fc_' + k + '.txt' 
#     plot_npeak_vs_fc(fc_thresholds, npeak_lst, plotName)


    #### Produce geneLists     
    
    file_bedmap =  '%s.bedmap.tsv'%k
    
    
#     cmd = 'bedmap --echo --echo-map-id-uniq --delim \'\t\' ' \
#       + '--range ' + gPar['TARGET_RANGE'] \
#       + ' ' + small_narrowPeak \
#       + ' ' + ANNOTATION_FILE + ' | tee ' + file_bedmap +'.tmp'
    
#     buf = StringIO.StringIO(pyutil.shellexec(cmd))
# #     buf = file_bed
# #     res_bedmap = pd.read_table(buf,header=None)
    
    
# #     if buf.read():
#     if buf.len:
#         buf.seek(0)
#         df = sutil.parseBedmap(fname = buf)
#     else:
#         header = sutil.bedHeader + ['hit']
#         df = pd.DataFrame(columns = header)
    
    
    cmd = '''
bedtools slop -b {RANGE} -i {ANNO} -g $GSIZE |bedtools sort > {ANNOBASE}.{RANGE}
bedtools closest -d -a {SNPK} -b {ANNOBASE}.{RANGE} | tee {FOUT}.tmp
'''.format(
        ANNO = ANNOTATION_FILE,
        ANNOBASE=ANNOTATION_FILE.split('/')[-1],
        SNPK = small_narrowPeak,
        RANGE=gPar['TARGET_RANGE'],
        FOUT = file_bedmap    
    ).strip()
    
    buf = StringIO.StringIO(pyutil.shellexec(cmd))
    if buf.len:
        buf.seek(0)
        df = sutil.parseBedClosest(fname = buf)
    else:
        assert 0,' Buffer is empty, check error msg'
    df['condition'] = k
    df = df[df['distance']==0]
    
#         raise e
    df = df.sort_values('FC',ascending = False,inplace = 0)
    
    #### deduplication on gene acc
    df = df.loc[~df.duplicated('hit')]
    res_bedmap = df
    df.to_csv(file_bedmap,sep='\t')
    
    
    genes = df
    
    outd['genes'] = None
    outd['nGene'] = len(df['hit'].unique())
    outd['file_bedmap'] = file_bedmap
#     outd['res_bedmap'] = res_bedmap

    fname = '%s/%s.gene.txt'% (SUMMARY_DIR,k)
    dfc = df.copy()[['hit','FC', 'acc',]]
    dfc.columns = ['geneAcc','maxFoldChange','peakAcc']
    dfc.to_csv(fname,sep= '\t')   
    
    
    outd['glst_filename'] =  fname
#     outd['goenrich_filename'] =  make_goenrichment_file(SUMMARY_DIR + '/' + k, genes)
    outd['goenrich_filename'] =  'NotImplemented'
    outd['plot_file'] = plotName 
    outd['peak_file'] = small_narrowPeak

    outd['key'] = k
    outd['extra'] = ''
    return outd             

def peakSummary(npkFS,gPar = None,dbg=0,
               FC    = 1.5,
               PVALUE=0.01,
               QVALUE=0.0005,
                **kwargs
               ):
    gPar = gPar or {
    "FC": FC, 
    "PVALUE": PVALUE, 
    "QVALUE": QVALUE, 
    "PAIRWISE_COMPARE": "Y", 
    "TARGET_RANGE": "3000", 
    "TITLE": "testRun"
}
    
#     f = functools.partial(process,gPar = gPar)
    f = lambda x: process(npkFile = x,gPar = gPar,dbg = dbg)
    condRes = res = map(f, npkFS)
    if dbg:
        with open('condRes.dbg','w') as f:
            print >>f,pyutil.ppJson(condRes)
    if dbg == 2:
        return condRes
    
    
    dfs = [pd.read_table(x['file_bedmap']).set_index('hit') for x in res]
#     for df in dfs:
#         print df.head(10)
#     print [type(df) for df in dfs]
    indAll = pd.concat( dfs,axis =1 ,join='outer').index
    print '[db1]',dfs[0].head()

    
    df = pd.concat( [ df.reindex( indAll) for df in dfs],
                   axis = 0) 
    df = df.reset_index()
#     df..reset_index()
#     df = df.set_index('hit')
    
    print '[db2]',df.head()

    cols = df.columns.to_series()
    cols[0]='index'
    df.columns = cols
    print '[db3]',df.head()
    df_fc = df.pivot_table(columns = 'condition', values='FC',index='index').fillna(0)

    sanitise = lambda x: x.split('.',1)[0]
    df_fc.index = map(sanitise, df_fc.index)

#     index = scount.vstack([dfs],as_index=1,how = 'outer')
    getPM =    lambda lst: ''.join(['+' if x!=0 else '-' for x in lst])
    vals = df_fc.apply(getPM,axis=1)
    df_fc.insert(0,'pmSummary',vals)
    print df_fc.head(10)
    df_fc.to_csv('FoldChange_table.csv')

    ##### write html summary report 

    print('... make html page %s' % (SUMMARY_FILE_NAME))

    f = open(SUMMARY_FILE_NAME, 'w')
    TITLE = 'test'
#     TITLE = gPar['TITLE']
    f.write('<html>')
    f.write('<head>')
    f.write('<title>%s</title>'% (TITLE) )
    f.write('<style> body {font-family:\"HelveticaNeue-Light\", \"Helvetica Neue Light\", \"Helvetica neue\"} </style>')
    f.write('</head>')
    f.write('<body>')
    f.write('<h2>%s</h2>' % (TITLE))

    #####################################################################
    f.write('<h3>Target genes and associated GO terms</h3>')

    table_string = '<table><tr><td>Sample</td><td>Target gene list</td><td>#target genes</td><td>GO enrichment</td><td>Peak selection thresholds</td><td>#peaks plot</td></tr>'
    rowFmt = '''
        <tr><td>{key}</td>
        <td><a href="{glst_filename}">target genes</a></td>
        <td align=right>{nGene:d}</td>
        <td><a href="{goenrich_filename}">enrichment</a></td>
        <td><a href="{peak_file}">{param}</a></td>
        <td><a href="{plot_file}">plot</a></td></tr>
        '''
    for d in res:
        table_string += rowFmt.format(**d)



    table_string += '</table>'
    f.write(table_string)

    #####################################################################
    f.write('<h3>Enriched GO terms associated to target genes in different conditions</h3>')
    f.write('<p>Most shared GO terms across conditions are on the top in the following table.</p>')
    filename = make_goenrichment_diff([d['goenrich_filename'] for d in res],
                                      SUMMARY_DIR)
    f.write('<a href=\"%s\">Each row is a GO term. Each column is a condition.</a>' % (filename))


    #####################################################################
#     if gPar['PAIRWISE_COMPARE'].lower().startswith('y'): 
#         f.write('<h3>Pairwise comparison between conditions</h3>')
#         f.write('<p>Each cell in the following table contains three numbers, X, Y and Z. X is the number of target genes that are in condition A but not in condition B. Z is the number of target genes that are in condition B but not in condition A.  Y is the number of target genes that are in both conditions.</p>')
#         make_comparison_table(f, gene_lists, agi2genename_dict, SUMMARY_DIR)



    #####################################################################
    
    file_bindingMat = 'FoldChange_table.csv' #% SUMMARY_DIR
   
    f.write('<h3>Binding to target genes in different conditions</h3>')
    f.write('<p><b>Note:</b> In the following, \'+\' means binding near a target gene in a particular experimental condition, and \'-\' means non-binding.</p>')
    f.write('<p>The columns are:<br/><br/>')
    colName = ['AGI_locus_name'] + [d['key'] for d in res] + ['gene_name (if available)']
    f.write( '<br/>'.join(colName))
    f.write('</p>')
    f.write('<p><a href=\"%s\">Download text version</a></p>' % (file_bindingMat))

    f.write('</body>')    
    f.write('</html>')    
    f.close()


    # clean up
    print('Done.')
    
#         delete_file('a.txt')
#         delete_file('a2.txt')
#         delete_file('b.txt')
#         delete_file('c.txt')
#         delete_file('d.txt')
#         delete_file('temp_gene_list.txt')
#         delete_file('temp_small.narrowPeak')
#         delete_file('predicctd')
#         delete_file('macs2-predictd.txt')
#         delete_file('predictd')
#         delete_files('*.xls')
#         delete_files('*.bed')
#         delete_files('*_peaks.narrowPeak')
# #         delete_files('AGI2-*.txt')
# #         delete_files('AGI2*.txt')
#         delete_files('small_*.narrowPeak.txt')
#         move_files('pipeline_*.sh', SUMMARY_DIR)

import os, argparse
import pymisca.util as pyutil
import pandas as pd
import synotil.util as sutil
import synotil.countMatrix as scount
import StringIO
if __name__=='__main__':
    #############################################################################################

    DEPENDENT_FILES_PATH        = '/media/pw_synology3/Software/chip-summary/'  # [path of chip-summary.py]
    DEFAULT_TARGET_RANGE        = '3000' # [change]  a string, not a number
    SUMMARY_FILE_NAME           = 'summary.html'
    SUMMARY_DIR                 = 'summary'
    PEAK_CALL_PIPELINE_TEMPLATE = os.path.join(DEPENDENT_FILES_PATH, 'depend/script/pipeline724-t.sh')
    PEAK_SELECT_SCRIPT          = ('/media/pw_synology3/Software/chip-summary/depend/script/select_peaks.py')
    GENELOCUS_TO_GENENAME_SCRIPT=('/media/pw_synology3/Software/chip-summary/depend/script/genelocus2genename.py')
    
    #### Slowest part to be refactored???
    EXTRACT_AGI_CODE_AND_FC     = ('/media/pw_synology3/Software/chip-summary/depend/script/extract_AGI_code_and_fold_change.py')
    GO_ENRICHMENT_SCRIPT        = os.path.join(DEPENDENT_FILES_PATH, 'depend/script/fe.sh')  # install goatools (GO enrichment) and edit fe.sh
    GO_ENRICHMENT_DIFF_SCRIPT   = os.path.join(DEPENDENT_FILES_PATH, 'depend/script/goterm-matrix.py')
    AGI_TO_GENE_NAMES           = ('/media/pw_synology3/Software/chip-summary/depend/data/AGI-to-gene-names.txt')
    ANNOTATION_FILE             = os.path.join(DEPENDENT_FILES_PATH, 'depend/data/genesTAIR10.bed') # for bedmap
    GENE_DESCRIPTION            = ('/media/pw_synology3/Software/chip-summary/depend/data/gene_description_20140101.txt')
    MAX_FOLD_CHANGE             = 10  # for number of peaks versus fold-change plot

    #############################################################################################
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', action='store_true')
    parser.add_argument('infiles', nargs='+', default=[],
                       help='*.narrowPeak files to be summarised')
    parser.add_argument('--ref','-r', nargs='?',
                        default='/media/pw_synology3/Software/chip-summary/depend/data/genesTAIR10.bed',
                       help='reference annotation')
    parser.add_argument('--output_dir','-o', nargs='?',
                        default='chipSum',
                       help='reference annotation')
    parser.add_argument('--FC',nargs='?',type=float,
                        default= 1.5 ,
                        help='')
    parser.add_argument('--PVALUE',nargs='?',type=float,
                        default= 0.01 ,
                        help='')
    parser.add_argument('--QVALUE',nargs='?',type=float,
                        default= 0.01 ,
                        help='')
    args =  parser.parse_args()
    ANNOTATION_FILE   = os.path.abspath(args.ref)
    infiles = map(os.path.abspath, args.infiles)

    VAR = 'GSIZE'
    assert VAR in os.environ, 'Cannot find bash environ: $%s'%VAR
    
    dbg = int(args.v)*10

    os.system('mkdir -p %s'%args.output_dir)
    os.chdir(args.output_dir)
    cmd= '''
mkdir -p summary; 
mkdir -p summary/npeaks_vs_fc_npk;
mkdir -p summary/npk;
cp -r {infiles} -t .
'''.format(
    infiles=' '.join(infiles))
    os.system(cmd)

    agi2genename_dict = make_AGI_to_gene_name_dict(AGI_TO_GENE_NAMES)
    
    peakSummary(infiles,
                dbg=dbg,
                **args.__dict__)
    
    SUMMARY_DIR = SUMMARY_DIR.rstrip('/')
    os.system('mkdir -p %s/peaks' % SUMMARY_DIR)
    # os.system('cp *_peaks.narrowPeak *.bed *.xls %s/peaks' % SUMMARY_DIR)

    #         res = {k:process(k) for k in condDict}

    sys.exit(0)
#     args =  parser.parse_args
#     parser.
#     parser
#     ## main #####################################################
#     if len(sys.argv) < 2:
#         print('Usage: python chip-summary.py summary.txt')
#         sys.exit()

#     # make sure path to chip-summary is set correctly
#     if not os.path.exists(PEAK_CALL_PIPELINE_TEMPLATE):
#         print('Set DEPENDENT_FILES_PATH in chip-summary.py the path of chip-summary.')
#         sys.exit()
#     else:
#         print('\nTips: set DEPENDENT_FILES_PATH in chip-summary.py to the full path of chip-summary (e.g., %s).\n      Then chip-summary.py is executable in other directories.\n') % (os.getcwd())


#     # make sure macs2 and bedmap are installed
#     program_installed('macs2')
#     program_installed('bedmap')

#     # make sure goatools is installed
#     goatools_installed()


#     f = sys.argv[1]
#      main(f)