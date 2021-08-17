import subprocess
import os, shutil, glob
import re
import sys
import argparse
from Bio import SeqIO
def count_gap(text):
    if text.strip() == "": # To take of care of all space input
        return 0
    #count = sum([1 if char=='N' else 0 for char in text ])
    count = text.count('N') # Your error is here, Only check for 1 space instead of 3 spaces
    count =count+ text.count('n')
    total_chars = len(text)

    return round(count / total_chars,3)*100
def ReverseComplement(Pattern):
    str=""
    for c in Pattern:
        if c=="A":
            str=str+"T"
        elif c=="T":
            str=str+"A"
        elif c=="G":
            str=str+"C"
        elif c=="C":
            str=str+"G"
        else:
            str=str+c
    return str[::-1]
def Complement(Pattern):
    str=""
    for c in Pattern:
        if c=="A":
            str=str+"T"
        elif c=="T":
            str=str+"A"
        elif c=="G":
            str=str+"C"
        elif c=="C":
            str=str+"G"
        else:
            str=str+c
    return str
def setupdb(args):
    """
    make blast database from fasta file in db folder,
    :param : fasta file (with folder'folder is the name of db and filename is 'sequences')
    :return:
    """
    #seqfile='/mnt/data/coronacheck/sarscov2.fasta'
    #gisaid_dir='/home/quang/Downloads/gisaid4000'
    input_dir=args.samples
    output_dir=args.output
    threadshol_n=args.threadshold
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    print('Setup database...')
    list_seq=[]

    for root, dirs, files in os.walk(input_dir):
        for _file in files:
            if _file.endswith(('.fasta')):
                #print(str(root)+'/'+_file)
                for seq in SeqIO.parse(str(root)+'/'+_file,'fasta'):
                    perc_gap=count_gap(seq.seq)
                    if perc_gap<=threadshol_n and len(seq.seq)>20000:
                       list_seq.append(seq)
    outfile=os.path.join(output_dir,'sequences_'+str(threadshol_n)+'N.fasta')
    SeqIO.write(list_seq,outfile,'fasta')
    print(len(list_seq))
    cmd="makeblastdb -in {path} -title {name} -dbtype {type} -logfile /dev/null".format(
                            path=outfile,
                            name='corona',
                            type='nucl'

    )
    print (cmd)
    os.system(cmd)
def setupdbfile(dbfile):
    """
    make blast database from fasta file in db folder,
    :param : fasta file (with folder'folder is the name of db and filename is 'sequences')
    :return:
    """

    cmd="makeblastdb -in {path} -title {name} -dbtype {type} -logfile /dev/null".format(
                            path=dbfile,
                            name='corona',
                            type='nucl'

    )
    print (cmd)
    os.system(cmd)
def export_file(sample,db,result,output,dict_cds):
    #    blast_fields={'qseqid':t[0], 'qstart':t[1], 'qend':t[2], 'qlen':t[3],\
    #     'sseqid':t[4], 'sstart':t[5], 'send':t[6], 'slen':t[7], 'sstrand':t[8],\
    #      'length':t[10], 'mismatch':t[11], 'qseq':t[12], 'sseq':t[13],\
    #       'stitle':t[14]}
    f=open(output,'w')
    f.write('PRIMER\tPRIMER START\tPRIMER END\tPRIMER LEN\tSAMPLE ID\tSAMPLE START\tSAMPLE END\tSAMPLE LENGTH\tSAMPLE STRAND\
    \tALIGNMENT LENGTH\t%COVER\tIDENTITY\tMISMATCH\tMISMATCH POSITION\tPRIMER ALIGNED SEQUENCE\tSAMPLE SEQUENCE\tSAMPLE TITLE\tGENE\n')
    #print(result)
    dict_sample_primer={}
    for s in result:
        #filter missmatch and not cover:
        isPrint=False
        if not s['sseqid'] in dict_sample_primer:
            dict_sample_primer[s['sseqid']]={}
        mismatch_c=''
        if not s['qseqid'] in dict_sample_primer[s['sseqid']]:
            dict_sample_primer[s['sseqid']][s['qseqid']]={}
            dict_sample_primer[s['sseqid']][s['qseqid']]['c']=''
            dict_sample_primer[s['sseqid']][s['qseqid']]['mm']=99999
            dict_sample_primer[s['sseqid']][s['qseqid']]['range']=''
        if int(s['mismatch'])>0:
            isPrint=True
            for i in range(len(s['sseq'])):
                if not s['qseq'][i]==s['sseq'][i]:
                    mismatch_c=mismatch_c+str(i+1)+':'+s['qseq'][i]+'-'+s['sseq'][i]+','

            dict_sample_primer[s['sseqid']][s['qseqid']]['mm']= int(s['mismatch'])
            dict_sample_primer[s['sseqid']][s['qseqid']]['c']='('+mismatch_c+')'

        if float(s['pident'])<100:
            isPrint=True
        if (int(s['qlen'])-int(s['length']))>3 :
             isPrint=False
        elif(int(s['qlen'])-int(s['length']))>0:
            isPrint=True
        else:
            if int(s['mismatch'])==0:
                isPrint=False
        #isPrint=True




        if (int(s['qlen'])-int(s['length'])) < dict_sample_primer[s['sseqid']][s['qseqid']]['mm']:
            if dict_sample_primer[s['sseqid']][s['qseqid']]['c']=='':
                dict_sample_primer[s['sseqid']][s['qseqid']]['mm']=int(s['qlen'])-int(s['length'])
            if (int(s['qlen'])-int(s['length']))>0:
                dict_sample_primer[s['sseqid']][s['qseqid']]['range']='('+str(s['qstart'])+'->'+str(s['qend'])+'/'+s['qlen']+','+str(s['sstart'])+'->'+str(s['send'])+')'


        if isPrint:
            #find genes:
            gene=''
            if not dict_cds==None:

                for key in dict_cds:
                    if s['sseqid'] in key:
                        for cds in dict_cds[key]:
                            if int(s['sstart']) >= dict_cds[key][cds]['cds_f'] and int(s['sstart']) <= dict_cds[key][cds]['cds_e']:
                                gene= dict_cds[key][cds]['gene']
            f.write(s['qseqid']+'\t'+str(s['qstart'])+'\t'+str(s['qend'])+'\t'+str(s['qlen'])+'\t'+\
            s['sseqid']+'\t'+str(s['sstart'])+'\t'+str(s['send'])+'\t'+str(s['slen'])+'\t'+\
            str(s['sstrand'])+'\t'+str(s['length'])+'\t'+str((int(int(s['length'])/int(s['qlen'])*100)))+'\t'+str(s['pident'])+'\t'+str(s['mismatch'])+'\t'+mismatch_c+'\t'+s['qseq']+\
            '\t'+s['sseq']+'\t'+s['stitle']+'\t'+gene+'\n')
    f.close()
    return dict_sample_primer

def blast(sample,db, identity=90, threads=1, mincov=90,dbtype='nucl'):
    """
    Call blastn with params
    :param query_file (in fasta), db (blast indexed db), number of threads and identity
    :return: list BLASTFields objects
    """
    #check db is indexed
    #dbfile=os.path.join(db_folder, 'sequences')




    # run blastn
    cmd ='blastn -query {query} -task blastn-short -max_target_seqs 1000000 -perc_identity {identity} -db {db} -outfmt \'6 qseqid qstart qend qlen sseqid sstart send slen sstrand length pident mismatch qseq sseq stitle\'  -num_threads {threads} >temp.tab'.format(

        query=sample,
        identity=identity,
        db=db,
        threads=threads

    )


    print(cmd)
    os.system(cmd)
    #parse result
    f=open('temp.tab')
    line = f.readline()
    result=[]
    while line:
        #result.append(line)
        t=line.strip().split('\t')
        blast_fields={'qseqid':t[0], 'qstart':t[1], 'qend':t[2], 'qlen':t[3],\
         'sseqid':t[4], 'sstart':t[5], 'send':t[6], 'slen':t[7], 'sstrand':t[8],\
           'length':t[9],'pident':t[10], 'mismatch':t[11], 'qseq':t[12], 'sseq':t[13],\
           'stitle':t[14]}
        result.append(blast_fields)
        line = f.readline()
    f.close()
    if os.path.exists('temp.tab'):
        os.remove('temp.tab')



    return result


def check_pos(list_primer_pos,ref_d):
    #check order
    for i in range(1,len(list_primer_pos)):
        if list_primer_pos[i]-list_primer_pos[i-1]<0:
            return False,'wrong order'
    #check length
    d=list_primer_pos[len(list_primer_pos)-1]-list_primer_pos[0]
    if d<ref_d-50 or d>ref_d+50:
        return False,'wrong distance'
    return True,'ok'
def readCDSFile(cds_file):
    dic_cds={}
    for seq in SeqIO.parse(cds_file,'fasta'):
        if seq.description.startswith('join'):
            cdses=seq.description.split('|')
            for i  in range(0,len(cdses)-1):
                cds=cdses[i].replace('join','').replace('(','').replace(')','')
                if not ':' in cds:
                    continue
                id=cds.split(':')[0].strip()
                c=cds.split(':')[1].strip()
                gene=cdses[len(cdses)-1]
                if not id in dic_cds:
                    dic_cds[id]={}
                dic_cds[id][c]={}
                dic_cds[id][c]['cds_f']=int(c.split('..')[0])
                dic_cds[id][c]['cds_e']=int(c.split('..')[1])
                dic_cds[id][c]['gene']=gene
        else:
            id=seq.name.split(':')[0].strip()
            cds=seq.name.split(':')[1].strip()
            gene=seq.description.split('|')[1]
            if not id in dic_cds:
                dic_cds[id]={}
            dic_cds[id][cds]={}

            dic_cds[id][cds]['cds_f']=int(cds.split('..')[0])
            dic_cds[id][cds]['cds_e']=int(cds.split('..')[1])
            dic_cds[id][cds]['gene']=gene


    return dic_cds
import subprocess
import os, shutil, glob
import re
import sys
import itertools
import csv
import argparse
import json
from Bio import SeqIO
def count_gap(text):
    if text.strip() == "": # To take of care of all space input
        return 0
    #count = sum([1 if char=='N' else 0 for char in text ])
    count = text.count('N') # Your error is here, Only check for 1 space instead of 3 spaces
    count =count+ text.count('n')
    total_chars = len(text)

    return round(count / total_chars,3)*100



def combileFastaFile(folder_in,file_out):
    list_seq=[]

    for root, dirs, files in os.walk(folder_in):
        for _file in files:
            if _file.endswith(('.fasta')):
                #print(str(root)+'/'+_file)
                for seq in SeqIO.parse(str(root)+'/'+_file,'fasta'):
                    perc_gap=count_gap(seq.seq)
                    if perc_gap<=0 and len(seq.seq)>20000:
                        list_seq.append(seq)
        SeqIO.write(list_seq,file_out,'fasta')
def filterFastaFile(fasta_in,file_out):
    list_seq=[]
    for seq in SeqIO.parse(fasta_in,'fasta'):
        perc_gap=count_gap(seq.seq)
        if perc_gap<=1 and len(seq.seq)>20000:

            list_seq.append(seq)
    SeqIO.write(list_seq,file_out,'fasta')
    return file_out
def getAMP(region,primerFi,primerRo):
    s_pos,mm,m=FittingAlignment(region,primerFi,1,1)
    e_pos,mm,m=FittingAlignment(region,primerRo,1,1)
    e_pos=e_pos+len(primerRo)
    return region[s_pos:e_pos],s_pos,e_pos
def setupdb(args):
    """
    make blast database from fasta file in db folder,
    :param : fasta file (with folder'folder is the name of db and filename is 'sequences')
    :return:
    """
    input_dir=args.samples
    output_dir=args.output
    threadshol_n=args.threadshold
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    print('Setup database...')
    list_seq=[]

    for root, dirs, files in os.walk(input_dir):
        for _file in files:
            if _file.endswith(('.fasta')):
                #print(str(root)+'/'+_file)
                for seq in SeqIO.parse(str(root)+'/'+_file,'fasta'):
                    list_seq.append(seq)
    outfile=os.path.join(output_dir,'sequences_'+str(threadshol_n)+'N.fasta')
    SeqIO.write(list_seq,outfile,'fasta')
    print(len(list_seq))
    cmd="makeblastdb -in {path} -title {name} -dbtype {type} -logfile /dev/null".format(
                            path=outfile,
                            name='corona',
                            type='nucl'

    )
    print (cmd)
    os.system(cmd)
def setupdbRef(ref_file):
    cmd="makeblastdb -in {path} -title {name} -dbtype {type} -logfile /dev/null".format(
                            path=ref_file,
                            name='corona',
                            type='nucl'

    )
    print (cmd)
    os.system(cmd)
import uuid
def blast(sample,db, identity=80, threads=10, mincov=70,dbtype='nucl'):
    """
    Call blastn with params
    :param query_file (in fasta), db (blast indexed db), number of threads and identity
    :return: list BLASTFields objects
    """
    #check db is indexed
    #dbfile=os.path.join(db_folder, 'sequences')
    # run blastn
    uid=str(uuid.uuid1().hex)
    temp_fasta=uid+'_temp.fasta'
    temp_tab=uid+'_temp.tab'
    f=open(temp_fasta,'w')
    f.write('>primer\n')
    f.write(sample)
    f.close()
    cmd ='blastn -query {query} -max_target_seqs 1000000 -perc_identity {identity} -db {db} -outfmt \'6 qseqid qstart qend qlen sseqid sstart send slen sstrand length pident mismatch qseq sseq stitle bitscore\'  -num_threads {threads} >{output}'.format(

        query=temp_fasta,
        identity=identity,
        db=db,
        threads=threads,
        output=temp_tab

    )
    print(cmd)
    os.system(cmd)
    #parse result
    f=open(temp_tab)
    line = f.readline()
    result=[]
    while line:
        #result.append(line)
        t=line.strip().split('\t')
        blast_fields={'qseqid':t[0], 'qstart':t[1], 'qend':t[2], 'qlen':t[3],\
         'sseqid':t[4], 'sstart':t[5], 'send':t[6], 'slen':t[7], 'sstrand':t[8],\
           'length':t[9],'pident':t[10], 'mismatch':t[11], 'qseq':t[12], 'sseq':t[13],\
           'stitle':t[14],'bitscore':t[15]}
        result.append(blast_fields)
        line = f.readline()
    f.close()
    if os.path.exists(temp_tab):
        os.remove(temp_tab)
    if os.path.exists(temp_fasta):
        os.remove(temp_fasta)
    return result
def scoreTable(a,b):
    if a==b:
        return 1
    else:
        return -1
def FittingAlignment(v,w,score,sigma):
    #print(v)
    #print(w)
    if v=='' or w=='':
        return 0,'',0
    s=[[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
    bk=[[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
    s[0][0]=0
    for i in range(1,len(v)+1):
        s[i][0]=0
    for j in range(1,len(w)+1):
        s[0][j]=s[0][j-1]-1
    for i in range(1,len(v)+1):
        for j in range(1,len(w)+1):
            if v[i-1]==w[j-1]:
                s[i][j]=max(s[i-1][j]-sigma,s[i][j-1]-sigma,s[i-1][j-1]+score)
            else:
                s[i][j]=max(s[i-1][j]-sigma,s[i][j-1]-sigma,s[i-1][j-1]-sigma)
            if s[i][j]==s[i-1][j]-sigma:
                bk[i][j]=1
            elif s[i][j]==s[i][j-1]-sigma:
                bk[i][j]=2
            else:
                bk[i][j]=3
    AlignV=''
    AlignW=''
    mi=0
    mj=len(w)
    for i in range(0,len(v)+1):
        if s[i][len(w)] >s[mi][len(w)]:
            mi=i
    i=mi
    j=mj
    while True:
        if j>0:
            if bk[i][j] ==3:

                AlignV=v[i-1]+AlignV
                AlignW=w[j-1]+AlignW
                i=i-1
                j=j-1
            elif bk[i][j] ==1:

                AlignW='-'+AlignW
                AlignV=v[i-1]+AlignV
                i =i-1
            else:

                AlignW=w[j-1]+AlignW
                AlignV='-'+AlignV
                j =j-1

        else:
            break

    mi=i
    mm=''
    m=0
    for i in range(len(AlignW)):
        if not AlignV[i]==AlignW[i]:
            mm=mm+str(i+1)+'('+AlignV[i]+'->'+AlignW[i]+')'+';'
        else:
            m=m+1
    #print(AlignW)
    #print(AlignV)
    count_gap=0
    for i in range(len(AlignV)):
        if not AlignV[i]=='-':
            break
        else:
            count_gap=count_gap+1
    mi=mi-count_gap
    return mi,mm,m/(len(w))
def GlobalAlignment(v,w,score,sigma):
    if v=='' or w=='':
        return '',0
    s=[[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
    bk=[[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
    s[0][0]=0
    for i in range(1,len(v)+1):
        s[i][0]=s[i-1][0]-sigma
    for j in range(1,len(w)+1):
        s[0][j]=s[0][j-1]-sigma
    for i in range(1,len(v)+1):
        for j in range(1,len(w)+1):
            # s[i][j]=max(s[i-1][j]-sigma,s[i][j-1]-sigma,s[i-1][j-1]+score[v[i-1]][w[j-1]])
            if not v[i-1]==w[j-1]:
                s[i][j]=max(s[i-1][j]-sigma,s[i][j-1]-sigma,s[i-1][j-1]-1)
            else:
                s[i][j]=s[i-1][j-1]+score
            if s[i][j]==s[i-1][j]-sigma:
                bk[i][j]=1
            elif s[i][j]==s[i][j-1]-sigma:
                bk[i][j]=2
            else:
                bk[i][j]=3
    AlignV=''
    AlignW=''
    i=len(v)
    j=len(w)
    while True:
        if i>0 and j>0:
            if bk[i][j] ==3:

                AlignV=v[i-1]+AlignV
                AlignW=w[j-1]+AlignW
                i=i-1
                j=j-1
            elif bk[i][j] ==1:

                AlignW='-'+AlignW
                AlignV=v[i-1]+AlignV
                i =i-1
            else:

                AlignW=w[j-1]+AlignW
                AlignV='-'+AlignV
                j =j-1
        elif i>0:
            AlignW='-'+AlignW
            AlignV=v[i-1]+AlignV
            i =i-1
        elif j>0:
            AlignW=w[j-1]+AlignW
            AlignV='-'+AlignV
            j =j-1
        else:
            break
    mm=''
    m=0
    for i in range(len(AlignV)):
        if not AlignV[i]==AlignW[i]:
            mm=mm+str(i+1)+'('+AlignV[i]+'->'+AlignW[i]+')'+';'
        else:
            m=m+1
    return mm,m/len(w)
def Align(primer,db):

    list_hit=[]




    for seq in SeqIO.parse(db,'fasta'):
        hit={}
        hit['stitle']=seq.description
        pos,mm,m=FittingAlignment(str(seq.seq).upper(),primer,1,1)
        hit['qstart']=1
        hit['sstart']=pos-len(primer)+1
        hit['qlen']=len(primer)
        list_hit.append(hit)
    return list_hit

def extensePrimerToHaplotype(ref_genome,primers):
    ref_sequence=''
    primers_extend=[]
    for seq in SeqIO.parse(ref_genome,'fasta'):
        ref_sequence=str(seq.seq).upper()
    for gp in primers:
        #blast Fo:
        ngp=gp
        blast_Fo=blast(gp['primer'][0]['seq'],ref_genome)
        #get best result
        bitscore=0
        for b in blast_Fo:
            if float(b['bitscore'])>bitscore:
                bitscore=float(b['bitscore'])
                #get haplotype by position, extend 20nu
                ht=ref_sequence[int(b['sstart'])-20:int(b['sstart'])+gp['length']+20]
                ngp['haplotype']=ht

        primers_extend.append(ngp)
    return primers_extend

def collectHaplotypeFromCorpus(genes,db):
    ref_sequence={}
    num_seq=0
    for seq in SeqIO.parse(db,'fasta'):
        num_seq=num_seq+1
        ref_sequence[seq.description]=str(seq.seq).upper()

    dict_haplotype={}
    dict_domain={}
    check_list={}
    haplotype_set=set()
    name_ht_dict={}
    for gene in genes:
        print('check genes' +gene.id)
        dict_haplotype[gene.id]={}
        dict_domain[gene.id]={}
        #print(len(str(gene.seq)))
        blast_haplotype=blast(str(gene.seq),db)
        check_list[gene.id]={}
        for b in blast_haplotype:
            ht=''


            position=int(b['sstart'])-1
            if int(b['qstart'])>1:
                position=int(b['sstart'])-int(b['qstart'])
            ht=b['sseq'].replace('-','')
            print(ht)
            if not b['stitle'] in dict_domain:
                dict_domain[b['stitle']]={}
            if not gene.id in dict_domain[b['stitle']]:
                dict_domain[b['stitle']][gene.id]=[]
            allele={'seq':ht,'start':position,'len':int(b['slen'])}
            print(b['pident'])
            print(b['sstrand'])
            print(b['mismatch'])


            if not ht in haplotype_set:

                newname=gene.id+'-'+str(len(dict_haplotype[gene.id].keys())+1)
                dict_haplotype[gene.id][newname]={}
                dict_haplotype[gene.id][newname]['seq']=ht
                p,mm,score=FittingAlignment(str(gene.seq),ht,1,1)
                print(mm)
                dict_haplotype[gene.id][newname]['mm']=mm
                dict_haplotype[gene.id][newname]['sample']=[]
                haplotype_set.add(ht)
                name_ht_dict[ht]=newname
                dict_haplotype[gene.id][newname]['sample'].append((b['stitle'],position))
            else:
                dict_haplotype[gene.id][name_ht_dict[ht]]['sample'].append((b['stitle'],position))
            #print(name_ht_dict)
            allele['ht']=name_ht_dict[ht]
            dict_domain[b['stitle']][gene.id].append(allele)
            print(b['stitle']+':'+allele['ht'])
    return dict_haplotype,dict_domain,num_seq

def Score(p, q):
    # your code here
    n=len(q)
    count=0
    list_mm=[]
    count_miss=0
    for i in range(0,n):
        if p[i]==q[i]:
            count=count+1
        elif p[i]=='Y' and (q[i]=='C' or q[i]=='T'):
            count=count+1
        elif q[i]=='Y' and (p[i]=='C' or p[i]=='T'):
            count=count+1
        elif p[i]=='R' and (q[i]=='A' or q[i]=='G'):
            count=count+1
        elif q[i]=='R' and (p[i]=='A' or p[i]=='G'):
            count=count+1
        elif p[i]=='S' and (q[i]=='G' or q[i]=='C'):
            count=count+1
        elif q[i]=='S' and (p[i]=='G' or p[i]=='C'):
            count=count+1
        else:
            count_miss=count_miss+1
            list_mm.append(str(i+1)+':'+p[i]+"->"+q[i])
    #if count_miss>3:
     #   list_mm=[]
      #  list_mm.append('not found')
    return count,list_mm
def getMMLong(p, q):
    # your code here
    n=len(q)
    count=0
    list_mm=[]
    count_miss=0
    for i in range(0,n):
        if p[i]==q[i]:
            count=count+1
        else:
            count_miss=count_miss+1
            list_mm.append(str(i+1)+':'+p[i]+"->"+q[i])
    #if count_miss>5:
    #   list_mm=[]
    #  list_mm.append('need alignment')
    return count,list_mm
def getMM(haplotype,gene):
    mm_s=''
    scores=[]
    k=len(gene)
    q=str(gene).upper()
    best=0
    mul_s=''
    list_m=[]
    for i in range(len(haplotype)-k+1):
        p=haplotype[i:i+k]
        score,list_mm=Score(q,p)
        if score>best:
            best=score
            mul_s=''
            for mm in list_mm:
                mul_s=mul_s+','+mm
    scores.append(best/k)
    if not mul_s=='':
        mm_s=mul_s
    return mm_s,scores
def readGeneFile(fastafile):
    list_seq=[]
    for seq in SeqIO.parse(fastafile,'fasta'):
        list_seq.append(seq)
    return list_seq
def checkAmpliconWithRef(amplicon,ref_region):
    #ret=blast(amplicon,ref_db)
    #get highest score hit
    # bestscore=0
    # best_pos_start=0
    # best_pos_end=0
    # for h in ret:
    #     if float(h['bitscore'])>bestscore:
    #         bestscore=float(h['bitscore'])
    #         best_pos_start=int(h['sstart'])-int(h['qstart'])
    #         best_pos_end=int(h['send'])+(int(h['qlen'])-int(h['qend']))

    # ref_amplicon= ref_seq[best_pos_start:best_pos_end]
    p,mm,m=FittingAlignment(ref_region,amplicon,1,1)
    return mm,m
def export_file(dict_domains,dict_haplotype,genes,output):
    #    blast_fields={'qseqid':t[0], 'qstart':t[1], 'qend':t[2], 'qlen':t[3],\
    #     'sseqid':t[4], 'sstart':t[5], 'send':t[6], 'slen':t[7], 'sstrand':t[8],\
    #      'length':t[10], 'mismatch':t[11], 'qseq':t[12], 'sseq':t[13],\
    #       'stitle':t[14]}
    file_out=os.path.join(output,'alleles.tsv')
    f=open(file_out,'w')
    f.write('SAMPLE\tGENEID\tALLELE\tPOSITION\tLENGHT\tSEQUENCE\tMISMATCH WITH REF\n')

    #statistic_haplotype={}
    statistic_sample={}

    for s in dict_domains:
        print('export result '+s)
        for g in dict_domains[s]:

            for a in dict_domains[s][g]:
                mm_s=dict_haplotype[g][a['ht']]['mm']
                #f.write(s+'\t'+g.id+'\t'+h+'\t'+dict_haplotype[gene][h]['seq']+'\t'+str(number_haplotype)+'\t'+str(total)+'\t'+str(freq)+'\t'+dict_domains[s][gene.id]+'\t'+dict_haplotype[gene][h]['mm']+'\n')
                f.write(s+'\t'+g+'\t'+a['ht']+'\t'+str(a['start'])+'\t'+str(a['len'])+'\t'+a['seq']+'\t'+mm_s+'\n')
    f.close()

def multithread(num,dict_domain,primers,db_file,output,ref):
    print('start thread '+str(num))
    outfile=os.path.join(output,'domain_primer_'+str(num)+'.tsv')
    export_domain_file(dict_domain,primers,db_file,outfile,ref)
    print('end thread '+str(num))
from multiprocessing import Process
def pipeline(args):
    db_file=args.db
    gene_file=args.genes
    genes=readGeneFile(gene_file)


    dict_haplotype,dict_domain,num_seq=collectHaplotypeFromCorpus(genes,db_file)

    export_file(dict_domain,dict_haplotype,genes,args.output)



def main(arguments=sys.argv[1:]):
    parser = argparse.ArgumentParser(
        prog='genen checker',
        description='Tool for checking gene in genome corpus')
    subparsers = parser.add_subparsers(title='sub command', help='sub command help')
    setup_cmd = subparsers.add_parser(
        'setupdb', description='Collect and filter from genome fasta files', help='Setup db file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    setup_cmd.set_defaults(func=setupdb)
    setup_cmd.add_argument('--samples', help='Directory of collection of samples in fasta format',type=str)
    setup_cmd.add_argument('--output', help='Directory of genome db file',type=str)
    setup_cmd.add_argument('--threadshold', help='Percentage of N character threadshold in sample', type=int)
    run_cmd = subparsers.add_parser(
            'run', description='Check haplotype status with specified gene', help='Check FN',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    run_cmd.set_defaults(func=pipeline)
    run_cmd.add_argument('--genes', help='Gene list in text file (Fasta format)',type=str)
    run_cmd.add_argument('--db', help='Genome db file',type=str)
    run_cmd.add_argument('--output', help='Output folder', type=str)
    run_cmd.add_argument('--threads', help='Number of thread', type=int)
    args = parser.parse_args(arguments)
    return args.func(args)
if __name__ == "__main__":
    main()
