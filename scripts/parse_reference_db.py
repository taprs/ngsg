#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Adapted from mathieu-genete/NGSgenotyp commit e1376e8 by taprs.
No lines removed, only commented out and added some new code.
Uses code from genotyp.py, load_dump.py, and checkDB.py scripts.
Checks fasta reference (1st argument) and writes a pickle file (2nd argument) for kmer filtering.
"""


import sys
import os
import argparse
import encodings
import subprocess
from Bio import SeqIO
from Bio import Seq
import codecs
import uuid
import gzip
import pickle

__version__ = "0.5"
AppName = "checkDB"

def run(infasta):
    outcheck=""
    outcheck+="=== {} v{} ===\n".format(AppName,__version__)
    outcheck+="input fasta: {}\n".format(infasta)
    fileOut,ftype,crlf=file_type(infasta)
    outencoded=convertToUnix(infasta)
    encodedfasta="".join(outencoded)
    fastadict=parse_fasta(encodedfasta)
    outcheck+="=> START CHECK\n"
    fastaOK,logout=analyse_fasta(fastadict)
    outcheck+=logout
    return fastaOK,outcheck
    
def analyse_fasta(fastadict):
    logout=""
    fastaOK=True
    NotAllowedChar=[',','.','/',' ',';','#','@',':']
    #NotAllowedChar=['.','/',' ',';','#','@']
    DNAalphabet=['A', 'C', 'B', 'D', 'G', 'H', 'K', 'M', 'N', 'S', 'R', 'T', 'W', 'V', 'Y', 'X']
    AllowedParameters=['Paralog','specie','grpRef','HaploId','gprId','allelePart','grpRefPart']
    comma_exception=['allelePart','grpRefPart']
    deprecatedParameters=['HaploId','gprId']
    warnerr={'E':'ERROR - ','W':'WARNING - '}

    #?:{'title':'?','err':[]}
    errorsout={1:{'title':'unsupported characters in sequence id:','err':[]},2:{'title':'unsupported characters in sequence string:','err':[]},3:{'title':'sequence id format','err':[]},4:{'title':'Gene ID duplicates','err':[]}}
    IdList=[]

    logout+="Sequences number: {}\n".format(len(fastadict))

    for idnbr,value in fastadict.items():
        
        parseID=value['description'].split("|")
        geneId=parseID[0]
        
        #===== search Id errors=====
        idCharERROR=False
        for c in NotAllowedChar:
            if c in value['description']:
                param_with_comma=[v.split('=')[0] for v in parseID if ',' in v and v.split('=')[0] not in comma_exception]
                if len(param_with_comma)==0:
                    continue
                cpos=[str(i) for i,ltr in enumerate(value['description']) if ltr==c]
                positionStr=list("."*len(value['description']))
                for p in cpos:
                    positionStr[int(p)]="^"
                lineheader="\n\t{we}id number: {idnb}\t".format(we=warnerr['E'],idnb=idnbr)
                errorsout[1]['err'].append("{head}{seqid} -- character: '{ch}' -- position(s): {pos}\n\t{lnh} {posstr}".format(head=lineheader,seqid=value['description'],ch=c,pos=",".join(cpos),posstr="".join(positionStr),lnh=" "*len(lineheader)))
                idCharERROR=True
                fastaOK=False
        
        if geneId not in IdList:
            IdList.append(parseID[0])
        else:
            errorsout[4]['err'].append("\t{}Duplicate found for id: '{}' -- id number: {}".format(warnerr['E'],geneId,idnbr))
            fastaOK=False

        #===== search parametters errors=====
        paramsSTR=parseID[1:]
        for p in paramsSTR:
            tmp=p.split('=')
            if len(tmp)==1 and tmp[0]=="":
                errorsout[3]['err'].append("\t{}id number: '{}' -- id: {} -- EMPTY parameter: '||'".format(warnerr['E'],idnbr,geneId))
                fastaOK=False
            elif len(tmp)==1 or len(tmp)>2:
                errorsout[3]['err'].append("\t{}id number: '{}' -- id: {} -- param error: {}".format(warnerr['E'],idnbr,geneId,p))
                fastaOK=False
            else:
                if tmp[0] not in AllowedParameters:
                    errorsout[3]['err'].append("\t{}id number: '{}' -- id: {} -- paramameter not allowed: {}".format(warnerr['E'],idnbr,geneId,tmp[0]))
                    fastaOK=False
                if tmp[0] in deprecatedParameters:
                    errorsout[3]['err'].append("\t{}id number: '{}' -- id: {} -- deprecated parameter: {}".format(warnerr['W'],idnbr,geneId,tmp[0]))
                if tmp[0]=='grpRef':
                    grpRefVal=tmp[1]
                    grpRefSplit=grpRefVal.split("-")
                    if len(grpRefSplit)==0 or len(grpRefSplit)>2:
                        errorsout[3]['err'].append("\t{}id number: '{}' -- id: {} -- grpRef value: {}".format(warnerr['E'],idnbr,geneId,grpRefVal))
                        fastaOK=False
                    elif len(grpRefSplit)==1 and len(grpRefVal)==5:
                        errorsout[3]['err'].append("\t{}id number: '{}' -- id: {} -- grpRef value deprecated use 'Hg-h': {}".format(warnerr['W'],idnbr,geneId,grpRefVal))

        #===== search sequences errors=====
        sequence=value['seq']
        SequenceNotAllowedAlphabet=[b for b in set(list(sequence.upper())) if b not in DNAalphabet]

        for base in SequenceNotAllowedAlphabet:
            bpos=[str(i) for i,ltr in enumerate(sequence) if ltr.upper()==base]
            errorsout[2]['err'].append("\t{}id number: '{}' -- id: {} -- base '{}' not in allowed alphabet -- position: {} bp".format(warnerr['E'],idnbr,geneId,base,",".join(bpos)))
            fastaOK=False

    for err in errorsout.values():
        if len(err['err'])>0:
            logout+="{}\n".format(err['title'])
            for chaine in err['err']:
                logout+="{}\n".format(chaine)

    return fastaOK,logout

def parse_fasta(outencoded):
    fastadict={}
    tmpname=os.path.join("/tmp",str(uuid.uuid4()))
    with open(tmpname,"w") as tmpfasta:
        for line in outencoded:
            tmpfasta.write(line)
    fasta=SeqIO.parse(tmpname,'fasta')
    seqNbr=0
    for rec in fasta:
        seqNbr+=1
        fastadict[seqNbr]={'id':str(rec.id),'description':str(rec.description),'seq':str(rec.seq)}
    os.remove(tmpname)
    return fastadict

def parse_GrpRef(grpRef,refname):
    grpID=None
    haploID=None
    if grpRef and grpRef[0].upper()=="H":
        if grpRef.count('-')==1:
            tmp=grpRef[1:].split('-')
            grpID=int(tmp[0])
            haploID=int(tmp[1])
        elif grpRef.count('-')==0:
            grpID=int(grpRef[1])
            haploID=int(grpRef[2:])
        else:
            stdout_print("WARNING grpRef parameter '{}' for reference '{}' not correctly formated".format(grpRef,refname),printInLog=True)
    return grpID,haploID
            
def file_type(infile):
    AllowedTypes=set(encodings.aliases.aliases.values())
    cmd="file -b {}".format(infile)
    job = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout,sdterr = job.communicate()
    stdout=stdout.decode("utf-8")
    ftype=None
    for t in AllowedTypes:
        if t.upper() in stdout.strip().upper():
            ftype=t
    crlf=("CRLF" in stdout.strip().upper())
    return stdout.strip(),ftype,crlf

def convertToUnix(infile):
    MsDOS_endline='\r\n'
    convertedFasta=[]
    with codecs.open(infile,'r') as file:
        for line in file.readlines():
            IS_MsDOS_EL=line[-2:]==MsDOS_endline
            if IS_MsDOS_EL:
                convertedFasta.append("{}\n".format(line[:-2]))
            else:
                convertedFasta.append("{}\n".format(line.strip()))
    return convertedFasta

def convert_value(value,vtype):
    if vtype=='bool':
        if value=="True": return True
        if value=="False": return False
        return bool(value)
    if value is not None:
        if vtype=='str':
            return str(value)
        elif vtype=='int':
            return int(value)
        elif vtype=='float':
            return float(value)
        elif vtype=='list':
            return list(value.split(','))

    return value
        
def pickle_dumpDic(outpickle,dic):
    pf = gzip.open(outpickle,'wb')
    pickle.dump(dic,pf,protocol=pickle.HIGHEST_PROTOCOL)
    pf.close()
    
def parse_referenceDB(infastafile, out_pkfile):
    # stdout_print(" -- parse reference file {} --".format(infastafile),printInLog=True)
    print(" -- parse reference file {} --".format(infastafile))
    dictfasta={}
    allelePart={}
    subAlleleList=[]
    fastafile_basename=os.path.basename(infastafile)
    
    #TO REMOVE YAML FILE
    # out_yamlfile=os.path.join(config['OutFolders']['reference_DB_fld'],"{}.yaml".format(fastafile_basename[:fastafile_basename.rfind('.')]))
    
    # out_pkfile=os.path.join(config['OutFolders']['reference_DB_fld'],"{}.pk".format(fastafile_basename[:fastafile_basename.rfind('.')]))
    #Parameters: Paralog, specie, grpRef
    Allowed_parameters={'Paralog':'bool', 'specie':'str', 'grpRef':'str','allelePart':'list','grpRefPart':'list'}

    reffasta=SeqIO.parse(infastafile,'fasta')
    for rec in reffasta:
        seqID_list=str(rec.id).split('|')
        seqID=seqID_list[0]
        dictfasta[seqID]={}
        dictfasta[seqID]['seq']=str(rec.seq)
        dictfasta[seqID]['len_seq']=len(str(rec.seq))
        dictfasta[seqID]['combined']=False
        for p,t in Allowed_parameters.items():
            dictfasta[seqID][p]=convert_value(None,t)
        
        if len(seqID_list)>1:
            for values in seqID_list[1:]:
                if values.count("=")==1:
                    tmp=values.split('=')
                    param_name=tmp[0]
                    param_val=tmp[1]

                    if param_name in Allowed_parameters.keys():
                        dictfasta[seqID][param_name]=convert_value(param_val,Allowed_parameters[param_name])
    for refname,refparam in dictfasta.items():
        grpID,haploID=parse_GrpRef(refparam['grpRef'],refname)
        if grpID is not None and haploID is not None:
            refparam['grpRef']="H{}-{}".format(grpID,haploID)
        dictfasta[refname]['grpIDPart']=None
        dictfasta[refname]['haploIDPart']=None
        
        if dictfasta[refname]['grpRefPart'] is not None and len(dictfasta[refname]['grpRefPart'])>0:
            subAlleleList.append(refname)
            dictfasta[refname]['grpID']=None
            dictfasta[refname]['haploID']=None
            dictfasta[refname]['grpIDPart']=[]
            dictfasta[refname]['haploIDPart']=[]
            for grpRefPart in dictfasta[refname]['grpRefPart']:
                grpIDPart,haploIDPart=parse_GrpRef(grpRefPart,refname)
                dictfasta[refname]['grpIDPart'].append(grpIDPart)
                dictfasta[refname]['haploIDPart'].append(haploIDPart)
        else:
            dictfasta[refname]['grpID']=grpID
            dictfasta[refname]['haploID']=haploID

    for subAllele in subAlleleList:
        for allele,group in zip(dictfasta[subAllele]['allelePart'],dictfasta[subAllele]['grpRefPart']):
            if allele not in allelePart.keys():
                allelePart[allele]={'subAllele':[],'len_seq':0,'subGrp':None,'subGrpID':None,'subHaploID':None}
            allelePart[allele]['subAllele'].append(subAllele)
            allelePart[allele]['subGrp']=group
            grpID,haploID=parse_GrpRef(group,allele)
            allelePart[allele]['subGrpID']=grpID
            allelePart[allele]['subHaploID']=haploID
            allelePart[allele]['len_seq']+=dictfasta[subAllele]['len_seq']
            
    for allele, alleles_val in allelePart.items():
        dictfasta[allele]={'combined':True,'Paralog': False, 'allelePart': None, 'grpID': alleles_val['subGrpID'], 'grpIDPart': None, 'grpRef': alleles_val['subGrp'], 'grpRefPart': None, 'haploID': alleles_val['subHaploID'], 'haploIDPart': None, 'len_seq': alleles_val['len_seq'], 'seq': None, 'specie': None}
        
        
    #TO REMOVE YAML FILE
    # ld.yaml_dumpDic(out_yamlfile,dictfasta)
    
    # ld.pickle_dumpDic(out_pkfile,dictfasta)
    pickle_dumpDic(out_pkfile,dictfasta)
        
    # return out_pkfile
            
ok,log = run(snakemake.input[0])
print(log, file=sys.stderr)
if not ok:
    sys.exit(1)
pkfile = parse_referenceDB(snakemake.input[0], snakemake.output[0])
