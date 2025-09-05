import sys
import os
import pickle

def get_paired_reads_from_list(reads_list):
    PairedList=[]
    pair=[None,None]
    for c,reads in enumerate(reads_list):
        pair[c%2]=reads
        if (c%2)==1:
            PairedList.append(pair)
            pair=[None,None]
    return PairedList

def get_fileNoExt(filename,extlist):
    filebname=os.path.basename(filename)
    for ext in extlist:
        ext_index=filebname.upper().rfind(ext.upper())
        if ext_index>0:
            filebname=filebname[:ext_index]
    return filebname

def load_readsList(readsList):
    readlist_error=False
    rslt = {}
    reads_dir=""
    Allowed_params={'format':'str','split':'bool','filtered':'bool'}
    default_value={'format':'single','split':False,'filtered':False}
    with open(readsList,'r') as rl:
        for read in rl:
            #nosplit value compatibility
            read=read.replace("nosplit=False","split=True")
            read=read.replace("nosplit=True","split=False")
            commentID = read.find('#')
            if commentID<0:
                tmp = read.strip()
            elif commentID>0:
                tmp = read[:commentID].strip()
            if commentID!=0:
                read_info = tmp.split(',')
                if len(read_info)>=2:
                    readlist_error=(reads_dir=="")
                    infos_dict={}
                    sample_name=read_info[0]
                    infos_dict['reads']=[]
                    infos_dict['params']={}
                    for p,t in Allowed_params.items():
                        infos_dict['params'][p]=convert_value(default_value[p],t)
                    for val in read_info[1:]:
                        if val.count("=")==1:
                            tmp_val=val.split('=')
                            if tmp_val[0] in Allowed_params.keys():
                                infos_dict['params'][tmp_val[0]]=convert_value(tmp_val[1],Allowed_params[tmp_val[0]])
                        else:
                            readfile=os.path.join(reads_dir,val)
                            if os.path.exists(readfile):
                                infos_dict['reads'].append(readfile)
                            else:
                                stdout_print("ERROR '{}' not exists".format(readfile),printInLog=True)
                                readlist_error=True
                                
                    rslt[sample_name] = infos_dict
                elif os.path.isdir(tmp):
                    reads_dir=tmp
    if readlist_error:
        exitProg("ERROR in read file '{}'".format(readsList))
        
    return rslt

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

readsConfig = load_readsList(snakemake.input[0])

with open(snakemake.output[0],'wb') as file:
    pickle.dump(readsConfig,file,protocol=pickle.HIGHEST_PROTOCOL)
