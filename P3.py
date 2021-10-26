import dbcan_profiler_function as dpf
import sys
import time
from dbcan_profiler_class import *

def Pipeline3(args):
    start_t = time.time()
    dpf.check_args(args)
    dpf.mk_dir(args)
    dpf.trim_reads(args)
    totalreadnum = dpf.read_count(args)
    sys.stderr.write('Prepare Run Time(second): ' + str(time.time() - start_t)+"\n")
    ### assemble read into metegenome first
    dpf.GoldStandard_assemble(args)
    dpf.copyrundbcan_result(args)
    dpf.CAZyCDS(args)
    dpf.Clear_fasta(args)
    dpf.extract_seq_with_seqid(args)
    if args.bowtie:
        dpf.Build_index(args)
        dpf.Run_Bowtie(args)
        paffile1,paffile2 = dpf.To_paf(args)
        #paffile1,paffile2 =["example0/data/bowtie.R1.paf","example0/data/bowtie.R2.paf"]
        paf1 = Paf(paffile1)
        paf2 = Paf(paffile2)
        cazyfpkm,readtable,cazy2seqid = dpf.Cal_FPKM(paf1,paf2,totalreadnum)
        dpf.FPKMToCsv(args,"Bowtie",cazyfpkm,readtable,cazy2seqid)
    dpf.rm_tmp(args)
