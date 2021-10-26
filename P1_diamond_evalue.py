import dbcan_profiler_function as dpf
import sys
import time
from dbcan_profiler_class import *

def Pipeline1(args):
    start_t = time.time()
    sys.stderr.write("[***]Checking input paramters!\n")
    dpf.check_args(args)
    sys.stderr.write("[***]Preparing output directory!\n")
    dpf.mk_dir(args)
    sys.stderr.write("[***]Start to trim reads!\n")
    dpf.trim_reads(args)
    totalreadnum = dpf.read_count(args)
    sys.stderr.write('Prepare Run Time(second): ' + str(time.time() - start_t)+"\n")
    if args.bowtie:
        dpf.Run_Bowtie(args)
        paffile1,paffile2 = dpf.To_paf(args)
        #paffile1,paffile2 =["example0/data/bowtie.R1.paf","example0/data/bowtie.R2.paf"]
        paf1 = Paf(paffile1)
        paf2 = Paf(paffile2)
        cazyfpkm,readtable,cazy2seqid = dpf.Cal_FPKM(paf1,paf2,totalreadnum)
        dpf.FPKMToCsv(args,"Bowtie",cazyfpkm,readtable,cazy2seqid)
    if args.diamond:
        dpf.Run_Diamond(args)
        paf1 = Paf(args.paf1)
        paf2 = Paf(args.paf2)
        cazyfpkm,readtable,cazy2seqid = dpf.Cal_FPKM(paf1,paf2,totalreadnum)
        dpf.FPKMToCsv(args,"Diamond",cazyfpkm,readtable,cazy2seqid)
    if args.bwa:
        dpf.Run_Bwa(args)
        paf1 = Paf(args.paf1)
        paf2 = Paf(args.paf2)
        cazyfpkm,readtable,cazy2seqid = dpf.Cal_FPKM(paf1,paf2,totalreadnum)
        dpf.FPKMToCsv(args,"Bwa",cazyfpkm,readtable,cazy2seqid)
    if args.minimap:
        dpf.Run_Minimap(args)
        paf1 = Paf(args.paf1)
        paf2 = Paf(args.paf2)
        cazyfpkm,readtable,cazy2seqid = dpf.Cal_FPKM(paf1,paf2,totalreadnum)
        dpf.FPKMToCsv(args,"Minimap",cazyfpkm,readtable,cazy2seqid)
