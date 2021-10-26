### P4 for previous work



import dbcan_profiler_function as dpf
import sys
import time
from dbcan_profiler_class import *

def Pipeline4(args):
    start_t = time.time()
    dpf.check_args(args)
    dpf.mk_dir(args)
    dpf.trim_reads(args)
    totalreadnum = dpf.read_count(args)
    sys.stderr.write('Prepare Run Time(second): ' + str(time.time() - start_t)+"\n")
    ### assemble read into metegenome first
    dpf.GoldStandard_assemble(args)
    dpf.copyrundbcan_result(args)
