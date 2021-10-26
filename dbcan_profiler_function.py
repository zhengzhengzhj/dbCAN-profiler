import os
import argparse
import sys
import Message as Me
import time
import re
from Bio import SeqIO
import shutil
from dbcan_profiler_class import *
###

def arg_parse():
    parser = argparse.ArgumentParser(description='dbCAN-profiler')
    parser.add_argument('-b', '--bowtie2', dest='bowtie', action='store_true', help='Run Bowtie2?', default=False)
    parser.add_argument('-m', '--minimap2', dest='minimap', action='store_true', help='Run minimap2?', default=False)
    parser.add_argument('-d', '--diamond', dest='diamond', action='store_true', help='Run diamond?', default=False)
    parser.add_argument('-a', '--bwa', dest='bwa', action='store_true', help='Run BWA?', default=False)
    parser.add_argument('-s', '--source', dest='source', action='store_true', help='Store source file?', default=False)
    parser.add_argument('-o', dest='output', type=str, help='Output location?', default="dbcan-profiler/")
    parser.add_argument('-r1', dest='input1', type=str, help='R1 File', required=True)
    parser.add_argument('-r2', dest='input2', type=str, help='R2 File', required=True)
    parser.add_argument('-c', dest='Name_CF', type=str, help='cazyfamily_Name',  default=False)
    parser.add_argument('-bi', dest='bowtieindex', type=str, help='Bowtie index for CAZy database (CDS).',  default=False)
    parser.add_argument('-mi', dest='minimapindex', type=str, help='Minimap2 index for CAZy database (CDS).',  default=False)
    parser.add_argument('-bwai', dest='bwaindex', type=str, help='BWA index for CAZy database (CDS).',  default=False)
    parser.add_argument('-di', dest='diamondindex', type=str, help='Diamond index for CAZy database (protein).',  default=False)
    parser.add_argument('-do', dest='domain', type=str, help='Domain ref or not.',  default=False)
    parser.add_argument('pipeline', help='which pipeline will be used (1,2,3).')
    parser.add_argument('-tr',"--trim_core",help='core number for trim reads! default:12 ',default="12",type=str)
    parser.add_argument('-bc',"--bowtie_core",help='core number for bowtie2 alignment! default:32 ',default="32",type=str)
    parser.add_argument('-gg',"--gold_genome",help='Golden standard genome.',type=str)
    parser.add_argument('-gt',"--gold_genome_type",help='Golden standard genome type. mouse or human.',type=str)
    parser.add_argument('-rm',"--rmtmp",help='remove the tmp files.',action='store_true',default=False)
    parser.add_argument('-st',"--sampletable",help='sample table.',type=str)
    parser.add_argument('-ab',"--abundance",help='abundance.',type=str)
    return parser.parse_args()

def rm_tmp(args):
    if args.rmtmp:
        sys.stderr.write("rm -rf " + args.output + "data/\n")
        #sys.stderr.write(os.getcwd()+"\n")
        shutil.rmtree(args.output + "data")
        #os.system(check_return("rm -rf " + args.output + "data/"))

def check_return(status):
    if status != 0:
        print("\033[1;31;40mError!\033[0m")
        exit(1)

def get_count_reads_fq(file):
    if file.endswith(".gz"):
        r = os.popen("zcat " + file + " | echo $((`wc -l`/4))")
    elif filename.endswith(".fq"):
        r = os.popen("cat " + file + " | echo $((`wc -l`/4))")
    text = r.read()
    r.close()
    return str(int(text))


def get_count_reads_fa(file):
    if file.endswith(".gz"):
        r = os.popen("zcat " + file + " | grep '>' " + " | wc -l")
    elif filename.endswith(".fa"):
        r = os.popen("grep '>' " + file + " | wc -l")
    text = r.read()
    r.close()
    return str(int(text))

### Prepare the output
def check_args(args):
    sys.stderr.write("[***]Checking input paramters!\n")
    if not any([args.bowtie, args.diamond, args.minimap, args.bwa]):
        sys.stderr.write(Me.HLError('Error') + " You must select a program to run(--bowtie, --diamond, --minimap or --bwa).\n")
        exit(1)

    if not args.output.endswith('/'):
        args.output = args.output + "/"
    sys.stderr.write("Writing to: " + args.output + "\n")

def mk_dir(args):
    sys.stderr.write("[***]Preparing output directory!\n")
    check_return(os.system("mkdir -p " + args.output))
    check_return(os.system("mkdir -p " + args.output + "data/"))
    check_return(os.system("mkdir -p " + args.output + "result/"))

## check file type for reads:
## fq,fq.gz
## fa,fa.gz
def check_read_type(filename):
    if filename.endswith("fq") or filename.endswith("fq.gz"):
        return "fq"
    elif filename.endswith("fa") or filename.endswith("fa.gz"):
        return "fa"
    else:
        sys.stderr.write( Me.HLError("Error") + " File type not supported, please provide .fa(fa.gz) or (fq)fq.gz datatype.\n")
        exit(1)

## trim reads using trim_galore

def trim_reads(args):
    sys.stderr.write("[***]Start to trim reads!\n")
    filetype = check_read_type(args.input1)
    if filetype == "fq":
        check_return(os.system(
            "cp " + args.input1 + " " + args.output + "data/input1.fq.gz && cp " + args.input2 + " " + args.output + "data/input2.fq.gz"))
        sys.stderr.write("Running trim_galore\n")
        if os.path.exists(args.output+"data/input1_val_1.fq.gz"): ### not need to trim 
            return 
        check_return(os.system(
            "trim_galore --gzip --no_report_file -j " + args.trim_core  + " --paired " + args.output + "data/input1.fq.gz " + args.output + "data/input2.fq.gz -o " + args.output + "data/"))
        check_return(os.system("rm -rf " + args.output + "data/input1.fq.gz " + args.output + "data/input2.fq.gz"))
    elif filetype == "fa":
        check_return(os.system(
            "cp " + args.input1 + " " + args.output + "data/R1.fa && cp " + args.input2 + " " + args.output + "data/R2.fa"))

def read_count(args):
    filetype = check_read_type(args.input1)
    if filetype == "fq":
        return get_count_reads_fq(args.output+"data/input1_val_1.fq.gz")
    elif filetype == "fa":
        return get_count_reads_fa(args.output + "data/R1.fa")


def To_paf(args):
    check_return(os.system(
        "bioconvert sam2paf --force " + args.output + "data/bowtie.R1.sam " + args.output + "data/bowtie.R1.paf"))
    check_return(os.system(
        "bioconvert sam2paf --force " + args.output + "data/bowtie.R2.sam " + args.output + "data/bowtie.R2.paf"))
    return [args.output + "data/bowtie.R1.paf", args.output + "data/bowtie.R2.paf"]
    
    ## test one time
    ## similary result to bioconvert
    #sam2paf(args.output + "data/bowtie.R1.sam" , args.output + "data/bowtie.R1.paf.own")
    #sam2paf(args.output + "data/bowtie.R2.sam" , args.output + "data/bowtie.R2.paf.own")

def Run_Bowtie(args):
    bowtie_start = time.time()
    sys.stderr.write("Running bowtie2\n")
    filetype = check_read_type(args.input1)
    ## map fq file
    if filetype == "fq":
        check_return(os.system(
            "bowtie2 -p "+ args.bowtie_core + " -x " +args.bowtieindex + " " + args.output + "data/input1_val_1.fq.gz -S " + args.output + "data/bowtie.R1.sam"))
        check_return(os.system(
            "bowtie2 -p " + args.bowtie_core + " -x " +args.bowtieindex + " " + args.output + "data/input2_val_2.fq.gz -S " + args.output + "data/bowtie.R2.sam"))
    ## map fa file
    elif filetype == "fa":
        check_return(os.system(
            "bowtie2 -p "+ args.bowtie_core + " -x " + args.bowtieindex + " -f " + args.output + "data/R1.fa -S " + args.output + "data/bowtie.R1.sam"))
        check_return(os.system(
            "bowtie2 -p "+ args.bowtie_core + " -x " + args.bowtieindex + " -f " + args.output + "data/R2.fa -S " + args.output + "data/bowtie.R2.sam"))
    return None

## .R2 alignment result filename with paf format
### output1: csv file for 

def merge_two_dicts(x, y):
    z = x.copy()
    z.update(y)
    return z

def getSeqlen(paf1,paf2):
    x = paf1.GetSeqLen()
    y = paf2.GetSeqLen()
    return merge_two_dicts(x,y)

def getCazySeqId(paf1,paf2):
    cazy2seqid = {}
    paf1.CAZy2SeqID(cazy2seqid)
    paf2.CAZy2SeqID(cazy2seqid)
    for cazy in cazy2seqid:
        cazy2seqid[cazy] = set(cazy2seqid[cazy])
    return cazy2seqid

def getSeqReadID(paf1,paf2):
    seqid2readid = {}
    paf1.SeqID2ReadID(seqid2readid)
    paf2.SeqID2ReadID(seqid2readid)
    return seqid2readid

def SeqReadCount(seqid2readid):
    ## 0.5 two reads should be one because of the input is pair end
    return{seqid:len(seqid2readid[seqid])*0.5 for seqid in seqid2readid}

###          reads per transcript
### FPKM  =  -------------------------------
###          total reads   transcript length 
###          ----------- X -----------------
###          10E6          1000

def SequenceFPKM(readtable,seq2len,totalreadnumber):
    seqfpkm = {}
    for seqid in readtable:
        tmp_total_read = float(totalreadnumber)/pow(10,6)
        tmp_trans_len  = float(seq2len[seqid])/1000
        read_count = float(readtable[seqid])
        tmp_fpkm = read_count/tmp_total_read/tmp_trans_len
        seqfpkm[seqid] = tmp_fpkm ### 0.5 because we used pair ends, two reads should be count as 1
    return seqfpkm

def CAZyFPKM(seqfpkm,cazy2seqid):
    cazyfpkm = {}
    for cazy in cazy2seqid:
        tmp_fpkm = 0.
        for seqid in cazy2seqid[cazy]:
            tmp_fpkm += float(seqfpkm[seqid])
        cazyfpkm[cazy] = tmp_fpkm
    return cazyfpkm


def outdict(aa):
    [print (f"{a}\t{aa[a]}") for a in aa]

### output2:dict CAZy family2FPKM, one-one

def Cal_FPKM(paf1,paf2,totalreadnumber):
    #ReadId1 = paf1.GetReadId()
    #ReadId2 = paf2.GetReadId()
    #Hitid1   = paf1.GetSeqId()
    #Hitid2   = paf2.GetSeqId()
    ## get sequence length from paf
    seq2len = getSeqlen(paf1,paf2)
    
    # get CAZy family to seq mapping table: CAZy ID 2 protein ID
    cazy2seqid = getCazySeqId(paf1,paf2)
    # outdict_list(cazy2seqid)

    ## get SeqID2ReadID to generate mapping table: protein ID 2 read ID
    seqid2readid = getSeqReadID(paf1,paf2)
    ### read table: protein ID 2 read count 
    readtable = SeqReadCount(seqid2readid)
    ## outdict(readtable)
    ## calculate fpkm for each protein seq
    seqfpkm = SequenceFPKM(readtable,seq2len,totalreadnumber)
    ## outdict(seqfpkm)
    cazyfpkm = CAZyFPKM(seqfpkm,cazy2seqid)
    #outdict(cazyfpkm)
    return cazyfpkm,readtable,cazy2seqid

def CAZyReadCount(cazyid,cazy2seqid,readtable):
    tmp_sum = 0
    for seqid in cazy2seqid[cazyid]:
        tmp_sum += readtable[seqid]
    return tmp_sum

def FPKMToCsv(args,tool,cazyfpkm,readtable,cazy2seqid):
    outfilename = args.output + f"/result/{tool}.CAZy.FPKM.csv"
    with open(outfilename,'w') as f:
        f.write(f"CAZy,SeqNum,ReadCount,FPKM\n")
        for cazyid in cazyfpkm:
            seqnum = len(cazy2seqid[cazyid])
            readcount = CAZyReadCount(cazyid,cazy2seqid,readtable)
            fpkm = cazyfpkm[cazyid]
            f.write(f"{cazyid},{seqnum},{readcount},{fpkm}\n")

def Megahit_assemble(args):
    args.R1 = args.output+"data/input1_val_1.fq.gz"
    args.R2 = args.output+"data/input2_val_2.fq.gz"
    args.megahit_out = args.output + "data/megahit/"
    #check_return(os.system("megahit -1 " + args.R1 + " -2 " + args.R2 + " -o " + args.megahit_out))
    args.megahit_fasta = "final.contigs.fa"

def GoldStandard_assemble(args):
    args.R1 = args.output+"data/input1_val_1.fq.gz"
    args.R2 = args.output+"data/input2_val_2.fq.gz"
    args.megahit_out = args.output + "data/"
    args.megahit_fasta = "golden.fasta"
    shutil.copy(args.gold_genome, args.megahit_out + args.megahit_fasta)
    
def rundbcan(args):
    args.contig = args.megahit_out + args.megahit_fasta
    args.dbcan_out  = args.output + "data/dbcan_out/"
    args.dbcan_db   = "/mnt/array2/pengxiang/annotation/pipeline2_Ge_merge/dbcan/db"
    #args.dbcan_db   = "/mnt/array2/jinfang/dbcan_da/dbcan_database/dbCAN_db"
    ## need to be addressed
    check_return(os.system(
        "run_dbcan.py " + args.contig + " meta --out_dir " + args.dbcan_out + " --db_dir " + args.dbcan_db + " --dia_cpu 12  --hmm_cpu 12 --hotpep_cpu 12 --tf_cpu 12 --tf_cpu 12 "))

def copyrundbcan_result(args):
    args.contig = args.megahit_out + args.megahit_fasta
    args.dbcan_out  = args.output + "data/dbcan_out/"
    args.dbcan_db   = "/mnt/array2/pengxiang/annotation/pipeline2_Ge_merge/dbcan/db"
    #args.dbcan_db   = "/mnt/array2/jinfang/dbcan_da/dbcan_database/dbCAN_db"
    ### for HMP rundbCAN result
    if args.gold_genome_type == "human":
        check_return(os.system(f"cp -r /mnt/array2/jinfang/dbcan_profiler/HMP_rundbcan/data/dbcan_out {args.dbcan_out}"))
    
    ### for mouse gut rundbCAN result
    elif args.gold_genome_type == "mouse":
        check_return(os.system(f"cp -r /mnt/array2/jinfang/dbcan_profiler/mousegut_rundbcan/data/dbcan_out {args.dbcan_out}"))
    ### for else, need to rerun 
    elif args.gold_genome_type == "previous":
        check_return(os.system(f"cp -r /mnt/array2/jinfang/dbcan_profiler/previous_work/previous/data/dbcan_out {args.dbcan_out}"))
    #else:
    #    check_return(os.system(
    #    "run_dbcan.py " + args.contig + " meta --out_dir " + args.dbcan_out + " --db_dir " + args.dbcan_db + " --dia_cpu 12  --hmm_cpu 12 --hotpep_cpu 12 --tf_cpu 12 --tf_cpu 12 "))
    

### filter CAZy CDS
def CAZyCDS(args):
    args.new_gff = args.output + "data/new.gff"
    args.tmpcds = args.output + "data/tmp.cds.fa"
    gff_filter(args.dbcan_out + "prodigal.gff",args.dbcan_out+"overview.txt",args.new_gff)
    check_return(os.system(
        "gff2bed < " + args.new_gff + " | sed \"s/cds-//\" | awk \'{if ($8==\"CDS\") print $0}\' | bedtools getfasta -fi " + args.contig + " -bed - -s -name > " + args.tmpcds ))
    return None

## 
## >k141_1003_1|HMMER=-|Hotpep=GH0(91)|DIAMOND=-(-)

def Clear_CAZy(tmpcazy):
    cazy = []
    for c in tmpcazy:
        if c[0] == "-":
            continue
        else:
            cs = c.split("(")[0].split("_")[0] ## for GH15_3(0-100)
            cazy.append(cs)
    return set(cazy)

def Clear_id(fastaid):
    tmpdes,hmmers,hotpep,diamond = fastaid.split("|")
    hmmers = hmmers[len("HMMER="):].split("+")
    hotpep = hotpep[len("Hotpep="):].split("+")
    diamond = diamond[len("DIAMOND="):].split("+")
    hmmers_cazy = Clear_CAZy(hmmers)
    hotpep_cazy = Clear_CAZy(hotpep)
    diamond_cazy = Clear_CAZy(diamond)
    CaZy = list((hmmers_cazy & hotpep_cazy)|(hotpep_cazy&diamond_cazy)|(diamond_cazy&hmmers_cazy))
    return tmpdes+"|"+"|".join(CaZy)

def Clear_fasta(args):
    seq2fasta = {record.id:record for record in SeqIO.parse(args.tmpcds,'fasta')}
    records = []
    for seqid in seq2fasta:
        newid = Clear_id(seqid)
        if newid.endswith("|"): ## no CAZy
            continue
        else:
            seq = seq2fasta[seqid]
            seq.id = newid
            records.append(seq)
    args.megahit_cds = args.output + "data/megahit.cds.fa"
    SeqIO.write(records,args.megahit_cds,'fasta')

def Build_index(args):
    args.bowtieindex = args.megahit_cds
    check_return(os.system("bowtie2-build " + args.megahit_cds + " " + args.megahit_cds))
    
    #if args.gold_genome_type == "human":
    #    check_return(os.system(f"cp -r /mnt/array2/jinfang/dbcan_profiler/HMP_rundbcan/data/megahit.cds.fa.*.bt2 {args.output}data"))
    #
    #### for mouse gut rundbCAN result
    #elif args.gold_genome_type == "mouse":
    #    check_return(os.system(f"cp -r /mnt/array2/jinfang/dbcan_profiler/mousegut_rundbcan/data/megahit.cds.fa.*.bt2 {args.output}data"))
    #
    #### for else
    #else:
    #    check_return(os.system("bowtie2-build " + args.megahit_cds + " " + args.megahit_cds))
    

### filter CAZy CDS
#def CAZyCDS(args):
#    args.new_gff = args.output + "data/new.gff"
#    check_return(os.system(f"cp {} {args.megahit_cds}"))

# Input(argv1): prodigal.gff file from dbCAN
# Input(argv2): overview.txt file from dbCAN
# Output(argv3): a filtered gff file

## gff_filter(gff_file=sys.argv[1], overview_file=sys.argv[2], output_file=sys.argv[3])

def gff_filter(gff_file, overview_file, output_file):
    with open(gff_file, 'r') as f:
        data = f.readlines()
    f.close()

    with open(overview_file, 'r') as f:
        raw_refs = f.readlines()
    f.close()

    ref_id = {}
    for raw_ref in raw_refs[1:]:
        id = raw_ref.split()[0]
        HMMER = raw_ref.split()[1]
        Hotpep = raw_ref.split()[2]
        DIAMOND = raw_ref.split()[3]
        ref_id[id] = {'HMMER': HMMER, 'Hotpep': Hotpep, 'DIAMOND': DIAMOND}

    result = []
    for row in data:
        if row.startswith('#'):
            continue
        else:
            brief = row.split('\n')[0]
            columns = brief.split()
            des = columns[-1]
            id = des.split(';partial')[0]
            id_suffix = id.split('_')[1]
            new_id = columns[0] + '_' + id_suffix

            if new_id in ref_id.keys():
                # print(new_id)
                reference = ref_id[new_id]
                output_id = 'ID=' + new_id + '|HMMER=' + reference['HMMER'] + '|Hotpep=' + reference[
                    'Hotpep'] + '|DIAMOND=' + reference['DIAMOND']
                new_brief = brief.replace(id, output_id, 1) + '\n'
                result.append(new_brief)

    with open(output_file, 'w') as file:
        file.writelines(result)
    file.close()

def check_result(new_id):
    with open('uniInput', 'r') as f:
        raw_checks = f.readlines()
    f.close()

    checks = []
    for raw_check in raw_checks:
        if raw_check.startswith('>'):
            raw_id = raw_check.split(' # ')[0]
            id = raw_id.split('>')[1]
            checks.append(id)

    # if new_id in checks:
    #     print('checked:' + new_id)
    # else:
    #     print(new_id)

def remove_not_hit(path, new_path):
    file1 = open(path, 'r')
    f = open(new_path, 'w')
    for line in file1.readlines():
        x = line.split("	")
        if x[2] != "*" and x[4] != 255:
            f.write(line)
    f.close()

def Build_diamond(args):
    #args.bowtieindex = args.megahit_cds
    #check_return(os.system("bowtie2-build " + args.megahit_cds + " " + args.megahit_cds))
    args.diamondindex = args.megahit_cds
    args.megahit_prot = args.output + "data/megahit.prot"
    check_return(os.system(f"seqkit translate {args.megahit_cds} > {args.megahit_prot}"))
    check_return(os.system(f"diamond makedb --in {args.megahit_prot} -d {args.diamondindex}"))

def Run_Diamond(args):
    diamond_start = time.time()
    filetype = check_read_type(args.input1)
    args.paf1 = args.output + "data/diamond.new.R1.paf"
    args.paf2 = args.output + "data/diamond.new.R2.paf"
    if filetype == "fq":
        sys.stderr.write("Running seqtk\n")
        check_return(os.system(
            "seqtk seq -a " + args.output + "data/input1_val_1.fq.gz  > " + args.output + "data/R1.fa"))
        check_return(os.system(
            "seqtk seq -a " + args.output + "data/input2_val_2.fq.gz  > " + args.output + "data/R2.fa"))
    sys.stderr.write("Running diamond\n")
    # CAZyDB.07312020.fa
    check_return(os.system(
        "diamond blastx -k 1 --strand both --evalue 1e-10 --query " + args.output + "data/R1.fa --db " +  args.diamondindex + " --threads 32 --out " + args.output + "data/diamond.R1.paf --outfmt 103"))
    check_return(os.system(
        "diamond blastx -k 1 --strand both --evalue 1e-10 --query " + args.output + "data/R2.fa --db "+ args.diamondindex +  " --threads 32 --out " + args.output + "data/diamond.R2.paf --outfmt 103"))
    # TODO: change Database

    sys.stderr.write("Running Removing unhit\n")
    remove_not_hit(args.output + "data/diamond.R1.paf", args.output + "data/diamond.new.R1.paf")
    remove_not_hit(args.output + "data/diamond.R2.paf", args.output + "data/diamond.new.R2.paf")
    sys.stderr.write("Running diamond analysis\n")
    diamond_end = time.time()
    sys.stderr.write('Run Diamond Time(second): ' + str(diamond_end - diamond_start)+'\n')

def Run_Bwa(args):
    bwa_start = time.time()
    filetype = check_read_type(args.input1)
    args.paf1 = args.output + "data/bwa.new.R1.paf"
    args.paf2 = args.output + "data/bwa.new.R2.paf"
    sys.stderr.write("Running BWA\n")
    if filetype == "fq":
        check_return(os.system(
            "bwa mem -t 32 " + args.bwaindex + " " + args.output + "data/input1_val_1.fq.gz >  " + args.output + "data/bwa_R1.sam"))
        check_return(os.system(
            "bwa mem -t 32 " + args.bwaindex +" " + args.output + "data/input2_val_2.fq.gz >  " + args.output + "data/bwa_R2.sam"))
    else:
        check_return(os.system(
            "bwa mem -t 32 " + args.bwaindex + " " + args.output + "data/R1.fa >  " + args.output + "data/bwa_R1.sam"))
        check_return(os.system(
            "bwa mem -t 32 "+ args.bwaindex +" "  + args.output + "data/R2.fa >  " + args.output + "data/bwa_R2.sam"))

    sys.stderr.write("Running bioconvert\n")
    check_return(os.system(
        "bioconvert sam2paf " + args.output + "data/bwa_R1.sam " + args.output + "data/bwa.R1.paf"))
    check_return(os.system(
        "bioconvert sam2paf " + args.output + "data/bwa_R2.sam " + args.output + "data/bwa.R2.paf"))
    sys.stderr.write("Removing unaligned reads\n")
    remove_not_hit(args.output + "data/bwa.R1.paf", args.output + "data/bwa.new.R1.paf")
    remove_not_hit(args.output + "data/bwa.R2.paf", args.output + "data/bwa.new.R2.paf")
    bwa_end = time.time()
    sys.stderr.write('Bwa Run Time(second): ' + str(bwa_end - bwa_start))

def Run_Minimap(args):
    minimap_start = time.time()
    filetype = check_read_type(args.input1)
    args.paf1 = args.output + "data/minimap.R1.paf"
    args.paf2 = args.output + "data/minimap.R2.paf"
    sys.stderr.write("Running Minimap2\n")
    if filetype == "fq":
        check_return(os.system(
            "minimap2 -t 32 -x sr " +args.minimapindex  + " "  + args.output + "data/input1_val_1.fq.gz  >  " + args.output + "data/minimap.R1.paf"))
        check_return(os.system(
            "minimap2 -t 32 -x sr " +args.minimapindex  + " "  + args.output + "data/input2_val_2.fq.gz  >  " + args.output + "data/minimap.R2.paf"))
    else:
        check_return(os.system(
            "minimap2 -t 32 -x sr " + args.minimapindex + " " + args.output + "data/R1.fa  >  " + args.output + "data/minimap.R1.paf"))
        check_return(os.system(
            "minimap2 -t 32 -x sr " + args.minimapindex + " " + args.output + "data/R2.fa  >  " + args.output + "data/minimap.R2.paf"))
    sys.stderr.write("Removing unaligned reads\n")
    #remove_not_hit(args.output + "data/minimap.R1.paf", args.output + "data/minimap.new.R1.paf")
    #remove_not_hit(args.output + "data/minimap.R2.paf", args.output + "data/minimap.new.R2.paf")
    minimap_end = time.time()
    sys.stderr.write('Minimap Run Time(second): ' + str(minimap_end - minimap_start))

def Read_abundance(filename):
    return [line.split()[0] for line in open(filename) if float(line.split()[1]) > 0]

def Get_abundance(filename):
    return dpc.Abundance(filename)

def expression_seqid(abundance,sampletable):
    seqids = []
    for otu in sampletable:
        if otu in abundance:
            seqids += sampletable[otu]
    return seqids

def cds_parse(args,seqids):
    records = [] ; tmpseqs = [seq.seqid for seq in seqids]
    for record in SeqIO.parse(args.megahit_cds,'fasta'):
        seqid = record.id.split("_")[0]
        if seqid in tmpseqs:
            records.append(record)
    SeqIO.write(records,args.megahit_cds,'fasta')

def extract_seq_with_seqid(args):
    sampletable =  Golden_smaple(args.sampletable).OTUdict()
    abundance = Read_abundance(args.abundance)
    expressedseqid = expression_seqid(abundance,sampletable)
    cds_parse(args,expressedseqid)

