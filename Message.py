### 
import sys

Golden_standard = {"HumanMapping":"/mnt/array2/jinfang/CAMISIM_data/HMP/short_read/mapping.table",
        
        } 


Pipeline_mess = {\
        "1":"Pipeline1 is an assembly free appraoch. The input reads will be directly mapped to CAZy protein or CDS sequence.",
        "2":"Pipeline2 is an assembly based appraoch. The reads will be assembled into metagenome by megahit, And then CAZy \
were preicted by dbCAN, next the reads will be mapped to CAZy!",
        "3":"Pipeline3 is for golden standard pipeline.",
        }

def HLError(mess):
    return f"\033[1;31;40m{mess}:\033[0m"

def Pipeline_info(args):
    sys.stderr.write(f"You are run pipeline {args.pipeline}!\n")
    sys.stderr.write(f"{Pipeline_mess[args.pipeline]}\n")
