import sys
import argparse
import os
import time
import Message as Me
import dbcan_profiler_function as dpf
### define the main funciton of dbCAN-profiler
### 

import P1 as P1
import P2 as P2
import P3 as P3

def main():
    args = dpf.arg_parse()
    Me.Pipeline_info(args)
    if args.pipeline ==  "1":
        P1.Pipeline1(args)
    if args.pipeline == "2":
        P2.Pipeline2(args)
    if args.pipeline == "3":
        P3.Pipeline3(args)

if __name__ == "__main__":
    main()
