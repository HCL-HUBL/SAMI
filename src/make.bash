#!/bin/bash

gcc -O2 -Wall collect-gaps.c htslib/libhts.a -lm -lz -lbz2 -llzma -lcurl -o collect-gaps
