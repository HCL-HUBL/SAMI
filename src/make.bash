#!/bin/bash

gcc -O2 -Wall harvest.c htslib/libhts.a -lm -lz -lbz2 -llzma -lcurl -o harvest
