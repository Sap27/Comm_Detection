#!/bin/bash
# modularity with resolution 0.1
python modularitynetwork.py ../data/subchallenge1/network.dat 0.3 > 3.final

#Core
(python coremoduleidentification.py 3.final ../data/subchallenge1/3_signal_anonym_directed_v3.txt)>3_signal_anonym_directed_v3.txt
sed -i 's/\(.*\)\t/\1/' 3_signal_anonym_directed_v3.txt
sed -i 's/ //g' 3_signal_anonym_directed_v3.txt
