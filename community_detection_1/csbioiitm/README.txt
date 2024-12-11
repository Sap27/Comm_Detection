# modularity with resolution 0.1
python modularitynetwork.py ../data/subchallenge1/3_signal_anonym_directed_v3.txt 0.1 > 3.final

#Core
(python coremoduleidentification.py 3.final ../data/subchallenge1/3_signal_anonym_directed_v3.txt)>3_signal_anonym_directed_v3.txt
