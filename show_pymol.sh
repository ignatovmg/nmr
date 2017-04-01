#!/bin/bash

pymol $(sort -n -k2 sandbox/RUN2res/chi_scores | cut -f1 | head -n500 | sed "s/\.dat/_aln\.pdb/g" | sed "s/\.\/sandbox\/refined/data\/RUN2_aligned/g") data/cilengitit_rosi/CilMD.pdb &
pymol $(awk '{$2=$2*1000; print}' sandbox/distance_filtered | sort -n -k2 | cut -f1 -d' ' | head -n500 | sed "s/RUN2/RUN2\_aligned/g" | sed "s/\.pdb/_aln\.pdb/g") data/cilengitit_rosi/CilMD.pdb

