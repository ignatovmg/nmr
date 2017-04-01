#!/bin/bash

TOPOLOGY=data/top_all22_prot_modified.rtf
INPUTPDB=data/testing/CIL_04844.pdb
RESIDUE=2
OUTPBD=${INPUTPDB/%.pdb/_mod.pdb}
OUTPSF=${INPUTPDB/%.pdb/_mod.psf}

SCRIPT="topology ${TOPOLOGY};
segment \"\" { pdb ${INPUTPDB}; first none; last none };
patch NMET :${RESIDUE};
coordpdb ${INPUTPDB};
guesscoord;
writepdb ${OUTPBD};
writepsf ${OUTPSF}"

SCRIPT="topology ${TOPOLOGY};
segment \"\" { pdb ${INPUTPDB} };
coordpdb ${INPUTPDB};
guesscoord;
writepdb ${OUTPBD};
writepsf ${OUTPSF}"

which vmd || ! echo "Install VMD" || exit 1

vmd <<< ${SCRIPT}
