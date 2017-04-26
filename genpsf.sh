#!/bin/bash

TOPOLOGY=data/formisha/dph/top_all36_prot_met.rtf #data/toppar/top_all22_prot_num.rtf #data/toppar/top_all22_prot_num.rtf #data/toppar/top_all36_prot.rtf #/home/mikhail/Projects/SAXS/pdb_prep/mol-prms/pdbamino.rtf #data/toppar/top_all22_prot.rtf #data/cilengitit_rosi/charmm_param.rtf #/home/mikhail/Projects/SAXS/pdb_prep/mol-prms/pdbamino.rtf #data/top_all22_prot_modified.rtf #data/cilengitit_rosi/charmm_param.rtf #data/top_all22_prot_modified.rtf
INPUTPDB=data/noe_testing/CilMD.pdb #data/cilengitit_rosi/b2m.pdb #data/cilengitit_rosi/CilMD.pdb #data/testing/CIL_04844.pdb
RESIDUE=5
OUTPBD=${INPUTPDB/%.pdb/_mod.pdb}
OUTPSF=${INPUTPDB/%.pdb/_mod.psf}

SCRIPT="topology ${TOPOLOGY};
#segment A0 { pdb ${INPUTPDB}; first none; last none; };
#patch NMET A0:${RESIDUE};
#coordpdb ${INPUTPDB} A0;
#guesscoord;
#writepdb ${OUTPBD};
#writepsf ${OUTPSF}"

SCRIPT="topology ${TOPOLOGY};
segment A0 { pdb ${INPUTPDB}; first none; last none; auto angles dihedrals };
coordpdb ${INPUTPDB} A0;
guesscoord;
writepdb ${OUTPBD};
writepsf charmm ${OUTPSF}"

SCRIPT="topology ${TOPOLOGY};
segment \"\" {pdb ${INPUTPDB}; first none; last none; auto angles dihedrals};
patch link :${RESIDUE} :1;
patch nmet :${RESIDUE};
regenerate angles dihedrals;
coordpdb ${INPUTPDB} ;
guesscoord;
regenerate angles dihedrals;
writepdb ${OUTPBD};
writepsf ${OUTPSF}"

which vmd || ! echo "Install VMD" || exit 1

vmd <<< ${SCRIPT}
