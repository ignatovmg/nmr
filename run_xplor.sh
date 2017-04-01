#!/bin/bash

#PDB="test\/CilMD_autopsf.pdb"
#PSF="test\/CilMD_autopsf.psf"
#OUT="test\/CilMD.spect"
#MIX="0.2"
#COR="1.0e-10"

#sed -e "/structure @/ s/@\S* end/@${PDB} end/" -e "/coor @/ s/@\S* /@${PSF} /" -e "/isotropic/ s/isotropic \S*/isotropic ${COR}/" -e "/taumix/ s/taumix D300 \S*/taumix D300 ${MIX}/" -e "/set print/ s/set print \S*/set print ${OUT}/"  src/spectrum.inp > tmp
#mv tmp src/spectrum.inp


xplor src/spectrum.inp


