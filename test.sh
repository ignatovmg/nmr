make
rm -f data/testing/CIL_04844_mod.spect
xplor src/spectrum.inp
cd data/testing
python convert_spect.py CIL_04844_mod.pdb CIL_04844_mod.spect CIL_04844_mod_volume
cd -
./run.sh
