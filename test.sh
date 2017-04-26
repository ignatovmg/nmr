make

rm -f data/noe_testing/some_cmin_xplor.spect
xplor src/spectrum.inp

cd data/noe_testing
python convert_spect.py some_cmin_xplor.pdb some_cmin_xplor.spect some_cmin_xplor_volume
cd -

./run.sh

