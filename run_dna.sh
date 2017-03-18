
pdbfile='/home/mikhail/Projects/nmr/data/Cyclo-dA/useful/c_dA_final_dnm.pdb'
subject='/home/mikhail/Projects/nmr/data/Cyclo-dA/useful/C_dA_80_pks.txt'
outpath='/home/mikhail/Projects/nmr/data/Cyclo-dA/useful/out'

python src/parse.py $pdbfile $subject $outpath
./run.sh
