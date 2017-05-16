#include <time.h>

#include "noe.h"

#include "mol2/benergy.h"
#include "mol2/gbsa.h"
#include "mol2/icharmm.h"
#include "mol2/minimize.h"
#include "mol2/nbenergy.h"


int main(int argc, char** argv)
{
	char* pdb = argv[1];
	char* psf = argv[2];
	char* prm = argv[3];
	char* rtf = argv[4];
	char* grp = argv[5];
	char* exp = argv[6];
	
	if (argc != 7)
	{
		fprintf(stderr, "argc wrong\n");
		exit(EXIT_FAILURE);
	}
	
	double omega = 0.00006; // e+9 Gz
	double t_cor = 0.1; // e-9 s
	double t_mix = 0.2; // s

	
	struct mol_atom_group *ag = mol_read_pdb(pdb);
	mol_atom_group_read_geometry(ag, psf, prm, rtf);
	
	// Read proton groups
	struct noe* spect = init_noe(grp);
	groups* grps = spect->grps;
	int size = spect->N;
	
	spect->rx_grad = init_gradient(grps->natoms, size);
	spect->in_grad = init_gradient(grps->natoms, size);
	spect->exp = read_exp(exp, size);
	spect->mask = get_mask(spect->exp, size);
	spect->omega = omega;
	spect->t_cor = t_cor;
	spect->t_mix = t_mix;


	// Intencities
	peaks(spect, ag);
	
	//my_pprint_matrix("sandbox/computed_matrix_raw.csv", spect->in, size);
	//my_pprint_matrix("sandbox/experimental_matrix.csv", spect->exp, size);	
	
	/*double mult = best_multiplier(spect->exp, spect->in, size);
	for (int i = 0; i < size * size; i++)
	{
		spect->in[i] *= mult;	
	}

	my_pprint_matrix("sandbox/computed_matrix.csv", spect->in, size);
	printf("Score: %.6f\n", fit_score(NULL, spect, 0.0));*/
	
	int natoms = grps->natoms;
	struct mol_vector3* fit_grad = calloc(ag->natoms, sizeof(struct mol_vector3));
	for (int i = 0; i < grps->natoms; i++)
	{
		fit_grad[grps->atoms[i]].X = 0.0;
		fit_grad[grps->atoms[i]].Y = 0.0;
		fit_grad[grps->atoms[i]].Z = 0.0;
	}
	
	clock_t start = clock();

	peaks(spect, ag);
	double fit = fit_score(fit_grad, spect, -1.0);
	
	printf("Analytical time: %f\n", (float)(clock() - start) / CLOCKS_PER_SEC);
	



	struct mol_vector3* fit_grad_num = calloc(ag->natoms, sizeof(struct mol_vector3));
	for (int i = 0; i < grps->natoms; i++)
	{
		fit_grad_num[grps->atoms[i]].X = 0.0;
		fit_grad_num[grps->atoms[i]].Y = 0.0;
		fit_grad_num[grps->atoms[i]].Z = 0.0;
	}
	start = clock();
	
	grad_numeric_in(spect, ag);
	fit = fit_score(fit_grad_num, spect, -1.0);
	
	printf("Numerical time: %f\n", (float)(clock() - start) / CLOCKS_PER_SEC);
	
	for (int i = 0; i < spect->grps->natoms; i++)
	{
		printf("%.6f %.6f\n", fit_grad[grps->atoms[i]].Z, fit_grad_num[grps->atoms[i]].Z);
	}
	
	
	/*double* fit_grad_num = calloc(3*natoms, sizeof(double));
	struct gradient* old_rx = spect->rx_grad;
	struct gradient* old_in = spect->in_grad;
	spect->rx_grad = NULL;
	spect->in_grad = NULL;

	double step = 0.0001;
	double pert_fit;
	
	for (int atm = 0; atm < natoms; atm++)
	{
		ag->coords[grps->atoms[atm]].X += step;
		peaks(spect, ag);
		ag->coords[grps->atoms[atm]].X -= step;
		pert_fit = fit_score(NULL, spect, 0);
		fit_grad_num[3*atm+0] = (pert_fit - fit) / step;
		
		ag->coords[grps->atoms[atm]].Y += step;
		peaks(spect, ag);
		ag->coords[grps->atoms[atm]].Y -= step;
		pert_fit = fit_score(NULL, spect, 0);
		//printf("%f\n", (pert_fit - fit)/step);
		fit_grad_num[3*atm+1] = (pert_fit - fit) / step;
		
		ag->coords[grps->atoms[atm]].Z += step;
		peaks(spect, ag);
		ag->coords[grps->atoms[atm]].Z -= step;
		pert_fit = fit_score(NULL, spect, 0);
		fit_grad_num[3*atm+2] = (pert_fit - fit) / step;
	}

	spect->rx_grad = old_rx;
	spect->in_grad = old_in;
	for (int i = 0; i < spect->grps->natoms; i++)
	{
		printf("%.6f %.6f\n", fit_grad[grps->atoms[i]].Z, fit_grad_num[3*i+2]);
	}*/
	
	/*peaks(spect, ag); peaks(spect, ag);
	struct gradient* an_grad = spect->in_grad; //init_gradient(grps->natoms, size);
	spect->in_grad = init_gradient(grps->natoms, size);
	
	start = clock();
	
	grad_numeric_in(spect, ag);
	printf("Numeric time: %f\n", (float)(clock() - start) / CLOCKS_PER_SEC);
	
	double dev = 0.0;
	int ctr = 0;
	for (int i = 0; i < an_grad->natoms*size*size; i++)
	{
		//printf("%.6f %.6f\n", an_grad->gradx[i], spect->in_grad->gradx[i]);
		
		if (an_grad->gradx[i] > 0.0000000001)
		{
			printf("%6i %.6f %.6f %.6f\n", i/(size*size), an_grad->gradx[i], spect->in_grad->gradx[i], an_grad->gradx[i]/spect->in_grad->gradx[i]);
			dev += pow(an_grad->gradx[i] - spect->in_grad->gradx[i], 2) / an_grad->gradx[i];
			ctr++;
		}
	}
	printf("\n");
	
	printf("num vs ana dev = %.6f\n", sqrt(dev / ctr));
	*/
	
	/*rx_mat(spect, ag);
	double* grad_rx = grad_numeric_rx(spect, ag);
	for (int i = 0; i < grps->natoms*size*size; i++)
	{
		double num = grad_rx[0*spect->grps->natoms*size*size + i];
		double anal = spect->rx_grad->gradx[i];
		if (num > 0.0000000001)
		printf("%4i %3i %3i %.6f %.6f %.6f\n", i/(size*size), (i%(size*size))/size, (i%(size*size))%size, num, anal, anal/num);
	}*/

	free_mol_atom_group(ag);
	free_noe(spect);

	return 0;
}

