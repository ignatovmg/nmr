#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "mol2/benergy.h"
#include "mol2/gbsa.h"
#include "mol2/icharmm.h"
#include "mol2/minimize.h"
#include "mol2/nbenergy.h"
#include "mol2/pdb.h"

#define __TOL__  5E-4
#define __DMIN__ 1E-6


typedef struct group_ // TODO: optimize size
{
	int N;
	int id;
	int group[3];
}group;

typedef struct groups_
{
	int N;
	int natoms;
	
	group* groups;
}groups;

struct gradient
{
	int N;
	int natoms;
	int size;
	double* gradx;
	double* grady;
	double* gradz;
};

groups* init_proton_groups(int N);
void free_proton_groups(groups* grps);
groups* read_proton_groups(char* path);
void print_proton_groups(groups* grps);

struct gradient* init_gradient(int natoms, int size);
void free_gradient();

int get_line_num(char* path);
void my_print_matrix(double* m, int size);
void my_fprint_matrix(FILE* f, double* m, int size);
void my_fprint_matrix_stacked(FILE* f, double* m, int size);
void my_pprint_matrix(char* path, double* m, int size);

//                                                             e+9 Gz        e-9 s  if grad = NULL then not computed
double* rx_mat(struct mol_atom_group *ag, groups* grps, double omega, double t_cor, struct gradient* grad);
gsl_matrix* matrix_exp(double* rx, int size, double t_mix, struct gradient* rx_grad, struct gradient* in_grad);
gsl_matrix* peaks(struct mol_atom_group *ag, groups* grps, double omega, double t_cor, double t_mix);
double fit_score(char* path, double* mat, int size);

#ifdef TEST
double best_multiplier(char* path, double* mat, int size);

double* grad_numeric(struct mol_atom_group *ag, groups* grps, double omega, double t_cor, double t_mix);
double* grad_numeric_rx(struct mol_atom_group *ag, groups* grps, double omega, double t_cor);
#endif


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
	
	double omega = 0.6; // e+9 Gz
	double t_cor = 0.1; // e-9 s
	double t_mix = 1.0; // s
	
	struct mol_atom_group *ag = mol_read_pdb(pdb);
	mol_atom_group_read_geometry(ag, psf, prm, rtf);
	
	// Read proton groups
	groups* grps = read_proton_groups(grp);
	int size = grps->N;

	// Relaxation matrix
	struct gradient* rx_grad = init_gradient(grps->natoms, size);
	double* rx = rx_mat(ag, grps, omega, t_cor, rx_grad);
	
#ifdef TEST
	my_pprint_matrix("sandbox/relaxation_matrix.csv", rx, size);
#endif	

	// Intencities
	struct gradient* in_grad = init_gradient(grps->natoms, size);
	gsl_matrix* in = matrix_exp(rx, size, t_mix, rx_grad, in_grad);
		
	// Normalization
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			//in->data[i*size+j] *= sqrt(grps->groups[i].N * grps->groups[j].N);
		}
	}
	
	for (int k = 0; k < grps->natoms; k++)
	{
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				//in_grad->gradx[k*size*size+i*size+j] *= sqrt(grps->groups[i].N * grps->groups[j].N);
				//in_grad->grady[k*size*size+i*size+j] *= sqrt(grps->groups[i].N * grps->groups[j].N);
				//in_grad->gradz[k*size*size+i*size+j] *= sqrt(grps->groups[i].N * grps->groups[j].N);
			}
		}
	}

	//print_proton_groups(grps);
	
#ifdef TEST
	my_pprint_matrix("sandbox/computed_matrix_raw.csv", in->data, size);
#endif	
	
#ifdef TEST
	/*double mult = best_multiplier(exp, in->data, size);
	for (int i = 0; i < size * size; i++)
	{
		in->data[i] *= mult;	
	}*/

	my_pprint_matrix("sandbox/computed_matrix.csv", in->data, size);
#endif
	
	//double fit = fit_score(exp, in->data, size);
	
#ifdef TEST
	//printf("fit: %.5f\n", fit);
#endif


#ifdef TEST	
	double* grad = grad_numeric(ag, grps, omega, t_cor, t_mix);
	//double* grad = grad_numeric_rx(ag, grps, omega, t_cor);
	for (int i = 0; i < 100; i++)
	{
		printf("%.6f %.6f\n", grad[3*size*size*(i / (size * size)) + (i % (size * size))], in_grad->gradx[i]);
	}
	
	/*printf("\n\n");
	for (int i = 0; i < 30; i++)
	{
		printf("%.6f ", grad[3*size*size + i]);
	}
	
	printf("\n");
	for (int i = 0; i < 30; i++)
	{
		printf("%.6f ", rx_grad->gradx[size*size + i]);
	}*/
#endif

	free_mol_atom_group(ag);
	free_proton_groups(grps);
	free(rx);

	return 0;
}



//                                                             e+9 Gz        e-9 s 
double* rx_mat(struct mol_atom_group *ag, groups* grps, double omega, double t_cor, struct gradient* grad)
{
	int size = grps->N;
	int loc, loc1;

	double* gradx;
	double* grady;
	double* gradz;
	
	if (grad != NULL)
	{
		gradx = grad->gradx;
		grady = grad->grady;
		gradz = grad->gradz;
		
		memset(gradx, 0, grad->N);
		memset(grady, 0, grad->N);
		memset(gradz, 0, grad->N);
	}
	
	double J_0 = t_cor;
	double J_1 = t_cor / (1.0 +  4.0 * ( M_PI * M_PI ) * omega * omega * t_cor * t_cor );  // *10^(-9) s
	double J_2 = t_cor / (1.0 + 16.0 * ( M_PI * M_PI ) * omega * omega * t_cor * t_cor );  // *10^(-9) s
	
	double gyr   = 2.6751965; // e+8
	double hbar  = 1.0545919; // e-34;
	double coef = pow(gyr, 4)*pow(hbar, 2);
	
	double  S0 =       J_0 * coef; // 1/s
	double  S1 = 1.5 * J_1 * coef; // 1/s
	double  S2 = 6.0 * J_2 * coef; // 1/s
	
	double* rx = calloc(size * size, sizeof(double));
	
	group *grpi, *grpj;	
	struct mol_vector3 *coordi, *coordj;
	double r2;
	double mult, mult1;
	
	int atomi = 0, atomj;
	float counter;
	for (int i = 0; i < size; i++)
	{
		grpi = &grps->groups[i];
		atomj = atomi;
		for (int j = i; j < size; j++)
		{
			grpj = &grps->groups[j];
			
			counter = 0.0;
			rx[size*i+j] = 0.0;
			
			if (i == j)
			{
				// Self relaxation for unresolved groups if i == j
				if (grpi->N > 1)
				{
					for (int i1 = 0; i1 < grpi->N; i1++)
					{
						for (int j1 = i1 + 1; j1 < grpj->N; j1++)
						{
							coordi = &ag->coords[grpi->group[i1]];
							coordj = &ag->coords[grpj->group[j1]];
					
							r2 = (coordi->X - coordj->X) * (coordi->X - coordj->X) + \
								 (coordi->Y - coordj->Y) * (coordi->Y - coordj->Y) + \
								 (coordi->Z - coordj->Z) * (coordi->Z - coordj->Z);
								 
							double r8 = r2 * r2 * r2 * r2;
							
							if (r2 > __DMIN__)
							{
								rx[size*i+j] += 1.0 / (r2 * r2 * r2);
								counter += 1.0;
								
								// accumulate gradient
								if (grad != NULL)
								{
									loc = (atomi + i1) * size * size + (i * size + j);
									gradx[loc] -= 6 * (coordi->X - coordj->X) / r8;
									grady[loc] -= 6 * (coordi->Y - coordj->Y) / r8;
									gradz[loc] -= 6 * (coordi->Z - coordj->Z) / r8;
									
									loc = (atomj + j1) * size * size + (i * size + j);
									gradx[loc] -= 6 * (coordj->X - coordi->X) / r8;
									grady[loc] -= 6 * (coordj->Y - coordi->Y) / r8;
									gradz[loc] -= 6 * (coordj->Z - coordi->Z) / r8;
								}
							}
						}
					}
					
					rx[size*i+j] /= counter;
					rx[size*i+j] *= 2 * (grpi->N - 1) * (S1 + S2);
					
					// times self relaxation
					if (grad != NULL)
					{
						for (int i1 = 0; i1 < grpi->N; i1++)
						{	
							loc = (atomi + i1) * size * size + (i * size + j);
							gradx[loc] *= 2 * (grpi->N - 1) * (S1 + S2) / counter;
							grady[loc] *= 2 * (grpi->N - 1) * (S1 + S2) / counter;
							gradz[loc] *= 2 * (grpi->N - 1) * (S1 + S2) / counter;
							
						}
					}
				}
			}
			else
			{
				// 6 average if i != j
				for (int i1 = 0; i1 < grpi->N; i1++)
				{
					for (int j1 = 0; j1 < grpj->N; j1++)
					{
						coordi = &ag->coords[grpi->group[i1]];
						coordj = &ag->coords[grpj->group[j1]];
					
						r2 = (coordi->X - coordj->X) * (coordi->X - coordj->X) + \
							 (coordi->Y - coordj->Y) * (coordi->Y - coordj->Y) + \
							 (coordi->Z - coordj->Z) * (coordi->Z - coordj->Z);
							
						double r8 = r2 * r2 * r2 * r2;	
							
						if (r2 > __DMIN__)
						{
							rx[size*i+j] += 1.0 / (r2 * r2 * r2);
							counter += 1.0;
							
							// accumulate gradient
							if (grad != NULL)
							{
								loc = (atomi + i1) * size * size + (i * size + j);
								gradx[loc] -= 6 * (coordi->X - coordj->X) / r8;
								grady[loc] -= 6 * (coordi->Y - coordj->Y) / r8;
								gradz[loc] -= 6 * (coordi->Z - coordj->Z) / r8;
								
								loc = (atomj + j1) * size * size + (i * size + j);
								gradx[loc] -= 6 * (coordj->X - coordi->X) / r8;
								grady[loc] -= 6 * (coordj->Y - coordi->Y) / r8;
								gradz[loc] -= 6 * (coordj->Z - coordi->Z) / r8;
							}
						}
					}
				}
				
				rx[size*i+j] /= counter;
				rx[size*j+i]  = rx[size*i+j];
				
				if (grad != NULL)
				{
					for (int i1 = 0; i1 < grpi->N; i1++)
					{
						loc = (atomi + i1) * size * size + (i * size + j);
						gradx[loc] /= counter;
						grady[loc] /= counter;
						gradz[loc] /= counter;
					}
					
					for (int j1 = 0; j1 < grpj->N; j1++)
					{
						loc = (atomj + j1) * size * size + (i * size + j);
						gradx[loc] /= counter;
						grady[loc] /= counter;
						gradz[loc] /= counter;
					}
				}
			}
			
			atomj += grpj->N;
		}
		
		atomi += grpi->N;
	}
	
	
	// Fill diagonal elements in R
	mult = S0 + 2 * S1 + S2;
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (i != j)
			{
				rx[i*size+i] += mult * grps->groups[j].N * rx[i*size+j];
			}
		}
	}
	
	// Gradient
	if (grad != NULL)
	{
		for (int k = 0; k < grps->natoms; k++)
		{
			loc = k * size * size;
			for (int i = 0; i < size; i++)
			{
				for (int j = 0; j < size; j++)
				{
					if (i != j)
					{
						mult1 = mult * grps->groups[j].N;
						
						if (i < j)
						{
							loc1 = (i * size + j);
						}
						else
						{
							loc1 = (j * size + i);
						}

						gradx[loc + (i * size + i)] += mult1 * gradx[loc + loc1];
						grady[loc + (i * size + i)] += mult1 * grady[loc + loc1];
						gradz[loc + (i * size + i)] += mult1 * gradz[loc + loc1];
					}
				}
			}
		}
	}
	
	// Fill the rest of R
	mult = S2 - S0;
	for (int i = 0; i < size - 1; i++)
	{
		for (int j = i + 1; j < size; j++)
		{
			rx[i*size+j] = mult * sqrt(grps->groups[i].N * grps->groups[j].N) * rx[i*size+j];
			rx[j*size+i] = rx[i*size+j];
		}
	}
	
	// Gradient
	if (grad != NULL)
	{
		for (int k = 0; k < grps->natoms; k++)
		{
			loc  = k * size * size;
			for (int i = 0; i < size - 1; i++)
			{
				for (int j = i + 1; j < size; j++)
				{
					mult1 = mult * sqrt(grps->groups[i].N * grps->groups[j].N);
					
					loc1 = i * size + j;
					int loc2 = j * size + i;
					gradx[loc + loc1] *= mult1;
					grady[loc + loc1] *= mult1;
					gradz[loc + loc1] *= mult1;
					
					gradx[loc + loc2] = gradx[loc + loc1];
					grady[loc + loc2] = grady[loc + loc1];
					gradz[loc + loc2] = gradz[loc + loc1];
				}
			}
		}
	}
	
	return rx;
}

gsl_matrix* matrix_exp(double* rx, int size, double t_mix, struct gradient* rx_grad, struct gradient* in_grad)
{
	//for (int i = 0; i < size*size; i++)
	//{
	//		rx[i] *= -t_mix; // TODO: change this
	//}
			
	gsl_matrix_view m = gsl_matrix_view_array(rx, size, size);	
			
	//print_matrix(&m.matrix);

	gsl_vector *eval = gsl_vector_alloc(size);
	gsl_matrix *evec = gsl_matrix_alloc(size, size);
	gsl_matrix *res  = gsl_matrix_alloc(size, size);
	gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(size);
	gsl_eigen_symmv(&m.matrix, eval, evec, w);
	
	//print_matrix(evec);
	double* expeval = calloc(size, sizeof(double));
	
	for (int i = 0; i < size; i++)
	{
		expeval[i] = exp(-t_mix * eval->data[i]);
		eval->data[i] *= -t_mix;
	}
	
	for (int i = 0; i < size; i++)
	{
		for (int j = i; j < size; j++)
		{
			int id = i*size+j;
			res->data[id] = 0.0;
			
			for (int k = 0; k < size; k++)
			{
				res->data[id] += expeval[k]*evec->data[i*size+k]*evec->data[j*size+k];
			}
				
			res->data[j*size+i] = res->data[i*size+j];
		}
	}

	if (rx_grad != NULL && in_grad != NULL)
	{
		double* helper = calloc(size * size, sizeof(double));
		double ab, ijx, ijy, ijz, coef;
		int loc;
		
		for (int i = 0; i < size; i++)
		{
			for (int j = i; j < size; j++)
			{
				for (int a = 0; a < size; a++)
				{
					for (int b = 0; b < size; b++)
					{

						ab = 0.0;
						for (int k = 0; k < size; k++)
						{
							for (int t = 0; t < size; t++)
							{
								if (k != t)
								{
									coef = -(expeval[k] - expeval[t]) / \
									       (eval->data[k] - eval->data[t]);
								}
								else
								{
									coef = - expeval[k];
								}

								ab += evec->data[i*size+k] * \
								      evec->data[a*size+k] * \
				 				      evec->data[b*size+t] * \
								      evec->data[j*size+t] * \
								      coef;

							}
						}
						
						helper[a*size+b] = ab;
					}
				}
				
				for (int atm = 0; atm < in_grad->natoms; atm++)
				{
					ijx = 0.0;
					ijy = 0.0;
					ijz = 0.0;
					
					// Trace
					for (int a = 0; a < size; a++)
					{
						for (int b = 0; b < size; b++)
						{
							ijx += rx_grad->gradx[atm*size*size + a*size + b] * helper[a*size+b];
							ijy += rx_grad->grady[atm*size*size + a*size + b] * helper[a*size+b];
							ijz += rx_grad->gradz[atm*size*size + a*size + b] * helper[a*size+b];
						}
					}
					
					loc = atm * size * size + i * size + j;
					in_grad->gradx[loc] = ijx;
					in_grad->grady[loc] = ijy;
					in_grad->gradz[loc] = ijz;
					
					loc = atm * size * size + j * size + i;
					in_grad->gradx[loc] = ijx;
					in_grad->grady[loc] = ijy;
					in_grad->gradz[loc] = ijz;
				}
			}
		}

	}		

	gsl_eigen_symmv_free(w);
	gsl_vector_free(eval);
	gsl_matrix_free(evec);
	free(expeval);
	
    return res;
}

gsl_matrix* peaks(struct mol_atom_group *ag, groups* grps, double omega, double t_cor, double t_mix)
{
	int size = grps->N;

	// Relaxation matrix
	double* rx = rx_mat(ag, grps, omega, t_cor, NULL);
	
	// Intencities
	gsl_matrix* in = matrix_exp(rx, size, t_mix, NULL, NULL);
	free(rx);
		
	// Normalization
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			//in->data[i*size+j] *= sqrt(grps->groups[i].N * grps->groups[j].N);
		}
	}
	
	return in;
}

double* grad_numeric(struct mol_atom_group *ag, groups* grps, double omega, double t_cor, double t_mix)
{
	gsl_matrix* pks1 = peaks(ag, grps, omega, t_cor, t_mix);
	gsl_matrix* pks2;
	
	int size = grps->N;
	double  step = 0.000001;
	double* grad = calloc(3 * grps->natoms * size * size, sizeof(double));
	
	int counter = 0;
	
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < grps->groups[i].N; j++)
		{
			int atom = grps->groups[i].group[j];
			
			// fill grad for X
			ag->coords[atom].X += step;
			pks2 = peaks(ag, grps, omega, t_cor, t_mix);
			ag->coords[atom].X -= step;
			
			for (int k = 0; k < size * size; k++)
			{
				grad[(3 * counter + 0) * size * size + k] = (pks2->data[k] - pks1->data[k]) / step;
			}
			
			gsl_matrix_free(pks2);
			
			// fill grad for Y
			ag->coords[atom].Y += step;
			pks2 = peaks(ag, grps, omega, t_cor, t_mix);
			ag->coords[atom].Y -= step;
			
			for (int k = 0; k < size * size; k++)
			{
				grad[(3 * counter + 1) * size * size + k] = (pks2->data[k] - pks1->data[k]) / step;
			}
			
			gsl_matrix_free(pks2);
			
			// fill grad for Z
			ag->coords[atom].Z += step;
			pks2 = peaks(ag, grps, omega, t_cor, t_mix);
			ag->coords[atom].Z -= step;
			
			for (int k = 0; k < size * size; k++)
			{
				grad[(3 * counter + 2) * size * size + k] = (pks2->data[k] - pks1->data[k]) / step;
			}
			
			gsl_matrix_free(pks2);
			
			counter++;
		}
	}
	
	gsl_matrix_free(pks1);
	
	return grad;
}

double* grad_numeric_rx(struct mol_atom_group *ag, groups* grps, double omega, double t_cor)
{
	double* rx1 = rx_mat(ag, grps, omega, t_cor, NULL);
	double* rx2;
	
	int size = grps->N;
	double  step = 0.0000001;
	double* grad = calloc(3 * grps->natoms * size * size, sizeof(double));
	
	int counter = 0;
	
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < grps->groups[i].N; j++)
		{
			int atom = grps->groups[i].group[j];
			
			// fill grad for X
			ag->coords[atom].X += step;
			rx2 = rx_mat(ag, grps, omega, t_cor, NULL);
			ag->coords[atom].X -= step;
			
			for (int k = 0; k < size * size; k++)
			{
				grad[(3 * counter + 0) * size * size + k] = (rx2[k] - rx1[k]) / step;
			}
			
			free(rx2);
			
			// fill grad for Y
			ag->coords[atom].Y += step;
			rx2 = rx_mat(ag, grps, omega, t_cor, NULL);
			ag->coords[atom].Y -= step;
			
			for (int k = 0; k < size * size; k++)
			{
				grad[(3 * counter + 1) * size * size + k] = (rx2[k] - rx1[k]) / step;
			}
			
			free(rx2);
			
			// fill grad for Z
			ag->coords[atom].Z += step;
			rx2 = rx_mat(ag, grps, omega, t_cor, NULL);
			ag->coords[atom].Z -= step;
			
			for (int k = 0; k < size * size; k++)
			{
				grad[(3 * counter + 2) * size * size + k] = (rx2[k] - rx1[k]) / step;
			}
			
			free(rx2);
			
			counter++;
		}
	}
	
	free(rx1);
	
	return grad;
}


#ifdef TEST
// For tests only
double best_multiplier(char* path, double* mat, int size)
{
	FILE* f = fopen(path, "r");
	if (f == NULL)
	{
		printf("%s not found\n", path);
		exit(EXIT_FAILURE);
	}
	
	int i, j;
	double val;
	double num = 0.0;
	double den = 0.0;
	while(fscanf(f, "%i\t%i\t%lf\n", &i, &j, &val) != EOF)
	{
		num += val*mat[i*size+j];
		//num += val;
		den += mat[i*size+j]*mat[i*size+j];
		//den += m[i*size+j];
	}
	rewind(f);
	double k = num / den;
	
	return k;
}
#endif

double fit_score(char* path, double* mat, int size)
{
	float power = 1.0 / 6.0;

	FILE* f = fopen(path, "r");
	if (f == NULL)
	{
		fprintf(stderr, "%s not found\n", path);
		exit(EXIT_FAILURE);
	}
	
	int i, j;
	double val;
	double num = 0.0;
	double den = 0.0;	
	
	int nchar;
	while((nchar = fscanf(f, "%i\t%i\t%lf\n", &i, &j, &val)) != EOF)
	{
		if (nchar != 3)
		{
			fprintf(stderr, "%s: Wrong input\n", __func__);
			exit(EXIT_FAILURE);
		}
		
		num += val*mat[i*size+j];
		//num += val;
		den += mat[i*size+j]*mat[i*size+j];
		//yden += mat[i*size+j];
	}
	rewind(f);
	
	double k = num / den;
	
#ifdef TEST
	double* reference = (double*)calloc(size*size, sizeof(double));
	
	while(fscanf(f, "%i\t%i\t%lf\n", &i, &j, &val) != EOF)
	{
		reference[i*size+j] = val;
		reference[j*size+i] = reference[i*size+j];
	}
	
	rewind(f);

	my_pprint_matrix("sandbox/experimental_matrix.csv", reference, size);
	free(reference);
#endif
	
	double trm;
	double fit = 0.0;
	double nrm = 0.0;
	int n = 0;
	
	while((nchar = fscanf(f, "%i\t%i\t%lf\n", &i, &j, &val)) != EOF)
	{
		if (nchar != 3)
		{
			fprintf(stderr, "%s: Wrong input\n", __func__);
			exit(EXIT_FAILURE);
		}
		
		trm = 0.0;
		if (fabs(val) > 0.0)
		{
			trm  = (val / fabs(val)) * pow(fabs(val), power);
		}
			
		if (fabs(k * mat[i*size+j]) > 0.0)
		{
			trm -= (k * mat[i*size+j] / fabs(k * mat[i*size+j])) * pow(fabs(k * mat[i*size+j]), power);
		}

		fit += trm * trm;
		nrm += pow(fabs(val), power);
		n++;
	}
	
	fclose(f);
	return fit / nrm;
}

groups* init_proton_groups(int N)
{
	groups* grps = calloc(1, sizeof(groups));

	grps->N = N;
	grps->groups = calloc(N, sizeof(group));
	//printf("N: %i\n", groups->N);
	
	return grps;
}

void free_proton_groups(groups* grps)
{
	free(grps->groups);
	free(grps);
}

groups* read_proton_groups(char* path)
{	
	groups* grps = init_proton_groups(get_line_num(path));

	FILE* grp_file = fopen(path, "r");
	if (grp_file == NULL)
	{
		fprintf(stderr, "%s not found\n", path);
		exit(EXIT_FAILURE);
	}
	
	char str[100];
	char *chunk;
	
	int counter = 0;
	int natoms  = 0;
	while(fgets(str, 100, grp_file) != NULL)
	{
		int grp_size = 0;
		
  		chunk = (char*)strtok(str, "\t");
  		grps->groups[counter].id = atoi(chunk);
  		
  		assert(grps->groups[counter].id == counter);
  		
  		chunk = (char*)strtok(NULL, "\t");
		while (chunk != NULL)
		{
			if (grp_size == 3)
			{
				fprintf(stderr, "Group cant me larger than 3\n");
				exit(EXIT_FAILURE);
			}
			
			grps->groups[counter].group[grp_size] = atoi(chunk) - 1;
			chunk = strtok(NULL, "\t");
			
			grp_size++;
			natoms++;
		}
		
		grps->groups[counter].N = grp_size;
		counter++;
	}
	
	grps->natoms = natoms;
	fclose(grp_file);

	return grps;
}

void print_proton_groups(groups* grps)
{
	for (int i = 0; i < grps->N; i++)
	{
		printf("group %i of size %i: ", grps->groups[i].id, grps->groups[i].N);
		for (int j = 0; j < grps->groups[i].N; j++)
			printf("%i\t", grps->groups[i].group[j]);
		printf("\n");
	}
}
                   
struct gradient* init_gradient(int natoms, int size)
{
	struct gradient* grad = malloc(sizeof(struct gradient));
	grad->N = natoms * size * size;
	grad->natoms = natoms;
	grad->size = size;
	
	grad->gradx = calloc(natoms * size * size, sizeof(double));
	grad->grady = calloc(natoms * size * size, sizeof(double));
	grad->gradz = calloc(natoms * size * size, sizeof(double));
	
	return grad;
}

void free_gradient(struct gradient* grad)
{
	free(grad->gradx);
	free(grad->grady);
	free(grad->gradz);
	
	free(grad);
}


int get_line_num(char* path)
{
	FILE* f = fopen(path, "r");
	if (f == NULL)
	{
		fprintf(stderr, "File not found\n");
		exit(EXIT_FAILURE);
	}
	char str[100];
	int i = 0;
	while(fgets(str, 100, f) != NULL) i++;
	
	fclose(f);

	return i;
}

void my_print_matrix(double* m, int size)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			printf("%.4f ", m[i*size+j]);
		}
		printf("\n");
	}
}

void my_fprint_matrix(FILE* f, double* m, int size)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size-1; j++)
			fprintf(f, "%.6f\t", m[i*size+j]);
		fprintf(f, "%.6f\n", m[i*size+size-1]);
	}
}

void my_fprint_matrix_stacked(FILE* f, double* m, int size)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			fprintf(f, "%i\t%i\t%.2f\n", i, j, m[i*size+j]);
		}
	}
}

void my_pprint_matrix(char* path, double* m, int size)
{
	FILE* file = fopen(path, "w");
	my_fprint_matrix(file, m, size);
	fclose(file);
}

