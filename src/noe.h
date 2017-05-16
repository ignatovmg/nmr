#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "mol2/pdb.h"

#define __DMIN__ 1E-8


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
	int* atoms;
	
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

struct noe
{
	int N;
	double omega;
	double t_cor;
	double t_mix;

	double* in;
	double* rx;
	double* exp;
	struct gradient* in_grad;
	struct gradient* rx_grad;
	int* mask;
	
	groups* grps;
};

groups* init_proton_groups(int N);
void free_proton_groups(groups* grps);
groups* read_proton_groups(char* path);
void print_proton_groups(groups* grps);

double* read_exp(char* exp, int size);

struct gradient* init_gradient(int natoms, int size);
void free_gradient();

struct noe* init_noe(char* grp);
void free_noe(struct noe* spect);

int get_line_num(char* path);
void my_print_matrix(double* m, int size);
void my_fprint_matrix(FILE* f, double* m, int size);
void my_fprint_matrix_stacked(FILE* f, double* m, int size);
void my_pprint_matrix(char* path, double* m, int size);

double* rx_mat(struct noe* spect, struct mol_atom_group *ag);
double* matrix_exp(struct noe* spect);
double* peaks(struct noe* spect, struct mol_atom_group *ag);

int* get_mask(double* exp, int size);
double best_multiplier(double* exp, double* spc, int size);
double fit_score(struct mol_vector3* grad, struct noe* spect, double gcoef);

void    grad_numeric_in(struct noe* spect, struct mol_atom_group *ag);
double* grad_numeric_rx(struct noe* spect, struct mol_atom_group *ag);








