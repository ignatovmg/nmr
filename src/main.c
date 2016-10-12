#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <string.h>
#include <assert.h>

typedef struct proton_group_
{
	int N;
	int id;
	int group[3];
}proton_group;

typedef struct proton_groups_
{
	int N;
	proton_group* groups;
}proton_groups;


int get_line_num(char* path)
{
	FILE* f = fopen(path, "r");
	if (f == NULL)
	{
		printf("File not found");
		exit(0);
	}
	char str[100];
	int i = 0;
	while(fgets(str, 100, f) != NULL)
		i++;
	fclose(f);

	return i;
}

int find_proton(int id, int* proton_ids, int size)
{
	int i = 0; 
	while(proton_ids[i] != id && i < size)
		 i++;
		 
	if (i >= size)
	{
		printf("proton %i not found", id);
		exit(0);
	}
	
	return i;
}

void init_proton_groups(proton_groups* groups, int N)
{
	groups->N = N;
	//printf("N: %i\n", groups->N);
	groups->groups = calloc(N, sizeof(proton_group));
}

void read_proton_groups(proton_groups* groups, char* path, int* proton_ids, int size)
{	
	init_proton_groups(groups, get_line_num(path));
	//printf("N: %i\n", groups->N);

	FILE* eq_groups_file = fopen(path, "r");
	if (eq_groups_file == NULL)
	{
		printf("%s not found", path);
		exit(0);
	}
	
	char str[100];
	char *chunk;
	
	int counter = 0;
	while(fgets(str, 100, eq_groups_file) != NULL)
	{
		int group_size = 0;
		
		//printf("groups %s", str);
  		chunk = (char*)strtok(str, "\t");
  		//printf("%s\n", chunk);
  		groups->groups[counter].id = atoi(chunk);
  		
  		assert(groups->groups[counter].id == counter);
  		
  		chunk = (char*)strtok(NULL, "\t");
		while (chunk != NULL)
		{
			//printf("%s\n", chunk);
			groups->groups[counter].group[group_size] = find_proton(atoi(chunk), proton_ids, size);
			chunk = strtok(NULL, "\t");
			group_size++;
		}
		groups->groups[counter].N = group_size;
		counter++;
	}
	fclose(eq_groups_file);
}

void print_proton_groups(proton_groups* groups)
{
	for (int i = 0; i < groups->N; i++)
	{
		printf("group %i of size %i: ", groups->groups[i].id, groups->groups[i].N);
		for (int j = 0; j < groups->groups[i].N; j++)
			printf("%i\t", groups->groups[i].group[j]);
		printf("\n");
	}
}

void collapse_groups(double* dst, gsl_matrix* src, proton_groups* grps)
{
	int src_size = src->size1;
	int dst_size = grps->N;

	for (int i = 0; i < dst_size; i++)
	{
		proton_group* grp_i = &grps->groups[i];

		for (int j = 0; j < dst_size; j++)
		{
			proton_group* grp_j = &grps->groups[j];
			dst[i*dst_size+j] = 0.0;

			for (int i_id = 0; i_id < grp_i->N; i_id++)
				for (int j_id = 0; j_id < grp_j->N; j_id++)
					dst[i*dst_size+j] += src->data[grp_i->group[i_id]*src_size + grp_j->group[j_id]];
		}
	}
} 

void free_proton_groups(proton_groups* groups)
{
	free(groups->groups);
	free(groups);
}

/*
 * singular decomposition of R (relaxation) matrix
 * R = E*S*E', followed by calculating intensity matrix
 * I = exp^(-R*mix_time) using formula 
 * I = E*exp^(-L*mix_time)*E', where L = 
 * = identity_matrix*eigenvalues_vector
**/
gsl_matrix* get_intensity(double* input, int size, double mix_time)
{
	int i, j, k;
	gsl_matrix_view m = gsl_matrix_view_array(input, size, size);
	//print_matrix(&m.matrix);

	gsl_vector *eval = gsl_vector_alloc(size);
	gsl_matrix *evec = gsl_matrix_alloc(size, size);
	gsl_matrix *res = gsl_matrix_alloc(size, size);
	gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(size);
	gsl_eigen_symmv(&m.matrix, eval, evec, w);
	
	//print_matrix(evec);
	
	for (i = 0; i < size; i++)
		eval->data[i] = exp(-mix_time*eval->data[i]);
	
	for (i = 0; i < size; i++)
		for (j = i; j < size; j++)
		{
			int id = i*size+j;
			res->data[id] = 0.0;
			for (k = 0; k < size; k++)
				res->data[id] += eval->data[k]*evec->data[i*size+k]*evec->data[j*size+k];
			res->data[j*size+i] = res->data[i*size+j];
		}
	//print_matrix(res);
	
	gsl_eigen_symmv_free(w);
	gsl_vector_free (eval);
	gsl_matrix_free (evec);
    return res;
}

void my_print_matrix(double* m, int size)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
			printf("%.4f ", m[i*size+j]);
		printf("\n");
	}
}

void my_fprint_matrix(FILE* f, double* m, int size)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size-1; j++)
			fprintf(f, "%.2f\t", m[i*size+j]);
		fprintf(f, "%.2f\n", m[i*size+size-1]);
	}
}

void my_fprint_matrix_stacked(FILE* f, double* m, int size)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
			fprintf(f, "%i\t%i\t%.2f\n", i, j, m[i*size+j]);
	}
}

double chi_score(char* path, double* m, int size)
{
	FILE* f = fopen(path, "r");
	if (f == NULL)
	{
		printf("%s not found", path);
		exit(0);
	}
	
	int i, j;
	double val;
	double num = 0.0;
	double den = 0.0;
	while(fscanf(f, "%i\t%i\t%lf\n", &i, &j, &val) != EOF)
	{
		num += val*m[i*size+j];
		den += m[i*size+j]*m[i*size+j];
	}
	rewind(f);
	double k = num/den;
	
	double chi = 0.0;
	int n = 0;
	while(fscanf(f, "%i\t%i\t%lf\n", &i, &j, &val) != EOF)
	{
		chi += (val-k*m[i*size+j])*(val-k*m[i*size+j]);
		n++;
	}
	
	fclose(f);
	//printf("k = %lf\n", k);
	return sqrt(chi);
}


int main(int argc, char** argv)
{
	int i, j, k;          
/*
 * zero_freq_coef = J(0)*gamma^4*h^2/(4pi^2*10)
 * sing_freq_coef = 3./2.*J(1)*gamma^4*h^2/(4pi^2*10)
 * doub_freq_coef = 6.*J(2*omega)*gamma^4*h^2/(4pi^2*10)
**/
	double zero_freq_coef = 5.67426; // TODO: Decide about the constants
	double doub_freq_coef = 3.0/2.0*5.67426/1.0047;  //
	double sing_freq_coef = 6.0*5.67426/1.0196;  //
	double mixing_time = 0.3;
	char* eq_groups_path = argv[1]; // path to equivalent proton groups
	char* complex_path = argv[2]; // path to peptide
	char* output_path = argv[3];
	
	printf("Computing complex %s\n", complex_path);
	//printf("argc=%i\n", argc);
	
	int size = get_line_num(complex_path);
	double rx_mat[size*size];
	memset(rx_mat, 0, sizeof(double)*size*size);
	
	FILE* f = fopen(complex_path, "r");
	if (f == NULL)
	{
		printf("%s not found", complex_path);
		exit(0);
	}

	double x,y,z;
	int proton_id;
	char s[1000];
	int* proton_ids = calloc(size, sizeof(int));
	k = 0;
	
	// Fill R cells with dij^6
	for (i = 0; i < size; i++)
	{
		double xo, yo, zo;
		for (j = 0; j <= i; j++)
			fscanf(f, "%i\t%s\t%lf\t%lf\t%lf\n", &proton_id, s, &xo, &yo ,&zo);
		//printf("%i --> \n", proton_id);
		proton_ids[i] = proton_id;
		
		while (fscanf(f, "%i\t%s\t%lf\t%lf\t%lf\n", &proton_id, s, &x, &y ,&z) != EOF)
		{
			k = i*size+j;
			rx_mat[k] = (x-xo)*(x-xo)+(y-yo)*(y-yo)+(z-zo)*(z-zo);
			rx_mat[k] = rx_mat[k]*rx_mat[k]*rx_mat[k];
			rx_mat[j*size+i] = rx_mat[k];
			//printf(" --> %i: %.3lf\n", proton_id, rx_mat[k]);
			j++;
		}
		rewind(f);
	}
	fclose(f);
	
	// Read proton groups
	proton_groups* eq_groups;
	eq_groups = calloc(1, sizeof(proton_groups));
	
	read_proton_groups(eq_groups, eq_groups_path, proton_ids, size);
	//print_proton_groups(eq_groups);
	
	// Fill diagonal elements in R
	double coef = zero_freq_coef+2.*sing_freq_coef+doub_freq_coef;
	for (i = 0; i < size; i++)
	{
		int id = i*size+i;
		for (j = 0; j < size; j++)
			if (i != j)
				rx_mat[id] += 1./rx_mat[i*size+j];
		rx_mat[id] *= coef;	
	}
	
	// Substitute excessive terms from diagonal elements
	for (int counter = 0; counter < eq_groups->N; counter++)
	{
		proton_group* grp_tmp = &eq_groups->groups[counter];
		for (i = 0; i < grp_tmp->N; i++)
			for (j = 0; j < grp_tmp->N; j++)
				if (i != j)
				{
					rx_mat[grp_tmp->group[i]*size+grp_tmp->group[i]] -= coef/rx_mat[grp_tmp->group[i]*size+grp_tmp->group[j]];
				}
	}
	
	// Fill the rest of R
	for (i = 0; i < size; i++)
		for (j = i; j < size; j++)
			if (i != j)
			{
				rx_mat[i*size+j] = (doub_freq_coef-zero_freq_coef)/rx_mat[i*size+j];
				rx_mat[j*size+i] = rx_mat[i*size+j];
			}
	
	// Set cross-relaxation for protons within each group to 0	
	for (int counter = 0; counter < eq_groups->N; counter++)
	{
		proton_group* grp_tmp = &eq_groups->groups[counter];
		for (i = 0; i < grp_tmp->N; i++)
			for (j = 0; j < grp_tmp->N; j++)
				if (i != j)
				{
					rx_mat[grp_tmp->group[i]*size+grp_tmp->group[j]] = 0.0;
				}
	}
	
	
	// Get proton pairwise intencity
	gsl_matrix* in = get_intensity(rx_mat, size, mixing_time);
	
	//double test_mat[] = {1, 2, 3, 4, 2, 4, 2, 7, 3, 2, 1, 3, 4, 7, 3, 2};
	//my_fprint_matrix(stdout, get_intensity(test_mat, 4, 1.00)->data, 4);
	//getchar();
	
	// Get group pairwise intencity
	double* grouped_in = calloc(eq_groups->N*eq_groups->N, sizeof(double));
	collapse_groups(grouped_in, in, eq_groups);
	
	FILE* comp_file = fopen("sandbox/computed_matrix", "a");
	my_fprint_matrix(comp_file, grouped_in, eq_groups->N);
	fclose(comp_file);
	
	// Add line to output file
	FILE* output = fopen(output_path, "a");
	fprintf(comp_file, "%s\t", complex_path);
	for (i = 0; i < argc-4; i++)
	{
		char* exp_path = argv[4+i];
		double chi = chi_score(exp_path, grouped_in, eq_groups->N);
		if (i < argc-5)
			fprintf(comp_file, "%.1lf\t", chi);
		else
			fprintf(comp_file, "%.1lf\n", chi);
	}
	fclose(output);
	
	// Free
	free(proton_ids);
	free_proton_groups(eq_groups);
	free(grouped_in);
	gsl_matrix_free(in);
    return 0;
}
