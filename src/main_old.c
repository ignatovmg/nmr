#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <string.h>
#include <assert.h>
#include <math.h>

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
		printf("File not found\n");
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
		printf("proton %i not found\n", id);
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
		printf("%s not found\n", path);
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

void collapse_groups(double* dst, double* src, int src_size, proton_groups* grps)
{
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
					dst[i*dst_size+j] += src[grp_i->group[i_id]*src_size + grp_j->group[j_id]];
		}
	}
} 

void collapse_groups_from_gsl(double* dst, gsl_matrix* src, proton_groups* grps)
{
	int src_size = src->size1;
	collapse_groups(dst, src->data, src_size, grps);
}

void free_proton_groups(proton_groups* groups)
{
	free(groups->groups);
	free(groups);
}

/*
 * singular decomposition of R relaxation matrix
 * R = E*S*E_T, followed by calculating the intensity matrix
 * I = exp^(-R*mix_time) using formula 
 * I = E*exp^(-L*mix_time)*E_T, where L = 
 * = identity_matrix*eigenvalues_vector
**/
gsl_matrix* get_intensity(double* input, int size, double mix_time)
{
	int i, j, k;
	
	for (i = 0; i < size*size; i++)
			input[i] *= -mix_time;
			
	gsl_matrix_view m = gsl_matrix_view_array(input, size, size);
			
	//print_matrix(&m.matrix);

	gsl_vector *eval = gsl_vector_alloc(size);
	gsl_matrix *evec = gsl_matrix_alloc(size, size);
	gsl_matrix *res = gsl_matrix_alloc(size, size);
	gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(size);
	gsl_eigen_symmv(&m.matrix, eval, evec, w);
	
	//print_matrix(evec);
	
	//for (i = 0; i < size; i++)
	//	eval->data[i] = exp(-mix_time*eval->data[i]);
	for (i = 0; i < size; i++)
		eval->data[i] = exp(eval->data[i]);
	
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
			fprintf(f, "%.6f\t", m[i*size+j]);
		fprintf(f, "%.6f\n", m[i*size+size-1]);
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

void my_pprint_matrix(char* path, double* m, int size)
{
	FILE* file = fopen(path, "w");
	my_fprint_matrix(file, m, size);
	fclose(file);
}

#ifdef TEST
// For tests only
double best_k(char* path, double* m, int size)
{
	FILE* f = fopen(path, "r");
	if (f == NULL)
	{
		printf("%s not found\n", path);
		exit(0);
	}
	
	int i, j;
	double val;
	double num = 0.0;
	double den = 0.0;
	while(fscanf(f, "%i\t%i\t%lf\n", &i, &j, &val) != EOF)
	{
		num += val*m[i*size+j];
		//num += val;
		den += m[i*size+j]*m[i*size+j];
		//den += m[i*size+j];
	}
	rewind(f);
	double k = num/den;
	
	return k;
}
#endif

double chi_score(char* path, double* m, int size, proton_groups* eq_groups)
{
	float power = 1.0/6.0;

	FILE* f = fopen(path, "r");
	if (f == NULL)
	{
		printf("%s not found\n", path);
		exit(1);
	}
	
	int i, j;
	double val;
	double num = 0.0;
	double den = 0.0;	
	while(fscanf(f, "%i\t%i\t%lf\n", &i, &j, &val) != EOF)
	{
		num += val*m[i*size+j];
		//num += val;
		den += m[i*size+j]*m[i*size+j];
		//yden += m[i*size+j];
	}
	rewind(f);
	double k = num/den;
	
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
	double chi = 0.0;
	double nrm = 0.0;
	int n = 0;
	while(fscanf(f, "%i\t%i\t%lf\n", &i, &j, &val) != EOF)
	{
		trm = 0.0;
		if (fabs(val) > 0.0)
			trm  = (val / fabs(val)) * pow(fabs(val), power);
		if (fabs(k*m[i*size+j]) > 0.0)
			trm -= (k*m[i*size+j] / fabs(k*m[i*size+j])) * pow(fabs(k*m[i*size+j]), power);

		chi += trm*trm;
		nrm += pow(fabs(val), power);
		n++;
	}
	
	fclose(f);
	return chi/nrm;
}


int main(int argc, char** argv)
{
#ifndef METHOD
	printf("How to treat unresolved groups? (0 or 1)\n");
	exit(1);
#endif

	int i, j, k;          
/*
 * zero_freq_coef = J(0)*gamma^4*h^2/(4pi^2*10)
 * sing_freq_coef = 3./2.*J(1*omega)*gamma^4*h^2/(4pi^2*10)
 * doub_freq_coef = 6.*J(2*omega)*gamma^4*h^2/(4pi^2*10)
**/
	//double zero_freq_coef = 5.67426; // TODO: Decide about the constants
	//double sing_freq_coef = 3.0/2.0*5.67426/1.0047;  //
	//double doub_freq_coef = 6.0*5.67426/1.0196;  //
	
	double omega       = 0.0006;    // *10^8 Gz
	double t_corr      = 0.1;       // *10^(-9) s
	double mixing_time = 0.2; // s
	
	double J_0 = t_corr; // *10^(-9) s
	//double J_1 = t_corr / (1.0 + pow(       omega * t_corr / 10.0, 2));  // *10^(-9) s
	//double J_2 = t_corr / (1.0 + pow( 2.0 * omega * t_corr / 10.0, 2));  // *10^(-9) s
	
	double J_1 = t_corr / (1.0 +  4.0 / ( M_PI * M_PI ) * omega * omega * t_corr * t_corr );  // *10^(-9) s
	double J_2 = t_corr / (1.0 + 16.0 / ( M_PI * M_PI ) * omega * omega * t_corr * t_corr );  // *10^(-9) s
	
	//double coeff = 5.838; // *10^7 <--- gamma^4*h^2/(4pi^2*10*10^(-48))
	double coeff = 56.95; // *10^7 <--- gamma^4*h^2/(4pi^2*10*10^(-48))
	//double coeff = 7.17*100;
	//double coeff = 34.986; 
	//coeff = 56.65;
	
	double zero_freq_coef =       J_0 * coeff; // 1/s
	double sing_freq_coef = 1.5 * J_1 * coeff; // 1/s
	double doub_freq_coef = 6.0 * J_2 * coeff; // 1/s
	
#ifndef METHOD
	printf("Ω_2 + 2*Ω_1 + Ω_0 = %.6f\n", zero_freq_coef + 2*sing_freq_coef + doub_freq_coef);
	printf("Ω_2 - Ω_0 = %.6f\n", doub_freq_coef - zero_freq_coef);
	printf("Ω_0 = %.6f\nΩ_1 = %.6f\nΩ_2 = %.6f\n", zero_freq_coef, sing_freq_coef, doub_freq_coef);
	printf("J_0 = %f\n", J_0);
	printf("J_1 = %f\n", J_1);
	printf("J_2 = %f\n", J_2);
#endif
	
	char* eq_groups_path = argv[1]; // path to equivalent proton groups
	char* complex_path   = argv[2]; // path to peptide
	char* output_path    = argv[3];
	
	//mprintf("Computing complex %s\n", complex_path);
	//printf("argc=%i\n", argc);
	
	int size = get_line_num(complex_path);
	double* rx_mat = calloc(size*size, sizeof(double));
	memset(rx_mat, 0, sizeof(double)*size*size);
	
	FILE* f = fopen(complex_path, "r");
	if (f == NULL)
	{
		printf("%s not found\n", complex_path);
		exit(1);
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
			rx_mat[k] = 0.0;
			
			if (i != j)
			{
				rx_mat[k] = (x-xo)*(x-xo)+(y-yo)*(y-yo)+(z-zo)*(z-zo);
				rx_mat[k] = 1.0/(rx_mat[k]*rx_mat[k]*rx_mat[k]);
				rx_mat[j*size+i] = rx_mat[k];
			}
			
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

	if (METHOD == 1)
	{
		double* rx_mat_tmp = calloc(eq_groups->N*eq_groups->N, sizeof(double));
		collapse_groups(rx_mat_tmp, rx_mat, size, eq_groups);
		
		free(rx_mat);
		rx_mat = rx_mat_tmp;
		
		size = eq_groups->N;
		
		for (i = 0; i < size; i++)
			for (j = 0; j < size; j++)
				if (i != j)
					rx_mat[i*size+j] /= (eq_groups->groups[i].N * eq_groups->groups[j].N);
				else
					// TODO: this value 85.45 can be used only with t_c = 1.0 & w = 0.0006 !!
					rx_mat[i*size+j] *=  85.45 / eq_groups->groups[i].N;  // 898 if t_c = 1.0 & w = 0.0006
				
	}
	
	// Print distance matrix
#ifdef TEST
	double* dist = calloc(size*size, sizeof(double));
	
	for (int i = 0; i < size*size; i++)
		dist[i] = pow(1.0/rx_mat[i], 1.0/6.0);
		
	my_pprint_matrix("sandbox/distance", dist, size);
	free(dist);
#endif
	
	
	// Fill diagonal elements in R
	double coef = zero_freq_coef + 2 * sing_freq_coef + doub_freq_coef;
	
	if (METHOD == 1)
	{
		for (i = 0; i < size; i++)
		{
			int id = i*size+i;    //    t_c = 1 & w = 6       t_c = 0.1 & w = 6
		  //rx_mat[id] *= 85.45;  //    23.0469;              58.38
		
			for (j = 0; j < size; j++)
				if (i != j)
					rx_mat[id] += rx_mat[i*size+j] * coef * eq_groups->groups[j].N;
		}
	}
	else
	{
		for (i = 0; i < size; i++)
		{
			int id = i*size+i;
			rx_mat[id] = 0.0;
			for (j = 0; j < size; j++)
				if (i != j)
					rx_mat[id] += rx_mat[i*size+j] * coef;
		}

		// Substract excessive terms from diagonal elements (not effective for method 1)
		for (int counter = 0; counter < eq_groups->N; counter++)
		{
			proton_group* grp_tmp = &eq_groups->groups[counter];
			for (i = 0; i < grp_tmp->N; i++)
				for (j = 0; j < grp_tmp->N; j++)
					if (i != j)
					{
						rx_mat[grp_tmp->group[i]*size+grp_tmp->group[i]] -= coef * rx_mat[grp_tmp->group[i] * size + grp_tmp->group[j]];
					}
		}
	}
	
	// Fill the rest of R
	coef = doub_freq_coef - zero_freq_coef;
	
	if (METHOD == 1)
	{
		for (i = 0; i < size; i++)
			for (j = i; j < size; j++)
				if (i != j)
				{
					rx_mat[i*size+j] = rx_mat[i*size+j] * coef * (eq_groups->groups[i].N * eq_groups->groups[j].N);
					rx_mat[j*size+i] = rx_mat[i*size+j];
				}
	}
	else
	{
		// Set cross-relaxation for protons within each group to 0	
		for (int counter = 0; counter < eq_groups->N; counter++)
		{
			proton_group* grp_tmp = &eq_groups->groups[counter];
			for (i = 0; i < grp_tmp->N; i++)
				for (j = 0; j < grp_tmp->N; j++)
					if (i != j)
						rx_mat[grp_tmp->group[i]*size+grp_tmp->group[j]] = 0.0;
		}
	}

		
#ifdef TEST
	my_pprint_matrix("sandbox/relaxation_matrix", rx_mat, size);
#endif

	// Get proton pairwise intencity
	gsl_matrix* in = get_intensity(rx_mat, size, mixing_time);
	
	if (METHOD == 1)
	{
		for (i = 0; i < size; i++)
			for (j = 0; j < size; j++)
			{
				if (i == j)
					in->data[i*size+j] *= (double)(eq_groups->groups[i].N);
				else
					in->data[i*size+j] *= 1;
			}
	}
	
#ifdef TEST
	my_pprint_matrix("sandbox/computed_matrix_raw", in->data, size);
#endif
	
	//double test_mat[] = {1, 2, 3, 4, 2, 4, 2, 7, 3, 2, 1, 3, 4, 7, 3, 2};
	//my_fprint_matrix(stdout, get_intensity(test_mat, 4, 1.00)->data, 4);
	//getchar();
	
	// Get group pairwise intencity
	double* grouped_in;
	
	if (METHOD == 1)
	{
		grouped_in = in->data;
	}
	else
	{
		grouped_in = calloc(eq_groups->N*eq_groups->N, sizeof(double));
		collapse_groups_from_gsl(grouped_in, in, eq_groups);
	}
	
#ifdef TEST
	double multip = best_k(argv[4], grouped_in, eq_groups->N);
	for (int i = 0; i < eq_groups->N*eq_groups->N; i++)
		grouped_in[i] = multip*grouped_in[i];	

	my_pprint_matrix("sandbox/computed_matrix", grouped_in, eq_groups->N);
#endif	
	
#ifdef TEST
	// TMP
	float max_val = 0.0;
	float mean = 0.0;
	int count = 0;
	for (i = 0; i < eq_groups->N-1; i++)
		for (j = i+1; j < eq_groups->N; j++)
		{
			if (fabs(grouped_in[i*eq_groups->N+j]) >= max_val)
				max_val = fabs(grouped_in[i]);
				
			mean += fabs(grouped_in[i*eq_groups->N+j]);
			count ++;
		}
			
	mean /= count;
	printf("mean : %.6f\n", mean);
#endif	
	
	// Add line to output file
	
	FILE* output = fopen(output_path, "a");
	fprintf(output, "%s\t", complex_path);
	double chi;
	for (i = 0; i < argc-4; i++)
	{
		char* exp_path = argv[4+i];
		chi = chi_score(exp_path, grouped_in, eq_groups->N, eq_groups);

		if (i < argc-5)
			fprintf(output, "%.6lf\t", chi);
		else
			fprintf(output, "%.6lf\n", chi);
	}
	fclose(output);

#ifdef TEST
	printf("chi: %.5f\n", chi);
#endif

	// Free
	if (METHOD != 1)
		free(grouped_in);
		
	free(proton_ids);
	free_proton_groups(eq_groups);
	free(rx_mat);
	gsl_matrix_free(in);
    return 0;
}
