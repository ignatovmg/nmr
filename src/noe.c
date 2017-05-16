#include "noe.h"

//                                                             e+9 Gz        e-9 s 
//double* rx_mat(double* rx, struct gradient* grad, struct mol_atom_group *ag, groups* grps, double omega, double t_cor);
double* rx_mat(struct noe* spect, struct mol_atom_group *ag)
{
	groups* grps = spect->grps;
	double* rx   = spect->rx;
	
	double omega = spect->omega;
	double t_cor = spect->t_cor;
	
	int size = grps->N;
	int loc, loc1;
	
	if (rx == NULL || grps == NULL)
	{
		fprintf(stderr, "NULL pointer in %s\n", __func__);
		exit(EXIT_FAILURE);
	}	
	
	double* gradx;
	double* grady;
	double* gradz;
	struct gradient* grad = spect->rx_grad;
	
	if (grad != NULL)
	{
		gradx = grad->gradx;
		grady = grad->grady;
		gradz = grad->gradz;
		
		memset(gradx, 0, grad->natoms*size*size*sizeof(double));
		memset(grady, 0, grad->natoms*size*size*sizeof(double));
		memset(gradz, 0, grad->natoms*size*size*sizeof(double));
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
	
	group *grpi, *grpj;	
	struct mol_vector3 *coordi, *coordj;
	double r2;
	double mult, mult1;
	
	int atomi = 0, atomj;
	double counter;
	for (int i = 0; i < size; i++)
	{
		grpi = &grps->groups[i];
		atomj = atomi;
		for (int j = i; j < size; j++)
		{
			grpj = &grps->groups[j];
			
			counter = 0.0;
			rx[size*i+j] = 0.0;
			rx[size*j+i] = 0.0;
			
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
							
							if (r2 > __DMIN__)
							{
								rx[size*i+j] += 1.0 / (r2 * r2 * r2);
								counter += 1.0;
								
								// accumulate gradient
								if (grad != NULL)
								{
									double r8 = r2 * r2 * r2 * r2;
									
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
							
						if (r2 > __DMIN__)
						{
							rx[size*i+j] += 1.0 / (r2 * r2 * r2);
							counter += 1.0;
							
							// accumulate gradient
							if (grad != NULL)
							{
								double r8 = r2 * r2 * r2 * r2;	
								
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

double* matrix_exp(struct noe* spect)
{
	double* in = spect->in;
	double* rx = spect->rx;
	int* mask = spect->mask;
	
	struct gradient* in_grad = spect->in_grad;
	struct gradient* rx_grad = spect->rx_grad;

	int size = spect->N; 
	double t_mix = spect->t_mix;

	if (in == NULL || rx == NULL || mask == NULL)
	{
		fprintf(stderr, "NULL pointer in %s\n", __func__);
		exit(EXIT_FAILURE);
	}	
		
	gsl_matrix_view m = gsl_matrix_view_array(rx, size, size);	
			
	//print_matrix(&m.matrix);

	gsl_vector *eval = gsl_vector_alloc(size);
	gsl_matrix *evec = gsl_matrix_alloc(size, size);
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
			if (mask[i*size+j] == 1)
			{
				in[i*size+j] = 0.0;
			
				for (int k = 0; k < size; k++)
				{
					in[i*size+j] += expeval[k]*evec->data[i*size+k]*evec->data[j*size+k];
				}
				
				in[j*size+i] = in[i*size+j];
			}
			else
			{
				in[i*size+j] = 0.0;
				in[j*size+i] = 0.0;
			}
		}
	}
	
	if (rx_grad != NULL && in_grad != NULL)
	{
		double* helperx = calloc(size * size, sizeof(double));
		double* helpery = calloc(size * size, sizeof(double));
		double* helperz = calloc(size * size, sizeof(double));
		
		double rux, ruy, ruz, ijx, ijy, ijz, coef;
		int loc;
		
		for (int atm = 0; atm < in_grad->natoms; atm++)
		{
			// http://ac.els-cdn.com/0022236489903600/1-s2.0-0022236489903600-main.pdf?_tid=1fc4f52e-2b89-11e7-a148-00000aacb35d&acdnat=1493325540_544dadaafac0507df3318e832cea0248
			for (int r = 0; r < size; r++)
			{
				for (int u = 0; u < size; u++)
				{
					rux = 0.0;
					ruy = 0.0;
					ruz = 0.0;
					
					for (int s = 0; s < size; s++)
					{
						for (int t = 0; t < size; t++)
						{
							if (r != u)
							{
								coef = -(expeval[r] - expeval[u]) / \
								       (eval->data[r] - eval->data[u]) * t_mix;
							}
							else
							{
								coef = - expeval[r] * t_mix;
							}
							
							coef *= evec->data[s*size+r] * evec->data[t*size+u];
							
							rux += coef * rx_grad->gradx[atm*size*size + s*size + t];
							ruy += coef * rx_grad->grady[atm*size*size + s*size + t];
							ruz += coef * rx_grad->gradz[atm*size*size + s*size + t];
						}
					}
					
					helperx[r*size+u] = rux;
					helpery[r*size+u] = ruy;
					helperz[r*size+u] = ruz;
				}
			}
			
			
			for (int i = 0; i < size; i++)
			{
				for (int j = i; j < size; j++)
				{
					ijx = 0.0;
					ijy = 0.0;
					ijz = 0.0;
					
					if (mask[i*size+j] == 1)
					{
						for (int r = 0; r < size; r++)
						{
							for (int u = 0; u < size; u++)
							{
								ijx += helperx[r*size+u] * \
									   evec->data[i*size+r] * \
									   evec->data[j*size+u];
									   
								ijy += helpery[r*size+u] * \
									   evec->data[i*size+r] * \
									   evec->data[j*size+u];
							
								ijz += helperz[r*size+u] * \
									   evec->data[i*size+r] * \
									   evec->data[j*size+u];
							}
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
		
		free(helperx);
		free(helpery);
		free(helperz);
	}
		
	gsl_eigen_symmv_free(w);
	gsl_vector_free(eval);	
	gsl_matrix_free(evec);
	free(expeval);	

    return in;
}

double* peaks(struct noe* spect, struct mol_atom_group *ag)
{
	double* in = spect->in; 
	double* rx = spect->rx; 
	groups* grps = spect->grps; 
	
	int size = grps->N;
	
	if (in == NULL || rx == NULL || grps == NULL)
	{
		fprintf(stderr, "NULL pointer in %s\n", __func__);
		exit(EXIT_FAILURE);
	}	

	// Relaxation matrix
	rx = rx_mat(spect, ag);
			
	// Intencities
	in = matrix_exp(spect);	
		
	// Normalization
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			in[i*size+j] *= sqrt(grps->groups[i].N * grps->groups[j].N);
		}
	}
	
	if (spect->in_grad != NULL)
	{
		for (int k = 0; k < grps->natoms; k++)
		{
			for (int i = 0; i < size; i++)
			{
				for (int j = 0; j < size; j++)
				{
					spect->in_grad->gradx[k*size*size+i*size+j] *= sqrt(grps->groups[i].N * grps->groups[j].N);
					spect->in_grad->grady[k*size*size+i*size+j] *= sqrt(grps->groups[i].N * grps->groups[j].N);
					spect->in_grad->gradz[k*size*size+i*size+j] *= sqrt(grps->groups[i].N * grps->groups[j].N);
				}
			}
		}
	}
	
	return in;
}

/*gsl_matrix* rx_perturb(struct mol_atom_group *ag, groups* grps, double* rx_dst, double* rx_src, double step, int grp_id, int atom_id, double omega, double t_cor)
{
	int size = grps->N;

	double J_0 = t_cor;
	double J_1 = t_cor / (1.0 +  4.0 * ( M_PI * M_PI ) * omega * omega * t_cor * t_cor );  // *10^(-9) s
	double J_2 = t_cor / (1.0 + 16.0 * ( M_PI * M_PI ) * omega * omega * t_cor * t_cor );  // *10^(-9) s
	
	double gyr  = 2.6751965; // e+8
	double hbar = 1.0545919; // e-34;
	double coef = pow(gyr, 4)*pow(hbar, 2);
	
	double  S0 =       J_0 * coef; // 1/s
	double  S1 = 1.5 * J_1 * coef; // 1/s
	double  S2 = 6.0 * J_2 * coef; // 1/s

	// Relaxation matrix
	struct mol_vector3 *crd1 = &ag->coords[grps->groups[grp_id]->group[atom_id]];
	
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			rx_dst[i*size+j] = rx_src[i*size+j];
			
			if (j == grp_id)
			{
				if (i == j)
				{
					if (grps->groups[grp_id].N > 1)
					{
						for (int k = 0; k < grps->groups[grp_id].N; k++)
						{
							struct mol_vector3 *crd2 = &ag->coords[grps->groups[grp_id]->group[k]];
					 
							r2 = (crd2->X - crd1->X) * (crd2->X - crd1->X) + \
								 (crd2->Y - crd1->Y) * (crd2->Y - crd1->Y) + \
								 (crd2->Z - crd1->Z) * (crd2->Z - crd1->Z);
								 
							if (r2 > __DMIN__)
							{
								
							}
							 
						}
						
						rx[size*i+j] *= 2 * grpi->N * (S1 + S2);
					}
					
				}
			}
		}
	}
	
	// Intencities
	gsl_matrix* in = matrix_exp(rx, size, t_mix, NULL, NULL);
	//free(rx);
		
	// Normalization
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			//in->data[i*size+j] *= sqrt(grps->groups[i].N * grps->groups[j].N);
		}
	}
	
	return in;
}*/

double* r2_mat(struct mol_atom_group *ag, groups* grps)
{
	int size = grps->natoms;
	double* r2 = calloc(size*size, sizeof(int));
	
	struct mol_vector3 *crdi, *crdj;
	int grpi = 0, atomi = 0, grpj, atomj;
	
	for (int i = 0; i < size; i++)
	{
		crdi = &ag->coords[grps->groups[grpi].group[atomi]];
		if (atomi == grps->groups[grpi].N)
		{
			atomi = 0;
			grpi++;
		}
		
		grpj  = grpi;
		atomj = atomi;
		for (int j = i; j < size; j++)
		{
			crdj = &ag->coords[grps->groups[grpj].group[atomj]];
			if (atomj == grps->groups[grpj].N)
			{
				atomj = 0;
				grpj++;
			}
			
			r2[i*size+j] = (crdj->X - crdi->X) * (crdj->X - crdi->X) + \
						   (crdj->Y - crdi->Y) * (crdj->Y - crdi->Y) + \
						   (crdj->Z - crdi->Z) * (crdj->Z - crdi->Z);
						   
			r2[j*size+i] = r2[i*size+j];
			atomj++;
		}
		
		atomi++;
	}
	
	return r2;
}

void grad_numeric_in(struct noe* spect, struct mol_atom_group *ag)
{
	groups* grps = spect->grps;
	int size = grps->N;
	int msize = size * size;

	// save gradients pointers
	struct gradient* in_grad = spect->in_grad;
	struct gradient* rx_grad = spect->rx_grad;
	
	// make NULL so they will not be computed analytically
	spect->in_grad = NULL;
	spect->rx_grad = NULL;

	double* pks1 = peaks(spect, ag);
	
	spect->in = calloc(size*size, sizeof(double));
	double* pks2;

	double  step = 0.000001;
	int counter = 0;
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < grps->groups[i].N; j++)
		{
			int atom = grps->groups[i].group[j];
			
			// fill grad for X
			ag->coords[atom].X += step;
			pks2 = peaks(spect, ag);
			ag->coords[atom].X -= step;
			
			for (int k = 0; k < size * size; k++)
			{
				in_grad->gradx[counter * size * size + k] = (pks2[k] - pks1[k]) / step;
			}
			
			// fill grad for Y
			ag->coords[atom].Y += step;
			pks2 = peaks(spect, ag);
			ag->coords[atom].Y -= step;
			
			for (int k = 0; k < size * size; k++)
			{
				in_grad->grady[counter * size * size + k] = (pks2[k] - pks1[k]) / step;
			}

			// fill grad for Z
			ag->coords[atom].Z += step;
			pks2 = peaks(spect, ag);
			ag->coords[atom].Z -= step;
			
			for (int k = 0; k < size * size; k++)
			{
				in_grad->gradz[counter * size * size + k] = (pks2[k] - pks1[k]) / step;
			}
			
			counter++;//= msize;
		}
	}

	free(pks2);
	spect->in = pks1;
	spect->in_grad = in_grad;
	spect->rx_grad = rx_grad;
}

double* grad_numeric_rx(struct noe* spect, struct mol_atom_group *ag)
{
	groups* grps = spect->grps;
	int size = grps->N;
	
	struct gradient* rx_grad = spect->rx_grad;
	struct gradient* in_grad = spect->in_grad;
	spect->in_grad = NULL;
	spect->rx_grad = NULL;
	
	rx_mat(spect, ag);
	double* rx1 = calloc(size*size, sizeof(double));
	memcpy(rx1, spect->rx,  size*size*sizeof(double));
	
	double* rx2;
	
	double  step = 0.0000000001;
	double* grad = calloc(3 * grps->natoms * size * size, sizeof(double));
	
	int counter = 0;
	
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < grps->groups[i].N; j++)
		{
			int atom = grps->groups[i].group[j];
			
			// fill grad for X
			ag->coords[atom].X += step;
			rx2 = rx_mat(spect, ag);
			ag->coords[atom].X -= step;
			
			for (int k = 0; k < size * size; k++)
			{
				grad[0 * grps->natoms * size * size + counter * size * size + k] = (rx2[k] - rx1[k]) / step;
			}
			
			// fill grad for Y
			ag->coords[atom].Y += step;
			rx2 = rx_mat(spect, ag);
			ag->coords[atom].Y -= step;
			
			for (int k = 0; k < size * size; k++)
			{
				grad[1 * grps->natoms * size * size + counter * size * size + k] = (rx2[k] - rx1[k]) / step;
			}
			
			// fill grad for Z
			ag->coords[atom].Z += step;
			rx2 = rx_mat(spect, ag);
			ag->coords[atom].Z -= step;
			
			for (int k = 0; k < size * size; k++)
			{
				grad[2 * grps->natoms * size * size + counter * size * size + k] = (rx2[k] - rx1[k]) / step;
			}
			
			counter++;
		}
	}
	
	free(rx1);
	//free(rx2);
	
	spect->in_grad = in_grad;
	spect->rx_grad = rx_grad;
	
	return grad;
}

double* read_exp(char* path, int size)
{
	FILE* f = fopen(path, "r");
	if (f == NULL)
	{
		fprintf(stderr, "%s not found\n", path);
		exit(EXIT_FAILURE);
	}
	
	double* exp = calloc(size*size, sizeof(double));
	
	int i, j, nchar;
	double val;
	while((nchar = fscanf(f, "%i\t%i\t%lf\n", &i, &j, &val)) != EOF)
	{
		if (nchar != 3)
		{
			fprintf(stderr, "%s: Wrong input\n", __func__);
			exit(EXIT_FAILURE);
		}
		
		if (fabs(val) > __DMIN__)
		{
			exp[i*size+j] = val;
			exp[j*size+i] = val;
		}
	}
	
	fclose(f);
	
	return exp;
}

double best_multiplier(double* exp, double* spc, int size)
{
	double num = 0.0;
	double den = 0.0;
	
	for (int i = 0; i < size; i++)
	{
		for (int j = i; j < size; j++)
		{
			if (fabs(exp[i*size+j]) > __DMIN__)
			{
				num += exp[i*size+j]*spc[i*size+j];
				//num += val;
				den += spc[i*size+j]*spc[i*size+j];
				//den += m[i*size+j];
			}
		}
	}

	double k = num / den;
	
	return k;
}

double fit_score(struct mol_vector3* grad, struct noe* spect, double gcoef)
{
	int natoms = spect->grps->natoms;
	int size = spect->N;

	double* exp = spect->exp;
	double* cmp = spect->in;
	double* gradx;
	double* grady;
	double* gradz;
	int* atoms;
	struct mol_vector3 grd;
	
	if (exp == NULL || cmp == NULL)
	{
		fprintf(stderr, "NULL pointer in %s\n", __func__);
		exit(EXIT_FAILURE);
	}

	float power = 1.0 / 6.0;
	double k = best_multiplier(exp, cmp, size);
	
	double trm, val, sgn, cmp_val, cmp_val6, coef;
	double fit = 0.0;
	double nrm = 0.0;
	int n = 0;
	
	//double* test = malloc(3*natoms*sizeof(double));
	//memset(test, 0, 3*natoms*sizeof(double));
	
	if (grad != NULL)
	{
		gradx = spect->in_grad->gradx;
		grady = spect->in_grad->grady;
		gradz = spect->in_grad->gradz;
	
		atoms = spect->grps->atoms;
	}
	
	for (int i = 0; i < size; i++)
	{
		for (int j = i; j < size; j++)
		{
			nrm += pow(fabs(exp[i*size+j]), power);
		}
	}
	
	for (int i = 0; i < size; i++)
	{
		for (int j = i; j < size; j++)
		{
			val = exp[i*size+j];
			if (fabs(val) > __DMIN__)
			{
				trm = (val / fabs(val)) * pow(fabs(val), power);
			
				cmp_val  = fabs(k * cmp[i*size+j]);
				cmp_val6 = pow(cmp_val, power);
				sgn = (k * cmp[i*size+j] / cmp_val);
				
				trm -= sgn * cmp_val6;
				fit += trm * trm;
				n++;
			
				if (grad != NULL)
				{
					coef = 2 * power * k * trm * cmp_val6 / cmp_val;
					for (int atm = 0; atm < natoms; atm++)
					{
						// SUBTRACTING gradient
						grd.X = gcoef * (- coef * gradx[atm*size*size + i*size + j]) / nrm;
						grd.Y = gcoef * (- coef * grady[atm*size*size + i*size + j]) / nrm;
						grd.Z = gcoef * (- coef * gradz[atm*size*size + i*size + j]) / nrm;
					
						//test[3*atm+0] -= grd.X;
						//test[3*atm+1] -= grd.Y;
						//test[3*atm+2] -= grd.Z;
						
						MOL_VEC_SUB(grad[atoms[atm]], grad[atoms[atm]], grd);
					}
				}
			}
		}
	}
	
	/*for (int atm = 0; atm < natoms; atm++)
	{
		printf("%.6f ", test[3*atm+0]);
	}
	printf("\n\n");
	free(test);*/
	
	return fit / nrm;
}


int* get_mask(double* exp, int size)
{
	int* mask = calloc(size * size, sizeof(int));

	for (int i = 0; i < size; i++)
	{
		for (int j = i; j < size; j++)
		{
			if (fabs(exp[i*size+j]) > __DMIN__)
			{
				mask[i*size+j] = 1;
				mask[j*size+i] = 1;
			} 
			else
			{
				mask[i*size+j] = 0;
				mask[j*size+i] = 0;
			}
		}
	}
	
	return mask;
}

struct noe* init_noe(char* grp)
{
	struct noe* spect = malloc(sizeof(struct noe));
	
	// Read proton groups
	spect->grps = read_proton_groups(grp);
	int size = spect->grps->N;

	spect->N = size;
	spect->in = calloc(size*size, sizeof(double)); 
	spect->rx = calloc(size*size, sizeof(double));
	
	return spect;
}

void free_noe(struct noe* spect)
{
	free_proton_groups(spect->grps);
	free(spect->in);
	free(spect->rx);
	
	if (spect->in_grad == NULL)
	{
		free_gradient(spect->in_grad);
	}
	
	if (spect->rx_grad == NULL)
	{
		free_gradient(spect->rx_grad);
	}
	
	if (spect->exp == NULL)
	{
		free(spect->exp);
	}
	
	if (spect->mask == NULL)
	{
		free(spect->mask);
	}
	
	free(spect);
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
	free(grps->atoms);
	
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
	
	int atm = 0;
	grps->atoms = calloc(natoms, sizeof(int));
	for (int i = 0; i < grps->N; i++)
	{
		for (int j = 0; j < grps->groups[i].N; j++)
		{
			grps->atoms[atm++] = grps->groups[i].group[j];
		}
	}

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
