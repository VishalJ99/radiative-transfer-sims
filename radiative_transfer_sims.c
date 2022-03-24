#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#define _USE_MATH_DEFINES // for C

// define funcs
void main();
double gen_uniform_deviate(int, int, long int);
double p_tau(double);
double inv_p_tau(double);

// define global vars
#define TAU_MAX 100000  

// params for the lcg random number func
#define GLOBAL_SEED 42
int a = 16807;
int c = 0;
long int m = 2147483647;

// define structs
typedef struct particle
{
   double x,y,z;
   float theta,phi;
} particle;

void main()
{
 // Q1a Generating P(τ)dτ = exp(−τ)dτ distribution
 
 // Rejection method:
 // -----------------------

 int N = 100000; 
 int x0 = 0; 
 int x1 = 10;
 int pmax = 1;
 double *X,*Y,*filtered_X;
 double uniform_deviate, norm_factor;
 X = (double*) malloc(N* sizeof(double));
 Y = (double*) malloc(N* sizeof(double));
 filtered_X = (double*) malloc(N* sizeof(double));
 
 // Generate N uniform deviates in the interval (x0,x1)
 norm_factor = x1-x0;
 for (int i = 0; i<N; i++)
 {
  uniform_deviate = gen_uniform_deviate(a,c,m);
  X[i] = (uniform_deviate * norm_factor) +  x0;
 }

 // Generate N uniform deviates y_i in the interval (0, pmax)
 for (int i = 0; i<N; i++)
 {
  uniform_deviate = gen_uniform_deviate(a,c,m);
  Y[i] = (uniform_deviate * pmax);
 }

 // if y_i < p_tau(x_i) accept x_i, else reject x_i 
 int j, num_rejected = 0;

 for (int i = 0; i<N; i++)
 {
  if (Y[i] < p_tau(X[i])) 
  {
   filtered_X[j] = X[i];
   j++;
  }
  
  else num_rejected++;
  
 }

 // write the first j elements of filtered_X to a text file
 FILE *pf;
 char rejection_outfile[50] = "rejection_method_deviates.txt";
 pf = fopen(rejection_outfile, "w+");
 for (int i = 0; i<j; i++)
 {
  fprintf(pf,"%lf\n",filtered_X[i]);
 }
 printf("Finished sampling distribution via rejection method\n");
 printf("%i/%i sample draws rejected\n",num_rejected,N);
 free(X);
 free(filtered_X);
 fclose(pf);
 
 // Direct mapping:
 // -----------------------
 // Take uniform deviates y_i and map them to deviates x_i using the law of ...
 double *direct_mapped_X;
 direct_mapped_X = (double*) malloc(N* sizeof(double));

 for (int i = 0; i<N; i++)
 {
  direct_mapped_X[i] = inv_p_tau(Y[i]);
 }
 
 // Write direct mapped deviates x_i to a text file
 char direct_map_outfile[] = "direct_mapping_method_deviates.txt";
 pf = fopen(direct_map_outfile, "w+");
 for (int i = 0; i<j; i++)
 {
  fprintf(pf,"%lf\n",direct_mapped_X[i]);
 }

 printf("Finished sampling distribution via direct mapping method\n");
 fclose(pf);

 // Q1b Monte Carlo scattering, isotropic 
	// -----------------------
	#define N_BINS 10
	double tau,step_size;
	struct particle photon; 
	float z_min = 0;
 float z_max = 1;
	float alpha = 10;
 float albedo = 1;
	int num_photons_simulated = 0;
	int num_photons_escaped = 0;
	int total_escaped = 1e6;
	double dtheta = 180. / N_BINS; // in degrees
	bool absorb_photon = false;
	int binned[N_BINS];
	int i_bin;
	double random_number;

	// clear out the theta bins first
	for (int i = 0; i < N_BINS; i++)
			binned[i] = 0;

	// run MC simulation
	while (num_photons_escaped != total_escaped)
	{
		// initialise photon at origin moving upwards
		photon.x = 0;
		photon.y = 0;
		photon.z = 0;

		photon.theta = gen_uniform_deviate(a,c,m) * M_PI / 2;
		photon.phi = gen_uniform_deviate(a,c,m) * 2 * M_PI; 
		
		// scatter photon till it escapes atmosphere or gets absorbed
		while ((photon.z >= 0) && (photon.z <= 1) && (absorb_photon != true))
			{	
				// draw tau value 
				uniform_deviate = gen_uniform_deviate(a,c,m);
				tau = inv_p_tau(uniform_deviate);
				
				// if tau is nan, set to a large value
				if (tau!=tau) tau = TAU_MAX;
				
				// take a step in the given orientation (theta,phi)
				step_size = tau / alpha;
				photon.x += step_size * sin(photon.theta) * sin(photon.phi);
				photon.y += step_size * sin(photon.theta) * cos(photon.phi);
				photon.z += step_size * sin(photon.theta);

				// check if photon is absorbed
				random_number = gen_uniform_deviate(a,c,m);
				if (random_number > albedo) absorb_photon = true;
				else
					{
						//draw theta and phi values and repeat
						photon.theta = M_PI * gen_uniform_deviate(a,c,m) - M_PI/2;
						photon.phi = gen_uniform_deviate(a,c,m) * 2 * M_PI; 
					}
			}
		// escaped in the right direction?
		if (photon.z > 1) 
		{
			// log escaped photon and bin theta val
			num_photons_escaped++;
			i_bin = (int) ((photon.theta * 180 / M_PI) + 90 ) / dtheta;
			binned[i_bin]++;
		}

		// log simulated photon
		num_photons_simulated++;
		printf("\r%i/%i photons made it through!",num_photons_escaped,total_escaped);
		fflush(stdout);
	}

	// simulation finished
	printf("\nFinished MC simulation of isotropic photon scattering! %i/%i photons made it through!\n",num_photons_escaped,num_photons_simulated);

	// save results to a txt file
 char isotropic_photon_scattering_outfile[] = "isotropic_photon_scattering.txt";
 pf = fopen(isotropic_photon_scattering_outfile, "w+");
 for (int i = 0; i<N_BINS; i++)
 {
  fprintf(pf,"%lf, %i\n",i*dtheta - 90,binned[i]);
 }

}

double gen_uniform_deviate(int a, int c, long int m)
{
 // Generates pseudo random numbers between 0 and 1 using a linear congruential generator
 static long long int seed = GLOBAL_SEED;
 seed = (a * seed + c) % m;
 return (double) seed / m;
}

double p_tau(double tau)
{
 // PDF for optical depths tau 
 return exp(-tau);
}

double inv_p_tau(double y)
{
 //inverted PDF for optical depths tau
 return -log(y);
}