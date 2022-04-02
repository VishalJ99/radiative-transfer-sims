#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#define _USE_MATH_DEFINES 


// define func prototypes
void main();
double gen_uniform_deviate(int, int, long int);
double p_tau(double);
double inv_p_tau(double);
double rayleigh_phase_func(double);
void init_Rz(double [3][3], double);
void init_Ry(double [3][3], double);
void matmul(int, int, int, double [*][*], double [*][*], double [*][*]);
void rotate_back_to_lab(double *, double *, double, double);

// define global consts
#define TAU_MAX 100000  
#define EPS 1.084202e-19
#define N_BINS_RAYLEIGH 50
#define N_BINS_ISOTROPIC 10
#define GLOBAL_SEED 42

// params for the lcg random number func
int a = 16807;
int c = 0;
long int m = 2147483647;

// define structs
typedef struct particle
{
   double x,y,z;
   double theta,phi;
} particle;

void main()
{
	/* Q1a Generating P(τ)dτ = exp(−τ)dτ distribution
	-------------------------------------------------*/
	int N = 1e6; // number of deviates to generate
	int x0 = 0; // x0 and x1 used to define domain pdf is sampled in
	int x1 = 10; 
	int pmax = 1;
	double *X_rejection, *X_direct, x, y;
	X_rejection = malloc(N*sizeof(double));
	X_direct = malloc(N*sizeof(double));
	FILE *pf;
	double norm_factor = x1-x0;
	int i = 0;
	int num_rejected = 0;
	
	/* Rejection method
	----------------------*/
	// start clock for method
	clock_t t;
    t = clock();
	
	while (i<N)
		{
			// Generate uniform deviate x in the interval (x0,x1)
			x = (gen_uniform_deviate(a,c,m) * norm_factor) +  x0;
			// Generate uniform deviate y in the interval (0, pmax)
			y = (gen_uniform_deviate(a,c,m) * pmax);
			// if y < p_tau(x) accept x, else reject x 
			while (y > p_tau(x)) 
				{
					y = (gen_uniform_deviate(a,c,m) * pmax);
					x = (gen_uniform_deviate(a,c,m) * norm_factor) +  x0;
					num_rejected++;
				}
			X_rejection[i] = x;
			i++;
		}
	
	// measure elapsed time
	t = clock() - t;
	double time_taken = ((double)t)/CLOCKS_PER_SEC; 
	
	// write results to txt file
	char rejection_outfile[50] = "rejection_method_deviates.txt";
	pf = fopen(rejection_outfile,"w+");
	for (int i = 0; i<N; i++) fprintf(pf,"%lf\n",X_rejection[i]);
	fclose(pf); 
	
	// output summary
	printf("Finished sampling %d points via rejection method in %lfs\n",N, time_taken);
	printf("%d sample draws rejected\n\n",num_rejected);

	/* Cumultive function mapping:
	-------------------------------------------------*/
	// Take uniform deviates y_i and map them to deviates x_i using the inverse cumulative distribution function
		
	// start clock for  method 
	t = clock();
	for (i = 0; i<N; i++)
		{
			// Generate uniform deviate y in the interval (0, pmax)
			y = (gen_uniform_deviate(a,c,m) * pmax);
			// map y to x
			X_direct[i] = inv_p_tau(y);
			// x = (gen_uniform_deviate(a,c,m) * norm_factor) +  x0;
			// write x to output file
			// fprintf(pf,"%lf\n",x);
		}

	// measure and print elapsed time
	t = clock() - t;
	time_taken = ((double)t)/CLOCKS_PER_SEC; 

	// write results to txt file
	char cumulative_map_outfile[] = "cumulative_mapping_method_deviates.txt";
	pf = fopen(cumulative_map_outfile,"w+");
	for (int i = 0; i<N; i++) fprintf(pf,"%lf\n",X_direct[i]);
	fclose(pf);

	// output summary
	printf("Finished sampling %d points via cumulative mapping method in %lfs\n",N, time_taken);

 	/* Q1b Monte Carlo isotropic scattering
	----------------------------------------- */
	double tau,step_size;
	struct particle photon; 
	double z_min = 0;
	double z_max = 1;
	double alpha = 10;
 	double albedo = 1;
	int num_photons_simulated = 0;
	int num_photons_escaped = 0;
	int total_escaped = 1e6;
	double dtheta = 90. / N_BINS_ISOTROPIC; // in degrees
	double dcostheta = 1. / N_BINS_ISOTROPIC;
	bool absorb_photon = false;
	int binned_theta_isotropic[N_BINS_ISOTROPIC];
	int binned_cos_theta_isotropic[N_BINS_ISOTROPIC];
	int theta_bin, costheta_bin;
	double random_number;
	
	// clear out the bins first
	for (int i = 0; i < N_BINS_ISOTROPIC; i++)
		{
			binned_theta_isotropic[i] = 0;
			binned_cos_theta_isotropic[i] = 0;
		}
	
	// file to store scattered angles and check if angles are isotropic
	int angle_counter = 0;
	FILE *pf2;
	pf2 = fopen("theta_vals.txt","w+");
	
	// run MC simulation
	printf("[INFO] Starting MC simulation of isotropic elastic scattering of %d photons\n",total_escaped);
	while (num_photons_escaped != total_escaped)
	{
		// initialise photon at origin moving upwards
		photon.x = 0;
		photon.y = 0;
		photon.z = 0;

		photon.theta = acos(gen_uniform_deviate(a,c,m));
		photon.phi = 2*M_PI*gen_uniform_deviate(a,c,m) - M_PI;
		
		int iteration = 1;
		// scatter photon till it escapes the atmosphere or gets absorbed
		while ((photon.z >= 0) && (photon.z <= 1) && (absorb_photon != true))
			{	
				// check if photon is absorbed
				random_number = gen_uniform_deviate(a,c,m);
				if (random_number > albedo) absorb_photon = true;
				
				//else scatter photon (unless first run)
				else if (iteration > 1 && absorb_photon != true)
					{
						photon.phi = 2*M_PI*gen_uniform_deviate(a,c,m) - M_PI;
						photon.theta = acos(2*gen_uniform_deviate(a,c,m)-1);
					}
				
				// draw tau value
				tau = inv_p_tau(gen_uniform_deviate(a,c,m));
				
				// if tau is nan, set to a large value
				if (tau!=tau) tau = TAU_MAX;
				
				// take a step in the given orientation if photon is not absorbed (theta,phi)
				if (absorb_photon != true)
				{
					step_size = tau / alpha;
					photon.x += step_size * sin(photon.theta) * cos(photon.phi);
					photon.y += step_size * sin(photon.theta) * sin(photon.phi);
					photon.z += step_size * cos(photon.theta);
				}

				iteration++;
			}

		// escaped in the right direction?
		if (photon.z > 1 && absorb_photon != true) 
		{
			// log escaped photon and bin theta val and costheta val
			num_photons_escaped++;
			theta_bin = (int) ((photon.theta * 180 / M_PI) / dtheta);
			binned_theta_isotropic[theta_bin]++;
			costheta_bin = (int) (cos(photon.theta) / dcostheta);
			binned_cos_theta_isotropic[costheta_bin]++;
			fprintf(pf2,"%lf\n",photon.theta);
			printf("\r%i/%i photons made it through!",num_photons_escaped,total_escaped);
		}

		// log simulated photon
		num_photons_simulated++;
		fflush(stdout);
	}


	// simulation finished
	printf("\nFinished MC simulation of isotropic photon scattering! %i/%i photons made it through!\n",num_photons_escaped,num_photons_simulated);
	
	fclose(pf2); // close angle storing file
	
	// save results to a txt file
	char isotropic_photon_scattering_outfile[] = "isotropic_photon_scattering.txt";
	pf = fopen(isotropic_photon_scattering_outfile, "w+");
	for (int i = 0; i<N_BINS_ISOTROPIC; i++)
	{
		fprintf(pf,"%lf,%i,%lf,%i\n",i*dtheta,binned_theta_isotropic[i],i*dcostheta,binned_cos_theta_isotropic[i]);
	}

	/* q1c Monte carlo photon rayleigh scattering 
	---------------------------------------------- */
	// define new vars and arrays for rayleigh scattering
	/*
	int binned_rayleigh[N_BINS_RAYLEIGH];
	double alphas[2] = {10,0.1};
	double rel_theta, rel_phi, y;
	// reset counters, dtheta and absorb_photon vars
	num_photons_simulated = 0;
	num_photons_escaped = 0;
	total_escaped = 1e5;
	absorb_photon = false;
	dtheta = 90. / N_BINS_RAYLEIGH; // in degrees

	// loop over photon types
	for (int num_photon_type = 0; num_photon_type < 2; num_photon_type++)
	{
		// reset the theta bins and counters 
		for (int i = 0; i < N_BINS_RAYLEIGH; i++) binned_rayleigh[i] = 0;
		num_photons_simulated = 0;
		num_photons_escaped = 0;
		angle_counter = 0;
		
		// open file to store angles
		FILE *pf2;
		char outfile_angles[1024]; 
		snprintf(outfile_angles, sizeof outfile_angles, "rayleigh_scattering_angles_%d.txt", num_photon_type);
		pf2 = fopen(outfile_angles,"w+");

		printf("[INFO] simulating photon type %d:\n",num_photon_type);
		
		// run MC simulation for photon type
		alpha = alphas[num_photon_type];

		while (num_photons_escaped != total_escaped)
		{
			// initialise photon at origin moving straight upwards
			photon.x = 0.;
			photon.y = 0.;
			photon.z = 0.;

			photon.theta = 0.;
			photon.phi = 2*M_PI*gen_uniform_deviate(a,c,m) - M_PI; 
			int iteration = 1;
			
			// scatter photon till it escapes atmosphere or gets absorbed
			while ((photon.z >= 0) && (photon.z <= 1) && (absorb_photon != true))
				{	
					// printf("photon.z %lf\n",photon.z);
					// check if photon is absorbed
					random_number = gen_uniform_deviate(a,c,m);
					if (random_number > albedo) absorb_photon = true;
					
					//else scatter photon (unless first run)
					else if (iteration > 1 && absorb_photon != true)
						{
							// draw rel phi at random 
							rel_phi = 2*M_PI*gen_uniform_deviate(a,c,m) - M_PI; 
							// draw rel theta via rejection method
							rel_theta = acos(2*gen_uniform_deviate(a,c,m) - 1);
							y = (3/8*M_PI)*gen_uniform_deviate(a,c,m);
							while (y>rayleigh_phase_func(rel_theta))
							{
								rel_theta = acos(2*gen_uniform_deviate(a,c,m) - 1);
								y = 2*gen_uniform_deviate(a,c,m);
							}
							
							// transform relative angles back to lab frame
							rotate_back_to_lab(&photon.theta, &photon.phi, rel_theta, rel_phi);
														

						}
					
					// draw tau value
					tau = inv_p_tau(gen_uniform_deviate(a,c,m));
					
					// if tau is nan, set to a large value
					if (tau!=tau) tau = TAU_MAX;
					
					// take a step in the given orientation if photon is not absorbed (theta,phi)
					if (absorb_photon != true)
					{
						// store the angles
						if (angle_counter<N)
							{
								fprintf(pf2,"%lf,%lf\n",photon.theta, photon.phi);
								angle_counter++;
							}
						step_size = tau / alpha;
						photon.x += step_size * sin(photon.theta) * cos(photon.phi);
						photon.y += step_size * sin(photon.theta) * sin(photon.phi);
						photon.z += step_size * cos(photon.theta);
					}

					iteration++;
				}

			// escaped in the right direction?
			if (photon.z > 1 && absorb_photon != true) 
				{
					// log escaped photon and bin theta val
					num_photons_escaped++;
					i_bin = (int) ((photon.theta * 180 / M_PI) ) / dtheta;
					// printf("theta:%lf, bin:%d\n",photon.theta*180/M_PI,i_bin);
					binned_rayleigh[i_bin]++;
				}	
			
			// log simulated photon
			num_photons_simulated++;
			// printf("\r%i/%i photons made it through!\n",num_photons_escaped,total_escaped);
			fflush(stdout);
		}
	
		// save results to a txt file
		char outfile[1024]; 
		snprintf(outfile, sizeof outfile, "rayleigh_scattering_%d.txt", num_photon_type);
		pf = fopen(outfile, "w+");
		for (int i = 0; i<N_BINS_RAYLEIGH; i++) fprintf(pf,"%lf, %i\n",i*dtheta +dtheta/2,binned_rayleigh[i]);
		fclose(pf);
		fclose(pf2);
	}
	*/
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

double rayleigh_phase_func(double rel_theta)
{
	// phase function for rayleigh scattering, 
	// argument: scattering angle of photon
	return (3/16*M_PI) * (1+pow(cos(rel_theta),2));
}

void rotate_back_to_lab(double *theta, double *phi, double rel_theta, double rel_phi)
{	
	// function to transform rel theta and phi angles defined in frame with z axis aligned with photons previous orientation back into lab frame.

	// outgoing photon defined in relative frame
	struct particle rel_photon; 
	rel_photon.x = sin(rel_theta) * cos(rel_phi);
	rel_photon.y = sin(rel_theta) * sin(rel_phi);
	rel_photon.z = cos(rel_theta);

	// define ACW rotation matrices Rz, Ry and compound rotation matrix R and position vectors p_rel, p_lab 
	double Rz[3][3], Ry[3][3], R[3][3];
	double p_rel[3][1] = {{rel_photon.x},{rel_photon.y},{rel_photon.z}}; 
	double p_lab[3][1]; 
	
	// align z axis with z'
	init_Ry(Ry,(*theta)); 
	if (fabs(*theta) > EPS) init_Rz(Rz,(*phi));
	else init_Rz(Rz,0); // if theta = 0, set phi to 0 (rotation about z unecessary if theta=0)

	// Combine rotation matrices and apply
	matmul(3,3,3,Rz,Ry,R);
	matmul(3,3,1,R,p_rel,p_lab);

	// update theta and phi vals for photon
	double r = pow(pow(p_lab[0][0],2)+pow(p_lab[1][0],2)+pow(p_lab[2][0],2),0.5);
	*theta = acos(p_lab[2][0]/r);
	*phi = atan2(p_lab[1][0],p_lab[0][0]); 
}
void init_Rz(double Rz[3][3], double theta)
{
	// define a rotation matrix to rotate points about z axis by theta degrees ACW
	Rz[0][0] = cos(theta);
	Rz[0][1] = -sin(theta);
	Rz[0][2] = 0;
	
	Rz[1][0] = sin(theta);
	Rz[1][1] = cos(theta);
	Rz[1][2] = 0;
	
	Rz[2][0] = 0;
	Rz[2][1] = 0;
	Rz[2][2] = 1;
}

void init_Ry(double Ry[3][3], double theta)
{
	// define a rotation matrix to rotate points about y axis by theta degrees ACW
	Ry[0][0] = cos(theta);
	Ry[0][1] = 0;
	Ry[0][2] = sin(theta);

	Ry[1][0] = 0;
	Ry[1][1] = 1;
	Ry[1][2] = 0;

	Ry[2][0] = -sin(theta);
	Ry[2][1] = 0;
	Ry[2][2] = cos(theta);
}

void matmul(int m, int n, int p, double A[m][n], double B[n][p], double C[m][p])
// Code to compute AB=C matrix multiplications
// code adapted from https://codereview.stackexchange.com/questions/179043/matrix-multiplication-using-functions-in-c
{
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < p; j++) {
			C[i][j] = 0;
			for (int k = 0; k < n; k++) {
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}