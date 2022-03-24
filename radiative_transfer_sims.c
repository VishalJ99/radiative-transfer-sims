/*
Question 1a: Random numbers from a non-flat distribution (40 marks)
For the radiative transfer problem there will be a number of non-flat probability
distributions to draw from when we take a Monte Carlo approach, so let us start
with developing a method to draw random numbers from a non-flat distribution.

Write code to reproduce a distribution of τ values obeying
P(τ)dτ = exp(−τ)dτ (1)

Use two methods: one ‘wasteful’ method (e.g. rejection or weighing method), and
one efficient method (e.g. direct mapping or cumulative function approach). As a
starting point, use the minimal standard generator we have discussed in class (and
of which you can find examples in the notes and class exercises). Do NOT use an
external random number library!

Show that your random draws are indeed approaching the distribution P(τ) by
using a figure illustrating one of your two approaches.

Demonstrate to me through your code output how your efficient method is indeed
more efficient, using a figure or numerical measure or combination of both of your
choice (as long as the point is made clear). If you wish to use an efficient method
different from direct mapping or the cumulative function approach, that is also fine,
as long as you can demonstrate its improvement over the inefficient methods. You
do not need to use C to produce your plots or any other visual, it is sufficient
to have the C program produce the list of values that is subsequently used by a
plotting program.
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// define funcs
void main();
long long int lcg_random_number(int, int, long int);
double p_tau(double);
double inv_p_tau(double);

// set global vars
// params for the lcg random number func
#define global_seed 42  
int a = 16807;
int c = 0;
long int m = 2147483647;
void main()
{
 // Generating P(τ)dτ = exp(−τ)dτ distribution
 
 // Rejection method:
 // -----------------------

 int N = 100000; 
 int x0 = 0; 
 int x1 = 100;
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
  uniform_deviate = (double) lcg_random_number(a,c,m) / m;
  X[i] = (uniform_deviate * norm_factor) +  x0;
 }

 // Generate N uniform deviates y_i in the interval (0, pmax)
 for (int i = 0; i<N; i++)
 {
  uniform_deviate = (double) lcg_random_number(a,c,m) / m;
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
 printf("%i/%i samples rejected",num_rejected,N);
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
 char direct_map_outfile[50] = "direct_mapping_method_deviates.txt";
 pf = fopen(direct_map_outfile, "w+");
 for (int i = 0; i<j; i++)
 {
  fprintf(pf,"%lf\n",direct_mapped_X[i]);
 }

 printf("Finished sampling distribution via direct mapping method\n");

}

long long int lcg_random_number(int a, int c, long int m)
{
 // Linear congruential generator 
 static long long int seed = global_seed;
 seed = (a * seed + c) % m;
 return seed;
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