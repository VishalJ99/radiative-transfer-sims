#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
// fix the angle for known case of theta = 0, phi = 0
// decide how angles are updated
// experiment with case theta = 0, phi = 0, check if that makes sense
void rotate_back_to_lab(double* , double*, double, double);
// TEST known cases
// 
// all zeros: new theta:-1.#IND00 new phi:0.000000
// rel_theta rel_phi = 0 returns old orientation, fine but phi is messed
// theta = 0 returns rel_theta rel_phi
// theta = 90 phi = 45, rel theta = 90 
// tests : phi = anything, theta = 0 : result for new theta and phi should be rel theta and  rel phi
// tests : theta = 54.7356 , phi = 45, rel theta = 90 - 54.7356, rel_phi = 45 
typedef struct particle
{
   double x,y,z;
   double theta,phi;
} particle;

void init_Rz(double Rz[3][3], double theta)
{
	 // define a rotation matrix to rotate points about z axis by theta degrees ACW
		// printf("\nsin theta val %lf\n",sin(theta));
		
		Rz[0][0] = cos(theta);
		Rz[0][1] = -sin(theta);
		Rz[0][2] = 0;
		
		Rz[1][0] = sin(theta);
		Rz[1][1] = cos(theta);
		Rz[1][2] = 0;
		
		Rz[2][0] = 0;
		Rz[2][1] = 0;
		Rz[2][2] = 1;
		
		// Rz[0][0] = 0;
		// Rz[0][1] = 1;
		// Rz[0][2] = 2;
		
		// Rz[1][0] = 3;
		// Rz[1][1] = 4;
		// Rz[1][2] = 5;
		
		// Rz[2][0] = 6;
		// Rz[2][1] = 7;
		// Rz[2][2] = 8;
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


void rotate_back_to_lab(double *theta, double *phi, double rel_theta, double rel_phi)
{
	#define EPS 1.084202e-19
	// function to transform rel theta and phi angles defined in frame with z axis aligned to photons previous orientation back into lab frame.

	// outgoing photon defined in relative frame
	struct particle rel_photon; 
	rel_photon.x = sin(rel_theta) * cos(rel_phi);
	rel_photon.y = sin(rel_theta) * sin(rel_phi);
	rel_photon.z = cos(rel_theta);

	// define ACW rotation matrices and position vector of photon in relative frame
	double Rz[3][3], Ry[3][3], R[3][3], r;
	double p[3][1] = {{rel_photon.x},{rel_photon.y},{rel_photon.z}};
	double pnew[3][1];
	
	printf("\n");
	printf("P before\n");
	for (int i = 0;i<3;i++){for (int j=0;j<1;j++){printf("%lf ",p[i][j]);} printf("\n");}
	// align z axis with z'
	init_Ry(Ry,(*theta)); 
	// if (fabs(*theta) > EPS) init_Rz(Rz,(*phi));
	// else init_Rz(Rz,0); // Rotation about z axis irrelevant since photon already in XZ plane.
	init_Rz(Rz,(*phi));
	// Combine rotation matrices and apply
	matmul(3,3,3,Rz,Ry,R);
	matmul(3,3,1,R,p,pnew);

	// Convert cartesian coordinates to spherical polar
	r = pow(pow(pnew[0][0],2)+pow(pnew[1][0],2)+pow(pnew[2][0],2),0.5);
	*theta = acos(pnew[2][0]/r);
	*phi = atan2(pnew[1][0],pnew[0][0]); // check this func vs arctan

	// update photons angles
	// check how to do this
  printf("Ry\n");
  for (int i = 0;i<3;i++){for (int j=0;j<3;j++){printf("%lf ",Ry[i][j]);} printf("\n");}
  printf("\n");
  printf("Rz\n");
  for (int i = 0;i<3;i++){for (int j=0;j<3;j++){printf("%lf ",Rz[i][j]);} printf("\n");}
  printf("\n");
  printf("R\n");
  for (int i = 0;i<3;i++){for (int j=0;j<3;j++){printf("%lf ",R[i][j]);} printf("\n");}
			printf("\n");
  printf("P after\n");
  for (int i = 0;i<3;i++){for (int j=0;j<1;j++){printf("%lf ",pnew[i][j]);} printf("\n");}
  // printf("theta:%lf\n",(*theta)*180/M_PI);

}
void main()
{
  // double;
  double rel_theta, rel_phi;
  struct particle photon = {.x = 0, .y = 0, .z = 0, .theta = 0 * (M_PI/180), .phi = 45 * (M_PI/180)};
  
  rel_theta = 54.7356 * (M_PI/180);
  rel_phi = 45 * (M_PI/180);
  
  printf("old theta:%lf old phi:%lf \n",photon.theta*180/M_PI,photon.phi*180/M_PI);
  printf("rel theta:%lf rel phi:%lf \n",rel_theta*180/M_PI,rel_phi*180/M_PI);
  rotate_back_to_lab(&photon.theta, &photon.phi, rel_theta, rel_phi);
  printf("\nnew theta:%lf new phi:%lf \n",photon.theta*180/M_PI,photon.phi*180/M_PI);


}