#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Definition of the functions
size_t read_Natoms(FILE* file){   // Reads the number of atoms from the input
	size_t Natoms;
	fscanf(file, "%zu", &Natoms); // %zu = Read the int without the sign
	return Natoms;
}

void read_molecule(FILE* file, // This functions reads the coordinates (x, y, z) and the masses
	size_t Natoms,
	double** coord,
	double* mass){

	for (size_t i = 0; i < Natoms; i++){
		fscanf(file, "%lf %lf %lf %lf",
				&coord[i][0], &coord[i][1], &coord[i][2], &mass[i]);
	}
}

void compute_distances(size_t Natoms, // This function calculates the distance between atoms, calculates the modulus of a vector
	double** coord,
	double** distance){

	for (size_t i = 0; i < Natoms; i++){
		for (size_t j = 0; j < Natoms; j++){
			double dx = coord[i][0] - coord[j][0];
			double dy = coord[i][1] - coord[j][1];
			double dz = coord[i][2] - coord[j][2];
			distance[i][j] = sqrt(dx*dx + dy*dy + dz*dz);
		}
	}
}

double** malloc_2d(size_t m, size_t n) { // Allocating function provide by the teachers
	double** a = malloc(m*sizeof(double*));
	if (a == NULL) {
		perror("Error allocating the memory");
	return NULL;
}
	a[0] = malloc(n*m*sizeof(double));
	if (a[0] == NULL) {
		free(a);
	return NULL;
}
	for (size_t i=1 ; i<m ; i++) {
		a[i] = a[i-1]+n;
	}
	return a;
}

void free_2d(double** a) { // Freeing function provide by the teachers
	free(a[0]);
	a[0] = NULL;
	free(a);
}

double V(double epsilon, // Lennard-Jones potential energy function
	double sigma,
	size_t Natoms,
	double** distance){
	// we have splited the equation y 4 terms
	double sigr; // calculates the value of sigma / r of every atom
	double t1; // t1 = (sigma/r)**12
	double t2; // t2 = (sigma/r)**6
	double etemp; // saves the operation 4*epsilon*(t1-t2) 
	double etot; // sums up all the potential energies of each iteration

	for (size_t i = 0; i < Natoms; i++){
		for (size_t j = 0; j < Natoms ; j++){
			if (distance[i][j] < 0.00001){;
				continue;
		}
			else{
				sigr = sigma / distance[i][j];
				t1 = pow(sigr, 12);
				t2 = pow(sigr, 6);
				etemp=4*epsilon*(t1-t2);
			//	printf("Etemp = %lf\n", etemp); Used to track if the loop is correct
				etot += etemp;
			//	printf("Etot = %lf\n", etot); Used to track if the loop is correct
			}
		}
	}
return etot;
}

double T(size_t Natoms,
	double** velocity,
	double* mass){

	double ektemp;
	double ektot;
	double t3;
	double t4;
	double veltemp;

	for (size_t i = 0; i < Natoms; i++){
		if (velocity[i] < 0.00001){
			continue;
		}
		else{
			veltemp = velocity[i];
			t3 = pow(veltemp,2);
		double	t4 = mass[i] * t3;
			ektemp = 1/2 * t4;
			printf("EKtemp = %lf\n", ektemp);
			ektot += ektemp;
			printf("EKtot = %lf\n", ektot);
		}
	}
}





int main(){
	FILE* file = fopen("inp.txt", "r"); //Opens the file
	if (file == NULL){
		perror("Error opening the file");
		return 1;
	}

	size_t Natoms = read_Natoms(file); // Read the number of atoms with the previous function
	
	double* mass = (double*)malloc(Natoms * sizeof(double)); // Here we allocate memory for mass, coordinates and distance variables
	double** coord = malloc_2d(Natoms, 3);
	double** distance = malloc_2d(Natoms,3);

	read_molecule(file, Natoms, coord, mass); // Obtaining the variables
	
	compute_distances(Natoms, coord, distance);

	printf("Nat = %zu\n", Natoms); // Printing Natoms and distance to check if they are correct
	for (size_t i = 0; i < Natoms; i++) {
	       for (size_t j = 0; j < Natoms; j++) {
        	   printf("%lf ", distance[i][j]);
	       }
	}


	double epsilon = 0.0661; //j/mol 
	double sigma = 0.3345; //nm
	double etot = V(epsilon, sigma, Natoms, distance);
	printf("V = %lf\n", etot);
	
	double** velocity = malloc_2d(3, Natoms);
	double ektot = T(Natoms, mass, velocity);
	printf("T = %lf\n", ektot);





















free_2d(coord); // Freeing the memory where were the variables
free_2d(distance);
free(mass);
free_2d(velocity);

}
