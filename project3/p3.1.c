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



void compute_distances(size_t Natoms, // This function calculates the distance between atoms
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

double** malloc_2d(size_t m, size_t n) {
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

void free_2d(double** a) {
	free(a[0]);
	a[0] = NULL;
	free(a);
}



int main(){
	FILE* file = fopen("inp.txt", "r");
	if (file == NULL){
		perror("Error opening the file");
		return 1;
	}

	size_t Natoms = read_Natoms(file);
	
	double* mass = (double*)malloc(Natoms * sizeof(double));
	double** coord = malloc_2d(Natoms, 3);
	double** distance = malloc_2d(Natoms,3);

	read_molecule(file, Natoms, coord, mass);
	
	compute_distances(Natoms, coord, distance);


printf("Nat = %zu\n", Natoms);
for (size_t i = 0; i < Natoms; i++) {
       for (size_t j = 0; j < Natoms; j++) {
           printf("%lf ", distance[i][j]);
       }
       printf("\n");
}


free_2d(coord);
free_2d(distance);
free(mass);

}
