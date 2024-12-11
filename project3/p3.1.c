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


int main(){
	FILE* file = fopen("inp.txt", "r");
	if (file == NULL){
		perror("Error opening the file");
		return 1;
	}

	size_t Natoms = read_Natoms(file);

	double** malloc_2d(size_t m, size_t n) {
		double** distance = malloc(m*sizeof(double*));
		if (distance == NULL) {
		return NULL;
}
	distance[0] = malloc(n*m*sizeof(double));
	if (distance[0] == NULL) {
		free(distance);
	return NULL;
}

	for (size_t i=1 ; i<m ; i++) {
		distance[i] = distance[i-1]+n;
}
	return distance;
}


void free_2d(double** distance) {
	free(distance[0]);
	distance[0] = NULL;
	free(distance);
}

double** coord = (double**)malloc(Natoms * sizeof(double*));
    for (size_t i = 0; i < Natoms; i++) {
        coord[i] = (double*)malloc(3 * sizeof(double)); // 3 coordenadas por átomo
    }

double* mass = (double*)malloc(Natoms * sizeof(double));

double** distance = (double**)malloc(Natoms * sizeof(double*));
    for (size_t i = 0; i < Natoms; i++) {
        distance[i] = (double*)malloc(Natoms * sizeof(double));
    }

read_molecule(file, Natoms, coord, mass);

compute_distances(Natoms, coord, distance);

printf("Nat = %u\n", Natoms);
for (size_t i = 0; i < Natoms; i++) {
       for (size_t j = 0; j < Natoms; j++) {
           printf("%lf ", distance[i][j]);
       }
       printf("\n");
}
}




