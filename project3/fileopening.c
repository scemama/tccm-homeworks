#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

size_t read_Natoms(FILE* input_file) { //it is the input name given by professor to read the input file
    char number[5];
    if (fgets(number, 5, input_file) == NULL) {
        printf("Error reading number of atoms.\n");
        return 0;
    }
    int N;
    sscanf(number, "%d", &N); // finds the number of atoms in the file
    printf("Number of Atoms = %d\n", N);

char cartesian[256]; // just a variable to read the rest of the input file. size cannot be less than 50 i // struct cart one[Natoms]; // the one array that holds every coordinate and atomic number of the file!!!

float** malloc_2d(size_t m, size_t n){
 float** a = malloc(m*sizeof(float*));
 if (a == NULL){
  return NULL;
 }
 a[0] = malloc(n*m*sizeof(float));
 if (a[0] == NULL) {
  free(a);
  return NULL;
 }
 for (size_t i=1;i<m;i++){
	a[i] = a[i-1]+n;
 }
 return a; 
}

float** coord = malloc_2d(N, 3);
float* mass = malloc(5 * N * sizeof(float));
for (int j = 0; j < N; j++) {
        if (fgets(cartesian, sizeof(cartesian), input_file) == NULL) {
            printf("Error reading coordinates and atomic number. Mostly there are more than 4 inputs in a line :/\n");
            break;
        }
        if (sscanf(cartesian, "%f %f %f %f", &coord[j][0], &coord[j][1], &coord[j][2],&mass[j]) != 4) {                                                                                                   printf("Error parsing coordinates and atomic number. Mostly there are more spaces than required :/ we usually ask only for 1 space between the numbers please \n");
            continue;
        }
        printf("%f %f %f with atomic number as %f\n", coord[j][0], coord[j][1], coord[j][2], mass[j]);
    }

printf("\n The distance matrix is given below:\n ");
float** distance = malloc_2d(N,N);
for (int i=0;i<N;i++){
	for (int j=0;j<N;j++){
             distance[i][j] = (coord[i][0]-coord[j][0])*(coord[i][0]-coord[j][0]);
	     distance[i][j] = distance[i][j]+(coord[i][1]-coord[j][1])*(coord[i][1]-coord[j][1]);
	     distance[i][j] = distance[i][j]+(coord[i][2]-coord[j][2])*(coord[i][2]-coord[j][2]);
	     distance[i][j]=sqrt(distance[i][j]);
	}
 printf("%f %f %f\n",distance[i][0],distance[i][1],distance[i][2]);
}





void free_2d(float** a){
free(a[0]);
a[0]=NULL;
free(a);
}
free_2d(coord);
free_2d(distance);
free(mass);
return N;
}

//int main() {
  //  FILE* input_file = fopen("inp.txt", "r"); //important command to open the file as a reader "r"
   // if (input_file == NULL) {
    //    printf("Error opening file.\n");
    //    return 1;
   // }

    //read_Natoms(input_file);

    //fclose(input_file); // important command to close the file
   // return 0;
//}
