#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

struct cart { // a data structure to define the coordinates and atomic number
    float x;
    float y;
    float z;
    float at;
};

size_t read_Natoms(FILE* input_file) { //it is the input name given by professor to read the input file
    char number[5];
    if (fgets(number, 5, input_file) == NULL) {
        printf("Error reading number of atoms.\n");
        return 0;
    }
    int N;
    sscanf(number, "%d", &N); // finds the number of atoms in the file
    printf("Number of Atoms = %d\n", N);

    char cartesian[256]; // just a variable to read the rest of the input file. size cannot be less than 50 i guess
    struct cart one[N]; // the one array that holds every coordinate and atomic number of the file!!!

    for (int j = 0; j < N; j++) {
        if (fgets(cartesian, sizeof(cartesian), input_file) == NULL) {
            printf("Error reading coordinates and atomic number. Mostly there are more than 4 inputs in a line :/\n");
            break;
        }
        if (sscanf(cartesian, "%f %f %f %f", &one[j].x, &one[j].y, &one[j].z, &one[j].at) != 4) {
            printf("Error parsing coordinates and atomic number. Mostly there are more spaces than required :/ we usually ask only for 1 space between the numbers please \n");
            continue;
        }
        printf("%f %f %f with atomic number as %f\n", one[j].x, one[j].y, one[j].z, one[j].at);
    }

    return N;
}

int main() {
    FILE* input_file = fopen("inp.txt", "r"); //important command to open the file as a reader "r"
    if (input_file == NULL) {
        printf("Error opening file.\n");
        return 1;
    }

    read_Natoms(input_file);

    fclose(input_file); // important command to close the file
    return 0;
}
