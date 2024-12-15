#include <trexio.h>
#include <stdio.h>
#include <stdlib.h>

int main(){
	
	//OPENING FILE

	trexio_exit_code rc;
	trexio_t* trexio_file = trexio_open("./data/h2o.h5",'r', TREXIO_AUTO, &rc);
	if (rc != TREXIO_SUCCESS) {
		printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
		exit(1);
	}


	//READING NUCLEUS REPULSIONS

	trexio_exit_code trexio_read_nucleus_repulsion(trexio_t* const trexio_file, 								double* const energy);
	double energy; // Variable where the energy is read
	rc = trexio_read_nucleus_repulsion(trexio_file, &energy);
	// Check the return code to be sure reading was OK
	if (rc != TREXIO_SUCCESS) {
  		printf("TREXIO Error reading nuclear repulsion energy:\n%s\n",
        	trexio_string_of_error(rc));
		exit(1);
	 }
	

	//READING NUMBER OF OCCUPIED MOLECULAR ORBITALS

	trexio_exit_code trexio_read_electron_up_num(trexio_t* const trexio_file,
                                             int32_t* const n_up);
	int n_up; //variable where number of occupied orbitals is read
	rc = trexio_read_electron_up_num(trexio_file, &n_up);
        // Check the return code to be sure reading was OK
        if (rc != TREXIO_SUCCESS) {
        	printf("TREXIO Error reading number of up-spin electrons:\n%s\n",
        	trexio_string_of_error(rc));
        	exit(1);
         }

	//READING NUMBER OF MOLECULAR ORBITALS

	trexio_exit_code trexio_read_mo_num(trexio_t* const trexio_file,
                                             int32_t* const mo_num);
        int mo_num; //variable where number of molecular orbitals is read
        rc = trexio_read_mo_num(trexio_file, &mo_num);
        if (rc != TREXIO_SUCCESS) {
        	printf("TREXIO Error reading numer of molecular orbitals:\n%s\n",
        	trexio_string_of_error(rc));
        	exit(1);
         }


	//READING ONE ELECTRON INTEGRALS
	

	 trexio_exit_code trexio_read_mo_1e_int_core_hamiltonian(trexio_t* const trexio_file, double* const mo_1e_int_core_hamiltonian);

	 double** data = malloc(mo_num * sizeof(double*));
         for (int i=1 ; i<mo_num ; i++) {
                data[i] = (double*)malloc(mo_num * sizeof(double));
         }

         rc = trexio_read_mo_1e_int_core_hamiltonian(trexio_file, *data);
         if (rc != TREXIO_SUCCESS) {
        	printf("TREXIO Error reading one electron integrals:\n%s\n",
        	trexio_string_of_error(rc));
        	exit(1);
         } 

	

	printf("%lf\n",energy);
	printf("%d\n",n_up);
	printf("%d\n",mo_num);
	
	
	for (int i = 0; i < 8; i++) {            // Loop through rows
        	for (int j = 0; j < 8; j++) {        // Loop through columns
            		printf("%f ", data[i][j]);
		
		printf("\n");
        	}
	}
		
	for (int i = 0; i < 8; i++) {
        	free(data[i]);  // Free each row
    	}
    	free(data); 

	
	rc = trexio_close(trexio_file);
	if (rc != TREXIO_SUCCESS) {
  	printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
	exit(1);
	 }
	trexio_file = NULL;

}


