#include <stdio.h>
#include <stdlib.h>
#include <string.h>  // for memset, memcmp
#include <math.h>
#include <CUnit/Automated.h>


int invert_matrixd(double *mat, double *inv, int rows);
void print_matrix_multd(double *a, double *b, int n);
void matrix_multd(double *a, double *b, double *c, int n);
int matrices_equald(double *a, double *b, int n);
void print_matrixd(double *a, int n);



//#define INVERSE_DEBUG 1


#define STDERR(x) fprintf(stderr, x)


/* CUnit stuff */


void basic_test()
{
	int i,j,k,tmp;
	
	int num_tests = 10000;
	//identity matrices
	double i4[16];
	double i8[64];
	double i16[256];
	
	//will hold result of matxxcopy*matxxinv
	double result4[16];
	double result8[64];
	double result16[256];
	
	double mat44[16];
	double mat44copy[16];
	double mat44inv[16];
	
	double mat88[64];
	double mat88copy[64];
	double mat88inv[64];
	
	double mat1616[256];
	double mat1616copy[256];
	double mat1616inv[256];
	
	//initialize identity matrices
	
	for (i=0; i<256; i++) {
		i16[i] = 0;
		if (i < 64)
			i8[i] = 0;
		if (i < 16)
			i4[i] = 0;
	}
	
	
	for (i=0; i<16; i++) {
		i16[i*16+i] = 1;
		if (i < 8)
			i8[i*8+i] = 1;
		if (i < 4)
			i4[i*4+i] = 1;
	}
	
	fprintf(stderr, "Running %d tests . . .\n\n", num_tests);
	//do a num_tests iterations of each type
	for (i=0; i<num_tests; i++) {
		for (j=0; j<16; j++)
			mat44[j] = rand();
		
		for (j=0; j<64; j++)
			mat88[j] = rand();
		
		for (j=0; j<256; j++)
			mat1616[j] = rand();
		
		memcpy(mat44copy, mat44, 16*sizeof(double));
		memcpy(mat88copy, mat88, 64*sizeof(double));
		memcpy(mat1616copy, mat1616, 256*sizeof(double));
		
		
		if (invert_matrixd(mat44, mat44inv, 4)) {
			STDERR("Bad 4x4 matrix.\n");
		} else {
			matrix_multd(mat44copy, mat44inv, result4, 4);
			CU_ASSERT(matrices_equald(i4, result4, 4));
		}
			
		if (invert_matrixd(mat88, mat88inv, 8)) {
			STDERR("Bad 8x8 matrix.\n");
		} else {
			matrix_multd(mat88copy, mat88inv, result8, 8);
			CU_ASSERT(matrices_equald(i8, result8, 8));
		}
			
		if (invert_matrixd(mat1616, mat1616inv, 16)) {
			STDERR("Bad 16x16 matrix.\n");
		} else {
			matrix_multd(mat1616copy, mat1616inv, result16, 16);
			CU_ASSERT(matrices_equald(i16, result16, 16));
		}	
	}
}



int invert_matrixd(double *mat, double *inv, int rows)
{
	int i,j,k,max;
	double tmp;
	
	//initialize inv to identity
	for (i=0; i<rows; i++)
		for(j=0; j<rows; j++)
			inv[i*rows+j] = 0;
	
	for (i=0; i<rows; i++)
		inv[i*rows+i] = 1;



	//inverse algorithm
	for (i=0; i<rows; i++) {

		//always swap in the row with the highest absolute value as the pivot
		max = i;
		tmp = mat[i*rows+i];
		for (j=i; j<rows; j++) {
			if (fabs(mat[j*rows+i]) > tmp) {
				max = j;
				tmp = mat[j*rows+i];
			}
		}
		
		//singular
		if (fabs(mat[max*rows+i]) < 0.000000000000001)
			return 1;
			
		if (max != i) {
			for (k=0; k<rows; k++) {
				tmp = mat[i*rows+k];
				mat[i*rows+k] = mat[max*rows+k];
				mat[max*rows+k] = tmp;
				
				tmp = inv[i*rows+k];
				inv[i*rows+k] = inv[max*rows+k];
				inv[max*rows+k] = tmp;
			}
		}
// 		}
		
		
		tmp = (double)1/mat[i*rows+i];
		for (j=0; j<rows; j++) {
			mat[i*rows+j] *= tmp; // mat[i*rows+j]*tmp;
			inv[i*rows+j] *= tmp; // gf_mul(inv[i*rows+j], tmp);
		}
		
		//should be 1 from above but for round of errors
		//maybe setting it explicitly will make it slightly more accurate
		//certainly can't hurt
		mat[i*rows+i] = 1;
		
		for (j=0; j<rows; j++) {
			if (j == i)
				continue;
			
			tmp = mat[j*rows+i];
			for (k=0; k<rows; k++) {
				inv[j*rows+k] += (-tmp)*inv[i*rows+k];  //gf_mul(tmp, inv[i*rows+k]);
				mat[j*rows+k] += (-tmp)*mat[i*rows+k]; //gf_mul(tmp, mat[i*rows+k]);
			}
		}
	}
	return 0;
}


int matrices_equald(double *a, double *b, int n)
{
	int i;
	for (i=0; i<n*n; i++)
		if (fabs(a[i]-b[i]) > 0.0000000001) {
			fprintf(stderr, "%.16lf\t%.16lf\n\n", a[i], b[i]);
			return 0;
		}
	return 1;
}



void matrix_multd(double *a, double *b, double *c, int n)
{
	int i,j,k;
	double d;

	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			d = 0;
			for(k=0; k<n; k++){
				d += a[n*i + k]*b[n*k + j];
			}
			c[i*n+j] = d;
		}
	}
}





void print_matrix_multd(double *a, double *b, int n)
{
	int i,j,k;
	double d;

	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			d = 0;
			for(k=0; k<n; k++){
				d += a[n*i + k]*b[n*k + j];
			}
			printf(" %.4lf", d);
		}
		printf("\n");
	}
	printf("\n");
}


void print_matrixd(double *a, int n)
{
	int i,j,k;
	double d;

	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			printf(" %.E", a[i*n+j]);			
		}
		printf("\n");
	}
	printf("\n");
}


CU_TestInfo tests[] = {
	{ "inverse_test",			basic_test },
};

CU_SuiteInfo suites[] = {
  { "inverse_tests", NULL, NULL, tests },
  CU_SUITE_INFO_NULL,
};


int main()
{

	if (CUE_SUCCESS != CU_initialize_registry())
	      return CU_get_error();

	CU_ErrorCode error = CU_register_suites(suites);

	if( error != CUE_SUCCESS )
		fprintf(stderr, "wtf!");


	CU_automated_run_tests();
	fprintf(stdout, "%d", CU_get_error());

	CU_cleanup_registry();

	return CU_get_error();

}
