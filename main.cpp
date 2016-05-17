# include <cstdlib>
# include <iostream>
# include <cmath>
# include <ctime>
# include <omp.h>
#include "utilities.h"
#include "Flat_DGHV.h"

// #include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <assert.h>

#include <fstream>
#include <NTL/ZZ.h>

#include "matrix.h"


// TODO : IMPLEMENTAREA SCHEMEI PE BAZA NTL
// DGHV securiy parameter
// dupa Fully homomorphic encryption over the Integers with Shorter Public Key by Coron et.al
// http://publications.uni.lu/bitstream/10993/12396/1/441.pdf
// valoare acceptata pentru x_0 si enc_0 reprezentati ca long
// #define SEC_PARAM 42		// toy 
// #define SEC_PARAM 52		// small
// #define SEC_PARAM 62		// medium 
// #define SEC_PARAM 72		// large

#define SEC_PARAM 10
#define NR_OPERATII 100

using namespace std;

ZZ *baza;

vector<int> compute_random_circuit()
{
    vector<int> circuit(NR_OPERATII);
    
    srand(time(NULL));
    for(int i=0; i<NR_OPERATII; i++)
    {
        circuit[i] = rand() % 2 ;
    }
    
    return circuit;
}

void test_Flat_DGHV_simple()
{
	cout << "Testare Flat_DGHV_simple ...\n\n";
	int lambda = SEC_PARAM;

    // const char *source_file = "../Data/DGHVparams42.txt";
	// Flat_DGHV fdghv(source_file);
	baza = new ZZ(2);
	Flat_DGHV fdghv(lambda, (*baza));

	Mat_ZZ C1;
	Mat_ZZ C2;
	int m1;
	int m2;
    double wtime = 0.0;
    
    cout << "\nTestare multiplicari ...\n\n";

	for (int test = 0; test < 20; test++)
	{
		vector<int> circuit = compute_random_circuit();

        srand(time(NULL));
        
		m1 = rand() % 2;
		C1 = fdghv.encrypt(m1);

#if defined(_OPENMP)
		wtime = omp_get_wtime ( );
#endif
		int i = 0;
		for (i = 0; i < NR_OPERATII; i++)
		{
			m2 = rand() % 2;
			C2 = fdghv.encrypt(m2);

            if( circuit[i] == 1 )
            {
                C1 = fdghv.hom_mult(C1, C2);
                m1 *= m2;
            }
            else
            {
                C1 = fdghv.hom_add(C1, C2);
                m1 = ( m1 + m2 ) % 2;
            }
                 

			if ( fdghv.decrypt(C1) != m1 )
			{

				cout << "Max mult depth FDGHV i = " << i << endl;
				break;
			}

		}
		cout << "Depth = " << i << endl;

#if defined(_OPENMP)
		wtime = omp_get_wtime ( ) - wtime;
#endif
        cout << "  TIME FARA OMP = " << wtime << "\n"; 

		// continue;
        
		m1 = 1;
		C1 = fdghv.encrypt(m1);
		m2 = 1;
		C2 = fdghv.encrypt(m2);

#if defined(_OPENMP)
        wtime = omp_get_wtime ( );
#endif

		for (int i = 0; i < NR_OPERATII; i++)
		{
			// m2 = rand() % 2;
			// C2 = fdghv.encrypt(m2);

			// C1 = fdghv.omp_hom_mult(C1, C2);
            // C1 = fdghv.omp_hom_mult_opt(C1, C2);
            
            if( circuit[i] == 1 )
            {
                C1 = fdghv.omp_hom_mult_opt(C1, C2);
                m1 *= m2;
            }
            else
            {
                C1 = fdghv.omp_hom_add(C1, C2);
                m1 = ( m1 + m2 ) % 2;
            }
                 

			if (fdghv.decrypt(C1) != m1 )
			{
				cout << "Max mult depth FDGHV i = " << i << endl;
				break;
			}
			else
			{
				// m1 = (m1 + m2) % 2;
				// m1 *= m2;
				// cout << "Succes\n";
			}

		}
        
        #if defined(_OPENMP)
		    wtime = omp_get_wtime ( ) - wtime;
        #endif
        cout << "  TIME CU OMP = " << wtime << "\n\n";

		// continue;

	}

	cout << "Final test Flat_DGHV_simple\n";
	delete baza;
	baza = NULL;

}

void test_FDGHV_base_w()
{
	LARGE_INTEGER frequency;
	LARGE_INTEGER start;
	LARGE_INTEGER end;
	double interval;

	cout << "\nTestare Flat_DGHV_simple with base different from 2 ...\n\n";
	int lambda = SEC_PARAM;

	int w_baza = pow(2, 15);

	baza = new ZZ(w_baza);
	Flat_DGHV fdghv(lambda, (*baza) );

	fdghv.create_context();

	int *matrix1 = nullptr, *matrix2 = nullptr;
	int *result = (int *)calloc(MATRIX_SIZE * MATRIX_SIZE, sizeof(int));

	get_matrix_int(&matrix1, MATRIX_SIZE, 0, 9);
	get_matrix_int(&matrix2, MATRIX_SIZE, 0, 9);
	fdghv.matrix_multiply("matrix_multiply_gpu", matrix1, matrix2, result);

	free(matrix1);
	free(matrix2);
	free(result);

	return;

	Mat_ZZ C1;
	Mat_ZZ C2;
	int m1;
	int m2;
	double wtime = 0.0;
	cout << "\nTestare operatii W = "<< *baza <<"  ...\n\n";

	for (int test = 0; test < 20; test++)
	{
		vector<int> circuit = compute_random_circuit();

		srand(time(NULL));

		QueryPerformanceFrequency(&frequency);
		QueryPerformanceCounter(&start);

		m1 = rand() % w_baza;
		C1 = fdghv.encrypt(m1);

		cout << "matrix size = ( " << C1.size() <<","<< C1[0].size() << " )" << endl;
		// cout << "Size_Ctxt = (" << C1.size() << ", " << C1[0].size() << ")" << endl;
		cout << "Elem = " << C1[0][0] << endl;

		break;

		int i;
		for (i = 0; i < NR_OPERATII/3; i++)
		{
			m2 = rand() % w_baza;
			C2 = fdghv.encrypt(m2);

			if (circuit[i] == 1)
			{
				C1 = fdghv.hom_mult_opt(C1, C2);
				m1 = (m1*m2) % w_baza;
			}
			else
			{
				C1 = fdghv.hom_add(C1, C2);
				m1 = (m1 + m2) % w_baza;
			}


			if (fdghv.decrypt(C1) != m1)
			{
				cout << "Max mult depth FDGHV i = " << i << endl << endl;
				break;
			}
		}

		m2 = 15;
		C2 = fdghv.encrypt(m2);

		for (; i < NR_OPERATII/2; i++)
		{
			C1 = fdghv.hom_mult_opt(C1, C2);
			m1 = (m1*m2) % w_baza;

			if (fdghv.decrypt(C1) != m1)
			{
				cout << "Max mult depth FDGHV i = " << i << endl << endl;
				break;
			}
		}

		m2 = rand() % w_baza;
		C2 = fdghv.encrypt(m2);

		for (; i < NR_OPERATII; i++)
		{
			C1 = fdghv.hom_add(C1, C2);
			m1 = (m1+m2) % w_baza;

			if (fdghv.decrypt(C1) != m1)
			{
				cout << "Max mult depth FDGHV i = " << i << endl << endl;
				break;
			}
		}

		cout << "Depth = " << i << endl;

		QueryPerformanceCounter(&end);
		interval = (double)(end.QuadPart - start.QuadPart) / frequency.QuadPart;

		printf("%f\n", interval);

		// break;
	}

	cout << "Final test Flat_DGHV_simple with base different from 2 \n";
	delete baza;
	baza = NULL;
}

int main()
{
	test_FDGHV_base_w();

	// test_Flat_DGHV_simple();

	return 0;
}
