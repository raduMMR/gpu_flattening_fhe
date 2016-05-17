#pragma once
#include <iostream>
#include <vector>
#include <CL/cl.h>
#include <NTL/ZZ.h>
#include "macros.h"

using namespace std;
using namespace NTL;

/*
@brief this class implements Flattening over Integers
	using a radix > 2 for bit decomposition

	the user must ensure that the value for radix
	is set appropriately with the security parameter
	so that it won't be overflows due to limited
	int representation

	this class enhance homomorphic operations using 
	gpu processing
*/
class Flat_GPU
{
	/*
	w_int - baza in care se face descompunerea
	*/
	int w;

	/*
	@field vectorul cheie secreta v = Powersof2(1)
	*/
	int *v;

	/*
	@field l = log x_0 + 1, numarul de biti necesari reprezentarii
	unui ciphertext al schemei DGHV
	dimensiunea unui ciphertext
	*/
	int l;

	/*
	cheia secreta pentru schema DGHV, care va fi folosita si la decriptarea
	schemei Flat_DGHV
	*/
	ZZ sk_DGHV;

	/*
	vectorul cheie publica pentru schema DGHV, pk[0] = x_0
	toate valorile dupa bitdecomp_1 se reduc modulo x_0
	*/
	vector<ZZ> pk_DGHV;

	/*
	matrice de dimensiune l x 1 care va fi completa la criptarea
	unui mesaj nou pe fiecare linie cu
	o criptare noua a lui 0, adica C_prim = (Enc(0)[0],
	Enc(0)[1],
	Enc(0)[2],
	...,
	Enc(0)[l-1] )
	*/
	vector<ZZ> C_prim;


	struct GPU_context
	{
		cl_context context;
		cl_command_queue command_queue;
		cl_program program;
		cl_kernel kernel;
		cl_int ret;
		cl_event ev;
		const char *source;
		size_t global_work_size[3];
		size_t local_work_size[3];
		char program_build_log[1024];
		cl_mem memory_obj1, memory_obj2, memory_result;
		int matrix_size;
		float serial_time;
		int kernel_local_chunck_size;
	};

	/*
	@brief GPU context
	*/
	GPU_context c;

	/*
	device used for GPU processing of homomorphic operations
	*/
	cl_device_id device;

	/*
	metodele private care vor fi folosite pentru tehnicile implementate de schema
	*/

	/*
	@brief metoda in care se calculeaza parametrii apartinand
	schemei DGHV : x_0, sk_DGHV, l - dimensiunea unei matrici ciphertext
	aceste valori depind se parametrul de securitate lambda
	*/
	void	compute_DGHV_settings(char *filename, int lambda);

	void	compute_DGHV_settings(int lambda);

	void	compute_FDGHV_settings();

	ZZ		encrypt_DGHV(int message)const;

	ZZ		decrypt_DGHV(ZZ &ctxt)const;

	void	gpu_bitdecomp(ZZ C_i, int* &result)const;					// BitDecomp(a) = {a_0, a_1, ..., a_n} , a = a_0 + 2*a_1+... +2^(n)*n

	void 	gpu_bitdecomp_1(int* &C_i, ZZ result)const;		// BitDecomp(a_0, a_1, ..., a_n) = a , cu a = a_0 + 2*a_1+... +2^(n)*n

	void	gpu_flatten(int* &C, int dim, int* &result)const;			// Flatten(C)=BitDecomp( BitDecomp_1(C) % x_0 )

	int		set_gpu_context();
	void	cleanup_gpu_context();

public:

	/*
	@brief constructor
	cu acest constructor se genereaza parametri si chei noi
	iar aceste valori nu se salveaza intr-un fisier
	*/
	Flat_GPU(int lambda, int baza = 1024);

	/*
	@brief criptare a unui mesaj intreg cu schema Flat_DGHV
	@param message - mesajul de criptat
	@return matricea C - ciphertext-ul rezultat in urma criptarii
	*/
	int*	encrypt(int message)const;

	/*
	@brief decriptare cu schema Flat_DGHV
	@param ref in C - ciphertext-ul care va fi decriptat
	@return valoarea intreaga obtinuta in urma decriptarii
	*/
	int		decrypt(int* &C)const;

	void matrix_multiply(const char *kernel_name, int* matrix1, int* matrix2,
		int* &result);

	/*
	@brief operatii homomorfice cu GPU
	*/
	
	void	gpu_hom_mult(int* &C1, int* &C2, int*& produs);

	/*
	de testat daca se merita sa fie implementata pe GPU
	*/
	void	gpu_hom_add(int* &matrix1, int* &matrix2, int*& result);

	~Flat_GPU()
	{
		cleanup_gpu_context();

		delete[] v;
	}
};

