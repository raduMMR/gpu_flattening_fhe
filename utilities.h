#pragma once

#include<NTL/ZZ.h>
#include <vector>
#include <iostream>
#include "Params.h"
using namespace NTL;
using namespace std;

/*
@brief esantioneaza eroarea r
*/
ZZ sample_r();
ZZ sample_r_for_encryption();

/*
@brief samples a big integer
	according to gamma and ro parameters
@return the sampled integer
*/
ZZ sample_integer(ZZ p);

/*
@brief genereaza cheia secreta si cheia publica
@param params - contextul cheii, parametrii gamma, ro, tau, eta, ro_prim
@param out sk = cheia secreta
@param out pk = vectorul cheii public <x_0, x_1, ..., x_tau>
*/
void generate_keys(ZZ &sk, vector<ZZ> &pk);

/*
@brief gaseste maximul dintr-un vector de intregi mari si il pune pe prima pozitie
@param vec vectorul din care se extrage maximul
@return vectorul structurat cu max pe prima pozitie
*/
vector<ZZ> ordonate_vector(vector<ZZ> vec, ZZ p);

/*
@brief cripteaza un intreg
@param sk - cheia publica
@param integer - mesajul de criptat
@return ciphertextul
*/
ZZ encrypt_integer(vector<ZZ> pk, ZZ integer);

/*
@brief decripteaza ciphertextul
@param ctxt ciphertextul din care se va recupera mesajul
@parama sk cheia secreta
@return mesajul intreg recuperat
*/
ZZ decrypt_ciphertext(ZZ ctxt, ZZ sk);

/*
@brief genereaza aleator o submultime din {1,...,tau}
*/
vector<int> choose_random_subset();

/*
@brief procedura de test pentru schema DGHV simetrica
*/
void test_symmentric_DGHV();



/*
@brief citire/scriere chei secrete si publice in si din fisier
*/
void read_DGHV_keys_from_file(vector<ZZ> &pk, ZZ &sk);
void write_DGHV_keys_in_file(vector<ZZ> &pk, ZZ &sk);

void test_file_IO_DGHV_params(char *filename, UL* params, vector<ZZ> &pk, ZZ &sk, int operatie);

/*
@brief citeste parametrii schemei DGHV din fisier
@param filename - numele fisierului sursa
@param params - parametrii schemei DGHV - lambda, gamma, eta, ro, tau, ro_prim
@params pk - cheia publica DGHV, vectorul de intregi x_i, cu x_0 cel mai mare dintre ei
@param sk - cheia secreta DGHV, intregul p
*/
void read_DGHV_params_from_file(const char *filename, UL* params, vector<ZZ> &pk, ZZ &sk);
void write_DGHV_params_in_file(const char *filename, UL* params, vector<ZZ> &pk, ZZ &sk);




/*****************            LARGE SETTINGS                      ****************/

vector<vector<ZZ> > ordonate_vector_LS(vector<vector<ZZ> > vec, ZZ p, UL index_pk);

void generate_keys_LS(ZZ &sk, vector<vector<ZZ> > &pk);

ZZ encrypt_integer_LS(vector<vector<ZZ> > pk, ZZ m);

