#include "utilities.h"
#include <assert.h>
#include <time.h>
#include <fstream>
#include <typeinfo>
# include <omp.h>

extern ZZ *baza;

ZZ sample_r()
{
	ZZ r;
	long ro_bound;
	// ZZ doi( (*baza) );
	ZZ doi(2);
	ZZ doi_la_ro;

	RandomBits(r, Params::getRo() + 1); // r apartine[0, 2^(ro+1) )
	conv(ro_bound, Params::getRo());
	doi_la_ro = power(doi, ro_bound);
	r = r - doi_la_ro; // aducem pe r in intervalul ( (-1)*2^ro, 2^ro)

	return r;
}

ZZ sample_r_for_encryption()
{
	ZZ r(0);
	long ro_bound;
	ZZ doi((*baza) );
	ZZ doi_la_ro;

	RandomBits(r, Params::getRoPrim() + 1); // r apartine[0, 2^(ro_prim+1) )
	conv(ro_bound, Params::getRoPrim());
	doi_la_ro = power(doi, ro_bound);
	r = r - doi_la_ro; // aducem pe r in intervalul ( -2^ro_prim, 2^ro_prim)

	return r;
}

ZZ sample_integer(ZZ p)
{
	ZZ x;
	ZZ q, r;

	long long_p;
	conv(long_p, p);

	RandomBits(q, Params::getGamma()/long_p);
	r = sample_r();

	// x = p*q+r;
	// x = p * q + 2 * r;
	x = p * q + (*baza) * r;

	return x;
}

void generate_keys(ZZ &sk, vector<ZZ> &pk)
{
	assert(Params::getEta() != 0);
	assert(Params::getGamma() != 0);
	assert(Params::getRo() != 0);
	assert(Params::getRoPrim() != 0);
	assert(Params::getTau() != 0);

	ZZ doi(2);
	ZZ doi_la_eta_1 = power(doi, Params::getEta() - 1);

	// generate secret key
	do
	{
		RandomBits(sk, Params::getEta()-1); // esantionam sk in [0, 2^(eta-1) ) 
		sk = sk + doi_la_eta_1; // sk apartine [ 2^(eta-1), 2^eta )

	} while (sk % (*baza) != 1);

	// ofstream out("gen.log", ios::out);

	UL tau = Params::getTau();
    UL i = 0;
    
#if defined(_OPENMP)
    int nt = omp_get_num_threads();
#endif
    

   do
   {
       pk.clear();

       try
       {
            // sample public key x_i's
                
            #if defined(_OPENMP)
            // #pragma omp parallel for shared(pk, sk, tau) private(i) \
                       schedule(static,1) num_threads(nt)
            for (i = 0; i < tau; i++)
            {
                   // ZZ x_i = sample_integer(sk);
                   // pk.push_back(x_i);
                   pk.push_back(sample_integer(sk));
            }
            #else
            for (i = 0; i < tau; i++)
            {
                   // ZZ x_i = sample_integer(sk);
                   // pk.push_back(x_i);
                   pk.push_back(sample_integer(sk));
            }
            #endif
        }
        catch (bad_alloc &ba)
        {
            cout << "bad alloc "<<ba.what()<<endl;
        }

        //out.close();
            

        // pune cel mai mare x pe prima pozitie
        pk = ordonate_vector(pk, sk);
        // x_0 - odd
        // r_p(x_0) - even

   } while ( ( pk[0] % (*baza) != 1 ) && (( pk[0] % sk ) % (*baza) != 0 ) );

	// cout << "pk.size = " << pk.size() << endl;
	// cout << "tau = " << Params::getTau() << endl;
	// assert(pk.size() == Params::getTau());


	/*cout << "ro = " << Params::getRo() << endl;
	cout << "gamma = " << Params::getGamma() << endl;
	cout << "p = " << sk << endl;
	cout << "x_0 = " << pk[0] << endl;
	cout << "r_p(x_0) = " << pk[0] % sk << endl;
	for (int i = 1; i < pk.size(); i++)
	{
		cout << "x_" << i << " = " << pk[i] << endl;
	}*/

}

vector<ZZ> ordonate_vector(vector<ZZ> vec, ZZ p)
{
	// cout << "vec.size = " << vec.size() << endl;
	// assert(vec.size() == Params::getTau());
	assert(vec.size() != 0);

	ZZ max = ZZ(-1); // set max egal cu -1, toti x_i sunt pozitivi

    for (int i = 0; i < vec.size(); i++)
	{
		if (vec[i] > max)
			max = vec[i];
	}

	vector<ZZ> ordered_vector;
	ordered_vector.push_back(max);

	bool first = true;

	for (int i = 0; i < vec.size(); i++)
	{
		if (vec[i] == max)
		{
			if (first == true)
			{
				first = false;
			}
			else
			{
				ZZ x_i = sample_integer(p);
				if (x_i >= max)
				{
					ordered_vector.push_back(x_i - max);
				}
				else
				{
					ordered_vector.push_back(x_i);
				}
			}
			
		}
		else
		{
			ordered_vector.push_back(vec[i]);
		}
			
	}

	// assert(ordered_vector.size() == Params::getTau());

	return ordered_vector;
}

ZZ encrypt_integer(vector<ZZ> pk, ZZ m)
{
	assert(pk.size() == Params::getTau());

	vector<int> random_subset = choose_random_subset();
	ZZ r = sample_r_for_encryption();
	ZZ c(0);
	
	// cout << "x0 = " << pk[0] << endl << endl;
	LL tau = Params::getTau();
    LL ran_size = random_subset.size();
    LL i;

    for (i = 0; i < ran_size; i++)
	{
		// cout << pk[random_subset[i]] << "\n";
		// assert(random_subset[i] < pk.size());
		if (random_subset[i] >=  tau)
		{
			cout << " index eronat = " << random_subset[i] << endl;
		}
		c += pk[random_subset[i]];
	}

	// c = c*2+2*r+m;
	c = c + (*baza) * r + m;

	c = c % pk[0]; // c = [ 2*SparseSubsetSum + 2*r + m] modulo x_0

	return c;
}

ZZ decrypt_ciphertext(ZZ ctxt, ZZ sk)
{
	ZZ miu;
	miu = ( ctxt % sk ) % (*baza);
	return miu;
}

vector<int> choose_random_subset()
{
	//pentru eficienta
	static bool generated = false;
	static vector<int> S;

	if (generated == true)
	{
		return S;
	}
	
	int resample = 0;
	int coeff = 1;

	int nr_de_coeff = rand() % (Params::getTau()-1) + 1; // nr apartine [1, tau]

	srand(time(NULL));
	for (int i = 0; i < nr_de_coeff; i++)
	{
		// coeff = rand() % (Params::getTau() - 1) + 1;
		do
		{	
			coeff = rand() % ( Params::getTau() - 1 ) + 1; // nr apartine [1, tau]
			resample = 0;
			for (int j = 0; j < S.size(); j++)
			{
				if (S[j] == coeff)
				{
					resample = 1;
					break;
				}
			}
		} while (resample == 1);
		
		S.push_back(coeff);
	}
	generated = true;

	return S;
}

void test_DGHV_Eval(ZZ sk, vector<ZZ> pk)
{
	ZZ m1(1);
	ZZ m2(2);
	ZZ clear(1);

	ZZ c1 = encrypt_integer(pk, m1);
	ZZ c2 = encrypt_integer(pk, m2);
	
	ZZ c_eval(1);
	ZZ miu;

	cout << "Test pentru inmultire\n";
	for (int i = 0; i < 100; i++)
	{
		c_eval = c_eval*c1*c2;
		clear = clear*m1*m2;

		miu = decrypt_ciphertext(c_eval, sk);

		if (miu != clear)
		{
			cout << "Eroare la decriptare\n";
			cout << "Iteratia i = " << i << endl;
			break;
		}
	}

	cout << "Test pentru adunare\n";
	c_eval = 0;
	clear = 0;
	for (int i = 0; i < 100; i++)
	{
		c_eval = c_eval+c1+c2;
		clear = clear+m1+m2;

		miu = decrypt_ciphertext(c_eval, sk);

		if (miu != clear)
		{
			cout << "Eroare la decriptare\n";
			cout << "Iteratia i = " << i << endl;
			break;
		}
	}

	cout << "Testele pentru operatii s-au incheiat\n";
}

void test_symmentric_DGHV()
{
	// lambda - paramentrul reprezentand securitatea
	int lambda = 80;

	int gamma = pow(lambda, 5);	// gamma = O(lambda^5)
	int eta = pow(lambda, 2);	// eta = O(lamda^2)
	int ro = lambda;
	int ro_prim = 2 * lambda;
	int tau = gamma + lambda;

	Params::set_params(gamma, eta, ro, tau, ro_prim);

	ZZ doi_la_eta_1 = power(ZZ(2), Params::getEta());

	ZZ p;
	ZZ c;
	ZZ q;
	ZZ r;
	int m = 0, miu = 1;

	// generate secret key
	do
	{
		RandomBits(p, Params::getEta() - 1); // esantionam sk in [0, 2^(eta-1) ) 
		p = p + doi_la_eta_1; // sk apartine [ 2^(eta-1), 2^eta )

	} while (p % 2 != 1);

	// encrypt
	long long_p;
	conv(long_p, p);
	RandomBits(q, Params::getGamma() / long_p);
	
	r = sample_r_for_encryption();

	for (int i = 0; i < 100; i++)
	{
		m = rand() % 2;
		c = p*q + 2 * r + m;

		miu = c%p % 2;

		if (miu != m)
		{
			cout << "Eroare la iteratia " << i << endl;
		}
	}

	cout << "Final test schema DGHV simetrica.\n";
}


// I/O
void read_DGHV_keys_from_file(vector<ZZ> &pk, ZZ &sk)
{
	ifstream file("dghv_param.txt", ios::in);

	file >> sk;

	// cout << "READ FC sk = " << sk << endl;

	long pk_size;
	file >> pk_size;
	pk.reserve(pk_size);
	ZZ elem;
	for (int i = 0; i < pk_size; i++)
	{
		file >> elem;
		pk.push_back(elem);
	}

	file.close();
}

void write_DGHV_keys_in_file(vector<ZZ> &pk, ZZ &sk)
{
	ofstream file("dghv_param.txt", ios::out);

	file << sk << endl;

	file << pk.size() << endl;
	for (int i = 0; i < pk.size(); i++)
	{
		file << pk[i] << " ";
	}
}

// operatie = 0 READ
// operatie = 1 WRITE
void test_file_IO_DGHV_params(char *filename, UL* params, vector<ZZ> &pk, ZZ &sk, int operatie)
{
	//if (operatie == 0)
	//{
		write_DGHV_params_in_file(filename, params, pk, sk);

		cout << "lambda = " << params[0] << endl;
		cout << "gamma = " << params[1] << endl;
		cout << "eta = " << params[2] << endl;
		cout << "ro = " << params[3] << endl;
		cout << "tau = " << params[4] << endl;
		cout << "ro_prim = " << params[5] << endl;

		cout << "sk = " << sk << endl;
		cout << "pk.size = " << pk.size() << endl;
		for (int i = 0; i < pk.size(); i++)
		{
			cout << pk[i] << " ";
		}

		cout << "\n\nScriere incheiata.\n";
	//}
	//else
	//{
		pk.clear();

		read_DGHV_params_from_file(filename, params, pk, sk);

		cout << "lambda = " << params[0] << endl;
		cout << "gamma = " << params[1] << endl;
		cout << "eta = " << params[2] << endl;
		cout << "ro = " << params[3] << endl;
		cout << "tau = " << params[4] << endl;
		cout << "ro_prim = " << params[5] << endl;

		cout << "sk = " << sk << endl;
		cout << "pk.size = " << pk.size() << endl;
		for (int i = 0; i < pk.size(); i++)
		{
			cout << pk[i] << " ";
		}

		cout << "\n\nCitire incheiata\n";
	//}
}

// I/O
void read_DGHV_params_from_file(const char *filename, UL* params, vector<ZZ> &pk, ZZ &sk)
{
	// params [0] , [1],   [2], [3],  [4], [5]
	//	  lambda,  gamma,  eta, ro ,  tau, ro_prim
	// ifstream file("dghv_param.txt", ios::in);
	ifstream file(filename, ios::in);

	for (int i = 0; i < 6; i++)
	{
		file >> params[i];
	}

	file >> sk;

	long pk_size;
	file >> pk_size;
	pk.reserve(pk_size);
	ZZ elem;
	for (int i = 0; i < pk_size; i++)
	{
		file >> elem;
		pk.push_back(elem);
	}

	file.close();
}

void write_DGHV_params_in_file(const char *filename, UL* params, vector<ZZ> &pk, ZZ &sk)
{
	// params [0] , [1],   [2], [3],  [4], [5]
	//	  lambda,  gamma,  eta, ro ,  tau, ro_prim
	// ofstream file("dghv_param.txt", ios::out);
	ofstream file(filename, ios::out);

	for (int i = 0; i < 6; i++)
	{
		file << params[i] << endl;
	}

	file << sk << endl;

	file << pk.size() << endl;
	for (int i = 0; i < pk.size(); i++)
	{
		file << pk[i] << " ";
	}

	file.flush();
	file.close();
}




/***************        LARGE SETTINGS						**********/
void generate_keys_LS(ZZ &sk, vector<vector<ZZ> > &pk)
{
	assert(Params::getEta() != 0);
	assert(Params::getGamma() != 0);
	assert(Params::getRo() != 0);
	assert(Params::getRoPrim() != 0);
	assert(Params::getTau() != 0);

	ZZ doi(2);
	ZZ doi_la_eta_1 = power(doi, Params::getEta() - 1);

	// generate secret key
	do
	{
		RandomBits(sk, Params::getEta() - 1); // esantionam sk in [0, 2^(eta-1) ) 
		sk = sk + doi_la_eta_1; // sk apartine [ 2^(eta-1), 2^eta )

	} while (sk % 2 != 1);

	ZZ tau(0);
	UL index_pk = 0;
	UL MAX_SIZE = 1000;
	vector<ZZ> empty;
	pk.push_back(empty);
	do
	{
		pk.clear();

		try
		{
			// sample public key x_i's
			tau = Params::get_tau_LS();
			for (ZZ i(0); i < tau; i++)
			{
				ZZ x_i = sample_integer(sk);
				pk[index_pk].push_back(x_i);
				
				if (index_pk == MAX_SIZE)
				{
					index_pk++;
					vector<ZZ> empty;
					pk.push_back(empty);
				}
			}
		}
		catch (bad_alloc &ba)
		{
			cout << "bad alloc " << ba.what() << endl;
		}

		pk = ordonate_vector_LS(pk, sk, index_pk);

	} while ((pk[0][0] % 2 != 1) && ((pk[0][0] % sk) % 2 != 0));

}

vector<vector<ZZ> > ordonate_vector_LS(vector<vector<ZZ> > vec, ZZ p, UL index_pk)
{
	assert(vec.size() != 0);

	ZZ max = ZZ(-1); // set max egal cu -1, toti x_i sunt pozitivi

	for (int i = 0; i < vec.size(); i++)
	{
		for (int j = 0; j < vec[i].size(); j++)
		{
			if (vec[i][j] > max)
				max = vec[i][j];
		}
	}

	vector<vector<ZZ> > ordered_vector;
	vector<ZZ> empty;
	ordered_vector.push_back(empty);
	ordered_vector[0].push_back(max);

	bool first = true;

	for (int i = 0; i < vec.size(); i++)
	{
		for (int j = 0; j < vec[i].size(); j++)
		{
			if (vec[i][j] == max)
			{
				if (first == true)
				{
					first = false;
				}
				else
				{
					ZZ x_i = sample_integer(p);
					if (x_i >= max)
					{
						ordered_vector[i].push_back(x_i - max);
					}
					else
					{
						ordered_vector[i].push_back(x_i);
					}
				}

			}
			else
			{
				ordered_vector[i].push_back(vec[i][j]);
			}
		}
		vector<ZZ> empty;
		ordered_vector.push_back(empty);

	}

	return ordered_vector;
}

ZZ encrypt_integer_LS(vector<vector<ZZ> > pk, ZZ m)
{
	assert(pk.size() == Params::getTau());

	//vector<int> random_subset = choose_random_subset();
	//ZZ r = sample_r_for_encryption();
	ZZ c(0);

	/*for (int i = 0; i < random_subset.size(); i++)
	{
		if (random_subset[i] >= Params::getTau())
		{
			cout << " index eronat = " << random_subset[i] << endl;
		}
		c += pk[random_subset[i]];
	}

	// c = c*2+2*r+m;
	c = c + 2 * r + m;

	c = c % pk[0]; // c = [ 2*SparseSubsetSum + 2*r + m] modulo x_0
	*/

	return c;
}