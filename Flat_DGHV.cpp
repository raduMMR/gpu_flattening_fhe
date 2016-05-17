#include "Flat_DGHV.h"
#include "Params.h"
#include "utilities.h"
#include <assert.h>


#define CALL_SAFE(func) \
	do { \
		cl_int res = func; \
		if (res != CL_SUCCESS) \
		{ \
			fprintf(stderr, "Error in File: %s, Function: %s, Line: %d, Error: %d\n", __FILE__, __FUNCTION__, __LINE__, res); \
			exit(EXIT_FAILURE); \
		} \
	} while(0)

#define NOT_NULL_RET(assertion) \
	do { \
		if (assertion == NULL || ret != CL_SUCCESS) \
		{ \
			fprintf(stderr, "Error in File: %s, Function: %s, Line: %d, Error: %d\n", __FILE__, __FUNCTION__, __LINE__, ret); \
			exit(EXIT_FAILURE); \
		} \
	} while(0)

#define NOT_NULL(assertion) \
	do { \
		if (assertion == NULL) \
		{ \
			fprintf(stderr, "Error in File: %s, Function: %s, Line: %d\n", __FILE__, __FUNCTION__, __LINE__); \
			exit(EXIT_FAILURE); \
		} \
	} while(0)

Flat_DGHV::Flat_DGHV(int lambda, ZZ baza)
{
	w = baza;
	compute_DGHV_settings(lambda);
	compute_FDGHV_settings();
}

Flat_DGHV::Flat_DGHV(const char *filename)
{
	UL params[6];

	read_DGHV_params_from_file(filename, params, pk_DGHV, sk_DGHV);

	// Params::set_params(gamma, eta, ro, tau, ro_prim);
	Params::set_params(params[1], params[2], params[3], params[4], params[5]);

	compute_FDGHV_settings();
}

Flat_DGHV::Flat_DGHV(char *filename, int lambda)
{
	// this->lambda = lambda;

	compute_DGHV_settings(filename, lambda);

	compute_FDGHV_settings();
}

void Flat_DGHV::compute_DGHV_settings(char *filename, int lambda)
{
	/*UL gamma = pow(lambda, 4);
	UL eta = lambda*lambda;
	UL ro = lambda;
	UL tau = gamma + lambda;
	UL ro_prim = 2 * lambda;

	Params::set_params(gamma, eta, ro, tau, ro_prim);

	generate_keys(sk_DGHV, pk_DGHV);
	assert(pk_DGHV.size() != 0);*/

	compute_DGHV_settings(lambda);

	UL params[6];
	params[0] = lambda;
	params[1] = Params::getGamma();
	params[2] = Params::getEta();
	params[3] = Params::getRo();
	params[4] = Params::getTau();
	params[5] = Params::getRoPrim();

	write_DGHV_params_in_file(filename, params, pk_DGHV, sk_DGHV);
}

void Flat_DGHV::compute_DGHV_settings(int lambda)
{
	UL gamma = pow(lambda, 4);
	UL eta = lambda*lambda;
	UL ro = lambda;
    cout << "\n\nSECURITATE COMPROMISA pt ca tau = lambda, cf def tau = O(gamma+lambda)\n\n";
    UL tau = lambda;
    // UL tau = gamma + lambda;
	UL ro_prim = 2 * lambda;

	Params::set_params(gamma, eta, ro, tau, ro_prim);

	generate_keys(sk_DGHV, pk_DGHV);
	assert(pk_DGHV.size() != 0);
}

void Flat_DGHV::compute_FDGHV_settings()
{
	l = 0;					// l = log x_0 + 1
	ZZ x_0 = pk_DGHV[0];
	while (x_0 != 0)
	{
		l++;
		x_0 = x_0 / w;
	}
	assert(l != 0);
	l += 1;

	v.reserve(l);			// v = Powersof2(1);
	ZZ two = w;
	ZZ pow_of_w(1);
	for (int i = 0; i < l; i++)
	{
		v.push_back(pow_of_w);
		pow_of_w *= w;
	}

	C_prim.reserve(l);
	for (int i = 0; i < l; i++)
	{
		// o criptare noua pentru fiecare linie a matricii C_prim
		Vec_ZZ empty;
		C_prim.push_back(empty);

		// pentru fiecare linie avem o criptare a lui zero diferita
		C_prim[i].push_back(encrypt_DGHV(0));
	}
}

ZZ	Flat_DGHV::encrypt_DGHV(int message)const
{
	return encrypt_integer(pk_DGHV, ZZ(message));
}

ZZ	Flat_DGHV::decrypt_DGHV(ZZ &ctxt)const
{
	// modifica schema DGHV => modulo t >> 2
	return ctxt % sk_DGHV % w ;
}

Vec_ZZ Flat_DGHV::bitdecomp(Vec_ZZ &C_i)const
{
	Vec_ZZ C_decomp;
	long length = C_i.size() * l;
	C_decomp.reserve(length);

	ZZ elem;
	for (int i = 0, j=-1; i < length; i++)
	{
		if (i % l == 0)
		{
			j++;
			elem = C_i[j];
		}

		C_decomp.push_back(ZZ(elem % w));
		elem = elem / w;
	}

	return C_decomp;
}

Vec_ZZ Flat_DGHV::bitdecomp_1(Vec_ZZ &C_i)const
{
	Vec_ZZ C_decomp_1;
	long length = C_i.size() / l;
	C_decomp_1.reserve(length);

	ZZ pow_of_w(1);
	// ZZ two(2);
	for (LL i = 0, j = -1; i < C_i.size(); i++)
	{
		if (i % l == 0)
		{
			pow_of_w = 1;
			j++;
			C_decomp_1.push_back(ZZ(0));		//C_decomp_1[j] = 0;
		}

		C_decomp_1[j] +=  C_i[i] * pow_of_w;
		pow_of_w *= w;
	}

	return C_decomp_1;
}

Mat_ZZ Flat_DGHV::flatten(Mat_ZZ &C)const
{
	Mat_ZZ Flat_C(l);
	Vec_ZZ bd_1;

	for (LL i = 0; i < C.size(); i++)
	{
		bd_1 = bitdecomp_1(C[i]);
		bd_1[0] = bd_1[0] % pk_DGHV[0];
		Flat_C[i] = bitdecomp(bd_1);
	}

	return Flat_C;
}

Mat_ZZ Flat_DGHV::encrypt(int message)const
{
	Mat_ZZ C(l);

	for (int i = 0; i < l; i++)
	{
		// Vec_ZZ empty;
		// C.push_back(empty);
		Vec_ZZ linie = C_prim[i];
		C[i] = bitdecomp(linie);
		C[i][i] += message;
	}
	flatten(C);

	return C;
}

ZZ Flat_DGHV::decrypt(Mat_ZZ &C)const
{
	ZZ message(0);
	for (int i = 0; i < l; i++)
	{
		message += C[0][i] * v[i];
	}
	return decrypt_DGHV(message);
}

Mat_ZZ Flat_DGHV::hom_add(Mat_ZZ &C1, Mat_ZZ &C2)const
{
	assert(C1.size() == C2.size());
	assert(C1.size() == l);
	assert(C1[0].size() == l);

	Mat_ZZ C_add(l);
	
	for (int i = 0; i < l; i++)
	{
		for (int j = 0; j < l; j++)
		{
			C_add[i].push_back(C1[i][j] + C2[i][j]);
		}
	}

	return flatten(C_add);
}

Mat_ZZ Flat_DGHV::hom_mult(Mat_ZZ &C1, Mat_ZZ &C2)const
{
	assert(C1.size() == C2.size());
	assert(C1.size() == l);
	assert(C1[0].size() == l);

	Mat_ZZ C_mult(l);

	ZZ z;
	for (int i = 0; i < l; i++)
	{
		// Vec_ZZ empty;
		// C_mult.push_back(empty);

		C_mult[i].reserve(l);
		
		for (int j = 0; j < l; j++)
		{
			C_mult[i].push_back(ZZ(0));
			for (int k = 0; k < l; k++)
			{
				C_mult[i][j] += C1[i][k] * C2[k][j];
			}
		}
	}

	return flatten(C_mult);
}

Mat_ZZ Flat_DGHV::hom_mult_opt(Mat_ZZ &C1, Mat_ZZ &C2)const
{
	assert(C1.size() == C2.size());
	assert(C1.size() == l);
	assert(C1[0].size() == l);

	Mat_ZZ C_mult(l);
	Mat_ZZ C1_prim(l);
    
	for (int i = 0; i < l; i++)
	{
        C1_prim[i] = bitdecomp_1(C1[i]);
    }

	for (int i = 0; i < l; i++)
	{
		ZZ elem(0);
		
		for (int j = 0; j < l; j++)
		{
			elem += C2[i][j] * C1_prim[j][0];
		}

		Vec_ZZ linie;
		linie.push_back(elem);
        C_mult[i] = bitdecomp(linie);;
	}

	return (C_mult);
}


/*****************				OMP						*******/
Mat_ZZ Flat_DGHV::omp_hom_add(Mat_ZZ &C1, Mat_ZZ &C2)const
{
#if defined(_OPENMP)

	assert(C1.size() == C2.size());
	assert(C1.size() == l);
	assert(C1[0].size() == l);

	Mat_ZZ C_add(l);
	LL i, j;
    
#pragma omp parallel for  \
    default(none) shared(C_add, C1, C2) private(i)
	for (i = 0; i < l; i++)
	{
		for (j = 0; j < l; j++)
		{
			C_add[i].push_back(C1[i][j] + C2[i][j]);
		}
    }
    
    return flatten(C_add);
    
#else
	return hom_add(C1, C2);
#endif

}

Mat_ZZ Flat_DGHV::omp_hom_mult(Mat_ZZ &C1, Mat_ZZ &C2)const
{
#if defined(_OPENMP)
	assert(C1.size() == C2.size());
	assert(C1.size() == l);
	assert(C1[0].size() == l);

	Mat_ZZ C_mult(l);

	LL i, j, k;
	ZZ elem;

#pragma omp parallel for \
	default(none) shared(C_mult, C1, C2) private(i, j, k, elem) 
		for (i = 0; i<l; i++)
		{
            //#pragma omp parallel for collapse(2) \
	            default(none) shared(C_mult, C1, C2) \
                private(i, j, k, elem) 
			for (j = 0; j<l; j++)
			{
				elem = 0;
				for (k = 0; k<l; k++)
				{
					elem += C1[i][k] * C2[k][j];
                    //C_mult[i][j] += C1[i][k] * C2[k][j];
                }	
                C_mult[i].push_back(elem);			
			}
		}

	return flatten(C_mult);
#else
	return hom_mult(C1, C2);
#endif

	
}

Mat_ZZ Flat_DGHV::omp_hom_mult_opt(Mat_ZZ &C1, Mat_ZZ &C2)const
{
#if defined(_OPENMP)

	assert(C1.size() == C2.size());
	assert(C1.size() == l);
	assert(C1[0].size() == l);

	Mat_ZZ C_mult(l);
    Mat_ZZ C1_prim(l);
	LL i, j, k;
    
	for (i = 0; i < l; i++)
	{
		C1_prim[i] = bitdecomp_1(C1[i]);
        C_mult[i].push_back( ZZ(0) );
	}
    
    //#pragma omp parallel for  \
        default(none) shared(C_mult, C1_prim, C2) private(i,j)
    for (i = 0; i < l; i++)
	{
		for (j = 0; j < l; j++)
		{
            C_mult[i][0] += C2[i][j] * C1_prim[j][0];
		}
	//}
    
    //for (i = 0; i < l; i++)
	//{
        C_mult[i] = bitdecomp(C_mult[i]);
    }

	return C_mult;

#else
	return hom_mult_opt(C1, C2);
#endif
}

Mat_ZZ Flat_DGHV::omp_encrypt(int message)const
{
#if defined(_OPENMP)
	Mat_ZZ C(l);
	LL i;
	Vec_ZZ linie;

#pragma omp parallel for shared(C) private(i, linie)
	for (i = 0; i < l; i++)
	{
        linie = C_prim[i];
		C[i] = bitdecomp(linie);
		C[i][i] += message;
	}

	return flatten(C);

#else
	return encrypt(message);
#endif
}

int Flat_DGHV::omp_decrypt(Mat_ZZ &C)const
{
#if defined(_OPENMP)

	ZZ message(0);
	LL i;
#pragma omp parallel for shared(message) private(i)
	for (i = 0; i < l; i++)
	{
		message += C[0][i] * v[i];
	}

	int val;
	conv(val, decrypt_DGHV(message));
	return val;

#else
	ZZ x = decrypt(C);
	int val;
	conv(val, x);
	return val;
#endif
	
}


/*Vec_ZZ Flat_DGHV::omp_bitdecomp(Vec_ZZ &C_i)const
{
    
#if defined(_OPENMP)
	long length = C_i.size() * l;
	Vec_ZZ C_decomp(length);
	ZZ elem;
	LL i;
    LL j = -1;
        
#pragma omp parallel for shared(C_decomp, C_i, j) private(elem, i)
	for (i = 0; i < length; i++)
	{
		
        #pragma omp critical
        {
            if (i % l == 0)
            {
                j++;
                elem = C_i[j];
            }
            C_decomp[i] = ZZ(elem % 2);
        }
		
		elem = elem / 2;
	}
    
	return C_decomp;
#else
	return bitdecomp(C_i);
#endif

}

Vec_ZZ Flat_DGHV::omp_bitdecomp_1(Vec_ZZ &C_i)const
{
#if defined(_OPENMP)
	long length = C_i.size() / l;
	Vec_ZZ C_decomp_1(length);
	ZZ pow_of_two(1);
	ZZ two(2);
	UL i, j;
#pragma omp parallel for shared(C_decomp_1) private(i,j)
	for (i = 0, j = -1; i < C_i.size(); i++)
	{
		if (i % l == 0)
		{
			pow_of_two = 1;
			j++;
			C_decomp_1[i] = ZZ(0);		//C_decomp_1[j] = 0;
		}

		C_decomp_1[j] += C_i[i] * pow_of_two;
		pow_of_two *= two;
	}

	return C_decomp_1;
#else
	return omp_bitdecomp_1(C_i);
#endif

}*/



/**************************    BATCHING       ***********************/

/*
@brief for batching we considered the Scale Invariance over the Integers scheme
for each Enc_i(0) we used a different scale 
*/
Mat_ZZ Flat_DGHV::batch_encrypt(vector<int> message)const
{
	Mat_ZZ C(l);

	for (int i = 0; i < l; i++)
	{
		Vec_ZZ linie = C_prim[i];
		C[i] = bitdecomp(linie);
		C[i][i] += message[i];
	}
	flatten(C);

	return C;
}

Vec_ZZ Flat_DGHV::batch_decrypt(Mat_ZZ &C)const
{
	Vec_ZZ messages(C.size());

	for (int i = 0; i < l; i++)
	{
		ZZ message(0);
		for (int j = 0; j < l; j++)
		{
			// message += C[i][j] * batch_v[i][j];
			message += C[i][j] * v[j];
		}
		messages[i] = decrypt_DGHV(message);
	}

	return messages;
}

/***************************************************************************/
/***************               pentru calcule pe  GPU      *****************/
void Flat_DGHV::Mat_ZZ_to_mat_int(Mat_ZZ& C, int* &mat)const
{
	int n = C.size();
	int m = C[0].size();

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			conv(mat[i*m+j],C[i][j]);
		}
	}
}

void Flat_DGHV::mat_int_to_Mat_ZZ(int* &mat, Mat_ZZ& C)const
{
	int n = C.size();
	int m = C[0].size();

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			conv(C[i][j], mat[i*m+j]);
		}
	}
}

int get_matrix_int_(int  **output, int size_n, int min_val, int max_val)
{
	int ret = 0;
	*output = (int *)calloc(size_n * size_n, sizeof(int));
	// srand(time(NULL));
	int i, j;
	int range = max_val - min_val + 1;

	for (i = 0; i < size_n; i++)
	{
		for (j = 0; j < size_n; j++)
		{
			(*output)[i * size_n + j] = (rand() % range) + min_val;
		}
	}

	return ret;
}

char *load_kernel_(const char *function)
{
	char *source;
	char file_name[256];
	FILE *file;
	long size;

	strcpy_s(file_name, function);
	strcat_s(file_name, ".cl");

	fopen_s(&file, file_name, "r");
	if (file == NULL)
	{
		fprintf(stderr, "Could not open %s!\n", file_name);
		exit(EXIT_FAILURE);
	}
	fseek(file, 0, SEEK_END);
	size = ftell(file);
	source = (char *)calloc(size, sizeof(char));
	fseek(file, 0, SEEK_SET);
	fread(source, size, sizeof(char), file);
	fclose(file);

	return source;
}

void Flat_DGHV::matrix_multiply(const char *kernel_name, int* matrix1, int* matrix2,
	int* &result)
{
	printf("Start Matrix Multiply!\n");

	/*cl_context context;
	cl_command_queue command_queue;
	cl_program program;
	cl_kernel kernel;
	cl_int ret;
	cl_ev ev;
	const char *source = load_kernel(__FUNCTION__);
	size_t global_work_size[3] = { GLOBAL_WORK_SIZE, 0, 0 }, local_work_size[3] = { LOCAL_WORK_SIZE, 0, 0 };
	char program_build_log[1024] = { 0, };
	cl_mem memory_obj1, memory_obj2, memory_result;
	// int *matrix1, *matrix2, *result = (int *)calloc(MATRIX_SIZE * MATRIX_SIZE, sizeof(int));
	int matrix_size = MATRIX_SIZE;
	float serial_time = 0;
	int kernel_local_chunck_size = MATRIX_SIZE / GLOBAL_WORK_SIZE;*/

	NOT_NULL((c.memory_obj1 = clCreateBuffer(c.context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(int) * MATRIX_SIZE * MATRIX_SIZE, matrix1, NULL)));
	NOT_NULL((c.memory_obj2 = clCreateBuffer(c.context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(int) * MATRIX_SIZE * MATRIX_SIZE, matrix2, NULL)));
	NOT_NULL((c.memory_result = clCreateBuffer(c.context, CL_MEM_READ_WRITE, sizeof(int) * MATRIX_SIZE * MATRIX_SIZE, NULL, NULL)));

	CALL_SAFE((clSetKernelArg(c.kernel, 0, sizeof(cl_mem), (void *)&c.memory_obj1)));
	CALL_SAFE((clSetKernelArg(c.kernel, 1, sizeof(cl_mem), (void *)&c.memory_obj2)));
	CALL_SAFE((clSetKernelArg(c.kernel, 2, sizeof(cl_mem), (void *)&c.memory_result)));
	CALL_SAFE((clSetKernelArg(c.kernel, 3, sizeof(int), (void *)&c.matrix_size)));

	CALL_SAFE((clSetKernelArg(c.kernel, 4, sizeof(int), NULL)));
	CALL_SAFE((clSetKernelArg(c.kernel, 5, sizeof(int) * MATRIX_SIZE, NULL)));

	CALL_SAFE(clEnqueueNDRangeKernel(c.command_queue, c.kernel, 1, NULL, c.global_work_size, c.local_work_size, 0, NULL, &c.ev));

	CALL_SAFE(clEnqueueReadBuffer(c.command_queue, c.memory_result, CL_TRUE, 0, sizeof(int) * MATRIX_SIZE * MATRIX_SIZE, result, 0, NULL, NULL));

	// if (matrix_mul_and_compare(matrix1, matrix2, MATRIX_SIZE, result, &c.serial_time) == 0)
	// {

	cl_ulong kernel_start_time, kernel_end_time;

	CALL_SAFE(clGetEventProfilingInfo(c.ev, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &kernel_start_time, NULL));
	CALL_SAFE(clGetEventProfilingInfo(c.ev, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &kernel_end_time, NULL));

	printf("Kernel Execution Time: %.6f seconds\n", (kernel_end_time - kernel_start_time) / 1000000000.0);

	/*printf("Serial Execution Time: %.6f seconds\n", c.serial_time);
	}
	free((char *)c.source);
	free(matrix1);
	free(matrix2);
	free(result);
	CALL_SAFE(clReleaseKernel(c.kernel));
	CALL_SAFE(clReleaseProgram(c.program));
	CALL_SAFE(clReleaseCommandQueue(c.command_queue));
	CALL_SAFE(clReleaseContext(c.context));*/

	printf("End Matrix Multiply!\n");
}

int Flat_DGHV::create_context()
{
	// int *matrix1 = nullptr, *matrix2 = nullptr, *result = (int *)calloc(MATRIX_SIZE * MATRIX_SIZE, sizeof(int));
	// GPU_context *context = new GPU_context();
	// GPU_context c;
	// find number of supported platforms

	cl_uint num_platforms;
	CALL_SAFE(clGetPlatformIDs(0, NULL, &num_platforms));
	printf("Found %u platforms\n", num_platforms);

	// get pointers to supported platforms
	cl_platform_id *platform_ids = (cl_platform_id *)malloc(num_platforms * sizeof(cl_platform_id));
	CALL_SAFE(clGetPlatformIDs(num_platforms, platform_ids, NULL));

	for (cl_uint i = 0; i < num_platforms; i++)
	{
		// display platform information
		char platform_name[256], version[256];
		CALL_SAFE(clGetPlatformInfo(platform_ids[i], CL_PLATFORM_NAME, sizeof(platform_name), platform_name, NULL));
		printf("\t%s\n", platform_name);

		CALL_SAFE(clGetPlatformInfo(platform_ids[i], CL_PLATFORM_VERSION, sizeof(version), version, NULL));
		printf("\t%s\n", version);

		// get number of devices
		cl_uint num_devices;
		CALL_SAFE(clGetDeviceIDs(platform_ids[i], CL_DEVICE_TYPE_ALL, 0, NULL, &num_devices));
		printf("\n\tFound %u devices\n", num_devices);

		// get pointers to supported devices
		cl_device_id *device_ids = (cl_device_id *)malloc(num_devices * sizeof(cl_device_id));
		CALL_SAFE(clGetDeviceIDs(platform_ids[i], CL_DEVICE_TYPE_ALL, num_devices, device_ids, NULL));

		for (cl_uint j = 0; j < num_devices; j++)
		{
			char device_name[256];
			cl_device_type device_type;

			CALL_SAFE(clGetDeviceInfo(device_ids[j], CL_DEVICE_NAME, sizeof(device_name), device_name, NULL));
			CALL_SAFE(clGetDeviceInfo(device_ids[j], CL_DEVICE_TYPE, sizeof(device_type), &device_type, NULL));

			switch (device_type)
			{
			case CL_DEVICE_TYPE_CPU:
				// printf("\t\t(CPU)\n");
				// matrix_multiply(device_ids[j], "matrix_multiply_cpu");
				break;
			case CL_DEVICE_TYPE_GPU:
				printf("\t\t%s (GPU)\n", device_name);

				c.source = load_kernel_("matrix_multiply");
				c.global_work_size[0] = GLOBAL_WORK_SIZE;
				c.global_work_size[1] = c.global_work_size[2] = 0;
				c.local_work_size[0] = LOCAL_WORK_SIZE;
				c.local_work_size[1] = c.local_work_size[2] = 0;
				c.program_build_log[0] = 0;
				c.matrix_size = MATRIX_SIZE;
				c.serial_time = 0;
				c.kernel_local_chunck_size = MATRIX_SIZE / GLOBAL_WORK_SIZE;

				cl_int ret;
				NOT_NULL_RET((c.context = clCreateContext(NULL, 1, &device_ids[j], NULL, NULL, &ret)));
				NOT_NULL_RET((c.command_queue = clCreateCommandQueue(c.context, device_ids[j], CL_QUEUE_PROFILING_ENABLE, &ret)));

				NOT_NULL_RET((c.program = clCreateProgramWithSource(c.context, 1, &c.source, NULL, &ret)));
				CALL_SAFE(clBuildProgram(c.program, 1, &device_ids[j], NULL, NULL, NULL));
				CALL_SAFE(clGetProgramBuildInfo(c.program, device_ids[j], CL_PROGRAM_BUILD_LOG, sizeof(c.program_build_log), c.program_build_log, NULL));
				if (strlen(c.program_build_log) > 0)
				{
					printf("LOG:\n%s\n", c.program_build_log);
				}

				NOT_NULL_RET((c.kernel = clCreateKernel(c.program, "matrix_multiply_gpu", &ret)));

				/*NOT_NULL((c.memory_obj1 = clCreateBuffer(c.context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(int) * MATRIX_SIZE * MATRIX_SIZE, matrix1, NULL)));
				NOT_NULL((c.memory_obj2 = clCreateBuffer(c.context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(int) * MATRIX_SIZE * MATRIX_SIZE, matrix2, NULL)));
				NOT_NULL((c.memory_result = clCreateBuffer(c.context, CL_MEM_READ_WRITE, sizeof(int) * MATRIX_SIZE * MATRIX_SIZE, NULL, NULL)));
				CALL_SAFE((clSetKernelArg(c.kernel, 0, sizeof(cl_mem), (void *)&c.memory_obj1)));
				CALL_SAFE((clSetKernelArg(c.kernel, 1, sizeof(cl_mem), (void *)&c.memory_obj2)));
				CALL_SAFE((clSetKernelArg(c.kernel, 2, sizeof(cl_mem), (void *)&c.memory_result)));
				CALL_SAFE((clSetKernelArg(c.kernel, 3, sizeof(int), (void *)&c.matrix_size)));
				CALL_SAFE((clSetKernelArg(c.kernel, 4, sizeof(int), NULL)));
				CALL_SAFE((clSetKernelArg(c.kernel, 5, sizeof(int) * MATRIX_SIZE, NULL)));
				CALL_SAFE(clEnqueueNDRangeKernel(c.command_queue, c.kernel, 1, NULL, c.global_work_size, c.local_work_size, 0, NULL, &c.ev));
				CALL_SAFE(clEnqueueReadBuffer(c.command_queue, c.memory_result, CL_TRUE, 0, sizeof(int) * MATRIX_SIZE * MATRIX_SIZE, result, 0, NULL, NULL));
				matrix_multiply(device_ids[j], "matrix_multiply_gpu", matrix1, matrix2, result);*/

				device = device_ids[j];
				break;
			default:
				break;
			}
		}

		free(device_ids);
	}

	free(platform_ids);
	// delete context;
	return 0;
}


void Flat_DGHV::cleanup_gpu_context()
{
	free((char *)c.source);
	CALL_SAFE(clReleaseKernel(c.kernel));
	CALL_SAFE(clReleaseProgram(c.program));
	CALL_SAFE(clReleaseCommandQueue(c.command_queue));
	CALL_SAFE(clReleaseContext(c.context));
}