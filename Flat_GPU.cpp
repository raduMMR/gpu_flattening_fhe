#include "Flat_GPU.h"
#include "Params.h"
#include "utilities.h"
#include <assert.h>

Flat_GPU::Flat_GPU(int lambda, int baza)
{
	w = baza;
	compute_DGHV_settings(lambda);
	compute_FDGHV_settings();
	assert(set_gpu_context() != 0);
}

void Flat_GPU::compute_DGHV_settings(char *filename, int lambda)
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

void Flat_GPU::compute_DGHV_settings(int lambda)
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

void Flat_GPU::compute_FDGHV_settings()
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

	v = new int[l];		// v = Powersof2(1);
	int pow_of_w = 1;
	for (int i = 0; i < l; i++)
	{
		v[i] = pow_of_w;
		pow_of_w *= w;
	}
	
	C_prim.reserve(l);
	for (int i = 0; i < l; i++)
	{
		// pentru fiecare linie avem o criptare a lui zero diferita
		C_prim.push_back(encrypt_DGHV(0));
	}
}

ZZ	Flat_GPU::encrypt_DGHV(int message)const
{
	return encrypt_integer(pk_DGHV, ZZ(message));
}

ZZ	Flat_GPU::decrypt_DGHV(ZZ &ctxt)const
{
	// modifica schema DGHV => modulo t >> 2
	ZZ baza;
	conv(baza, w);
	return ctxt % sk_DGHV % baza;
}


/***************************************************************************/
/***************               pentru calcule pe  GPU      *****************/
void Flat_GPU::matrix_multiply(const char *kernel_name, int* matrix1, int* matrix2, int* &result)
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

int Flat_GPU::set_gpu_context()
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

void Flat_GPU::cleanup_gpu_context()
{
	free((char *)c.source);
	CALL_SAFE(clReleaseKernel(c.kernel));
	CALL_SAFE(clReleaseProgram(c.program));
	CALL_SAFE(clReleaseCommandQueue(c.command_queue));
	CALL_SAFE(clReleaseContext(c.context));
}

int* Flat_GPU::encrypt(int message)const
{
	int* C = new int[l*l];
	int* ptr_C = C;

	for (int i = 0; i < l; i++)
	{
		gpu_bitdecomp( C_prim[i], ptr_C);
		ptr_C[i] += message;
		ptr_C += l;
	}

	gpu_flatten(C, l*l, C);

	return C;
}

int	Flat_GPU::decrypt(int* &C)const
{
	ZZ message(0);
	for (int i = 0; i < l; i++)
	{
		message += C[i] * v[i];
	}

	int miu = 0;
	conv(miu, decrypt_DGHV(message));

	return miu;
}

void Flat_GPU::gpu_hom_mult(int* &matrix1, int* &matrix2, int*& result)
{
	matrix_multiply("matrix_multiply_gpu", matrix1, matrix2, result);

	// flattening
	gpu_flatten(result, l*l, result);
}

void Flat_GPU::gpu_hom_add(int* &matrix1, int* &matrix2, int*& result)
{
	// matrix_multiply("matrix_multiply_gpu", matrix1, matrix2, result);

	// flattening
	gpu_flatten(result, l*l, result);
}

void Flat_GPU::gpu_bitdecomp(ZZ C_i, int* &result)const
{
	ZZ elem(0);
	for (int i = 0; i < l; i++)
	{
		result[i] = elem % w;
		elem /= w;
	}
}

void Flat_GPU::gpu_bitdecomp_1(int* &C_i, ZZ result)const
{
	ZZ pow_of_2(1);
	ZZ two(2);
	for (int i = 0; i < l; i++)
	{
		result += pow_of_2*C_i[i];
		pow_of_2 *= two;
	}
}

// vezi unde se poate folosi gpu-ul pentru a creste viteza
void Flat_GPU::gpu_flatten(int* &C, int dim, int* &result)const
{
	vector<ZZ> bd_1(l);	// !!! declara-l static sau aloca-l in constructor
	ZZ x_0 = pk_DGHV[0];
	int* ptr_result = result;

	for (int i = 0; i < l; i++)
	{
		// bitdecomp_1 pentru fiecare linie a matricii C
		gpu_bitdecomp_1(C, bd_1[i]);
		// reducere modulo x_0
		bd_1[i] = bd_1[i] % x_0;
		// bitdecomp pentru fiecare linie a matricii C
		gpu_bitdecomp(bd_1[i], ptr_result);
		ptr_result += l;
	}

}


