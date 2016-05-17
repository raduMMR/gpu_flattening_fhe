#pragma once
#include <CL/cl.h>

#define MATRIX_SIZE			1024

// cum sa aleg eficient aceste valori ??????????????
#define GLOBAL_WORK_SIZE	1024 // 256
#define LOCAL_WORK_SIZE		256

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