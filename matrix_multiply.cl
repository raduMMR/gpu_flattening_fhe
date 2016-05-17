#pragma OPENCL EXTENSION cl_amd_printf : enable

#define element(matrix, n, x, y)	matrix[x * n + y]
#define MAX_MATRIX_SIZE				1024

__kernel void matrix_multiply_gpu(__global const int *matrix1, __global const int *matrix2, __global int *result, const int size, 
	__local int *matrix1_local, __local int *matrix2_local)
{
	size_t local_id = get_local_id(0);
	size_t local_size = get_local_size(0);
    size_t global_id = get_global_id(0);
	
	size_t i, j;
	int tmp;
	__private int line[MAX_MATRIX_SIZE];
	
	for (i = 0; i < size; i++)
	{
		line[i] = element(matrix1, size, global_id, i); 
	}

	for (i = 0; i < size; i++)
	{
		tmp = 0;

		for (j = local_id; j < size; j += local_size)
		{
			matrix2_local[j] = element(matrix2, size, j, i);
		}
		
		barrier(CLK_LOCAL_MEM_FENCE);
		
		for (j = 0; j < size; j++)
		{
			tmp += line[j] * matrix2_local[j];
			// tmp += line[j] * element(matrix2, size, j, i);
			// tmp += element(matrix1, size, global_id, j) * element(matrix2, size, j, i);		
		}
		element(result, size, global_id, i) = tmp;
	}
}

__kernel void matrix_multiply_cpu(__global const int *matrix1, __global const int *matrix2, __global int *result, const int size,
	__local int *matrix1_local, __local int *matrix2_local)
{
	size_t local_id = get_local_id(0);
    size_t global_id = get_global_id(0);
	size_t global_size = get_global_size(0);
	size_t start_line = size / global_size * global_id, stop_line = start_line + size / global_size;
	size_t i, j, k;

	for (i = start_line; i < stop_line; i++)
	{
		for (j = 0; j < size; j++)
		{
			for (k = 0; k < size; k++)
			{
				element(result, size, i, j) += element(matrix1, size, i, k) * element(matrix2, size, k, j);
			}
		}
	}
}