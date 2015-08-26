/*
   This program is a modified version of the Hill cipher. It generates an nxn
   encryption and decryption key, and reads in a text file as an argument 
   which is stored in an nxn array. The data in the array is expanded to form 
   an nxnx8 array of bits -- each character from the file is expanded to its 
   ascii/binary form and stored in a "bit plane." Each bit plane is then 
   encrypted/decrypted using the encryption/decryption key. Multiplication 
   algorithms use are Square-Matrix-Multiply and Strassen-Recursive.

*/


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdint.h>

#include "matrixmul_kernel.cu"

extern int** identi(int);
extern int** decryptKey(int**, int);

int** allocate(int size)
{
	int i;
	int** m= (int**)malloc(size * sizeof(int*));
	if(m == NULL)
	{
		printf("Out of memory");
		exit(0);
	}
	for(i=0; i<size; i++)
	{
		m[i] = (int*)malloc(size * sizeof(int));
		if(m[i] == NULL)
		{
			printf("out of memory");
			exit(0);
		}
	}
	return m;
}

void deallocate(int** m, int s)
{
	int i;
	for(i=0; i<s; i++)
	{
		free(m[i]);
	}
	free(m);
}

int** add(int** a, int** b, int s)
{
	int i,j;
	int** m = allocate(s);
	for(i=0; i< s; i++)
	{
		for(j=0; j< s; j++)
		{
			m[i][j] = (a[i][j] + b[i][j])%2;	
		}
	}
	return m;
}

int** subtract(int** a, int** b, int s)
{
	int i,j;
	int** m = allocate(s);
	for(i=0; i < s; i++)
	{
		for(j=0; j < s; j++)
		{
			m[i][j] = (a[i][j] - b[i][j])%2;
		}
	}
	return m;
}

int** multiply(int** a, int** b, int size)
{
	int i,j;
	int** c = allocate(size);
	
	if(size == 1)
	{
		c[0][0] = (a[0][0] * b[0][0])%2;
		return c;
	}

	if(size <= 2)
	{
		int a11,a12,a21,a22,b11,b12,b21,b22;	
		a11 = a[0][0];
		a12 = a[0][1];
		a21 = a[1][0];
		a22 = a[1][1];
		b11 = b[0][0];
		b12 = b[0][1];
		b21 = b[1][0];
		b22 = b[1][1];
		
		c[0][0] = (a11*b11 + a12*b21)%2;
        	c[0][1] = (a11*b12 + a12*b22)%2;
        	c[1][0] = (a21*b11 + a22*b21)%2;
		c[1][1] = (a21*b12 + a22*b22)%2;
        	return c;
	}

	size = size/2;

	int** A11 = allocate(size);
	int** A12 = allocate(size);
	int** A21 = allocate(size);
	int** A22 = allocate(size);
	int** B11 = allocate(size);
	int** B12 = allocate(size);
	int** B21 = allocate(size);
	int** B22 = allocate(size);

	for(i=0;i<size;i++)
	{
		for(j=0;j<size;j++)
		{
			A11[i][j] = a[i][j];	
			A12[i][j] = a[i][j+size];
			A21[i][j] = a[i+size][j];
			A22[i][j] = a[i + size][j + size];
			B11[i][j] = b[i][j];
			B12[i][j] = b[i][j + size];
			B21[i][j] = b[i + size][j];
			B22[i][j] = b[i + size][j + size];
		}
	}
	
	int** S1 = subtract(B12,B22,size);
	int** S2 = add(A11,A12, size);
	int** S3 = add(A21,A22, size);
	int** S4 = subtract(B21,B11, size);
	int** S5 = add(A11,A22, size);
	int** S6 = add(B11,B22, size);
	int** S7 = subtract(A12,A22, size);
	int** S8 = add(B21,B22, size);
	int** S9 = subtract(A11,A21, size);
	int** S10 = add(B11,B12, size);

	int** P1 = multiply(A11, S1, size);
	int** P2 = multiply(S2, B22, size);
	int** P3 = multiply(S3, B11, size);
	int** P4 = multiply(A22, S4, size);
	int** P5 = multiply(S5, S6, size);
	int** P6 = multiply(S7, S8, size);
	int** P7 = multiply(S9, S10,size);

	int** c11 = subtract(add(P5,P4,size), add(P2,P6,size), size);
	int** c12 = add(P1,P2,size);
	int** c21 = add(P3,P4,size);
	int** c22 = subtract(add(P5,P1,size), subtract(P3,P7,size), size);
	
	int** temp = add(P5,P4,size);
	int** temp2 = add(P2,P6, size);

	for(i=0; i< size; i++)
	{
		for(j=0; j< size; j++)
		{
			c[i][j] = abs(c11[i][j] % 2);			
			c[i][j+size] = abs(c12[i][j] % 2);
			c[i+size][j] = abs(c21[i][j] % 2);
			c[i+size][j+size] = abs(c22[i][j] % 2);
		}
	}

	deallocate(A11, size);
	deallocate(A12, size);
	deallocate(A21, size);
	deallocate(A22, size);
	deallocate(B11, size);
	deallocate(B12, size);
	deallocate(B21, size);
	deallocate(B22, size);
	deallocate(c11, size);
	deallocate(c12, size);
	deallocate(c21, size);
	deallocate(c22, size);
	deallocate(P1, size);
	deallocate(P2, size);
	deallocate(P3, size);
	deallocate(P4, size);
	deallocate(P5, size);
	deallocate(P6, size);
	deallocate(P7, size);
	deallocate(S1, size);
	deallocate(S2, size);
	deallocate(S3, size);
	deallocate(S4, size);
	deallocate(S5, size);
	deallocate(S6, size);
	deallocate(S7, size);
	deallocate(S8, size);
	deallocate(S9, size);
	deallocate(S10, size);
	deallocate(temp, size);
	deallocate(temp2, size);
	return c;
}

int** squareMatrixMultiply(int** a, int** b, int size)
{
	int i,j,k;
	int** c = allocate(size);
	for(i = 0; i < size; i++)
	{
		for(j = 0; j < size; j++)
		{
			c[i][j] = 0;	
			for(k = 0; k < size; k++)
			{
				c[i][j] = abs(c[i][j] + (a[i][k] * b[k][j]))%2;
			}
		}
	}
	return c;
}

int** generateKey(int n)
{
	int i,j;
	int** k = (int**)malloc(n * sizeof(int*));
	for(i = 0; i < n; i++)
	{
		k[i] = (int*)malloc(n * sizeof(int));
	}
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			k[i][j] = rand() % 2;
		}
	}
	return(k);
}

void strassenRecursive(int input, int*** nd, int yaxis, int*** em, int** key, int count)
{
	int i,j,n,k,m;
	srand(time(NULL));
	cudaEvent_t start;
	cudaEvent_t stop;
	float msecTotal;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	n = input;
	int** b = (int**)malloc(n * sizeof(int*));
	int** c;
	for(i=0; i<n; i++)
	{
		b[i] = (int*)malloc(n * sizeof(int));
	}
	int v = 0;
	cudaEventRecord(start,NULL);
	for(k = 0; k < 8; k++)
	{
		for(m = 0; m < yaxis/n; m++)
		{
			for(i=0; i < n; i++)
			{
				for(j = 0; j < n; j++)
				{
					b[i][j] = nd[i][j+(m*n)][k];
					v++;
				}
			}
			/* Multiply encryption matrix for each nxn matrix in the bitplane */
			if(count == 0)
				c = multiply(key,b,n);
			else if(count == 1)
				c = squareMatrixMultiply(key, b, n);
			for(i=0; i<n; i++)
			{
				for(j = 0; j < n; j++)
				{					/* Copy encrypted values to new 3D array */

					em[i][j+(m*n)][k] = c[i][j];
				}
			}
			deallocate(c, n);
		}
	}
	cudaEventRecord(stop,NULL);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&msecTotal, start, stop);

	printf("%.3fms\n", msecTotal);
	deallocate(b, n);
}

void retrieveFromBinary(int*** nd, int size, int yaxis, int ex)
{
	int i,j,k;
	char c;
	int count;
	for(i = 0; i < size; i++)
	{
		for(j = 0; j < yaxis; j++)
		{
			for(k = 7; k >= 0; k--)
			{
				c <<= 1;
	 			c += nd[i][j][k];
			}
			count++;
			if(ex == 0)
				printf("%c", c);
			else
				printf("%d", c);
		}
	}
	printf("\n\n");
}

void storeToBinary(int letter, int index, int height, int*** nd, int size)
{
	int i,j,k;
	for(i = 7; i >= 0; i--)
	{
		j = letter >> i;
		if(j & 1)
			k = 1;
		else
			k = 0;
		nd[index%size][height][i] = k;
	}
	
}

int compare(int** A, int size)
{
	int** identity; 
	identity = (int**)(intptr_t)identi(size);
	int i,j;
	for(i = 0; i < size; i++)	
	{
		for(j = 0; j < size; j++)
		{
			if(!(A[i][j] == identity[i][j]))				
			{
				deallocate(identity, size);
				return(0);
			}
		}
	}
	deallocate(identity, size);
	return(1);
}

void hillCipher(int size, char* fname, int** key, int** dkey, int count)
{
	FILE* file = fopen(fname, "r");
	char *s;
	long bufsize;
	int yaxis;
	size_t length;
        if(fseek(file, 0L, SEEK_END) == 0)
        {
                bufsize = ftell(file);
                if(bufsize == -1)
                {
                        printf("Error in buff");
                        exit(0);
                }
                s = (char*)malloc(sizeof(char) * (bufsize + 1));
                if(fseek(file, 0L, SEEK_SET) != 0)
                {
                        printf("Error in seek");
                }
                length = fread(s, sizeof(char), bufsize, file);
                if(length == 0)
                {
                        fputs("Error", stderr);
                }
                else
                {
                        s[++length] = '\0';
                }
        }
	int x,y;
	int z = 0;
	int*** dm = (int***)malloc(size * sizeof(int**));
	int*** enMatrix = (int***)malloc(size * sizeof(int**));
	int*** dMatrix = (int***)malloc(size * sizeof(int**));
	/* determine the number of rows */
	yaxis = length/size;
	while(yaxis%size != 0) //fill in rows with 0's until divide evenly
		yaxis += 1;
	for(x = 0; x < size; x++)
	{
		dm[x] = (int**)malloc(yaxis * sizeof(int*));
		enMatrix[x] = (int**)malloc(yaxis * sizeof(int*));
		dMatrix[x] = (int**)malloc(yaxis * sizeof(int*));
		for(y = 0; y < yaxis; y++)
		{
			dm[x][y] = (int*)malloc(8 * sizeof(int));
			enMatrix[x][y] = (int*)malloc(8 * sizeof(int));
			dMatrix[x][y] = (int*)malloc(8 * sizeof(int));
			if(z < (int)length)
			{
				storeToBinary((int)s[z], x, y, dm, size);
				printf("%c", s[z]);
			}
			else
			{
				storeToBinary(0, x, y, dm, size);
			}
			z++;
		}
	}
	printf("Encrypt\n");
	strassenRecursive(size, dm, yaxis, enMatrix, key, count);
	retrieveFromBinary(enMatrix, size, yaxis, 1);
	printf("Decrypt\n");
	strassenRecursive(size, enMatrix, yaxis, dMatrix, dkey, count);
	retrieveFromBinary(dMatrix, size, yaxis, 0);
	free(dm);
	free(enMatrix);
	free(dMatrix);
}

int main(int argc, char *argv[])
{
	srand(time(NULL));
	int k, det, count,ka;
	int i,j,n,m;
	int cuda = 0;
//Cuda variables
	unsigned int mem_size_A;
	unsigned int mem_size_B;
	unsigned int mem_size_C;
	cudaEvent_t start;
	cudaEvent_t stop;
	float msecTotal;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	int** key;
	int** test;
	char* f = argv[1];
	int input[7] = {8,16,32,64,128,256,512};
	int** dkey;
	for(count = 0; count < 1; count++)
	{
		if(count == 0){
			printf("Strassen Times\n");
		}
		else if(count == 1){	
			printf("Square Matrix Multiply\n");
		}
		else
			printf("Cuda Times (encryption, decryption) are\n");
	for(k = 0; k < 1; k++){
		det = 0;
		while(det == 0){
			key = generateKey(input[k]);
			dkey = (int**)(intptr_t)decryptKey(key, input[k]);
			test = squareMatrixMultiply(key, dkey, input[k]);
			det = compare(test, input[k]);
			if(det == 0){
				deallocate(key, input[k]);
				deallocate(dkey, input[k]);
				deallocate(test, input[k]);
			}
		}
		if(count < 2){
	//		cudaEventRecord(start, NULL);
			hillCipher(input[k], f, key, dkey,count);
	//		cudaEventRecord(stop, NULL);
	//		cudaEventSynchronize(stop);
	//		cudaEventElapsedTime(&msecTotal, start, stop);
	//		printf("%.3fms\n", msecTotal);
			deallocate(key, input[k]);
			deallocate(dkey, input[k]);
			deallocate(test, input[k]);
		}
		else{
			FILE* file = fopen(f, "r");
			char *s;
			long bufsize;
			int yaxis;
			size_t length;
			if(fseek(file, 0L, SEEK_END) == 0)
			{
				bufsize = ftell(file);
				if(bufsize == -1)
				{
				        printf("Error in buff");
				        exit(0);
				}
				s = (char*)malloc(sizeof(char) * (bufsize + 1));
				if(fseek(file, 0L, SEEK_SET) != 0)
				{
				        printf("Error in seek");
				}
				length = fread(s, sizeof(char), bufsize, file);
				if(length == 0)
				{
				        fputs("Error", stderr);
				}
				else
				{
				        s[++length] = '\0';
				}
			}
			int x,y;
			int z = 0;
			int*** dm = (int***)malloc(input[k] * sizeof(int**));
			int*** enMatrix = (int***)malloc(input[k] * sizeof(int**));
			int*** dMatrix = (int***)malloc(input[k] * sizeof(int**));
			yaxis = length/input[k];
			mem_size_A = sizeof(float) * input[k] * input[k];
			fclose(file);
			while(yaxis%input[k] != 0) //fill in rows with 0's until divide evenly
				yaxis += 1;
			mem_size_B = sizeof(float) * input[k] * yaxis * 8;
			mem_size_C = sizeof(float) * input[k] * yaxis * 8;
			for(x = 0; x < input[k]; x++)
			{
				dm[x] = (int**)malloc(yaxis * sizeof(int*));
				enMatrix[x] = (int**)malloc(yaxis * sizeof(int*));
				dMatrix[x] = (int**)malloc(yaxis * sizeof(int*));
				for(y = 0; y < yaxis; y++)
				{
					dm[x][y] = (int*)malloc(8 * sizeof(int));
					enMatrix[x][y] = (int*)malloc(8 * sizeof(int));
					dMatrix[x][y] = (int*)malloc(8 * sizeof(int));
					if(z < (int)length)
					{
						storeToBinary((int)s[z], x, y, dm, input[k]);
					}
					else
					{
						storeToBinary(0, x, y, dm, input[k]);
					}
					z++;
				}
			}
//FIRST STRASSEN
				n = input[k];
			//	int** b = (int**)malloc(n * sizeof(int*));
			//	int** c;
				float* d_A, *d_B, *d_C;
				cudaMalloc((void**) &d_A, mem_size_A);
				cudaMalloc((void**) &d_B, mem_size_B);
				cudaMalloc((void**) &d_C, mem_size_C);
				float* h_A = (float*)malloc(mem_size_A);
				float* h_B = (float*)malloc(mem_size_B);
				float* h_C = (float*)malloc(mem_size_C);
				for(i = 0; i < n; i++)
				{
				//	b[i] = (int*)malloc(n * sizeof(int));
				}
//Encryption
				for(ka = 0; ka < 8; ka++)
				{
					for(m = 0; m < yaxis; m++)
					{
						for(i=0; i < n; i++)
						{
							h_B[i+(m*n)+(ka*(n*yaxis))] = dm[i][m][ka];
						}
						//	c = multiply(key,b,n);
						for(i=0; i<n; i++)
						{
							for(j = 0; j < n; j++){

						//		enMatrix[i][j+(m*n)][ka] = c[i][j];
							}
						}
					//	deallocate(c, n);
					}
				}
				for(i = 0; i < n; i++)
				{
					for(j = 0; j < n; j++)
					{
						h_A[j+(i*n)] = key[i][j];
					}
				}
				cudaEventRecord(start,NULL);
				matrixMul<<<input[k],16>>>(d_C, d_A, d_B, input[k], input[k]);
				cudaThreadSynchronize();
				cudaMemcpy(h_C, d_C, mem_size_C, cudaMemcpyDeviceToHost);
				cudaEventRecord(stop, NULL);
				cudaEventSynchronize(stop);
				cudaEventElapsedTime(&msecTotal, start, stop);
				printf("Encrypt: %.3fms\n", msecTotal);
				//retrieveFromBinary(enMatrix, input[k], yaxis);
//SECOND STRASSEN
				n = input[k];
			//	int** ba = (int**)malloc(n * sizeof(int*));
			//	int** ca;
				for(i = 0; i < n; i++)
				{
			//		ba[i] = (int*)malloc(n * sizeof(int));
				}
				for(ka = 0; ka < 8; ka++)
				{
					for(m = 0; m < yaxis; m++)
					{
						for(i=0; i < n; i++)
						{
							h_B[i+(m*n)+(ka*(n*yaxis))] = enMatrix[i][m][ka];
						}
					}
				}
				for(i = 0; i < n; i++)
				{
					for(j = 0; j < n; j++)
					{
						h_A[j+(i*n)] = dkey[i][j];
					}
				}
				cudaEventRecord(start,NULL);
				matrixMul<<<input[k],16>>>(d_C, d_A, d_B, input[k], input[k]);
				cudaThreadSynchronize();
				cudaMemcpy(h_C, d_C, mem_size_C, cudaMemcpyDeviceToHost);
//				deallocate(ba, n);*/
			//	retrieveFromBinary(dMatrix, input[k], yaxis);
				cudaEventRecord(stop,NULL);
				cudaEventSynchronize(stop);
				cudaEventElapsedTime(&msecTotal, start, stop);
				printf("Decrypt: %.3fms\n", msecTotal);
				free(dm);
				free(enMatrix);
				free(dMatrix);
				free(h_A);
				free(h_B);
				free(h_C);
				cudaFree(d_A);
				cudaFree(d_B);
				cudaFree(d_C);
		}
	}
	}
return 0;
}
