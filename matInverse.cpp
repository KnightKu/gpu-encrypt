#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int** identi(int s){
	int i,j;
	int** a = (int**)malloc(s * sizeof(int*));
	for(i = 0; i < s; i++)
	{
		a[i] = (int*)malloc(s * sizeof(int));
		for(j = 0; j < s; j++)
		{
			if(i == j)			
				a[i][j] = 1;
			else
				a[i][j] = 0;
		}
	}
	return(a);
}

int** cat(int** A, int** ident, int size)
{
	int i, j, k, m;
	int** n = (int**)malloc(size * sizeof(int*));
	for(i = 0; i < size; i++)
	{
		n[i] = (int*)malloc(2*size * sizeof(int));
		for(j = 0; j < size; j++)
		{
			n[i][j] = A[i][j];
		}
		for(j = size; j < 2*size; j++)
		{
			n[i][j] = ident[i][j-size];	
		}
	}
	return(n);
}

int** zeros(int size)
{
        int i,j;
        int** a = (int**)malloc(size * sizeof(int*));
        for(i = 0; i < size; i++)
        {
                a[i] = (int*)malloc(size * sizeof(int));
                for(j = 0; j < size; j++)
                {
                       a[i][j] = 0;
                }
        }
        return(a);	
}

int** decryptKey(int** key, int size){
	int ex,wy;
	int n =  size;	
	int** id = identi(size);
	int** X = cat(key, id, n); //conjoin 2 matrices next to each other
//	deallocate(id);
	int** B;
	int i,b,x, ii, j, a;
	for(i = 0; i < n; i++)
	{
		x = X[i][i];
		ii = i;
		while(x == 0)
		{
			if(ii == (n-1))
			{
				B = zeros(n);
				return(B);
			}
			ii = ii + 1;
			x = ((X[ii][i] + X[i][i]) % 2);//added -1 to account for n
		}
		if(ii != i)
		{
			for(b = 0; b < 2*n; b++)	
			{
				X[i][b] = ((X[i][b] + X[ii][b]) % 2);
			}
		}
		for(j = 0; j < n; j++)
		{
			if(j != i)
			{
				x = X[j][i];
				for(b = 0; b < 2*n; b++)
				{
					X[j][b] = ((X[j][b] - x * X[i][b]) % 2);
				}
			}
		}
	}
	B = (int**)malloc(size * sizeof(int*));
	for(b = 0; b < n; b++)
	{
		B[b] = (int*)malloc(size * sizeof(int));
		for(a = n; a < 2*n; a++)
		{
			B[b][a-n] = X[b][a];
		}
	}
	for(a = 0; a < size; a++)
	{
		free(X[a]);
	}
	free(X);
	return(B);
}
