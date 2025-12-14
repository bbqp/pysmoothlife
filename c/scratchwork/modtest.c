#include <stdio.h>

int main(void)
{
	int i;
	int M = 64;
	int cspan = 5;
	int start = -M - cspan;
	int end = M + cspan;
	
	for (i = start; i <= end; i++) {
		printf("i = %4d, i / M = %4d, (i / M) * M = %4d, i - (i / M) * M = %4d, i %% M = %4d\n", i, i / M, (i / M) * M, i - (i / M) * M, i % M);
	}

	return 0;
}