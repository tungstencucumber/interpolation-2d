#include <stdio.h>

int main()
{
	FILE *out = fopen("untitled_0.025.smsh","w");
	double h = 0.025, w = 1.0;
	int N = w/h+1;
	fprintf(out, "%d\n", N*N);
	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++)
			fprintf(out, "%d   %lf %lf\n", i*N+j, j*h, i*h);
	fprintf(out, "\n\n%d\n", (N-1)*(N-1)*2);
	for (int i=0; i<N-1; i++)
		for (int j=0; j<N-1; j++)
		{
			fprintf(out, "%d   %d %d %d\n", 2*(i*(N-1)+j)  , i*N+j, (i+1)*N+j, i*N+j+1);	
			fprintf(out, "%d   %d %d %d\n", 2*(i*(N-1)+j)+1, (i+1)*N+j, (i+1)*N+j+1, i*N+j+1);	
		}
	fclose(out);
	return 0;
}

