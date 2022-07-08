#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main(int argc, char* argv[])
{
	int i,j, n;
	double **MI, **fxy, *miu, *sigma, **zx, **zy;
	FILE *fp;
	char *inputfile,*outputfile;


	inputfile = malloc(strlen(argv[1]));
	outputfile = malloc(strlen(argv[1]) + 5);
	strcpy(inputfile,argv[1]);
	strcpy(outputfile,"f_");
	strcat(outputfile,inputfile);
    printf("Input file: %s Output file: %s\n", inputfile, outputfile);

    n=atoi(argv[2]);

	MI    = malloc(n*sizeof(double*));
	MI[0] = malloc(n*n*sizeof(double));
	for(i=1; i<n; i++) MI[i]=MI[i-1]+n;

	fxy	   = malloc(n*sizeof(double*));
    fxy[0] = malloc(n*n*sizeof(double));
    for(i=1; i<n; i++) fxy[i]=fxy[i-1]+n;

    miu = malloc(n*sizeof(double));
    sigma=malloc(n*sizeof(double));

    zx    = malloc(n*sizeof(double*));
    zx[0] =malloc(n*n*sizeof(double));
    for(i=1; i<n; i++) zx[i]=zx[i-1]+n;

    zy    = malloc(n*sizeof(double*));
    zy[0] =malloc(n*n*sizeof(double));
    for(i=1; i<n; i++) zy[i]=zy[i-1]+n;

	fp=fopen(inputfile, "r");

	for(i=0; i<n; i++)
		for(j=0; j<n; j++)
			fscanf(fp, "%lf", &MI[i][j]);
	fclose(fp);

	for(i=0; i<n; i++){
		MI[i][i]=0;
		miu[i]=0;
		for(j=0; j<n; j++) 
			miu[i] += MI[i][j];
		miu[i]=miu[i]/n;
		//printf("miu%d is %g\n", i+1, miu[i]);
	}

	for(i=0; i<n; i++){
		sigma[i]=0;
	    for(j=0; j<n; j++)
	    	sigma[i] += (MI[i][j]-miu[i])*(MI[i][j]-miu[i]);
	    sigma[i]=sqrt(sigma[i]/(n-1));
	   // printf("sigma%d is %g\n", i+1, sigma[i]);
	}

	for(i=0; i<n; i++)
		for(j=i; j<n; j++){
			zx[i][j]=(MI[i][j]-miu[i])/sigma[i];
			if(zx[i][j]<0) zx[i][j]=0;
			//zx[i][j]=zx[i][j]*zx[i][j];
			zx[j][i]=zx[i][j];
			zy[i][j]=(MI[i][j]-miu[j])/sigma[j];
			if(zy[i][j]<0) zy[i][j]=0;
			//zy[i][j]=zy[i][j]*zy[i][j];
			zy[j][i]=zy[i][j];
		}

	for(i=0; i<n; i++)
		for(j=i; j<n; j++){
			fxy[i][j]=sqrt(zx[i][j]*zx[i][j]+zy[i][j]*zy[i][j]);
			fxy[j][i]=fxy[i][j];
		}

	fp=fopen(outputfile, "w");
	for(i=0; i<n; i++){
		for(j=0; j<n; j++)
			fprintf(fp, "%8.6f ", fxy[i][j]);
		fprintf(fp, "\n");
	}
	fclose(fp);

//
        for(i=0; i<10; i++){
           for(j=0; j<10; j++)
              printf("%8.6f ", fxy[i][j]);
           printf("\n");
         }

//
    free(MI[0]);
    free(MI);
    free(fxy[0]);
    free(fxy);
    free(zx[0]);
    free(zx);
    free(zy[0]);
    free(zy);
    free(miu);
    free(sigma);
    free(inputfile);
    free(outputfile);
	return 0;
}
