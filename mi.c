/*
This program calculates mutual information between gene pairs from their expression profiles using 
B-Spline method (Daub et al.)
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define toler 1E-10

double **B_spline(int k, int M, int N, double *x){
	int i,j,l, *t;
	double max=0, min=100, width, *z, **B, *BB, *w;

	z=malloc(N*sizeof(double));
	t=malloc((M+k)*sizeof(int));
    
    /*open uniform knot vector*/
	for(i=0; i<M+k; i++) {
		if (i<k) t[i]=0;
		else if (i<=M-1) t[i]=i-k+1;
		else t[i]=M-1-k+2;
	}


	for(i=0; i<N; i++) {
		if (x[i]>max) max=x[i];
		if (x[i]<min) min=x[i];
	}
	width=max-min;
	for(i=0; i<N; i++) {
		z[i]=(x[i]-min)*(M-k+1)/width;

	}
    
    B=malloc(N*sizeof(double *));
    for(i=0; i<N; i++){
    	B[i]=malloc((M+k-1)*sizeof(double));
    	for(j=0; j<M+k-1; j++) B[i][j]=0;
    	if (fabs(z[i]-(M-k+1))<toler) B[i][M-1]=1;
    	else {
    		j=0;
    		while(t[j+1]<=z[i]) {
    			j++;
    		}
    		B[i][j]=1;
    	}

    	/*use B values of order l to calculated B values of order l+1, until we get B for order k*/
    	for(l=1; l<k; l++) {
    		BB=malloc((M+k-1-l)*sizeof(double));
    		w=malloc((M+k-1-l+1)*sizeof(double));
    		for(j=0; j<M+k-1-l+1; j++) {
    			if (t[j]==t[j+l]) 
    				w[j]=0;
    			else 
    				w[j]=(z[i]-t[j])/(t[j+l]-t[j]);
    		}
    		for(j=0; j<M+k-1-l; j++)
    			BB[j]=w[j]*B[i][j]+(1-w[j+1])*B[i][j+1];
    		free(B[i]);
    		free(w);
    		B[i]=malloc((M+k-1-l)*sizeof(double));
    		for(j=0; j<M+k-1-l; j++){
    			B[i][j]=BB[j];
            }
    		free(BB);
    	}
        //for(j=0; j<M; j++)
            //printf("%9.6f ", B[i][j]);
        //printf("\n");
    }

    free(z);
    free(t);
	
	return(B);	
}

int main(int argc, char *argv[]){
	int i,j,l, ii, jj, n, N, M=10, k=3;
	char a[20], distr[30];
	double explvl, z, **expr, ***B, H0, H1, H;
    double **pm, ****pj, **MI;
	FILE *fp, *fp2;

    char *inputfile,*outputfile;


    inputfile = malloc(strlen(argv[1]));
    outputfile = malloc(strlen(argv[1]) + 5);
    strcpy(inputfile,argv[1]);
    strcpy(outputfile,"MI_");
    strcat(outputfile,inputfile);
    printf("Input file: %s Output file: %s\n", inputfile, outputfile);

    n=atoi(argv[2]); 
    N=atoi(argv[3]); 

	expr=malloc(n*sizeof(double *));
	expr[0]=malloc(n*N*sizeof(double));
	for (i=1; i<n; i++) expr[i]=expr[i-1]+N;
	pm=malloc(n*sizeof(double*));  //marginal probability
    pm[0]=malloc(n*M*sizeof(double));
    for (i=1; i<n; i++) pm[i]=pm[i-1]+M;
    B=malloc(n*sizeof(double **));
    pj=malloc(n*sizeof(double ***));  //joint probability
    for (i=0; i<n; i++) {
    	pj[i]=malloc(n*sizeof(double **));
    	for (j=0; j<n; j++) {
    		pj[i][j]=malloc(M*sizeof(double *));
    		for (l=0; l<M; l++) 
    			pj[i][j][l]=malloc(M*sizeof(double));
    	}
    }
    MI=malloc(n*sizeof(double *));
    MI[0]=malloc(n*n*sizeof(double));
    for(i=1; i<n; i++) MI[i]=MI[i-1]+n;
    
        
	fp=fopen(inputfile, "r");
    /*calculate marginal probability*/
    for(i=0; i<n; i++) {
    	fscanf(fp,"%s", a);
	    printf("%s\n", a);
	    for(j=0; j<N; j++) {
	    	fscanf(fp, "%lf", &expr[i][j]);
	    }
		
	    B[i]=B_spline(k,M,N,expr[i]);
	    for (j=0; j<M; j++){
	    	pm[i][j]=0;
	    	for(l=0; l<N; l++)
	    		pm[i][j] += B[i][l][j];
	    	pm[i][j]=pm[i][j]/N;
	    }
	    
	}
    fclose(fp);
    
    /*calculate joint probability*/
    for (i=0; i<n; i++) 
        for (j=i; j<n; j++)
            for (ii=0; ii<M; ii++)
                for (jj=0; jj<M; jj++) {
                    pj[i][j][ii][jj]=0;
                    for (l=0; l<N; l++) 
                        pj[i][j][ii][jj] += B[i][l][ii]*B[j][l][jj];
                    pj[i][j][ii][jj]=pj[i][j][ii][jj]/N;
                }

    
    /*calculate and export MI matrix*/
    fp=fopen(outputfile, "w");
    for (i=0; i<n; i++) {
        for (j=i; j<n; j++){
            MI[i][j]=0;
            for(ii=0; ii<M; ii++)
                for(jj=0; jj<M; jj++)
                    if(pj[i][j][ii][jj]>toler) {
                    MI[i][j] += pj[i][j][ii][jj]*log(pj[i][j][ii][jj]/(pm[i][ii]*pm[j][jj]))/log(2);
                    MI[j][i]=MI[i][j];
                }
        }
    }

    for(i=0; i<n; i++) {
        for(j=0; j<n; j++) 
            fprintf(fp, "%.6f ", MI[i][j]);
        fprintf(fp, "\n");
    }
    fclose(fp);     


	free(expr[0]);
	free(expr);
	free(pm[0]);
	free(pm);
    free(MI[0]);
    free(MI);
    
    for(i=0; i<n; i++) {
        free(B[i]);
    }
    free(B);
    for(i=0; i<n; i++) {
        for(j=0; j<n; j++) {
            for(l=0; l<M; l++)
                free(pj[i][j][l]);
            free(pj[i][j]);
        }
        free(pj[i]);
    }
    free(pj); 

    free(inputfile);
    free(outputfile);
	return 0;
}

