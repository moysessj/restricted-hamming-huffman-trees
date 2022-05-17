/*
authors: Min C. Lin, Fabiano S. Oliveira, Paulo E. D. Pinto, Moysés S. Sampaio Jr and Jayme L. Szwarcfiter
year: 2022

This program is part of the research paper titled "Restricted Hamming-Huffman trees"
published in "RAIRO - Operations Research"

Evaluation of the costs used in section "5.1. uniform [3]-HHT optimality hypothesis"

*/

#include<stdio.h>
#include<stdlib.h>


const int N = 13;

const int  L = (1<<(N-1));
const int MINF = -1000000000;
const int INF = 2000000000;
const int DELTA = 4;//Theorem 4

void buildPTriangle(int pTriangle[N+1][N+1],int d);
int phi(int pTriangle[N+1][N+1],int l, int n);


void buildRhoF(int rhoF[L+1][N+1][L+1]);
void buildCF(int cF[L+1][L+1][DELTA+1]);

void printTables(int cF[L+1][L+1][DELTA+1]);
void printTablesTex(int cF[L+1][L+1][DELTA+1]);

int divCeil(int x, int y);
int log2(int n);
int h(int l);

struct cel{int l1; int h1;}; 
cel anscF[L+1][L+1][DELTA+1];

int rhoF[L+1][N+1][L+1], cF[L+1][L+1][DELTA+1];
int pTriangle[N+1][N+1];

int main(){
	buildPTriangle(pTriangle, N);	

	printf("Processing...\n\n");
	buildCF(cF);
	printTables(cF);
	//printTablesTex(cF);

	return 0;
}


void buildPTriangle(int pTriangle[N+1][N+1], int d){
	for (int n = 0; n < d; n++){
		for (int k = 0; k <= n; k++){
			if (n == k || !k)	pTriangle[n][k] = 1;
			else			pTriangle[n][k] = pTriangle[n - 1][k - 1] + pTriangle[n- 1][k];
			printf("%d ",  pTriangle[n][k]);
		}
		puts("");
	}
}

//phiVal = G(l, n-1)
int phi(int pTriangle[N+1][N+1], int l, int n){
	if(!l)	return 0;

	int ll = l;
	int ln = n-1;
	int m = n-1;

	int phiVal;

	// first component of the (n-1)-bounded canonical representation
	while(ll && ll >= pTriangle[ln][m]){
		ll -= pTriangle[ln][m];
		m--;	
	}
	
	phiVal = l - ll;
	phiVal +=  pTriangle[ln][m];

	// second component of the (n-1)-bounded canonical representation, starting at binom(a_k,k)
	while(ll){
		while(ll < pTriangle[ln][m])	ln--;
		
		ll -= pTriangle[ln][m];
		phiVal += pTriangle[ln][m-1];
		m--;
	}

	return phiVal;
}

void buildRhoF(int rhoF[L+1][N+1][L+1]){
	for(int r = 0; r<=L ;r++){
		printf("Evaluating RhoF, r %d of %d\n",r,L);
		for(int nv = 0; nv<=N; nv++){
			int maxL  = r*(1<<(nv-1)); 

			for(int l = 0; l<=L; l++){
				if(!l)								rhoF[r][nv][l] = (1<<nv)*r;
				else if(!r || !nv || l > maxL)		rhoF[r][nv][l] = MINF;
				else{
					int maxLaux =  (1<<(nv-1))>=L?L:(1<<(nv-1));
 					int maxFreeNodes = MINF;

					rhoF[r][nv][l]	= MINF;	
					for(int laux = 0; laux<=maxLaux;laux++){
						int freeNodes =  (1<<nv) - phi(pTriangle,laux, nv) - laux + rhoF[r-1][nv][l-laux];					
						if(freeNodes > maxFreeNodes)	
							maxFreeNodes = freeNodes;				
					}
					rhoF[r][nv][l]	= maxFreeNodes;
				}
			}
		}
	}
	puts("");	
}

void buildCF(int cF[L+1][L+1][DELTA+1]){
	buildRhoF(rhoF);
	
	for(int l = 1; l<=L; l++)
		for(int j = 0; j<=L ;j++){
			cF[l][j][0] = INF;
			cF[l][0][j] = INF;
		}		
	
	for(int r = 1; r<=L ;r++)
		for(int k = 1; k<=DELTA; k++)
			cF[0][r][k] = 0;

	for(int r = 1; r<=L; r++)
		for(int l = 1; l<=L; l++){
			cF[l][r][1] = h(divCeil(l,r))*l;
			anscF[l][r][1].l1 = l;
			anscF[l][r][1].h1 = h(divCeil(l,r));
		}

	for(int k = 2; k<=DELTA; k++){
		printf("Evaluating cF, k %d of %d\n",k,DELTA);
		for(int r = 1; r<=L; r++)
			for(int l = 1; l<=L; l++){
				int maxL1  = l-k+1; 
				int minCost = INF;
				int	ansl1 = 0, ansh1 = 0;

				for(int l1 = 1; l1<=maxL1; l1++){
					int h1lb = h(divCeil(l1,r));
					int h1ub =  h(divCeil(l,r));

					for(int h1 = h1lb; h1 <= h1ub; h1++){
						if(rhoF[r][h1][l1] >= 0){
							int freeNodes = rhoF[r][h1][l1]>=L ? L : rhoF[r][h1][l1];
							int cost = h1*l + cF[l-l1][freeNodes][k-1];
	
							if(cost < minCost){
								minCost = cost;
								ansl1 = l1;	
								ansh1 = h1;						
							}
						}
					}
				}
				cF[l][r][k] = minCost;
				anscF[l][r][k].l1 = ansl1;
				anscF[l][r][k].h1 = ansh1;
			}
	}
	puts("");	
}					

void printTables(int cF[L+1][L+1][DELTA+1]){
	int minK, minCost;
	int twohhtCost;	
	
	puts("---------------------------------------------\t\t---------------------------------------------");
	printf("|%4s|%10s|%8s|%11s|%6s|\t\t|%4s|%10s|%8s|%11s|%6s|\n", "l", "2-HHT cost", "K-HHT: k", "K-HHT: cost", "diff", "l", "2-HHT cost", "K-HHT: k", "K-HHT: cost", "diff");
	puts("---------------------------------------------\t\t---------------------------------------------");

	int desloc = ((L/102)%2)==0?(L/102)+1:L/102; 	
	for(int l = 3; l <= L; l+=desloc){
		minCost = INF;
		int laux = l/10;

		for(int k = 1;k <= DELTA; k++)
			if(cF[l][1][k]< minCost){
				minCost = cF[l][1][k];
				minK= k;
			}

		if(minK==1)
			twohhtCost = cF[l][1][1];
		else
			twohhtCost = cF[l][1][2];

		printf("|%4d|%10.3lf|%8d|%11.3lf|%6.3lf|",l,twohhtCost/((double) l),minK,minCost/((double) l),(twohhtCost*100.0)/minCost -100);
		
		if((l%2))
			printf("\t\t");
		else
			printf("\n");	

		if((l+desloc)>L && l!=L)
			l = L-desloc;	
	}
	
	puts("---------------------------------------------\t\t---------------------------------------------");
	
	puts("\n\n\nAll 2-HHT costs (to be used in a graph): \n");
	
	puts("-----------------");
	printf("|%4s|%10s|\n", "l", " cost");
	puts("-----------------");
	for(int l = 1; l <= L; l++){
		minCost = INF;
		for(int k = 1;k <= DELTA; k++)
		if(cF[l][1][k]< minCost){
			minCost = cF[l][1][k];
			minK= k;
		}
		if(minK==1)
			twohhtCost = cF[l][1][1];
		else
			twohhtCost = cF[l][1][2];
	
		printf(" %4d %10.3lf\n",l,twohhtCost/((double) l));
	}

	puts("\n\n\nAll K-HHT costs (to be used in a graph): ");
	
	puts("-----------------");
	printf("|%4s|%10s|\n", "l", " cost");
	puts("-----------------");
	for(int l = 1; l <= L; l++){
		minCost = INF;
		for(int k = 1;k <= DELTA; k++)
		if(cF[l][1][k]< minCost){
			minCost = cF[l][1][k];
			minK= k;
		}
		printf(" %4d %10.3lf\n",l,minCost/((double) l));
	}
}

void printTablesTex(int cF[L+1][L+1][DELTA+1]){
	int minK, minCost;
	int twohhtCost;	

	printf("\\begin{table}[]\n"); 
	printf("\t \\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|}\n"); 
	printf("\t\t \\cline{1-5} \\cline{7-11}\n");
	printf("\t\t \\multirow{2}{*}{$\\ell$} & \\multirow{2}{*}{\\begin{tabular}[c]{@{}c@{}} $2$-HHT \\\\ cost \\end{tabular}} & \\multicolumn{2}{c|}{\\begin{tabular}[c]{@{}c@{}}lower bound \\\\ $k$-HHT\\end{tabular}} &  \\multirow{2}{*}{\\%% diff} &  & \\multirow{2}{*}{$\\ell$} & \\multirow{2}{*}{\\begin{tabular}[c]{@{}c@{}}$2$-HHT \\\\ cost\\end{tabular}}  &  \\multicolumn{2}{c|}{\\begin{tabular}[c]{@{}c@{}}lower bound \\\\ $k$-HHT\\end{tabular}}  &  \\multirow{2}{*}{\\%% diff} \\\\ \\cline{3-4} \\cline{9-10}\n"); 
    printf("\t\t & & $k$ & cost & &  & & & \\multicolumn{1}{c|}{$k$} & Cost & \\\\ \\cline{1-5} \\cline{7-11}  \n");
	
	int desloc = ((L/102)%2)==0?(L/102)+1:L/102; 
	for(int l = 3; l <= L; l+=desloc){
		minCost = INF;
		int laux = l/10;

		for(int k = 1;k <= DELTA; k++)
		if(cF[l][1][k]< minCost){
			minCost = cF[l][1][k];
			minK= k;
		}
		
		if((l%2))
			printf("\t\t ");

		if(minK==1)
			twohhtCost = cF[l][1][1];
		else
			twohhtCost = cF[l][1][2];

		printf("%d & %.3lf & %d & %.3lf & %.3lf",l,twohhtCost/((double) l),minK,minCost/((double) l),(twohhtCost*100.0)/minCost -100);
		
		if((l%2))
			printf(" & & ");
		else
			printf("\\\\ \\cline{1-5} \\cline{7-11} \n");	
		if((l+desloc)>L && l!=L)
			l = L-desloc;	
	}
	printf("\t \\end{tabular}\n");
	printf("\\end{table}");	
	printf("\n\n");
}


int divCeil(int x, int y){
	return  (x + y - 1) / y;
}

int log2(int n){
	int log = -1;
	int ln = n;

	while(ln){
		ln>>=1;	log++;
	}

	return log;
}

int h(int l){
	int logL = log2(l);

	int log2Ceil = l != (1<<logL) ? logL+1:logL;
	
	return log2Ceil+1;
}

