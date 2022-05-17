/*
authors: Min C. Lin, Fabiano S. Oliveira, Paulo E. D. Pinto, Moysés S. Sampaio Jr and Jayme L. Szwarcfiter
year: 2022

This program is part of the research paper titled "Restricted Hamming-Huffman trees"
published in "RAIRO - Operations Research"

Evaluation of the costs in section "5.2. [2]-HHTs efficiency"

*/


#include<stdio.h>
#include<stdlib.h>

const int N = 23;
const int L = 1<<N;

typedef struct node* arvh;

struct node{
	int symbol;
	double frequency;
	arvh left, right;
};

struct tableCell{
	int l;
	double huffmanCost, twoHHTCost, diff;
};

void buildPTriangle(int pTriangle[N+1][N+1],int d);
int phi(int pTriangle[N+1][N+1],int l, int n);


arvh newNode(int simbolo,arvh le,arvh ld,double freq);
void buildHeap(int k,arvh* heap);
void downHeapify(int k,int m,arvh* heap);
arvh buildHuffmanTree(int l,arvh* heap);

double twoHHTCost(int l);
double huffmanCost(arvh root,int h);

int log2(int n);
int divCeil(int x, int y);
int h(int l);
int min(int a,int b);

void printTables(struct tableCell uniform[L+1], struct tableCell zipf[L+1], int length);

int pTriangle[N+1][N+1];

double sumFreqs[L+1], freqs[L+1], cost2HHT, costHFF;

arvh heap[L+1], root;

struct tableCell uniform[L+1], zipf[L+1];

int main(){
	int length = 1;
	int spaceBetweenRows = 10;
	int consecutiveRows = 0;
	double sumZipfFreqs = 0;
	
	buildPTriangle(pTriangle, N);
	
	printf("Processing...\n\n");
	for(int l = 10; spaceBetweenRows<=100000; consecutiveRows++, l+=spaceBetweenRows, length++){
		sumZipfFreqs = 0;

		if(consecutiveRows == 10){
			spaceBetweenRows *=10; 
			consecutiveRows = 0;
		}

		for(int i = 1;i<= l;i++){
			freqs[i] = 1.0;
			heap[i]= newNode(i, NULL, NULL, freqs[i]);
			sumFreqs[i]= sumFreqs[i-1] + freqs[i];
		}

		root = buildHuffmanTree(l,heap);
		costHFF = huffmanCost(root,0)/l;
		cost2HHT = twoHHTCost(l)/l;
		
		uniform[length].l = l;		uniform[length].huffmanCost = costHFF;		uniform[length].twoHHTCost = cost2HHT;
		uniform[length].diff = (cost2HHT/costHFF)*100-100;		

		for(int i = 1;i<= l;i++){
			freqs[i] = 1.0/i;
			sumZipfFreqs +=freqs[i];
			heap[i]= newNode(i,NULL,NULL,freqs[i]);
			sumFreqs[i]= sumFreqs[i-1]+freqs[i];
		}


		root = buildHuffmanTree(l,heap);
		costHFF = huffmanCost(root,0)/sumZipfFreqs;
		cost2HHT = twoHHTCost(l)/sumZipfFreqs;
		
		zipf[length].l = l;		zipf[length].huffmanCost = costHFF;		zipf[length].twoHHTCost =  cost2HHT;
		zipf[length].diff = (cost2HHT/costHFF)*100-100;
	}

	printTables(uniform, zipf, length-1);

	return 0;
}

void buildPTriangle(int pTriangle[N+1][N+1], int d){
	for (int n = 0; n < d; n++){
		for (int k = 0; k <= n; k++){
			if (n == k || !k)	pTriangle[n][k] = 1;
			else			pTriangle[n][k] = pTriangle[n - 1][k - 1] + pTriangle[n- 1][k];
		}
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

arvh newNode(int symbol, arvh left, arvh right, double frequency){
    arvh node;
    
    node = (arvh)malloc(sizeof(*node));
    
    node->symbol = symbol;
    node->left = left;
    node->right = right;
    node->frequency = frequency;
    
    return node;
}

void buildHeap(int k, arvh* heap){
    int x;
    
    for(x = k/2 ;x >= 1; x--)
        downHeapify(x,k,heap);
}

void downHeapify(int k, int m, arvh* heap){
    int x=0,j;
    arvh t;
    
    t = heap[k];  
    if(k<=(m/2))
    	x = 1;
    while(x){
        j = k*2;
        if( (j < m ) && (heap[j]->frequency > heap[j+1]->frequency ) )
            j++;
        if( t->frequency <= heap[j]->frequency)
            x = 0;
        else{
            heap[k] = heap[j]; 
            k = j;
            if(k <= (m/2))
            	x = 1;
            else
            	x = 0;
        }
    }
    heap[k] = t;
}


arvh buildHuffmanTree(int l, arvh* heap){
	arvh p,q,j,aux,root;
	int x;
	
	if(l != 0){
		buildHeap(l, heap);
   		for(x=1; x<l; x++){
        	p = heap[1];

        	aux = heap[l-x+1];
        	heap[l-x+1] = heap[1];
        	heap[1] = aux;

        	downHeapify(1, l-x, heap);
        	q = heap[1];

        	j = newNode(0,p,q,q->frequency + p->frequency);
        	heap[1] = j;
        	downHeapify(1,l-x,heap);
    	}
    	root = heap[1];
    	return root;
	}
	else	
		return NULL;
}


double twoHHTCost(int l){
	double minCost = 0;
	double treeCost;
	int h2;
	int minH1,minH2,minL1;

	for(int i = 1;i <= l; i++)
		minCost += freqs[i];

	minCost *= h(l); //1-HHT cost
   
	for(int h1 = 1; h1 < h(l); h1++)
		for(int l1 = 1; l1 <= min(l-1, 1<<(h1-1)); l1++){
			treeCost = 0;
			int rho = (1<<h1) - l1- phi(pTriangle, l1, h1);

			if(rho > 0){
				treeCost+=h1*sumFreqs[l];

				h2 = h(divCeil(l-l1, rho));
				treeCost+= h2 * (sumFreqs[l] - sumFreqs[l1]);				
				if(treeCost<minCost){
					minL1 = l1;	minH1 = h1; minH2 = h2;
					minCost = treeCost;
				}				
			}
		}

	return minCost;
}

double huffmanCost(arvh root, int level){

	if(root == NULL)	return 0;
	
	if(root->right == NULL)
		return root->frequency*level;
	else{
		return huffmanCost(root->right, level+1) + huffmanCost(root->left, level+1);	
	}
}

int divCeil(int x, int y){
	return  (x + y - 1) / y;
}

int log2(int l){
	int log = -1;
	int ll = l;

	while(ll){
		ll>>=1;	log++;
	}

	return log;
}

int h(int l){
	int logL = log2(l);

	int log2Ceil = l != (1<<logL) ? logL+1:logL;
	
	return log2Ceil+1;
}

int min(int a,int b){
	if(a<b)	return a;
	return b;
}

void printTables(struct tableCell uniform[L+1], struct tableCell zipf[L+1], int length){
	//header
	puts("----------------------------------------------\t\t----------------------------------------------");
	printf("|%45s|\t\t|%45s|\n", "Uniform Probabilities", "Zipf Distribution");
	puts("----------------------------------------------\t\t----------------------------------------------");
	printf("|%10s|%10s|%12s|%10s|\t\t|%10s|%10s|%10s|%10s|\n", "l", "2-HHT cost", "Huffman cost", "diff", "l", "2-HHT cost", "Huffman cost", "diff");
	puts("----------------------------------------------\t\t----------------------------------------------");
	
	for(int i = 1 ; i <= length; i++){
		printf("|%10d|%10.3lf|%12.3lf|%10.3lf|\t\t", uniform[i].l, uniform[i].twoHHTCost, uniform[i].huffmanCost, uniform[i].diff);
		printf("|%10d|%10.3lf|%12.3lf|%10.3lf|\n", zipf[i].l , zipf[i].twoHHTCost, zipf[i].huffmanCost, zipf[i].diff);		
	}
	puts("----------------------------------------------\t\t----------------------------------------------");	
}

