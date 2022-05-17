/*
authors: Min C. Lin, Fabiano S. Oliveira, Paulo E. D. Pinto, Moysés S. Sampaio Jr and Jayme L. Szwarcfiter
year: 2022

This program is part of the research paper titled "Restricted Hamming-Huffman trees"
published in "RAIRO - Operations Research"

Evaluation of the error detection capabilities in section "5.2. [2]-HHTs efficiency"

*/


#include<stdio.h>
#include<stdlib.h>
#include <string>
#include <map>
#include <queue>  
#include <iostream>  
#include <time.h>
#include<string.h>

using namespace std;

const int L = 500000;
const int N = 23;

typedef struct node* arvh;

struct node{
	int symbol;
	double frequency;
	arvh left, right;
};

struct cod{
	int tam;
	unsigned int cod;
};//struct de codificação de cada símbolo


void buildPTriangle(int pTriangle[N+1][N+1],int d);
int phi(int pTriangle[N+1][N+1], int l, int n);

arvh newNode(int simbolo, arvh left, arvh right, double freq);
void buildHeap(int k, arvh* heap);
void downHeapify(int k, int m, arvh* heap);

arvh buildHuffmanTree(int l, arvh* heap);
void getHTSymbolCodes(arvh root, string codification);
bool decodeMessage(arvh treeRoot, string message);

int findWorst2HHT(int lg);
arvh create2HHT(int l);
arvh build2HHT(int l, int level, string codification, int h1, int h2);


double twoHHTCost(int l);
double huffmanCost(arvh root,int h);


void evaluatePte(int l, arvh twoHHTRoot,arvh htRoot);


int log2(int n);
int divCeil(int x, int y);
int h(int l);
int min(int a,int b);

map<string,int> codeCtS;//2HHT: codification to symbol
map<int,string> codeStC;//2HHT: symbol to codification
map<string,int> errorLeaves;// 2HHT:codification to erro leaves

map<int,string> codeStCHT;//HT: symbol to codification

int pTriangle[N+1][N+1];

double sumFreqs[L+1], freqs[L+1], cost2hht, cost2hht2, costHff;

arvh heap[L+1], twoHHTRoot, htRoot;


int main(){
	double zipfFreqSum = 0;
	
	buildPTriangle(pTriangle, N);
		
	int l = 496023;// findWorst2HHT(500000)= 496023 with ratio = 0.999464
	
	srand (time(NULL));
	
		
	codeCtS.clear();
	codeStCHT.clear();
	codeStC.clear();
	errorLeaves.clear();
	zipfFreqSum = 0;

	for(int i = 1; i<= l;i++){
		freqs[i] = 1.0/i;
		zipfFreqSum += freqs[i];
		heap[i]= newNode(i, NULL, NULL, freqs[i]);
		sumFreqs[i]= sumFreqs[i-1]+freqs[i];
	}

	htRoot = buildHuffmanTree(l,heap);
	getHTSymbolCodes(htRoot, "");
	
	twoHHTRoot = create2HHT(l);

	cost2hht = twoHHTCost(l)/zipfFreqSum;
	cost2hht2 = huffmanCost(twoHHTRoot,0)/zipfFreqSum;
	costHff = huffmanCost(htRoot,0)/zipfFreqSum;
	
	if((int) cost2hht*100000 != (int) cost2hht2*100000)
		printf("Error -> l = %d; %lf !=  %lf ?> %.3lf\n",l, cost2hht, cost2hht2,costHff);
		
	evaluatePte(l, twoHHTRoot,htRoot);

	return 0;
}


int findWorst2HHT(int lg){
	int h2;
	int worstL = 0;	
	double zipfFreqSum;
	double worstRatio = 0.0;
	
	printf("Processing 1 <= l <= %d that minimizes the 2-HHT error detection capabilities (for zipf distribuition)\n\n",lg);
	for(int l = 10; l <= lg; l++){
		int minH1 = h(l), minH2=0, minL1=l, errorNodes = (1<<minH1)-minL1;
		double minCost = 0;
		double treeCost;
		zipfFreqSum = 0;

		for(int i = 1;i<= l;i++){
			freqs[i] = 1.0/i;
			sumFreqs[i]= sumFreqs[i-1]+freqs[i];
			zipfFreqSum +=freqs[i];
		}
		for(int i = 1;i <= l; i++)
			minCost += freqs[i];
		minCost *= h(l);
	  
		for(int h1 = 1; h1 < h(l); h1++){
			for(int l1 = 1; l1 <= min(l-1, (1<<(h1-1))); l1++){
				treeCost = 0;
				int rho = (1<<h1) - l1- phi(pTriangle, l1, h1);
				if(rho > 0){
					treeCost += h1*sumFreqs[l];
					h2 = h(divCeil(l-l1, rho));
					treeCost += h2*(sumFreqs[l] - sumFreqs[l1]);		
					if(treeCost < minCost){
						minL1 = l1;	minH1 = h1; minH2 = h2; 
						errorNodes = phi(pTriangle,l1,h1) + (1<<h2)*rho - (l - l1);
						minCost = treeCost;
					}					
				}
			}
		}
		
		if(((double) l)/errorNodes>worstRatio){
			worstRatio = ((double) l)/errorNodes;
			worstL=l;
		}
		
		if(!(l % (lg/10)))
			printf("for l = %d: worst found l = %d with ratio = %lf\n",l, worstL, worstRatio);
	}
	
	printf("Finished processing\n");
	printf("Worst found l = %d with ratio = %lf\n", worstL,worstRatio);
	
	return worstL;
}

arvh create2HHT(int l){
	double minCost = 0;
	double treeCost;
	int h2;
	int minH1=h(l), minH2=0, minL1=l;

	for(int i = 1;i <= l; i++)
		minCost += freqs[i];
		
	minCost *= h(l);//upper bound of the cost  
	
	//find first and second level with symbols of an optimal tree 
	for(int h1 = 1; h1 < h(l); h1++)
		for(int l1 = 1; l1 <= min(l-1, (1<<(h1-1))); l1++){
			treeCost = 0;
			int rho = (1<<h1) - l1- phi(pTriangle, l1, h1);
			
			if(rho > 0){
				treeCost += h1*sumFreqs[l];
				h2 = h(divCeil(l-l1, rho));
				treeCost += h2*(sumFreqs[l] - sumFreqs[l1]);
								
				if(treeCost<minCost){
					minL1 = l1;	minH1 = h1; minH2 = h2;
					minCost = treeCost;
				}				
			}
		}
	
	string codification = "";
	for(int i =1; i <= minH1;i++)	codification+="0";
	
	map<string,int> qmap;
	queue<string> qeo; 
	qeo.push(codification);
	qmap[codification]=1;
	
	if(minH2!=0){
		while(qeo.size()<minL1){//generate first level enconding starting by enconding 0...0
			int qSize = qeo.size();		

			for(int i=1; i <= qSize; qeo.pop(), i++){
				codification = qeo.front();
				for(int j =0; j < minH1;j++){//puts codification's neighboors into the queue
					string aux = codification;
					aux[j] = aux[j] == '0' ? '1' : '0';//flips j-th bit
					if(!qmap.count(aux)){
						qeo.push(aux);
						qmap[aux] = 1;				
					}
				}
				qmap.erase(codification);
			}
		}
		
		for(int i=1; i <= minL1; qeo.pop(), i++){
			codification = qeo.front();
			codeCtS[codification]=i;
			codeStC[i] = codification;
			for(int i=0; i < minH1;i++){
				string aux = codification;
				aux[i] = aux[i]=='0'?'1':'0';
				errorLeaves[aux] = -1;				
			}
		}
	}
	return build2HHT(l,0,"",minH1,minH1+minH2);
}

arvh build2HHT(int l, int level, string codification, int h1,int last_level){
	arvh node;
	
	if(level < last_level){
		node = newNode(0,NULL,NULL,0);
		
		if(codeCtS.count(codification)){	//is a symbol leaf in level h1
			node->symbol = codeCtS[codification]; 
			node->frequency = freqs[codeCtS[codification]];
		}else if(errorLeaves.count(codification)){	//is an error leaf in level h1
			node->symbol = errorLeaves[codification]; 
			node->frequency = 0;
		}else{	//is an internal node
			node->right = build2HHT(l, level+1, codification+"1", h1, last_level);
			node->left = build2HHT(l, level+1, codification+"0", h1, last_level);
			node->frequency = node->right->frequency + node->left->frequency;
		}	
	}else{	/*for level last_level all nodes are either symbol or error leaves, 
			depending on the parity of the  codification*/
		node = newNode(-1,NULL,NULL,0);
		errorLeaves[codification] = -1;	//by default set new leaf as an error leaf		

		int countBitsOne = 0;
		for(int i=0;i < codification.size();i++)
			if(codification[i] == '1')
				countBitsOne++;

		bool evenCode = countBitsOne & 1 ? false : true;
		if(evenCode){
			if(codeCtS.size()<l){// creates necessary symbol leaves
				codeCtS[codification] = codeCtS.size()+1;
				codeStC[codeCtS[codification]] = codification;
				
				node->symbol = codeCtS[codification]; 
				node->frequency = freqs[codeCtS[codification]];
					
				errorLeaves.erase(codification);		
			}
		}
				
	}
	return node;
}

void evaluatePte(int l, arvh twoHHTRoot, arvh htRoot){
	int bValues[10]={10,25,50,100,250,500,1000,2500,5000};
	int randomBit;

	printf("\nProcessing -> b in [%d", bValues[0]);
	for(int i = 0; i <= 8; i++)
		printf(", %d", bValues[i]);
	printf("]\n\n");
	
	puts("--------------------------------------");
	printf("|%10s|%12s|%12s|\n", "b", "Huffman Tree","2-HHT");
	puts("--------------------------------------");

	for(int bIndex = 0; bIndex <= 8; bIndex++){
		int b = bValues[bIndex];
		int printInterval = 250000000/b;
		int countErr2HHT, countErrHT, totalTests;
		string msgHT, msgHTError, msg2HHT, msg2HHTError;
		
		totalTests = countErr2HHT = countErrHT = 0;		
		msgHT = msg2HHT ="";
		
		for(int i = 1;i <= b;i++){
			int randomSymbol = rand() % l + 1;
			msg2HHT+=codeStC[randomSymbol];
			msgHT+=codeStCHT[randomSymbol];
		}//random encoded messages with b symbols

		for(int i = 1; i <= 1000000;i++){
			for(int errorBits = 1; errorBits <= min(20,b); errorBits++){
				totalTests++;
				msg2HHTError = msg2HHT;
				msgHTError = msgHT;
				
				for(int j = 1; j <= errorBits; j++){
					randomBit = rand() % msg2HHTError.size();
					msg2HHTError[randomBit] = msg2HHTError[randomBit]=='0'?'1':'0';	

					randomBit = rand() % msgHTError.size();
					msgHTError[randomBit] = msgHTError[randomBit]=='0'?'1':'0';				
				}
				
				if(!decodeMessage(twoHHTRoot,msg2HHTError))
					countErr2HHT++;
				if(!decodeMessage(htRoot,msgHTError))
					countErrHT++;
			}
		}
		
		printf("|%10d|%12.4lf|%12.4lf|\n", b,(countErrHT/((double)totalTests))*100,(countErr2HHT/((double)totalTests))*100);
	}
	
	puts("--------------------------------------");
}

bool decodeMessage(arvh treeRoot, string message){
	arvh node = treeRoot;

	for(int i = 0; i < message.size();i++){
		node = message[i] == '1'?node->right: node->left;
		
		if(node->right==NULL){
			if(node->symbol == -1)
				return false;
			node = treeRoot;	
		}	
	}
	
	if(node != treeRoot)
		return false;
	return true;
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

void getHTSymbolCodes(arvh root, string codification){
	if(root!= NULL){
		if(root->symbol > 0)
			codeStCHT[root->symbol]  = codification;
		getHTSymbolCodes(root->right, codification+"1");
		getHTSymbolCodes(root->left, codification+"0");		
	}
}

double twoHHTCost(int l){
	double minCost = 0;
	double cost;
	int h2;
	int minH1,minH2,minL1;

	for(int i = 1;i <= l; i++)
		minCost += freqs[i];

	minCost *= h(l); //1-HHT cost
   
	for(int h1 = 1; h1 <= h(l); h1++)
		for(int l1 = 1; l1 <= min(l-1, 1<<(h1-1)); l1++){
			cost = 0;
			int rho = (1<<h1) - l1- phi(pTriangle, l1, h1);

			if(rho > 0){
				cost+=h1*sumFreqs[l];

				h2 = h(divCeil(l-l1, rho));
				cost+= h2 * (sumFreqs[l] - sumFreqs[l1]);				
				if(cost<minCost){
					minL1 = l1;	minH1 = h1; minH2 = h2;
					minCost = cost;
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
