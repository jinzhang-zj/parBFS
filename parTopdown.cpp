#include<stdio.h>
#include<fstream>
#include<string>
#include<iostream>
#include<sstream>
#include<omp.h>
#include<mpi.h>
#include<vector>
using namespace std;

class Node {
public:
    int* neighbor;
    int num;
    void setNum(int x) {
        num = x;
    }
    void setNeighbor(int* x, int N) {
        neighbor = new int[N];
        for(int i = 0; i < N; i++)
            neighbor[i] = x[i];
    }
};

//create graph from filename;
//numNode returns # of nodes in current partition
//totalNode returns total # of nodes
//function returns a Node array with length numNode
Node* createGraph(int offset, string filename, int *numNode, int *totalNode) {
    ifstream input;
    int size, total, id, num;
    input.open(filename.c_str());
    input >> total;
    input >> size;
    *totalNode = total;
    *numNode = size;
    Node* list = new Node[size];
    for(int i = 0; i < size; i++) {
        input >> id;
        input >> num;
        list[id - offset].setNum(num);
        if(num > 0) {
            list[id - offset].neighbor = new int[num];
            for(int j = 0; j < num; j++) {
                input >> list[id - offset].neighbor[j];
            }
        }
        else list[id-offset].neighbor = NULL;
    }
    return list;
};

int determineSocket(int idx, int numNode){
	return 0;
};

//layer: the array with length totalNode which contains the layer info for all nodes, unvisited node has layer -1
//nextlayer: new layer to search (number)
//parent: array with length numNode which contains the parent info for all local nodes
//this function returns an array with length newFrontierNum which contains new local frontier
//parent is changed. layer will be kept the same
int* topDown(int nextlayer, int* layer, Node* graph, int numNode, int totalNode, int *parent, int *newFrontierNum, int frontierSize, int* frontier, int tasks, int rank, int p){
	printf("entering topdown!\n");
	printf("nextlayer=%d, frontier size=%d\n", nextlayer, frontierSize);

	int* newFrontier;
	*newFrontierNum = 0;
	
	int *l = new int[p];
	for (int i=0; i<p; i++){
		l[i]=0;
	}

	int **sl = new int*[tasks];
	int* rankcum = new int[tasks];

	for (int i=0; i<tasks; i++){
		sl[i] = new int[p];
		for (int j=0; j<p; j++)
			sl[i][j]=0;
		rankcum[i] = 0;
	}

	int rbSize = 0;
	int* sendbuff;
	int* recvbuff;


	#pragma omp parallel
	{
		int tid = omp_get_thread_num();

		// split current frontier
		int from, to;
		from = tid * frontierSize / p;
		to = (tid+1) * frontierSize / p;
		
		vector<int> newsubFrontier;
		vector< vector<int> >subSQ (tasks) ;

		for (int i=from; i<to; i++){
			int u = frontier[i];
			for (int j=0; j<graph[u].num; j++){
				if (layer[graph[u].neighbor[j]] == -1){
					int v = graph[u].neighbor[j];
					int s = determineSocket(v,numNode);
					if (s == rank){
						int vlayer;
						#pragma omp critical
						{	
							vlayer = layer[v];
							layer[v] = nextlayer;
						}
						if (vlayer == -1){
							//printf ("process %d add %d from %d\n", tid, v, u);
							parent[v] = u;
							newsubFrontier.push_back(v);
							l[tid]++;
						}
					}else{
						subSQ.at(s).push_back(u);
						subSQ.at(s).push_back(v);
						sl[s][tid]+=2;
					}
				}
			}
		}
		#pragma omp barrier
		
		#pragma omp master
		{
			for (int i=0; i<tasks && i!=rank; i++){
				for (int j=1; j<p; j++)
					sl[i][j] += sl[i][j-1];
				if (!i)
					rankcum[i] = rankcum[i-1] + sl[i-1][p-1];
			}
		}
		#pragma omp barrier

		// merge the send buffer
		for (int i=0; i<tasks && i!=rank; i++){
			int pcum = 0;
			if (tid) pcum = sl[i][tid-1];

			for (int j=0; j< sl[i][tid] - sl[i][tid-1]; j++){
				sendbuff[rankcum[i] + pcum + j] = subSQ.at(i).at(j);
			}
		}
		#pragma omp barrier

		//printf("process %d done merging\n", tid);
//============================================MPI part start==================================================		
		// send/recieve queue
		// sendbuff: array of tuples sent to socket i
		// MPI_all_to_all_personalized
		
		#pragma omp master 
		{
		int* sendSize = new int[tasks];
		int* recvSize = new int[tasks];
		//send the size of data to each socket
				



		//alltoall send and receive

		}
		#pragma omp barrier
//============================================MPI part end====================================================

		// gather buffer from other socket and process the buffer
		// assume all the data is gathered in Nodetuple* recvbuff
		// recvbuff: array store all tuples
		// rbSize: total number of tuples
		/*
		from = tid * rbSize / p;
		to = (tid+1) * rbSize / p;

		for (int i=from; i<to; i+=2){
			int u = recvbuff[i];
			int v = recvbuff[i+1];
			if (layer[v] == -1){
				int vlayer;
				#pragma omp critical
				{
					vlayer = layer[v];
					layer[v] = nextlayer;
				}
				if (vlayer == -1){
					parent[v] = u;
					newsubFrontier.push_back(v);
					l[tid]++;
				}
			}
		}
		#pragma omp barrier
		*/

		#pragma omp master
		{
			for (int i=1; i<p; i++)	l[i] += l[i-1];
			if ( l[p-1]>0)
				newFrontier = new int[l[p-1]];
			else	
				newFrontier = NULL;
			*newFrontierNum = l[p-1];
		}	
		#pragma omp barrier
		//printf("process %d has newsubFrontier size: %d\n", tid, newsubFrontier.size());
		//for (int i=0; i<newsubFrontier.size(); i++)
		//	cout << newsubFrontier.at(i);
		//cout << endl;

		// merge next frontier
		int start=0;
		if (tid == 0)
			start = 0;
		else
			start = l[tid-1];
		for (int i=0; i< l[tid] - start; i++){
			newFrontier[i + start] = newsubFrontier.at(i);
		}
	}

	if (! *newFrontierNum)	return NULL;
	
	return newFrontier;
};

int main(){
	int rank, size;
	int num_threads = 4;
	/*
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	*/
	int offset;
	rank=0;
	string file= "data";
	Node *graph;
	int numNode, totalNode;

	if(rank==0)	offset=0;

	graph = createGraph(offset, file, &numNode, &totalNode);
	printf ("graph completed!\n");

	int* layer = new int[totalNode];
	int* parent = new int[numNode];
	
	omp_set_num_threads(num_threads);

	#pragma omp parallel for
	for(int i=0; i<totalNode; i++){
		layer[i] = -1;
		parent[i] = -1;
	}

	// root's layer is 0
	layer[0] = 0;
	int nextlayer = 1;

	int* frontier = new int[1];
	frontier[0] = 0;
	int frontiersize = 1;
	int tasks = 1;
	int nextfrontierNum;

	while(frontiersize){
		frontier = topDown(nextlayer, layer, graph, numNode, totalNode, parent, &nextfrontierNum, frontiersize, frontier, tasks, rank, num_threads);
		//printf("while_loop %d\n",frontier[0]);
		frontiersize = nextfrontierNum;
		nextlayer++;
		nextfrontierNum=0;
	}

	for (int i=0; i<numNode; i++){
		printf( "node=%d, parent=%d, layer=%d\n", i, parent[i], layer[i]);
	}


	return 0;
}

