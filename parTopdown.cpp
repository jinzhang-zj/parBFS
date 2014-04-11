#include<stdio.h>
#include<fstream>
#include<string>
#include<iostream>
#include<sstream>
#include<omp.h>
#include<mpi.h>
#include<vector>

using namespace std;

typedef vector<int> Vec;
typedef vector< vector<pair<int,int> > > Matpair;

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
	return idx/numNode;
};

void scan(int *count,int* displ,int size){
    displ[0]=0;
    for(int i = 0; i<size-1; i++){
        displ[i+1] = displ[i]+count[i];
    }
}

void alltoallPersonalized(int* sendPar, int* sendNode, int* sendCount, int size,//input
                int* &recPar, int* &recNode, int &recS )
{
    int recCount[size];
    for (int i=0;i<size;i++)
        recCount[i] = 0;
    MPI_Alltoall(sendCount,1,MPI_INT,recCount,1, MPI_INT, MPI_COMM_WORLD);

    int sendDispl[size];
    int recDispl[size];
    scan(sendCount, sendDispl,size);
    scan(recCount,recDispl,size);
    recS = recDispl[size-1]+recCount[size-1];
    recPar = new int[recS];
    recNode = new int[recS];

    MPI_Alltoallv(sendPar,sendCount,sendDispl,MPI_INT,recPar,recCount,recDispl,MPI_INT,MPI_COMM_WORLD);
    MPI_Alltoallv(sendNode,sendCount,sendDispl,MPI_INT,recNode,recCount,recDispl,MPI_INT,MPI_COMM_WORLD);
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
	int* sendpare;
	int* sendnode;
	int* recvpare;
	int* recvnode;

	#pragma omp parallel
	{
		int tid = omp_get_thread_num();

		// split current frontier
		int from, to;
		from = tid * frontierSize / p;
		to = (tid+1) * frontierSize / p;
		
		vector<int> newsubFrontier;
		Matpair subSQ (tasks) ;

		for (int i=from; i<to; i++){
			int u = frontier[i];
			for (int j=0; j<graph[u].num; j++){
				int v = graph[u].neighbor[j];
				int s = determineSocket(v,numNode);
				if (s == rank){
					if (layer[v] == -1){
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
					}
				}else{
					subSQ.at(s).push_back( make_pair(u,v)  );
					sl[s][tid]+=1;
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

			for (int j=0; j< sl[i][tid] - pcum; j++){
				sendpare[rankcum[i] + pcum + j] = subSQ[i][j].first;
				sendnode[rankcum[i] + pcum + j] = subSQ[i][j].second;
			}
		}
		#pragma omp barrier

		//printf("process %d done merging\n", tid);
		//printf("mpi communication\n", tid);

		// send/recieve buffer
		// sendbuff: array of int sent to socket each socket
		// MPI_all_to_all_personalized
		
		#pragma omp master 
		{
		int* sendSize = new int[tasks];
		int rbSize;

		for (int i=0; i<tasks; i++){
			sendSize[i] = sl[i][p-1];
		}


		alltoallPersonalized(sendpare,sendnode,sendSize,tasks,recvpare,recvnode,rbSize);


		}
		#pragma omp barrier
		//printf("mpi communication done!\n", tid);

		// gather buffer from other socket and process the buffer
		// assume all the data is gathered in int* recvbuff
		// recvbuff: array store all ints
		// rbSize: total number of ints
		
		from = tid * rbSize / p;
		to = (tid+1) * rbSize / p;

		for (int i=from; i<to; i++){
			int u = recvpare[i];
			int v = recvnode[i+1];
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

int main(int argc, char* argv[]){
	int rank, size;
	int num_threads = 4;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
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

	MPI_Finalize();

	for (int i=0; i<numNode; i++){
		printf( "node=%d, parent=%d, layer=%d\n", i, parent[i], layer[i]);
	}

	return 0;
}

