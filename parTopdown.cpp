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
class Graph
{
    int V;
    int N;
public:
    Vec *adj;
    Graph(int V);
    Graph(int rank, string filename);
    void addEdge(int v,int w);
    int localNode() const {return V;}
    int totalNode() const {return N;}
};

Graph::Graph(int V)
{
    this->V = V;
    adj = new Vec[V];
}

void Graph::addEdge(int v,int w)
{
    adj[v].push_back(w);
}

Graph::Graph(int rank, string filename)
{
    ifstream input;
    input.open(filename.c_str());
    input >> this->N;
    input >> this->V;
    adj = new Vec[V];
    int offset = rank*V;
    int v,w, numngh;
    for (int i = 0; i < V; i++)
    {
        input >> v;
        input >> numngh;
        if (numngh>0){
            for (int j = 0; j<numngh; j++)
            {
                input >> w;
                addEdge(v-offset,w);
            }
        }
    }
}

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
    //printf("entering communication\n");
    int recCount[size];
    for (int i=0;i<size;i++)
        recCount[i] = 0;
    //printf("set up the barrier\n");
    MPI_Barrier(MPI_COMM_WORLD);

    //printf("getting the received buff size\n");
    MPI_Alltoall(sendCount,1,MPI_INT,recCount,1, MPI_INT, MPI_COMM_WORLD);
    //printf("getting the received buff size done\n");
    


    int sendDispl[size];
    int recDispl[size];
    scan(sendCount, sendDispl,size);
    scan(recCount,recDispl,size);
    recS = recDispl[size-1]+recCount[size-1];
    recPar = new int[recS];
    recNode = new int[recS];

    MPI_Alltoallv(sendPar,sendCount,sendDispl,MPI_INT,recPar,recCount,recDispl,MPI_INT,MPI_COMM_WORLD);
    MPI_Alltoallv(sendNode,sendCount,sendDispl,MPI_INT,recNode,recCount,recDispl,MPI_INT,MPI_COMM_WORLD);
}

//layer: the array with length totalNode which contains the layer info for all nodes, unvisited node has layer -1
//nextlayer: new layer to search (number)
//parent: array with length numNode which contains the parent info for all local nodes
//this function returns an array with length newFrontierNum which contains new local frontier
//parent is changed. layer will be kept the same
int* topDown(int nextlayer, int* layer, Graph &graph, int numNode, int totalNode, int *parent, int *newFrontierNum, int frontierSize, int* frontier, int tasks, int rank, int p){
	//printf("machine %d entering topdown!\n",rank);
	//printf("nextlayer=%d, frontier size=%d\n", nextlayer, frontierSize);
	int offset=rank*numNode;
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
		//printf("threads %d on rank %d is working\n",tid,rank);
		
	
		// split current frontier
		int from, to;
		from = tid * frontierSize / p;
		to = (tid+1) * frontierSize / p;
		
		//printf("rank %d thread %d process frontier from %d to %d\n", rank, tid, from, to);

		vector<int> newsubFrontier;
		Matpair subSQ (tasks) ;

		for (int i=from; i<to; i++){
			int u = frontier[i];
			for (int j=0; j<graph.adj[u-offset].size(); j++){
				int v = graph.adj[u-offset][j];
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
							//printf ("machine %d process %d add %d from %d\n", rank, tid, v, u);
							parent[v-offset] = u;
							newsubFrontier.push_back(v);
							l[tid]++;
						}
					}
				}else{
					//printf ("machine %d process %d packaging %d for %d\n", rank, tid, v, s);
					subSQ.at(s).push_back( make_pair(u,v)  );
					sl[s][tid]++;
				}
			}
		}
		#pragma omp barrier
		//printf("rank %d done packing\n", rank);
		
		#pragma omp master
		{
		//printf("estimating space...\n");
			for (int i=0; i<tasks; i++){
				//printf("rank %d is accumulating size for rank %d\n", rank, i);
				for (int j=1; j<p; j++){
					sl[i][j] += sl[i][j-1];
				}
				if (i)
					rankcum[i] = rankcum[i-1] + sl[i-1][p-1];
			}

			sendpare = new int[rankcum[tasks-1]+sl[tasks-1][p-1]];
			sendnode = new int[rankcum[tasks-1]+sl[tasks-1][p-1]];
		}
		#pragma omp barrier
		//printf("rank %d done space estimation\n", rank);

		// merge the send buffer
		for (int i=0; i<tasks; i++){
			int pcum = 0;
			if (tid) pcum = sl[i][tid-1];

			for (int j=0; j< sl[i][tid] - pcum; j++){
				sendpare[rankcum[i] + pcum + j] = subSQ[i][j].first;
				sendnode[rankcum[i] + pcum + j] = subSQ[i][j].second;
			}
		}
		//printf("rank %d thread %d done merging\n", rank, tid);
		#pragma omp barrier

		//printf("mpi communication\n", tid);

		// send/recieve buffer
		// sendbuff: array of int sent to socket each socket
		// MPI_all_to_all_personalized
		
		#pragma omp master 
		{
		int* sendSize = new int[tasks];

		for (int i=0; i<tasks; i++){
			sendSize[i] = sl[i][p-1];
			//printf("rank %d send %d elements to rank %d\n", rank, sendSize[i], i);
		}

		//printf("rank %d is communicating \n ",rank);
		alltoallPersonalized(sendpare,sendnode,sendSize,tasks,recvpare,recvnode,rbSize);
		//printf("rank %d received %d data in total\n", rank, rbSize);

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
			int v = recvnode[i];
			if (layer[v] == -1){
				int vlayer;
				#pragma omp critical
				{
					vlayer = layer[v];
					layer[v] = nextlayer;
				}
				if (vlayer == -1){
					
					parent[v-offset] = u;
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
	
	printf("%d total machine, I'm machine %d\n", size, rank);
	
	int offset;
	
    	ostringstream convert;
   	convert << rank+1;
    	string file = "data"+ convert.str();
    	int numNode, totalNode;
    	Graph graph(rank,file);
	//printf ("graph completed!\n");

	numNode= graph.localNode();
	totalNode = graph.totalNode();


	offset = rank*numNode;

	int* layer = new int[totalNode];
	int* parent = new int[numNode];
	
	omp_set_num_threads(num_threads);

	#pragma omp parallel for
	for(int i=0; i<totalNode; i++){
		layer[i] = -1;
	}

	#pragma omp parallel for
	for (int i=0;i<numNode;i++)
		parent[i] = -1;
	
	int root = 0;
	layer[root] = 0;
	// root's layer is 0
	int frontiersize = 0;
	int* frontier;
	
	if (rank == root/numNode){
		parent[root - offset] = -2;
		frontier = new int[1];
		frontiersize = 1;
		frontier[0] = root;
	}
	
	int nextlayer = 1;
	
	int nextfrontierNum;
	int globalfinish=0;
	
	
	while(1){
		MPI_Allreduce(&frontiersize,&globalfinish,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
		if (globalfinish==0)
			break;
		
		frontier = topDown(nextlayer, layer, graph, numNode, totalNode, parent, &nextfrontierNum, frontiersize, frontier, size, rank, num_threads);
		frontiersize = nextfrontierNum;
		
		nextlayer++;
		nextfrontierNum=0;
		break;
	}

	MPI_Finalize();

	for (int i=0; i<numNode; i++){
		printf( "node=%d, parent=%d, layer=%d\n", i+offset, parent[i], layer[i+offset]);
	}

	return 0;
}

