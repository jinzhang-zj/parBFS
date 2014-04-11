#include<stdio.h>
#include<fstream>
#include<string>
#include<iostream>
#include<sstream>
//#include<omp.h>
#include <mpi.h>
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

void nextFrBuf(const Graph &g, Vec Fr, int rank, int size,  Matpair &buf, int offset )
{
    if (Fr.empty())
        return;
    int r, frNode;
    buf.clear();
    Vec::iterator i,j;
    for(i = Fr.begin(); i != Fr.end(); i++){
        frNode = *i - offset;
        for (j = g.adj[frNode].begin(); j != g.adj[frNode].end(); j++){
            r = *j/g.localNode();
            
            buf[r].push_back(make_pair(*i,*j));
        }
    }
}

void nextFr(int* recPar, int* recNode,int recS ,int offset, int* parents, Vec &Fr )
{
    for (int i = 0; i < recS; i++)
        if (parents[recNode[i]-offset] == -1){
            Fr.push_back(recNode[i]);
            parents[recNode[i]-offset] = recPar[i];
        }
}

void commCount(Matpair &buf,int* count,  int size){
    for(int i = 0; i< size; i++){
        count[i] = buf[i].size();
    }
}

int vecpairSize(Matpair &buf, int size)
{   
    int s = 0;
    for (int i=0; i<size; i++){
        s += buf[i].size();
    }
    return s;
}


void vecpair2arrays(Matpair &buf,int* par, int* node, int size  ){  
    int count[size];
    int curPosi = 0;
    commCount(buf,count,size);
    for (int i=0;i<size;i++){
        for (int j=0;j < count[i]; j++){
            par[curPosi+j] = buf[i][j].first;
            node[curPosi+j] = buf[i][j].second;
        }
        curPosi += count[i];
    }
}

void scan(int *count,int* displ,int size){
    displ[0]=0;
    for(int i = 0; i<size-1; i++){
        displ[i+1] = displ[i]+count[i];
    }
}

void alltoallPersonalized(Matpair &buf, int* parents,  Vec &nxFr, int size,const int &offset)
{
    //how many pairs will be collected from other rank
    nxFr.clear();
    int sendCount[size];
    int recCount[size];
    for (int i=0;i<size;i++){
        sendCount[i] = 0;
        recCount[i] = 0;
    }
    commCount(buf,sendCount,size);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Alltoall(sendCount,1,MPI_INT,recCount,1, MPI_INT, MPI_COMM_WORLD);
    int sendDispl[size];
    int recDispl[size];
    scan(sendCount, sendDispl,size);
    scan(recCount,recDispl,size);

    int sendS = sendDispl[size-1]+sendCount[size-1];
    int recS = recDispl[size-1]+recCount[size-1];

    int sendPar[sendS],recPar[recS],sendNode[sendS],recNode[recS];
    vecpair2arrays(buf,sendPar,sendNode,size);
    MPI_Alltoallv(sendPar,sendCount,sendDispl,MPI_INT,recPar,recCount,recDispl,MPI_INT,MPI_COMM_WORLD);
    MPI_Alltoallv(sendNode,sendCount,sendDispl,MPI_INT,recNode,recCount,recDispl,MPI_INT,MPI_COMM_WORLD);

    nextFr(recPar,recNode,recS,offset,parents,nxFr);

}



int main(int argc, char* argv[])
{
    int rank,size,offset;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  
    ostringstream convert;
    convert << rank+1;
    string file = "data"+ convert.str();
    int numLocal, totalNode;
    Graph g(rank,file);

    Matpair buf(4);
    Vec F;
    numLocal = g.localNode();
    offset = rank*numLocal;
    int root = 512;
    int finish,localfinish;
    int parents[numLocal];
    for (int i=0;i<numLocal;i++){
        parents[i] = -1;
    }
    if (rank == root/numLocal){
        F.push_back(root);
        parents[root-offset] = -2;}
    localfinish = !F.empty(); 

    Vec nxFr;
   
    int dep=0;
    int depth[numLocal];
    for (int i = 0; i < numLocal; i++){
	    depth[i]=0;
    }

    while(1){
    	localfinish=!F.empty();
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&localfinish,&finish,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    	if (finish==0)
		break;
        if (!F.empty()){
	   for (Vec::iterator i = F.begin();i != F.end(); i++){
		depth[*i-offset] = dep+1;
		cout << "node "<< *i <<" depth: "<<dep+1<<endl;
	   }
	}

        nextFrBuf(g,F,rank,size,buf,offset);
        alltoallPersonalized(buf,parents,nxFr,size,offset);
	F=nxFr;
	dep++;
    } 

    MPI_Finalize();
    return 0;
}


