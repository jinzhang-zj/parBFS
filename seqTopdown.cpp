#include<stdio.h>
#include<fstream>
#include<string>
#include<iostream>
#include<sstream>
#include<queue>

using namespace std;

class Node {
	public:
	int* neighbor;
	int num;
	int visited;

	void setNum(int x) {
		num = x;
	}
	void setNeighbor(int* x, int N) {
	neighbor = new int[N];
	for(int i = 0; i < N; i++)
		neighbor[i] = x[i];
	}
	int setVisited(int x){
		visited = x;
	}
	
};

//create graph from filename;
//numNode returns # of nodes in current partition
//totalNode returns total # of nodes
//function returns a Node array with length numNode

Node* createGraph(string filename, int *numNode, int *totalNode) {
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
		list[id].setNum(num);
		list[id].setVisited(-1);
		if(num > 0) {
			list[id].neighbor = new int[num];
			for(int j = 0; j < num; j++) {
				input >> list[id].neighbor[j];
			}
		}
		else list[id].neighbor = NULL;
	}
	return list;
};

// implement sequential topdown BFS algorithm
int topDown(Node* graph, queue<int> &frontier, queue<int> &nextFrontier, int* parent, int* layer, int currentlayer){
	while(! frontier.empty()){
		int v = frontier.front();
		layer[v] = currentlayer;
		frontier.pop();
		for (int i=0; i< graph[v].num; i++){
			int u = graph[v].neighbor[i];
			if (graph[u].visited == -1){
				graph[u].setVisited(1);
				parent[u] = v;
				nextFrontier.push(u);
			}
		}
	}
	return 0;
};

int main(int argc, char* argv[])
{
	string file = "data";
	int numNode, totalNode;
	Node *graph = createGraph(file, &numNode, &totalNode);
	queue <int> frontier, nextFrontier;
	int * parent = new int[totalNode];
	int * layer = new int[totalNode];
	for (int i=0; i< totalNode; i++){ 
		parent[i]= -1;
		layer[i] = 0;
	}
	int source = 0;
	int currentlayer = 0;
	frontier.push(source);

	while (!frontier.empty()){
		topDown(graph, frontier, nextFrontier, parent, layer, currentlayer);	
		queue<int> empty;
		swap(frontier,nextFrontier);
		swap(nextFrontier, empty);
		currentlayer++;
	}

	for (int i=0;i<totalNode;i++)
		printf("node=%d, parent=%d, layer=%d\n",i,parent[i],layer[i]);

	return 0;
}
