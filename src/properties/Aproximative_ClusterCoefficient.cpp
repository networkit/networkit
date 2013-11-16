#include "Aproximative_ClusterCoefficient.h"
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <cstdlib>
#include <ctime>
#include <iostream>

namespace NetworKit {

Aproximative_ClusterCoefficient::Aproximative_ClusterCoefficient() {

}

Aproximative_ClusterCoefficient::~Aproximative_ClusterCoefficient() {

}
 

int FindIndex(double a, int* summe, int Gr)
{
	int Links =1;
	int Rechts = Gr;
	while(Links <= Rechts)
	{
			int Mitte =Links+((Rechts-Links)/2);
			if( (summe[Mitte-1]<a) && (a<= summe[Mitte]))
			return Mitte;
			if (summe[Mitte]< a)
			Links=Mitte+1;
			if( summe[Mitte-1]>=a)
			Rechts= Mitte-1;
			
	}	
	
  return 0;
}

 float Aproximative_ClusterCoefficient::run(const Graph& G, int k) {
	srand(time(NULL));
	unsigned int n=G.numberOfNodes();
	//node degree[n+1];
	int Summe[n+1];
	int Knoten;
	int Knoten2;
	int l=0;
	long sum =0;
	Summe[0]=0;
	int q=0;
	int Knotenarray[n+1];
	for(int i=0;i<n;i++)
	{
		// degree[i+1]=G.degree(i);
		if(G.degree(i)>=2)
		{
			sum += (G.degree(i)*(G.degree(i)-1))/2;
			q++;
			Summe[q] = sum;
			Knotenarray[q]=i;
		}
		std::cout<<"Summe["<< q <<"]="<< sum;
	}
	std::cout<< "sum="<< sum;
	for(int i=1;i<=k;i++)
	
	{
		int size =0;
		int a = (rand()% sum)+1;
		std::cout<<"a="<< a;

		int index = FindIndex(a, &Summe[0], q);
		Knoten=Knotenarray[index];
		std::cout<< "Knoten=" << Knoten;
		node neighbors[n];
		G.forEdgesOf(Knoten, [&](node Knoten, node w)
		{	
			size++;
			neighbors[size]=w;
			
		
		});
		
		int u1 = (rand() %size)+1;
		int Knoten1 = neighbors[u1];
		int u2;
		do{	
			u2=(rand() %size)+1;
		}while(u1==u2);
		std::cout << "u1=" << u1 << ", u2=" << u2 << "\n";
		Knoten2 = neighbors[u2];
		if (G.hasEdge(Knoten1,Knoten2))
			l++;
			std::cout << "l=" << l;
	}
	float g=(float)l/k;
	std::cout << "l=" << l << ", k=" << k << "\n";
	std::cout << "g=" << g;


	return g;


}

}


 














