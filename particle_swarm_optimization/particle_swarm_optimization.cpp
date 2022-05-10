#include "particle_swarm_optimization.h"

using namespace std;

int main()
{
	srand(time(nullptr));
	particle_swarm_optimization PSO(100,2,2);

	PSO.algorithm(200);

	system("pause");
	return 0;
}