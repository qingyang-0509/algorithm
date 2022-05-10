//
// Created by 30205 on 2022/5/10.
//

#ifndef ALGORITHM_PARTICLE_SWARM_OPTIMIZATION_PARTICLE_SWARM_OPTIMIZATION_H_
#define ALGORITHM_PARTICLE_SWARM_OPTIMIZATION_PARTICLE_SWARM_OPTIMIZATION_H_


#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>

using namespace std;

double double_random(double low, double high)
{
	return low + (double)rand() / RAND_MAX * (high - low);
}

class particle_swarm_optimization
{
	friend double double_random(double low, double high);
public:
	particle_swarm_optimization(int n, int a, int b)
		: N(n), c_1(a), c_2(b)
	{
		x = new double[N];
		f = new double[N];
		v = new double[N];
		local_best = new double[N];

		for (int i = 0; i < N; i++)
		{
			x[i] = double_random(0, 2);
			fitness_function();
			v[i] = 0.1;
			local_best[i] = f[i];
		}
		for (int i = 0; i < N; i++)
		{
			if (global_best < f[i])
			{
				global_best = f[i];
			}
		}
	}

	void fitness_function();

	void algorithm(int iteration_steps);

private:
	double* x;
	double* f;
	double* v;
	int N;
	double* local_best;
	double global_best;

	int c_1;
	int c_2;
};

void particle_swarm_optimization::fitness_function()
{
	for (int i = 0; i < N; i++)
	{
		f[i] = (x[i] - 0.7) * (x[i] - 1);
	}
}

void particle_swarm_optimization::algorithm(int iteration_steps)
{
	for (int i = 0; i < iteration_steps; i++)
	{
		double w = 0.4;

		for (int j = 0; j < N; j++)
		{
			v[j] = w * v[i] + c_1 * double_random(0, 1) * (local_best[j] - x[j]) + c_2 * double_random(0, 1) *
				(global_best - x[j]);
			if (v[j] > 0.1)
			{
				v[j] = 0.1;
			}
			x[j] += v[j];

			if (x[j] > 2.0)
			{
				x[j] = 2.0;
			}
			if (x[j] < 0.0)
			{
				x[j] = 0.0;
			}
		}

		fitness_function();

		for (int k = 0; k < N; k++)
		{
			if (local_best[k] < f[k])
			{
				local_best[k] = f[k];
			}
			if (global_best < local_best[k])
			{
				global_best = local_best[k];
			}
		}
		cout << (i + 1) << " " << global_best << endl;
	}
}
#endif //ALGORITHM_PARTICLE_SWARM_OPTIMIZATION_PARTICLE_SWARM_OPTIMIZATION_H_
