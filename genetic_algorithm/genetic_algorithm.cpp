#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

#define N 10					//种群大小（解个数）
#define C 10					//城市个数
#define T 10					//染色体基因个数（T=K+LV+1）
#define cross_rate 0.75			//交叉率
#define muta_rate 0.1			//变异率
#define I 200					//迭代次数

//10个城市的坐标：
double citys_position[C][2] =
	{
		{1304,2312},{3639,1315},{4177,2244},{3712,1399},{3488,1535},
		{3326,1556},{3238,1229},{4196,1004},{4312,7900},{4386,570}
	};

class Genetic_Algorithm
{
private:
	int i, j;
	int population[N][T];			//N-T,初始种群
	int it;							//迭代控制变量
	//int init_population[N][T];		//N-T,初始种群
	double F1[N];					//N，存放初始种群每个解方案（正）目标值
	double F2[N];					//N，存放变异种群每个解方案（正）目标值
	double R[N];					//N,存放每个解方案（倒）目标值
	double p[N];					//N，存放轮盘赌参考概率（每个解被选概率）
	int rand_population[N][T];		//N-T,轮盘赌种群
	int cross_population[N][T];		//N-T,交叉子代种群
	int muta_population[N][T];		//N-T,变异子代种群
	int mix_population[N][T];		//N-T,混合种群
public:
	/*《函数GA申明：遗传算法主体》*/
	void GA();
	/*《函数Initial_Population申明：生成初始种群》*/
	void Initial_Population();
	/*《函数Fitness申明：对传入解计算适应度值》*/
	double Fitness(int* input_solution);
	/*《distance函数声明：计算两两节点的距离，传入两个城市各自的坐标信息》*/
	double distance(double* city1, double* city2);
	/*《函数Selection申明：轮盘赌选出N个个体参与交叉》*/
	void Selection();
	/*《函数Crossover申明：概率交叉》*/
	void Crossover();
	/*《函数Order_Crossover申明：顺序交叉OX，并将交叉结果直接传给cross_population》*/
	void Order_Crossover(int* father, int* mother, int Z1, int Z2);
	/*《函数Mutation申明：概率变异》*/
	void Mutation();
	/*《函数Best_Solution定义：输出本次迭代搜索到的最好解，包括初始种群和变异种群》*/
	void Best_Solution();
	/*《函数Mixing_population声明：为下一迭代准备新的种群》
	构成：70%来自变异种群，20%来自初始种群，10%来自新产生个体。*/
	void Mixing_population();
};

int main()
{
	srand((unsigned)time(NULL));

	Genetic_Algorithm GA1;
	GA1.GA();

	system("pause");
	return 0;
}
/*《函数GA定义：遗传算法主体》*/
void Genetic_Algorithm::GA()
{
	if (N % 2 != 0)
	{
		cout << "种群N设置为2的倍数！" << endl;
	}
	/*①《调用函数Initial_Population：生成初始种群》*/
	Initial_Population();
	/*②《调用函数Fitness：对初始种群传入解计算适应度值》*/
	for (i = 0; i < N; i++)
	{
		F1[i] = Fitness(population[i]);
	}
	/*迭代*/
	for (it = 0; it < I; it++)
	{
		/*③《调用函数Selection：轮盘赌选出N个个体参与交叉》*/
		Selection();
		/*④《调用函数Crossover：概率交叉》*/
		Crossover();
		/*⑤《调用函数Mutation：概率变异》*/
		Mutation();
		/*⑥《调用函数Fitness：对变异种群传入解计算适应度值》*/
		for (i = 0; i < N; i++)
		{
			F2[i] = Fitness(muta_population[i]);
		}
		/*⑦《调用函数Best_Solution：输出本次迭代搜索到的最好解，包括初始种群和变异种群》*/
		Best_Solution();
		/*⑧《调用函数Mixing_population：为下一迭代准备新的种群》
		构成：70%来自变异种群，20%来自初始种群，10%来自新产生个体。*/
		Mixing_population();
	}
}
/*《函数Initial_Population定义：生成初始种群》*/
void Genetic_Algorithm::Initial_Population()
{
	//①一个城市序列
	vector<int> temp_city;
	for (int i = 0; i < C; i++)
	{
		temp_city.push_back(i + 1);
	}
	//②打乱后生成初始种群
	for (i = 0; i < N; i++)
	{
		random_shuffle(temp_city.begin(), temp_city.end());
		for (int j = 0; j < temp_city.size(); j++)
		{
			population[i][j] = temp_city[j];
		}
	}
	cout << "init_population:" << endl;
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < T; j++)
		{
			cout << population[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}
/*《Fitness函数定义：对传入解计算适应度值》*/
double Genetic_Algorithm::Fitness(int* input_solution)
{
	//①初始化路径长度
	double cost = 0;
	//②从第一个城市遍历到最后一个城市
	for (int j = 0; j < C - 1; j++)
	{
		/*《调用distance函数：计算两两节点的距离，传入两个城市各自的坐标信息》*/
		cost += distance(citys_position[input_solution[j] - 1], citys_position[input_solution[j + 1] - 1]);
	}
	//③从第一个城市回到第一个城市
	cost += distance(citys_position[input_solution[C - 1] - 1], citys_position[input_solution[0] - 1]);

	return cost;
}
/*《distance函数声明：计算两两节点的距离，传入两个城市各自的坐标信息》*/
double Genetic_Algorithm::distance(double* city1, double* city2)
{
	/*计算相邻两城市间距离*/
	double x1 = city1[0];			//城市1横坐标
	double y1 = city1[1];
	double x2 = city2[0];
	double y2 = city2[1];			//城市2纵坐标
	double dist = pow((pow((x1 - x2), 2) + pow((y1 - y2), 2)), 0.5);

	return dist;					//返回距离值
}
/*《函数Selection定义：轮盘赌选出N个个体参与交叉》*/
void Genetic_Algorithm::Selection()
{
	double sum = 0;
	double SUM = 0;
	//①计算目标值倒数
	for (i = 0; i < N; i++)
	{
		R[i] = 1 / (F1[i]);
		SUM += R[i];
	}
	//②计算倒目标值的累计概率
	for (i = 0; i < N; i++)
	{
		R[i] = 1 / (F1[i]);
		sum += R[i];
		p[i] = sum / SUM;
	}
	//③轮盘赌选出参与交叉的种群
	for (i = 0; i < N; i++)
	{
		double r = { 0.0 };
		r = (double)rand() / RAND_MAX;
		for (int k = 0; k < N; k++)
		{
			if (r <= p[k])
			{
				for (j = 0; j < T; j++)
				{
					rand_population[i][j] = population[k][j];
				}
				break;
			}
		}
	}
}
/*《函数Crossover定义：概率交叉》*/
void Genetic_Algorithm::Crossover()
{
	int chromoN1, chromoN2;
	int Z1 = 0;
	int Z2 = 1;
	int father[T];				//T,选择交叉父体
	int mother[T];				//T,选择交叉母体
	//顺序交叉N/2次
	for (i = 0; i < N / 2; i++)
	{
		//①选择双亲
		chromoN1 = rand() % N;
		chromoN2 = rand() % N;
		while (chromoN1 == chromoN2)
		{
			chromoN2 = rand() % N;
		}
		for (j = 0; j < T; j++)
		{
			father[j] = rand_population[chromoN1][j];	//父
			mother[j] = rand_population[chromoN2][j];	//母
		}
		//②几率交叉
		double r = (double)rand() / RAND_MAX;
		if (r <= cross_rate)
		{
			/*《调用函数Order_Crossover：顺序交叉OX，并将交叉结果直接传给cross_population》*/
			Order_Crossover(father, mother, Z1, Z2);
			Z1 += 2;//为了cross_population[Z1][]的行在不重复地递增
			Z2 += 2;//为了cross_population[Z2][]的行在不重复地递增
		}
		else
		{
			//不交叉
			for (j = 0; j < T; j++)
			{
				cross_population[Z1][j] = father[j];
				cross_population[Z2][j] = mother[j];
			}
			Z1 += 2;//为了cross_population[Z1][]的行在不重复地递增
			Z2 += 2;//为了cross_population[Z2][]的行在不重复地递增
		}
	}

}
/*《函数Order_Crossover定义：顺序交叉OX并将交叉结果直接传给cross_population》*/
void Genetic_Algorithm::Order_Crossover(int* father, int* mother, int Z1, int Z2)
{
	int cut_point1, cut_point2;
	int child_indiv1[T] = { 0 };	//T，交叉子代1
	int child_indiv2[T] = { 0 };	//T，交叉子代2
	//①随机生成交叉点
	cut_point1 = rand() % (T - 1);	//cutpoint[5]={0,1..T-3}
	cut_point2 = rand() % T;		//cutpoint[5]={0,1,2..T-2}

	while (cut_point1 >= cut_point2)
	{
		cut_point2 = rand() % T;	//cut_point2>cut_point1
	}
	/*
	cout << "cut_point1:" << cut_point1 << "  cut_point2:" << cut_point2 << endl;
	cout << "father:";
	for (j = 0; j < T; j++)
	{
		cout << father[j] << " ";
	}
	cout << endl;
	cout << "mother:";
	for (j = 0; j < T; j++)
	{
		cout << mother[j] << " ";
	}
	cout << endl;
	*/
	//②先取出交叉段中间部分的基因
	for (int x = cut_point1; x <= cut_point2; x++)
	{
		child_indiv1[x] = father[x];		 //子代1取父代的中段基因
		child_indiv2[x] = mother[x];		 //子代2取母代的中段基因
	}
	/*如果第一交换点不是在开头*/
	if (cut_point1 != 0)
	{
		//③先对子代1取首尾段基因
		int index1 = 0;							 //子代1的索引
		for (int y = 0; y < T; y++)				 //搜索母代基因
		{
			bool bt = true;
			for (int z = cut_point1; z <= cut_point2; z++)
			{
				//搜索子代中段基因对比
				if (mother[y] == child_indiv1[z])
				{
					bt = false;					 //改变bool值，使其不执行后续代码
					//搜到重复则跳出当前最近for循环，执行以下，但是bt=false，故下面程序不会被执行！
					break;
				}
			}
			//经过以上对比，若bool值没有改变，即该基因是新基因，则赋给子代
			if (bt == true)
			{
				child_indiv1[index1] = mother[y];//母代赋值给子代1
				index1 += 1;					 //索引子代下一个基因
				//如果首段子代基因赋值完毕，则跳到尾段基因赋值部分
				if (index1 == cut_point1)
				{
					//如果第二交换点不在末尾，那么还有第三段基因
					if (cut_point2 != T - 1)
					{
						index1 = cut_point2 + 1;
					}
						//第二交换点正好在末尾，那么退出
					else
					{
						break;
					}
				}
			}
		}
		//④再对子代2取首尾段基因
		int index2 = 0;							 //子代2的索引
		for (int y = 0; y < T; y++)				 //搜索父代基因
		{
			bool bt = true;
			for (int z = cut_point1; z <= cut_point2; z++)
			{
				//搜索子代中段基因对比
				if (father[y] == child_indiv2[z])
				{
					bt = false;					 //改变bool值，使其不执行后续代码
					//搜到重复则跳出当前最近for循环，执行以下，但是bt=false，故下面程序不会被执行！
					break;
				}
			}
			//经过以上对比，若bool值没有改变，即该基因是新基因，则赋给子代
			if (bt == true)
			{
				child_indiv2[index2] = father[y];//父代赋值给子代2
				index2 += 1;					 //索引子代下一个基因
				//如果首段子代基因赋值完毕，则跳到尾段基因赋值部分
				if (index2 == cut_point1)
				{
					//如果第二交换点不在末尾，那么还有第三段基因
					if (cut_point2 != T - 1)
					{
						index2 = cut_point2 + 1;
					}
						//第二交换点正好在末尾，直接退出
					else
					{
						break;
					}
				}
			}
		}
		//⑤将交叉子代结果传给cross_population
		for (j = 0; j < T; j++)
		{
			cross_population[Z1][j] = child_indiv1[j];
			cross_population[Z2][j] = child_indiv2[j];
		}
		/*
		cout << "child1:";
		for (j = 0; j < T; j++)
		{
			cout << child_indiv1[j] << "  ";
		}
		cout << endl;
		cout << "child2:";
		for (j = 0; j < T; j++)
		{
			cout << child_indiv2[j] << "  ";
		}
		cout << endl;
		*/
	}
		/*如果首交换点在开头*/
	else
	{
		/*如果第二交换点正好也在尾部，那么就不用交叉了*/
		if (cut_point2 == T - 1)
		{
			//⑤将父代们直接传给cross_population
			for (j = 0; j < T; j++)
			{
				cross_population[Z1][j] = father[j];
				cross_population[Z2][j] = mother[j];
			}
		}
			/*如果第二交换点不在尾部*/
		else
		{
			//③先对子代1取第二段基因
			int index1 = cut_point2 + 1;			 //子代1的索引
			for (int y = 0; y < T; y++)				 //搜索母代基因
			{
				bool bt = true;
				for (int z = cut_point1; z <= cut_point2; z++)
				{
					//搜索子代中段基因对比
					if (mother[y] == child_indiv1[z])
					{
						bt = false;					 //改变bool值，使其不执行后续代码
						//搜到重复则跳出当前最近for循环，执行以下，但是bt=false，故下面程序不会被执行！
						break;
					}
				}
				//经过以上对比，若bool值没有改变，即该基因是新基因，则赋给子代
				if (bt == true)
				{
					child_indiv1[index1] = mother[y];//母代赋值给子代1
					index1 += 1;					 //索引子代下一个基因
					/*如果下一个索引没超界*/
					if (index1 != T)
					{
						index1 += 0;
					}
						/*如果下一个索引超界了，则退出*/
					else
					{
						break;
					}
				}
			}
			//④再对子代2取首尾段基因
			int index2 = cut_point2 + 1;			 //子代2的索引
			for (int y = 0; y < T; y++)				 //搜索父代基因
			{
				bool bt = true;
				for (int z = cut_point1; z <= cut_point2; z++)
				{
					//搜索子代中段基因对比
					if (father[y] == child_indiv2[z])
					{
						bt = false;					 //改变bool值，使其不执行后续代码
						//搜到重复则跳出当前最近for循环，执行以下，但是bt=false，故下面程序不会被执行！
						break;
					}
				}
				//经过以上对比，若bool值没有改变，即该基因是新基因，则赋给子代
				if (bt == true)
				{
					child_indiv2[index2] = father[y];//父代赋值给子代2
					index2 += 1;					 //索引子代下一个基因
					/*如果下一个索引没超界*/
					if (index2 != T)
					{
						index2 += 0;
					}
						/*如果下一个索引超界了，则退出*/
					else
					{
						break;
					}
				}
			}
			//⑤将交叉子代结果传给cross_population
			for (j = 0; j < T; j++)
			{
				cross_population[Z1][j] = child_indiv1[j];
				cross_population[Z2][j] = child_indiv2[j];
			}
			/*
			cout << "child1:";
			for (j = 0; j < T; j++)
			{
				cout << child_indiv1[j] << "  ";
			}
			cout << endl;
			cout << "child2:";
			for (j = 0; j < T; j++)
			{
				cout << child_indiv2[j] << "  ";
			}
			cout << endl;
			*/
		}
	}

}
/*《函数Mutation定义：概率变异》*/
void Genetic_Algorithm::Mutation()
{
	//①初始化变异种群为交叉种群
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < T; j++)
		{
			muta_population[i][j] = cross_population[i][j];
		}
	}
	//②变异算子
	int point1, point2;					//变异点
	for (i = 0; i < N; i++)
	{
		double r = { 0.0 };
		r = (double)rand() / RAND_MAX;	//随机一个变异率
		//变异
		if (r <= muta_rate)
		{
			point1 = rand() % T;
			point2 = rand() % T;
			while (point1 == point2)
			{
				point2 = rand() % T;	//避免交换点重复
			}
			swap(muta_population[i][point1], muta_population[i][point2]);
		}
	}
}
/*《函数Best_Solution定义：输出本次迭代搜索到的最好解，包括初始种群和变异种群》*/
void Genetic_Algorithm::Best_Solution()
{
	int big_population[N + N][T];			//融合population和muta_population
	double big_fitness[N + N];				//融合population和muta_population适应度值
	//①融合population和muta_population
	for (i = 0; i < (N + N); i++)
	{
		if (i < N)
		{
			for (j = 0; j < T; j++)
			{
				big_population[i][j] = population[i][j];
			}
			big_fitness[i] = F1[i];
		}
		else
		{
			for (j = 0; j < T; j++)
			{
				big_population[i][j] = muta_population[i - N][j];
			}
			big_fitness[i] = F2[i - N];
		}
	}
	//②选出big_population大种群中最好方案
	int min_index = 0;
	int min_value = big_fitness[0];
	int k;
	for (k = 0; k < (N + N); k++)
	{
		if (big_fitness[k] < min_value)
		{
			min_value = big_fitness[k];
			min_index = k;
		}
	}
	if (min_index < N)
	{
		cout << "第" << it << "代最优方案来自本次初始种群：" << endl;
	}
	else
	{
		cout << "第" << it << "代最优方案来自本次变异种群：" << endl;
	}
	for (j = 0; j < T; j++)
	{
		cout << big_population[min_index][j] << "—>";
	}
	cout << big_population[min_index][0] << endl;
	cout << "适应值为：" << big_fitness[min_index] << endl;

}
/*《函数Mixing_population声明：为下一迭代准备新的种群》
	构成：70%来自变异种群，20%来自初始种群，10%来自新产生个体。*/
void Genetic_Algorithm::Mixing_population()
{
	int muta_num = round(0.7*N);					//变异种群选出个数
	int init_num = round(0.2*N);					//初始种群选出个数
	int gnew_num = N - muta_num - init_num;			//新生成个体数量
	int mix_index = 0;								//混合种群的索引
	double F1_copy[N];								//复制F1[]
	copy(F1, F1 + N, F1_copy);
	double F2_copy[N];
	copy(F2, F2 + N, F2_copy);						//复制F2[]
	double F3[N];									//存放mix_popualtion[][]适应度值

	/*（1）从初始种群选，从优到劣*/
	//①从小到大对适应值排序
	double sort_F1[N];								//存放F1[]排序后值
	int sort_F1_preindex[N];						//存放temp_F1[]排序后的原来索引
	for (i = 0; i < N; i++)
	{
		sort_F1[i] = F1[i];//初始化
	}
	//从小到大排序
	sort(sort_F1, sort_F1 + N);
	//②对排序后的适应值寻找它们之前对应的下标
	for (i = 0; i < N; i++)
	{
		//find()函数返回的是指针地址，所以减去首地址可得下标索引。
		sort_F1_preindex[i] = (find(F1_copy, F1_copy + N, sort_F1[i]) - F1_copy);
		//对已经搜过的值排除
		F1_copy[sort_F1_preindex[i]] = 0;
	}
	/*
	for (i = 0; i < N; i++)
	{
		cout << sort_F1_preindex[i] << "  ";
	}
	cout << endl;
	*/
	//③选出排名靠前的所需个体存进mix_population
	for (i = 0; i < init_num; i++)
	{
		//复制个体
		for (j = 0; j < T; j++)
		{
			mix_population[mix_index][j] = population[sort_F1_preindex[i]][j];
		}
		//复制个体的适应度值
		F3[mix_index] = F1[sort_F1_preindex[i]];
		//混合种群的索引加一
		mix_index += 1;
	}
	/*（2）从变异种群选，从优到劣*/
	//①从小到大对适应值排序
	double sort_F2[N];								//复制变异种群适应度值F2[]
	int sort_F2_preindex[N];						//存放temp_F2[]排序后的原来索引
	for (i = 0; i < N; i++)
	{
		sort_F2[i] = F2[i];
	}
	sort(sort_F2, sort_F2 + N);						//从小到大排序
	//②对排序后的适应值寻找它们之前对应的下标
	for (i = 0; i < N; i++)
	{
		sort_F2_preindex[i] = (find(F2_copy, F2_copy + N, sort_F2[i]) - F2_copy);
		//对已经搜过的值排除
		F2_copy[sort_F2_preindex[i]] = 0;
	}
	/*
	for (i = 0; i < N; i++)
	{
		cout << sort_F2_preindex[i] << "  ";
	}
	cout << endl;
	*/
	//③选出排名靠前的所需个体存进mix_population
	for (i = 0; i < muta_num; i++)
	{
		//复制个体
		for (j = 0; j < T; j++)
		{
			mix_population[mix_index][j] = muta_population[sort_F2_preindex[i]][j];
		}
		//复制个体的适应度值
		F3[mix_index] = F2[sort_F2_preindex[i]];
		//混合种群的索引加一
		mix_index += 1;
	}
	/*（3）新生成个体，并计算适应度值*/
	//①一个城市序列
	int temp_city[C];
	for (int i = 0; i < C; i++)
	{
		temp_city[i] = (i + 1);
	}
	//②打乱后生成新个体
	for (i = 0; i < gnew_num; i++)
	{
		random_shuffle(temp_city, temp_city + C);
		for (j = 0; j < C; j++)
		{
			mix_population[mix_index][j] = temp_city[j];
		}
		F3[mix_index] = Fitness(temp_city);
		mix_index += 1;
	}
	/*（4）将mix_population[][]复制给population[][]
		   将F3[]复制给F1[]*/
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < T; j++)
		{
			population[i][j] = mix_population[i][j];
		}
		F1[i] = F3[i];
	}
	cout << "mix_population:" << endl;
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < T; j++)
		{
			cout << population[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}