#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

#define N 10					//��Ⱥ��С���������
#define C 10					//���и���
#define T 10					//Ⱦɫ����������T=K+LV+1��
#define cross_rate 0.75			//������
#define muta_rate 0.1			//������
#define I 200					//��������

//10�����е����꣺
double citys_position[C][2] =
	{
		{1304,2312},{3639,1315},{4177,2244},{3712,1399},{3488,1535},
		{3326,1556},{3238,1229},{4196,1004},{4312,7900},{4386,570}
	};

class Genetic_Algorithm
{
private:
	int i, j;
	int population[N][T];			//N-T,��ʼ��Ⱥ
	int it;							//�������Ʊ���
	//int init_population[N][T];		//N-T,��ʼ��Ⱥ
	double F1[N];					//N����ų�ʼ��Ⱥÿ���ⷽ��������Ŀ��ֵ
	double F2[N];					//N����ű�����Ⱥÿ���ⷽ��������Ŀ��ֵ
	double R[N];					//N,���ÿ���ⷽ��������Ŀ��ֵ
	double p[N];					//N��������̶Ĳο����ʣ�ÿ���ⱻѡ���ʣ�
	int rand_population[N][T];		//N-T,���̶���Ⱥ
	int cross_population[N][T];		//N-T,�����Ӵ���Ⱥ
	int muta_population[N][T];		//N-T,�����Ӵ���Ⱥ
	int mix_population[N][T];		//N-T,�����Ⱥ
public:
	/*������GA�������Ŵ��㷨���塷*/
	void GA();
	/*������Initial_Population���������ɳ�ʼ��Ⱥ��*/
	void Initial_Population();
	/*������Fitness�������Դ���������Ӧ��ֵ��*/
	double Fitness(int* input_solution);
	/*��distance�������������������ڵ�ľ��룬�����������и��Ե�������Ϣ��*/
	double distance(double* city1, double* city2);
	/*������Selection���������̶�ѡ��N��������뽻�桷*/
	void Selection();
	/*������Crossover���������ʽ��桷*/
	void Crossover();
	/*������Order_Crossover������˳�򽻲�OX������������ֱ�Ӵ���cross_population��*/
	void Order_Crossover(int* father, int* mother, int Z1, int Z2);
	/*������Mutation���������ʱ��졷*/
	void Mutation();
	/*������Best_Solution���壺������ε�������������ý⣬������ʼ��Ⱥ�ͱ�����Ⱥ��*/
	void Best_Solution();
	/*������Mixing_population������Ϊ��һ����׼���µ���Ⱥ��
	���ɣ�70%���Ա�����Ⱥ��20%���Գ�ʼ��Ⱥ��10%�����²������塣*/
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
/*������GA���壺�Ŵ��㷨���塷*/
void Genetic_Algorithm::GA()
{
	if (N % 2 != 0)
	{
		cout << "��ȺN����Ϊ2�ı�����" << endl;
	}
	/*�١����ú���Initial_Population�����ɳ�ʼ��Ⱥ��*/
	Initial_Population();
	/*�ڡ����ú���Fitness���Գ�ʼ��Ⱥ����������Ӧ��ֵ��*/
	for (i = 0; i < N; i++)
	{
		F1[i] = Fitness(population[i]);
	}
	/*����*/
	for (it = 0; it < I; it++)
	{
		/*�ۡ����ú���Selection�����̶�ѡ��N��������뽻�桷*/
		Selection();
		/*�ܡ����ú���Crossover�����ʽ��桷*/
		Crossover();
		/*�ݡ����ú���Mutation�����ʱ��졷*/
		Mutation();
		/*�ޡ����ú���Fitness���Ա�����Ⱥ����������Ӧ��ֵ��*/
		for (i = 0; i < N; i++)
		{
			F2[i] = Fitness(muta_population[i]);
		}
		/*�ߡ����ú���Best_Solution��������ε�������������ý⣬������ʼ��Ⱥ�ͱ�����Ⱥ��*/
		Best_Solution();
		/*�ࡶ���ú���Mixing_population��Ϊ��һ����׼���µ���Ⱥ��
		���ɣ�70%���Ա�����Ⱥ��20%���Գ�ʼ��Ⱥ��10%�����²������塣*/
		Mixing_population();
	}
}
/*������Initial_Population���壺���ɳ�ʼ��Ⱥ��*/
void Genetic_Algorithm::Initial_Population()
{
	//��һ����������
	vector<int> temp_city;
	for (int i = 0; i < C; i++)
	{
		temp_city.push_back(i + 1);
	}
	//�ڴ��Һ����ɳ�ʼ��Ⱥ
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
/*��Fitness�������壺�Դ���������Ӧ��ֵ��*/
double Genetic_Algorithm::Fitness(int* input_solution)
{
	//�ٳ�ʼ��·������
	double cost = 0;
	//�ڴӵ�һ�����б��������һ������
	for (int j = 0; j < C - 1; j++)
	{
		/*������distance���������������ڵ�ľ��룬�����������и��Ե�������Ϣ��*/
		cost += distance(citys_position[input_solution[j] - 1], citys_position[input_solution[j + 1] - 1]);
	}
	//�۴ӵ�һ�����лص���һ������
	cost += distance(citys_position[input_solution[C - 1] - 1], citys_position[input_solution[0] - 1]);

	return cost;
}
/*��distance�������������������ڵ�ľ��룬�����������и��Ե�������Ϣ��*/
double Genetic_Algorithm::distance(double* city1, double* city2)
{
	/*�������������м����*/
	double x1 = city1[0];			//����1������
	double y1 = city1[1];
	double x2 = city2[0];
	double y2 = city2[1];			//����2������
	double dist = pow((pow((x1 - x2), 2) + pow((y1 - y2), 2)), 0.5);

	return dist;					//���ؾ���ֵ
}
/*������Selection���壺���̶�ѡ��N��������뽻�桷*/
void Genetic_Algorithm::Selection()
{
	double sum = 0;
	double SUM = 0;
	//�ټ���Ŀ��ֵ����
	for (i = 0; i < N; i++)
	{
		R[i] = 1 / (F1[i]);
		SUM += R[i];
	}
	//�ڼ��㵹Ŀ��ֵ���ۼƸ���
	for (i = 0; i < N; i++)
	{
		R[i] = 1 / (F1[i]);
		sum += R[i];
		p[i] = sum / SUM;
	}
	//�����̶�ѡ�����뽻�����Ⱥ
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
/*������Crossover���壺���ʽ��桷*/
void Genetic_Algorithm::Crossover()
{
	int chromoN1, chromoN2;
	int Z1 = 0;
	int Z2 = 1;
	int father[T];				//T,ѡ�񽻲游��
	int mother[T];				//T,ѡ�񽻲�ĸ��
	//˳�򽻲�N/2��
	for (i = 0; i < N / 2; i++)
	{
		//��ѡ��˫��
		chromoN1 = rand() % N;
		chromoN2 = rand() % N;
		while (chromoN1 == chromoN2)
		{
			chromoN2 = rand() % N;
		}
		for (j = 0; j < T; j++)
		{
			father[j] = rand_population[chromoN1][j];	//��
			mother[j] = rand_population[chromoN2][j];	//ĸ
		}
		//�ڼ��ʽ���
		double r = (double)rand() / RAND_MAX;
		if (r <= cross_rate)
		{
			/*�����ú���Order_Crossover��˳�򽻲�OX������������ֱ�Ӵ���cross_population��*/
			Order_Crossover(father, mother, Z1, Z2);
			Z1 += 2;//Ϊ��cross_population[Z1][]�����ڲ��ظ��ص���
			Z2 += 2;//Ϊ��cross_population[Z2][]�����ڲ��ظ��ص���
		}
		else
		{
			//������
			for (j = 0; j < T; j++)
			{
				cross_population[Z1][j] = father[j];
				cross_population[Z2][j] = mother[j];
			}
			Z1 += 2;//Ϊ��cross_population[Z1][]�����ڲ��ظ��ص���
			Z2 += 2;//Ϊ��cross_population[Z2][]�����ڲ��ظ��ص���
		}
	}

}
/*������Order_Crossover���壺˳�򽻲�OX����������ֱ�Ӵ���cross_population��*/
void Genetic_Algorithm::Order_Crossover(int* father, int* mother, int Z1, int Z2)
{
	int cut_point1, cut_point2;
	int child_indiv1[T] = { 0 };	//T�������Ӵ�1
	int child_indiv2[T] = { 0 };	//T�������Ӵ�2
	//��������ɽ����
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
	//����ȡ��������м䲿�ֵĻ���
	for (int x = cut_point1; x <= cut_point2; x++)
	{
		child_indiv1[x] = father[x];		 //�Ӵ�1ȡ�������жλ���
		child_indiv2[x] = mother[x];		 //�Ӵ�2ȡĸ�����жλ���
	}
	/*�����һ�����㲻���ڿ�ͷ*/
	if (cut_point1 != 0)
	{
		//���ȶ��Ӵ�1ȡ��β�λ���
		int index1 = 0;							 //�Ӵ�1������
		for (int y = 0; y < T; y++)				 //����ĸ������
		{
			bool bt = true;
			for (int z = cut_point1; z <= cut_point2; z++)
			{
				//�����Ӵ��жλ���Ա�
				if (mother[y] == child_indiv1[z])
				{
					bt = false;					 //�ı�boolֵ��ʹ�䲻ִ�к�������
					//�ѵ��ظ���������ǰ���forѭ����ִ�����£�����bt=false����������򲻻ᱻִ�У�
					break;
				}
			}
			//�������϶Աȣ���boolֵû�иı䣬���û������»����򸳸��Ӵ�
			if (bt == true)
			{
				child_indiv1[index1] = mother[y];//ĸ����ֵ���Ӵ�1
				index1 += 1;					 //�����Ӵ���һ������
				//����׶��Ӵ�����ֵ��ϣ�������β�λ���ֵ����
				if (index1 == cut_point1)
				{
					//����ڶ������㲻��ĩβ����ô���е����λ���
					if (cut_point2 != T - 1)
					{
						index1 = cut_point2 + 1;
					}
						//�ڶ�������������ĩβ����ô�˳�
					else
					{
						break;
					}
				}
			}
		}
		//���ٶ��Ӵ�2ȡ��β�λ���
		int index2 = 0;							 //�Ӵ�2������
		for (int y = 0; y < T; y++)				 //������������
		{
			bool bt = true;
			for (int z = cut_point1; z <= cut_point2; z++)
			{
				//�����Ӵ��жλ���Ա�
				if (father[y] == child_indiv2[z])
				{
					bt = false;					 //�ı�boolֵ��ʹ�䲻ִ�к�������
					//�ѵ��ظ���������ǰ���forѭ����ִ�����£�����bt=false����������򲻻ᱻִ�У�
					break;
				}
			}
			//�������϶Աȣ���boolֵû�иı䣬���û������»����򸳸��Ӵ�
			if (bt == true)
			{
				child_indiv2[index2] = father[y];//������ֵ���Ӵ�2
				index2 += 1;					 //�����Ӵ���һ������
				//����׶��Ӵ�����ֵ��ϣ�������β�λ���ֵ����
				if (index2 == cut_point1)
				{
					//����ڶ������㲻��ĩβ����ô���е����λ���
					if (cut_point2 != T - 1)
					{
						index2 = cut_point2 + 1;
					}
						//�ڶ�������������ĩβ��ֱ���˳�
					else
					{
						break;
					}
				}
			}
		}
		//�ݽ������Ӵ��������cross_population
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
		/*����׽������ڿ�ͷ*/
	else
	{
		/*����ڶ�����������Ҳ��β������ô�Ͳ��ý�����*/
		if (cut_point2 == T - 1)
		{
			//�ݽ�������ֱ�Ӵ���cross_population
			for (j = 0; j < T; j++)
			{
				cross_population[Z1][j] = father[j];
				cross_population[Z2][j] = mother[j];
			}
		}
			/*����ڶ������㲻��β��*/
		else
		{
			//���ȶ��Ӵ�1ȡ�ڶ��λ���
			int index1 = cut_point2 + 1;			 //�Ӵ�1������
			for (int y = 0; y < T; y++)				 //����ĸ������
			{
				bool bt = true;
				for (int z = cut_point1; z <= cut_point2; z++)
				{
					//�����Ӵ��жλ���Ա�
					if (mother[y] == child_indiv1[z])
					{
						bt = false;					 //�ı�boolֵ��ʹ�䲻ִ�к�������
						//�ѵ��ظ���������ǰ���forѭ����ִ�����£�����bt=false����������򲻻ᱻִ�У�
						break;
					}
				}
				//�������϶Աȣ���boolֵû�иı䣬���û������»����򸳸��Ӵ�
				if (bt == true)
				{
					child_indiv1[index1] = mother[y];//ĸ����ֵ���Ӵ�1
					index1 += 1;					 //�����Ӵ���һ������
					/*�����һ������û����*/
					if (index1 != T)
					{
						index1 += 0;
					}
						/*�����һ�����������ˣ����˳�*/
					else
					{
						break;
					}
				}
			}
			//���ٶ��Ӵ�2ȡ��β�λ���
			int index2 = cut_point2 + 1;			 //�Ӵ�2������
			for (int y = 0; y < T; y++)				 //������������
			{
				bool bt = true;
				for (int z = cut_point1; z <= cut_point2; z++)
				{
					//�����Ӵ��жλ���Ա�
					if (father[y] == child_indiv2[z])
					{
						bt = false;					 //�ı�boolֵ��ʹ�䲻ִ�к�������
						//�ѵ��ظ���������ǰ���forѭ����ִ�����£�����bt=false����������򲻻ᱻִ�У�
						break;
					}
				}
				//�������϶Աȣ���boolֵû�иı䣬���û������»����򸳸��Ӵ�
				if (bt == true)
				{
					child_indiv2[index2] = father[y];//������ֵ���Ӵ�2
					index2 += 1;					 //�����Ӵ���һ������
					/*�����һ������û����*/
					if (index2 != T)
					{
						index2 += 0;
					}
						/*�����һ�����������ˣ����˳�*/
					else
					{
						break;
					}
				}
			}
			//�ݽ������Ӵ��������cross_population
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
/*������Mutation���壺���ʱ��졷*/
void Genetic_Algorithm::Mutation()
{
	//�ٳ�ʼ��������ȺΪ������Ⱥ
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < T; j++)
		{
			muta_population[i][j] = cross_population[i][j];
		}
	}
	//�ڱ�������
	int point1, point2;					//�����
	for (i = 0; i < N; i++)
	{
		double r = { 0.0 };
		r = (double)rand() / RAND_MAX;	//���һ��������
		//����
		if (r <= muta_rate)
		{
			point1 = rand() % T;
			point2 = rand() % T;
			while (point1 == point2)
			{
				point2 = rand() % T;	//���⽻�����ظ�
			}
			swap(muta_population[i][point1], muta_population[i][point2]);
		}
	}
}
/*������Best_Solution���壺������ε�������������ý⣬������ʼ��Ⱥ�ͱ�����Ⱥ��*/
void Genetic_Algorithm::Best_Solution()
{
	int big_population[N + N][T];			//�ں�population��muta_population
	double big_fitness[N + N];				//�ں�population��muta_population��Ӧ��ֵ
	//���ں�population��muta_population
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
	//��ѡ��big_population����Ⱥ����÷���
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
		cout << "��" << it << "�����ŷ������Ա��γ�ʼ��Ⱥ��" << endl;
	}
	else
	{
		cout << "��" << it << "�����ŷ������Ա��α�����Ⱥ��" << endl;
	}
	for (j = 0; j < T; j++)
	{
		cout << big_population[min_index][j] << "��>";
	}
	cout << big_population[min_index][0] << endl;
	cout << "��ӦֵΪ��" << big_fitness[min_index] << endl;

}
/*������Mixing_population������Ϊ��һ����׼���µ���Ⱥ��
	���ɣ�70%���Ա�����Ⱥ��20%���Գ�ʼ��Ⱥ��10%�����²������塣*/
void Genetic_Algorithm::Mixing_population()
{
	int muta_num = round(0.7*N);					//������Ⱥѡ������
	int init_num = round(0.2*N);					//��ʼ��Ⱥѡ������
	int gnew_num = N - muta_num - init_num;			//�����ɸ�������
	int mix_index = 0;								//�����Ⱥ������
	double F1_copy[N];								//����F1[]
	copy(F1, F1 + N, F1_copy);
	double F2_copy[N];
	copy(F2, F2 + N, F2_copy);						//����F2[]
	double F3[N];									//���mix_popualtion[][]��Ӧ��ֵ

	/*��1���ӳ�ʼ��Ⱥѡ�����ŵ���*/
	//�ٴ�С�������Ӧֵ����
	double sort_F1[N];								//���F1[]�����ֵ
	int sort_F1_preindex[N];						//���temp_F1[]������ԭ������
	for (i = 0; i < N; i++)
	{
		sort_F1[i] = F1[i];//��ʼ��
	}
	//��С��������
	sort(sort_F1, sort_F1 + N);
	//�ڶ���������ӦֵѰ������֮ǰ��Ӧ���±�
	for (i = 0; i < N; i++)
	{
		//find()�������ص���ָ���ַ�����Լ�ȥ�׵�ַ�ɵ��±�������
		sort_F1_preindex[i] = (find(F1_copy, F1_copy + N, sort_F1[i]) - F1_copy);
		//���Ѿ��ѹ���ֵ�ų�
		F1_copy[sort_F1_preindex[i]] = 0;
	}
	/*
	for (i = 0; i < N; i++)
	{
		cout << sort_F1_preindex[i] << "  ";
	}
	cout << endl;
	*/
	//��ѡ��������ǰ�����������mix_population
	for (i = 0; i < init_num; i++)
	{
		//���Ƹ���
		for (j = 0; j < T; j++)
		{
			mix_population[mix_index][j] = population[sort_F1_preindex[i]][j];
		}
		//���Ƹ������Ӧ��ֵ
		F3[mix_index] = F1[sort_F1_preindex[i]];
		//�����Ⱥ��������һ
		mix_index += 1;
	}
	/*��2���ӱ�����Ⱥѡ�����ŵ���*/
	//�ٴ�С�������Ӧֵ����
	double sort_F2[N];								//���Ʊ�����Ⱥ��Ӧ��ֵF2[]
	int sort_F2_preindex[N];						//���temp_F2[]������ԭ������
	for (i = 0; i < N; i++)
	{
		sort_F2[i] = F2[i];
	}
	sort(sort_F2, sort_F2 + N);						//��С��������
	//�ڶ���������ӦֵѰ������֮ǰ��Ӧ���±�
	for (i = 0; i < N; i++)
	{
		sort_F2_preindex[i] = (find(F2_copy, F2_copy + N, sort_F2[i]) - F2_copy);
		//���Ѿ��ѹ���ֵ�ų�
		F2_copy[sort_F2_preindex[i]] = 0;
	}
	/*
	for (i = 0; i < N; i++)
	{
		cout << sort_F2_preindex[i] << "  ";
	}
	cout << endl;
	*/
	//��ѡ��������ǰ�����������mix_population
	for (i = 0; i < muta_num; i++)
	{
		//���Ƹ���
		for (j = 0; j < T; j++)
		{
			mix_population[mix_index][j] = muta_population[sort_F2_preindex[i]][j];
		}
		//���Ƹ������Ӧ��ֵ
		F3[mix_index] = F2[sort_F2_preindex[i]];
		//�����Ⱥ��������һ
		mix_index += 1;
	}
	/*��3�������ɸ��壬��������Ӧ��ֵ*/
	//��һ����������
	int temp_city[C];
	for (int i = 0; i < C; i++)
	{
		temp_city[i] = (i + 1);
	}
	//�ڴ��Һ������¸���
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
	/*��4����mix_population[][]���Ƹ�population[][]
		   ��F3[]���Ƹ�F1[]*/
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