//
// Created by 30205 on 2022/5/9.
//

#ifndef SORING_SORTING_H
#define SORING_SORTING_H

// sorting.h: 标准系统包含文件的包含文件
// 或项目特定的包含文件。

#pragma once
#include <iostream>

using namespace std;

/*随机整数*/
int int_random(int x)
{
	return (int)(rand() << 8 ^ rand()) % x;
}
/*随机小数*/
double double_random(double low, double high)
{
	return low + (double)rand() / RAND_MAX * (high - low);
}

class sorting
{
 public:
	sorting(int number) : N(number)
	{
		array = new int[N];
		for (int i = 0; i < N; i++)
		{
			array[i] = int_random(N);
		}
	}
	void print();
	void swap(int& t_1, int& t_2);
	void selection_sorting();
	void insertion_sorting();
	void shell_sorting();
 private:
	int* array;
	int N;
};

void sorting::print()
{
	for (int i = 0; i < N; i++)
	{
		cout << array[i] << " ";
	}
	cout << endl;
}
void sorting::swap(int& t_1, int& t_2)
{
	int temp;
	temp = t_1;
	t_1 = t_2;
	t_2 = temp;
}
void sorting::selection_sorting()
{
	cout<< "selection_sorting" << endl;
	sorting::print();
	for (int i = 0; i < N - 1; i++)
	{
		int min_value_position = i;

		for (int j = i + 1; j < N; j++)
		{
			min_value_position = array[j] < array[min_value_position] ? j : min_value_position;
		}
		sorting::swap(array[i], array[min_value_position]);
	}
	sorting::print();
}

void sorting::insertion_sorting()
{
	cout<< "insertion_sorting" << endl;
	sorting::print();
	for (int i = 1; i < N - 1; i++)
	{
		for (int j = i; j > 0; j--)
		{
			if(array[j] < array[j - 1])
			{
				sorting::swap(array[j], array[j - 1]);
			}
		}
	}
	sorting::print();
}
#endif //SORING_SORTING_H
