#include "sorting.h"
#include <iostream>

using namespace std;

int main()
{
    sorting sort(10);

    sort.selection_sorting();
	cout << endl;
	sort.insertion_sorting();
	cout << endl;
	system("pause");

    return 0;
}
