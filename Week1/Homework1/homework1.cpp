// This program sums the N elements of a vector.

#include <iostream>
#include <vector>
using namespace std;
const int N = 40;

//Returns the result of summing all the elements of the vector
template <class summable>
summable sum(const vector<summable> vect, summable accum = 0)
{
    for(int i = 0; i < vect.size(); ++i)
        accum += vect[i];
    return accum;
}

int main()
{
    vector<int> data(N);
    for(int i = 0; i < data.size(); ++i)
        data[i] = i;
    cout << "sum is " << sum(data) << endl;
}
