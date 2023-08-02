#include <iostream>
#include <vector>

using std::vector;

int main()
{
    vector<vector<double>> a;
    vector<vector<double>> b;

    b.resize(100);
    for (int i = 0; i < 100; ++i) {
        b[i].resize(10);
    }


    a = vector<vector<double>> (10000, vector<double>(10000, 0));
    for (int i = 0; i < 100; ++i) {
        a = vector<vector<double>> (10000+i, vector<double>(10000, 0));;
        std::cout << i << '\n';
    }
    std::cout << "Finished\n";

}