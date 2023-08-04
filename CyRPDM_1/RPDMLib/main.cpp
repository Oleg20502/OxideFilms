#include <iostream>

#include "film.hpp"


int main()
{
    Film A;

    int b = A.init();
    std::cout << b << '\n';

    return 0;
}