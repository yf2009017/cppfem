#include <iostream>

//#include "mesh.h"
#include "typedefs.h"
#include "gaussquad.h"

using namespace Eigen;
using namespace std;

int main()
{
    Matrix <Float, 1, 3> m;
    m = gausst3(1);

    cout << m;

    return 0;
}

