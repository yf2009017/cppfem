#include <iostream>
//#include "typedefs.h"
#include "fem.h"

int main()
{
    Mesh temp;
    temp.readTriangleEle("../triangle/A.1.ele");
    temp.readTriangleNode("../triangle/A.1.node");

    return 0;
}

