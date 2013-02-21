#include <iostream>
#include "fem.h"

int main()
{
    Mesh temp;
    temp.readTriangleEle("triangle/A.1.ele");
    temp.readTriangleNode("triangle/A.1.node");
    temp.steadystate();
    temp.writeT("A.1.T");

    return 0;
}

