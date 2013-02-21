#include <iostream>
#include "fem.h"

int main()
{
    Mesh temp;
    temp.readTriangleEle("triangle/A.1.ele");
    temp.readTriangleNode("triangle/A.1.node");
    temp.init();

    temp.findKf();

    temp.ebc(0, 3.0);
    //temp.nbc(..);

    temp.steadystate();
    temp.writeT("A.1.T");

    return 0;
}

