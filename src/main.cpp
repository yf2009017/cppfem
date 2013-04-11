#include <iostream>
#include "fem.h"

int main()
{
    Mesh temp;
    temp.readTriangleEle("triangle/A.1.ele");
    temp.readTriangleNode("triangle/A.1.node");
    
    temp.init();

    std::cout << "findKf" << std::endl;
    temp.findKf();

    std::cout << "ebc" << std::endl;
    temp.ebc(0, 3.0);
    //temp.nbc(..);

    //temp.printK();
    //temp.printf();

    std::cout << "steadystate" << std::endl;
    temp.steadystate();
    temp.writeT("A.1.T");

    return 0;
}

