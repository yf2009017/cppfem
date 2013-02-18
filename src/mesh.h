#ifndef MESH_H_
#define MESH_H_

#include "../eigen/Eigen/Core"

using namespace Eigen;

class Mesh {

    private:

        // Node coordinates (one per row)
        Matrix <Float, Dynamic, 2> COORD;

        // Node indexes in element (one element per row)
        Matrix <Int, Dynamic, 3> TOPO;

    public:

        // Constructor
        Mesh();

        



#endif
