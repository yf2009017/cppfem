#ifndef FEM_H_
#define FEM_H_

#include <list>

#include "typedefs.h"
#include "../eigen/Eigen/Core"

using namespace Eigen;


class Node {

    private:

        // Node coordinate
        Matrix <Float, 1, 2> coord;

        // Node value
        Float T;



    public:

        // Constructor
        Node(Matrix <Float, 1, 2> coord, Float T = 0.0);

        // Get node value
        Float getT();

        // Set node value
        void setT(Float T);

};


class Element {

    private:

        // Element conductivity term value
        Float k;

        // Element source term value
        Float A;

        // Element node coordinate matrix
        Matrix <Float, 3, 2> Coord;

        // Element local gradient matrix
        Matrix <Float, 2, 3> dphi_l;

        // Element Jacobian matrix
        Matrix <Float, 2, 2> J;
        Float detJ;

        // Element global gradient matrix
        Matrix <Float, 2, 3> dphi_g;

        // Element load vector
        Matrix <Float, 3, 1> f_e;

        // Element stiffness matrix
        Matrix <Float, 3, 3> K_e;

        // List of node indexes
        std::list <Node> nodelist;

        // Add node
        void addnode(Node node);

        // Shape functions
        Matrix <Float, 1, 3> shplint3(
                Matrix <Float, Dynamic, 1> r,
                Matrix <Float, Dynamic, 1> s);

        // Update element coordinate matrix (Coord)
        void updateCoord();

        // Set local gradient matrix (dphi_l)
        void gradlt3();

        // Set Jacobian matrix (J)
        void jacobian();

        // Set global gradient matrix (dphi_g)
        void gradg();

        // Set the element load vector and stiffness matrix (K and f)
        void loadandstiffness(int N_ip = 1);




    public:

        // Constructor
        Element(Float k = 1.5, Float A = 1.0e-6);

        // Return node list
        std::list getnodes();

};


class Mesh {

    private:

        // Node coordinates (one per row)
        Matrix <Float, Dynamic, 2> COORD;

        // Node indexes in element (one element per row)
        Matrix <Int, Dynamic, 3> TOPO;

        // List of elements
        std::list <Element> elementlist;
        
        // Add element
        void addelement(Element element);
        
        // Iterate over elements to update the TOPO matrix
        // updateTOPO()


        // Iterate over elements (who iterate over nodes)
        // to update the COORD matrix
        // updateCOORD()


    public:

        // Constructor
        Mesh();

};

#endif
