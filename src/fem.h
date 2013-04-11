#ifndef FEM_H_
#define FEM_H_

#include <string>
#include <vector>

#include "typedefs.h"
#include "../eigen/Eigen/Core"

using namespace Eigen;


class Element {

    private:

        // Element index
        Int idx;

        // Element conductivity term value
        Float k;

        // Element source term value
        Float A;

        // Element node coordinate matrix
        Matrix <Float, 3, 2> Coord;

        // Element topology matrix
        Matrix <Int, 1, 3> Topo;

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

        // Shape functions
        Matrix <Float, 1, 3> shplint3(
                Matrix <Float, Dynamic, 1> r,
                Matrix <Float, Dynamic, 1> s);

        // Set local gradient matrix (dphi_l)
        void gradlt3();

        // Set Jacobian matrix (J)
        void jacobian();

        // Set global gradient matrix (dphi_g)
        void gradg();


    public:

        // Constructor
        Element(Int idx, Matrix <Float, 3, 2> Coord, Float k = 1.5, Float A = 1.0e-6);

        // Get element index
        Int getidx();

        // Find element stiffness matrix (K_e) and load vector (f_e)
        void findK_ef_e(int N_ip = 1);

        // Return element stiffness matrix (K_e)
        Matrix <Float, 3, 3> getK_e();

        // Return element load vector (f_e)
        Matrix <Float, 3, 1> getf_e();
};


class Mesh {

    private:

        // Number of elements
        Int N_e;

        // Number of nodes
        Int N;

        // Node coordinates (one per row)
        Matrix <Float, Dynamic, 2> COORD;

        // Node indexes in element (one element per row)
        Matrix <Int, Dynamic, 3> TOPO;

        // Stiffness matrix
        Matrix <Float, Dynamic, Dynamic> K;

        // Load vector
        Matrix <Float, Dynamic, 1> f;

        // Solution
        Matrix <Float, Dynamic, 1> T;

        // List of elements
        std::vector <Element> elementlist;
        
        
    public:

        // Constructor
        Mesh();
        
        // Initialize mesh, after COORD and TOPO are read
        void init();

        // Essential boundary condition (fixed val) on element i
        void ebc(const Int i, const Float val);

        // Natural boundary condition (fixed flux) on element pair i,j
        void nbc(const Int node_i, 
                const Int node_j,
                const Float flux,
                const Int N_ip = 1);

        // Find global stiffnes matrix (K) and load vector (f)
        void findKf();

        // Display stiffness matrix K
        void printK();

        // Display load vector f
        void printf();

        // Solve the system in the steady state
        void steadystate();

        // Triangle interfacing, should be called after "./triangle A.poly"
        // Read A.1.ele, containing topology data
        void readTriangleEle(const std::string filename);
        // Read A.1.poly, containing topology data
        void readTriangleNode(const std::string filename);

        // Write temperatures
        void writeT(const std::string filename);
};

#endif
