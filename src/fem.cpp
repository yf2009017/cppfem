#include <iostream>
#include <cstdio>

#include "../eigen/Eigen/Core"
#include "../eigen/Eigen/LU"
#include "../eigen/Eigen/Geometry"

#include "fem.h"
#include "gaussquad.h"
#include "typedefs.h"

using namespace Eigen;

////////// NODE ///////////

// Node class constructor
Node::Node(Int idx, Matrix <Float, 1, 2> coord, Float T)
    : idx(idx), coord(coord), T(T) { }

// Get node index
Int Node::getidx()
{ return idx; }

// Get node value
Float Node::getT()
{ return T; }

// Set node value
void Node::setT(Float T_new)
{ T = T_new; }

// Get node coordinate
Matrix <Float, 1, 2> Node::getcoord ()
{ return coord; }


////////// ELEMENT ///////////

// Element class constructor
Element::Element(Float k, Float A)
    : k(k), A(A) { }

// Add node to element
void Element::addnode(Node node)
{
    if (nodelist.size() > 2)
        throw "Error, element full of nodes";
    nodelist.push_back(node);
}

// Linear shape function for 3-node triangle (eq. 17)
Matrix <Float, 1, 3> Element::shplint3(
        Matrix <Float, Dynamic, 1> r,
        Matrix <Float, Dynamic, 1> s)
{
    Matrix <Float, Dynamic, 3> m;
    m.row(0) = 1.0 - r.array() - s.array();
    m.row(1) = r;
    m.row(2) = s;
    return m;
}

// Update element node coordinate matrix
void Element::updateCoord() 
{
    int i = 0;
    std::list<Node>::iterator n = nodelist.begin();
    while (n != nodelist.end()) {
        Coord.row(i++) = n->getcoord();
    }
}

// Update element topograpy matrix
void Element::updateTopo()
{
    int i = 0;
    std::list<Node>::iterator n = nodelist.begin();
    while (n != nodelist.end()) {
        Topo(0, i++) = n->getidx();
    }
}

// Sets local gradient matrix (eq. 23)
void Element::gradlt3()
{
    dphi_l << -1.0, 1.0, 0.0,
      -1.0, 0.0, 1.0;
}

// Sets the Jacobian matrix (eq. 26)
void Element::jacobian()
{
    J = dphi_l * Coord;
    detJ = J.determinant();
}

// Sets the global gradient matrix (eq. 22)
void Element::gradg()
{
    dphi_g = J.inverse() * dphi_l;
}

// Find element load vector and stiffness matrix (eqs. 38 and 39)
void Element::loadandstiffness(int N_ip)
{
    // Get Gauss quadrature data [r_i, s_i, w_i], one row per integration point
    Matrix <Float, Dynamic, 3> g_i = gausst3(N_ip);

    // Find local derivatives of the shape function
    gradlt3();

    // Get interpolation function values
    Matrix <Float, Dynamic, 3> phi = shplint3(g_i.row(0), g_i.row(1));

    // Find the Jacobian and it's determinant
    jacobian();

    // Find global gradients of the Jacobian
    gradg();

    K_e.fill(0.0);
    f_e.fill(0.0);

    // Loop over integration points
    for (int i = 0; i<N_ip; ++i) {
        K_e += g_i(i,2) * k * dphi_l.transpose() * dphi_l * detJ;
        f_e += g_i(i,2) * phi.transpose() * A * detJ;
    }
}

// Return node index vector
std::list <Node> Element::getnodes()
{ return nodelist; }
        
// Return element topology matrix
Matrix <Int, 1, 3> Element::getTopo()
{
    updateTopo();
    return Topo;
}

// Return element coordinate matrix
Matrix <Float, 3, 2> Element::getCoord()
{
    updateCoord();
    return Coord;
}



////////// MESH ///////////

// Mesh class constructor
Mesh::Mesh() { }

// Add element to mesh
void Mesh::addelement(Element element)
{ elementlist.push_back(element); }

// Update topology matrix by calling all elements, 
// who in turn call all their nodes
void Mesh::updateTOPO()
{
    int i = 0;
    std::list<Element>::iterator e = elementlist.begin();
    while (e != elementlist.end()) {
        TOPO.row(i++) = e->getTopo();
    }
}

// Update coordinate matrix by calling all elements, 
// who in turn call all their nodes
void Mesh::updateCOORD()
{
    int i = 0;
    std::list<Element>::iterator e = elementlist.begin();
    while (e != elementlist.end()) {
        Matrix <Float, 3, 2> Coord = e->getCoord();
        for (int j=0; j<3; ++j)
            COORD.row(i++ * 3 + j) = Coord.row(j);
    }
}

// Enforce essential (Dirichlet) boundary condition on nodes
// with index i to value val.
// Using algorithm displayed in fig. 5, p. 12
void Mesh::ebc(const Int i, const Float val)
{
    f -= val * K.col(i);
    f(i) = val;
    K.row(i).fill(0.0);
    K.col(i).fill(0.0);
    K(i,i) = 1.0;
}

// Enforce natural (Neumann) boundary condition to node pair
void Mesh::nbc(
        const Int node_i, 
        const Int node_j,
        const Float val,
        const Int N_ip)
{
    //Node n_i = 
    //Node n_j = 
    //Matrix <Float, 1, 2> coord_i = n_i.getcoord();
    //Matrix <Float, 1, 2> coord_j = n_j.getcoord();
    Matrix <Float, 1, 2> coord_i = COORD.row(node_i);
    Matrix <Float, 1, 2> coord_j = COORD.row(node_j);
    Matrix <Float, 1, 2> x_ij = coord_j - coord_i;
    Matrix <Float, 2, 1> f_new;
    Float dl = sqrt(x_ij(0)*x_ij(0) + x_ij(1)*x_ij(1));

    Matrix <Float, Dynamic, 2> gi = gauss1d(N_ip);
    //Matrix <Float, 1, 2> p;
    Matrix <Float, 2, 1> p;

    for (int ip = 0; ip<N_ip; ++ip) {
        p << 1.0+gi(ip,0), 1.0-gi(ip,0);
        f_new = gi(ip,1) * 0.5 * p * val * dl/2.0;
    }

    f(node_i) += f_new(0);
    f(node_j) += f_new(1);
}


// Solve the steady state problem of the system
// See http://eigen.tuxfamily.org/dox/TutorialLinearAlgebra.html 
// for information on the different decompositions/solvers
void Mesh::steadystate()
{ 
    findKf();
    T = K.fullPivLu().solve(f);     // Speed: -, Accuracy: +++
    //T = K.householderQr().solve(f); // Speed: ++, Accuracy: +
}

// Read triangle topology file (*.ele)
void Mesh::readTriangleEle(const std::string filename)
{
    FILE* fp;
    if ((fp = fopen(filename.c_str(), "r")) == NULL) {
        std::cerr << "Error opening " << filename << std::endl;
        throw "Read error in readTriangleEle";
    }

    int nd, boundarymarkers, n1, n2, n3;
    long int tmp;
    fscanf(fp, "%ld  %d  %d", &N_e, &nd, &boundarymarkers);
    TOPO.resize(N_e, 3);

    for (int i=0; i<N_e; ++i) {
        fscanf(fp, "%4ld    %4d  %4d  %4d", &tmp, &n1, &n2, &n3);
        TOPO(i,0) = n1;
        TOPO(i,1) = n2;
        TOPO(i,2) = n3;
    }

    fclose(fp);
}


// Read triangle coordinate file (*.node)
void Mesh::readTriangleNode(const std::string filename)
{
    FILE* fp;
    if ((fp = fopen(filename.c_str(), "r")) == NULL) {
        std::cerr << "Error opening " << filename << std::endl;
        throw "Read error in readTriangleNode";
    }

    Float x, y;
    int nd, nobound;
    long int tmp;
    fscanf(fp, "%ld  %d  %d  %d", &N, &nd, &tmp, &nobound);
    COORD.resize(N, 2);

    for (int i=0; i<N_e; ++i) {
        fscanf(fp, "%4d    %.17g  %.17g", &tmp, &x, &y);
        COORD(i,0) = x;
        COORD(i,1) = y;
    }

    fclose(fp);
}
        
// Write temperatures to file
void Mesh::writeT(const std::string filename)
{
    FILE* fp;
    if ((fp = fopen(filename.c_str(), "w")) == NULL) {
        std::cerr << "Error opening " << filename << std::endl;
        throw "Read error in writeT";
    }

    std::cout << T;

    fclose(fp);
}

