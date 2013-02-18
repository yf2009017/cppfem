#include "../eigen/Eigen/Core"
#include "../eigen/Eigen/LU"
#include "../eigen/Eigen/Geometry"

#include "fem.h"
#include "gaussquad.h"
#include "typedefs.h"

using namespace Eigen;


// Node class constructor
Node::Node(Matrix <Float, 1, 2> coord, Float T)
    : coord(coord), T(T) { }

// Get node value
Float Node::getT()
{
    return T;
}

// Set node value
Float Node::setT(Float T_new)
{
    T = T_new;
}



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

//void Element::updateCoord() { }

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

// Return node list
std::list <Node> getnodes()
{
    return nodelist;
}





// Mesh class constructor
Mesh::Mesh() { }

// Add element to mesh
void Mesh::addelement(Element element)
{
    elementlist.push_back(element);
}



