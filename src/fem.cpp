#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>

#include "../eigen/Eigen/Core"
#include "../eigen/Eigen/LU"
#include "../eigen/Eigen/Geometry"

#include "fem.h"
#include "gaussquad.h"
#include "typedefs.h"

using namespace Eigen;



////////// ELEMENT ///////////

// Element class constructor
Element::Element(Int idx, Matrix <Float, 3, 2> Coord, Float k, Float A)
    : idx(idx), Coord(Coord), k(k), A(A) { }

// Get element index
Int Element::getidx()
{ return idx; }

// Linear shape function for 3-node triangle (eq. 17)
Matrix <Float, 1, 3> Element::shplint3(
        Matrix <Float, Dynamic, 1> r,
        Matrix <Float, Dynamic, 1> s)
{
    Matrix <Float, Dynamic, 3> m;
    m.resize(r.size(), 3);
    for (int row=0; row<r.size(); ++row) {
        m(row, 0) = 1.0 - r(row) - s(row);
        m(row, 1) = r(row);
        m(row, 2) = s(row);
    }
    return m;
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
void Element::findK_ef_e(int N_ip)
{
    //std::cout << "\n---------------------------------------------" << std::endl;
    //std::cout << "element " << idx << " Coord = \n" << Coord << std::endl;

    // Get Gauss quadrature data [r_i, s_i, w_i], one row per integration point
    Matrix <Float, Dynamic, 3> g_i = gausst3(N_ip);
    //std::cout << "element " << idx << " g_i = \n" << g_i << std::endl;

    // Find local derivatives of the shape function
    gradlt3();
    //std::cout << "element " << idx << " dphi_l = \n" << dphi_l << std::endl;

    // Get interpolation function values
    Matrix <Float, Dynamic, 3> phi = shplint3(g_i.col(0), g_i.col(1));
    //std::cout << "element " << idx << " phi = \n" << phi << std::endl;

    // Find the Jacobian and it's determinant
    jacobian();
    //std::cout << "element " << idx << " J = \n" << J << std::endl;
    //std::cout << "element " << idx << " detJ = \n" << detJ << std::endl;

    // Find global gradients of the Jacobian
    gradg();
    //std::cout << "element " << idx << " dphi_g = \n" << dphi_g << std::endl;

    K_e.fill(0.0);
    f_e.fill(0.0);

    // Loop over integration points
    for (int i = 0; i<N_ip; ++i) {
        K_e += g_i(i,2) * k * dphi_l.transpose() * dphi_l * detJ;
        f_e += g_i(i,2) * phi.transpose() * A * detJ;
    }

    //std::cout << "element " << idx << " K_e = \n" << K_e << std::endl;
    //std::cout << "element " << idx << " f_e = \n" << f_e << std::endl;
}

// Return element stiffness matrix (K_e)
Matrix <Float, 3, 3> Element::getK_e()
{ return K_e; }

// Return element load vector (f_e)
Matrix <Float, 3, 1> Element::getf_e()
{ return f_e; }


////////// MESH ///////////

// Mesh class constructor
Mesh::Mesh() { }

// Initialize mesh with a predefined grid
void Mesh::init()
{
    Matrix <Int, 1, 3> Topo;
    Matrix <Float, 3, 2> Coord;
    // Create each element
    for (Int ie=0; ie<N_e; ++ie) {
        Element* e;
        Topo = TOPO.row(ie);
        Coord.row(0) = COORD.row(Topo(0));
        Coord.row(1) = COORD.row(Topo(1));
        Coord.row(2) = COORD.row(Topo(2));
        //std::cout << "element " << ie << " Coord = \n" << Coord << std::endl;
        e = new Element(ie, Coord);
        elementlist.push_back(*e);
    }
}

// Find global stiffness matrix (K) and load vector (f)
// by performing a volumetric integration of the elements
void Mesh::findKf()
{
    K.resize(N, N);
    f.resize(N, 1);
    K.setZero(N, N);
    f.setZero(N, 1);

    std::vector<Element>::iterator e;
    for (e = elementlist.begin(); e != elementlist.end(); ++e) {

        e->findK_ef_e();
        Matrix <Float, 3, 3> K_e = e->getK_e();
        Matrix <Float, 3, 1> f_e = e->getf_e();

        Matrix <Int, 1, 3> nodes = TOPO.row(e->getidx());
        
        for (int i=0; i<3; ++i)
            f(nodes(i)) += f_e(i);

        int i = 0;
        int r, c;
        for (int j=0; j<3; ++j) {
            for (int k=0; k<3; ++k) {
                r = nodes(j);
                c = nodes(k);
                K(r,c) += K_e(i++);
            }
        }

    }
    //std::cout << K << std::endl;
    //std::cout << f << std::endl;
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
    T = K.fullPivLu().solve(f);     // Speed: -, Accuracy: +++
    //T = K.householderQr().solve(f); // Speed: ++, Accuracy: +
}

// Read triangle topology file (*.ele)
void Mesh::readTriangleEle(const std::string filename)
{
    FILE* fp;
    if ((fp = fopen(filename.c_str(), "r")) == NULL) {
        std::cerr << "Error opening " << filename << std::endl;
    }

    int nd, boundarymarkers, n1, n2, n3;
    long int tmp;
    fscanf(fp, "%ld %d %d\n", &N_e, &nd, &boundarymarkers);
    TOPO.resize(N_e, 3);

    for (Int i=0; i<N_e; ++i) {
        fscanf(fp, "%ld %d %d %d\n", &tmp, &n1, &n2, &n3);
        TOPO(i,0) = n1-1;
        TOPO(i,1) = n2-1;
        TOPO(i,2) = n3-1;
    }

    //std::cout << "TOPO = \n " << TOPO << std::endl;
    fclose(fp);
}

// Read triangle coordinate file (*.node)
void Mesh::readTriangleNode(const std::string filename)
{
    FILE* fp;
    if ((fp = fopen(filename.c_str(), "r")) == NULL) {
        std::cerr << "Error opening " << filename << std::endl;
    }

    float fx, fy, tmpf;
    int nd, nobound, bound, tmp;
    long int tmpl;
    fscanf(fp, "%ld %d %ld %d", &N, &nd, &tmpl, &nobound);
    COORD.resize(N, 2);

    for (Int i=0; i<N; ++i) {
        fscanf(fp, "%4ld  %f %f  %f    %d\n", &tmpl, &fx, &fy, &tmpf, &bound);
        COORD(i,0) = fx;
        COORD(i,1) = fy;
    }

    fclose(fp);
}

// Print stiffness matrix K
void Mesh::printK()
{
    std::cout << "K = \n" << K << std::endl;
}
        
// Print load vector f
void Mesh::printf()
{
    std::cout << "f = \n" << f << std::endl;
}
        
// Write temperatures to file
void Mesh::writeT(const std::string filename)
{
    std::ofstream ofs(filename.c_str());
    if (!ofs) {
        std::cerr << "Error opening " << filename << std::endl;
    }

    ofs << T;

    if (ofs.is_open())
        ofs.close();
}

