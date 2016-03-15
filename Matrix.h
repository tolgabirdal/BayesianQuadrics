#ifndef Matrix_H
#define Matrix_H

#include "Eigen/Dense"
#include "Eigen/LU"
#include "Eigen/Eigenvalues"

#include <iostream>
#include <fstream>
#include <vector>
#include <list>

// These matrices are dynamic and slower than their
// static counterparts. The speed difference is in the allocation
// of memory on the heap.
namespace Cu
{ // The namespace is necessary to avoid conflicts with another library, so it is also quite short
    typedef Eigen::MatrixXd Matrix;
    typedef Eigen::VectorXd Vector;
    void print_matrix(const Cu::Matrix&, const std::string&);
    void print_matrix_octave(const Cu::Matrix&, const std::string&);
    void print_matrix_octave(const std::vector<std::vector<double>>& A, const std::string& name);
    void print_matrix_octave(const double* data, uint32_t M, uint32_t N, const std::string& name);
}

typedef std::pair< std::vector< std::vector<double> >, std::vector< std::vector<uint32_t> > > VerticesFaces;

Cu::Matrix ConvertToEigen(const double *in, const uint32_t M, const uint32_t N );

namespace MatrixFuncs
{
    void HomogeneousNormalise( Cu::Matrix& ); ///< Normalise a Cu::Cu::Matrix using homogeneous form
    std::vector<uint32_t> find( const Cu::Matrix& in, const std::function<bool(float)>& func  ); ///< Find indexes in the Cu::Cu::Matrix which satisfy the lambda function

    /**
     * @brief EigenValues
     * Return the eigenvectors and values of a matrix A. The eigenvectors are normalised by the determinant
     * of the matrix U, to ensure that the determinant is 1 and not -1. Either are valid solutions.
     * Only the real parts of the eigenvectors are returned.
     * @param A
     * @return
     */
    std::pair<Cu::Matrix, Cu::Matrix> EigenValues( const Cu::Matrix& A  );
}

namespace MeshFuncs
{
    void RemoveVerticesFaces(VerticesFaces& vf, std::vector<uint32_t> &VerticesToRemove );
        ///< Remove all faces and vertices from vf which contain the faces in FacesToRemove - note the second list is sorted
    void WritePly(const VerticesFaces& vf, const std::string filename);
}

namespace Funcs
{
    double normal_pdf(const double x, const double m, const double s);
}

#endif // Cu::Matrix_H
