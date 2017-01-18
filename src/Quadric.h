#ifndef QUADRIC_H
#define QUADRIC_H

#include "stdint.h"
#include <vector>

#include "Matrix.h"

/**
 * \brief Functions and classes relating to computational geometry
 */
namespace Geometry
{
    /**
    * \brief A class representing a quadric. 
    * This class was written for fitting a quadric to a noisy point cloud with a Bayesian prior. The problem
    * was that when there is missing data, using a standard fit can produce an unwanted surface, such as a 
    * hyperboloid. If there is some extra knowledge indicating the type of quadric that the fit should be,
    * one can tune the prior to encourage the surface to take on a particular type.
    * There are more details in the paper,
    *   <a href="http://link.springer.com/article/10.1007/s41095-016-0041-9">"Fitting quadrics with a Bayesian prior" - Daniel Beale et al. Journal of Computational Visual Media, 2016.</a>
    */
    class Quadric
    {
    public:
        /**
        * \brief A class containing the parameters of the quadric.
        */
        class Parameters
        {
        public:
            double m_A[3][3];
            double m_b[3];
            double m_c;
        };

        /**
        * \brief This struct provides an internal datastore for the quadric fitting algorithm, when fitting 
        * multiple quadrics points. It is pre-allocated to avoid uneccessary allocation on the heap (since 
        * the Eigen matrices are dynamic).
        */
        struct FastDataStore
        {
            Eigen::Matrix<double,13,13> Q, eye;
            Eigen::Matrix<double,13,1> z, sol;
            double Z[13];
            double a;
        };

        /**
        * \brief Construct an empty quadric
        */
        Quadric();

        /**
        * \brief Construct a quadric given a data matrix. 
        * The parameters will be estimated on construction.
        * @param data A pointer to the first element of a column major contiguous array of size 3N
        * @param N the number of columns, or datapoints.
        * @param sigma The hyper parameter
        */
        Quadric(const double* data, const uint32_t N, const double sigma);

        /**
        * \brief Construct a quadric from the known parameters
        * @param The parameters of the quadric
        */
        Quadric(const Parameters& params);

        /**
        * \brief Return the parameters of the quadric
        * @return The parameters
        */
        Parameters GetParameters() const;

        /**
        * \brief Fit the quadric to a data matrix
        * @param data A pointer to the first element of a column major contiguous array of size 3N
        * @param N the number of columns, or datapoints.
        * @param sigma The hyper parameter
        */
        double Fit(const double* data, const uint32_t N, const double sigma);

        /**
        * \brief Fit the quadric to a data matrix given a pre allocated FastDataStore.
        * In order to execute this method, one would do the following,
        * \code 
        * FastDataStore *fds = new FastDataStore();
        * for(uint32_t k=0; k<NExecs; ++k)
        * {
        *   quadric[k].FitFast(data[k], N[k], sigma[k], fds);
        * }
        * \endcode
        * This FitFast method will run without any heap allocations. 
        * @param data A pointer to the first element of a column major contiguous array of size 3N
        * @param N the number of columns, or datapoints.
        * @param sigma The hyper parameter
        */
        void FitFast( const double* data, const uint32_t N, const double sigma, FastDataStore* fds );

        /**
        * \brief Return a mesh representing the quadric.
        * The method takes a canonical surface such as an Ellipsoid and transforms it 
        * using the parameters of the quadric. Only an Ellipsoid has been implemented
        * for demonstration purposes.
        * @param N Generate N^2 vertices on the mesh
        * @param data A pointer to the colum major, contigouous, data matrix (with 3NPts elements).
        * @param NPts The number of points in the data matrix
        */
        VerticesFaces GetMeshRep(const uint32_t N, const double *data, const uint32_t NPts) const; 
  
        /**
        * \brief Return a bounded mesh representing the quadric.
        * The method takes a canonical surface such as an Ellipsoid and transforms it 
        * using the parameters of the quadric. Only an Ellipsoid has been implemented
        * for demonstration purposes. The mesh is then truncated so that the vertices on the 
        * input point cloud do not fall outside the mesh.
        * @param N Generate N^2 vertices on the mesh
        * @param data A pointer to the colum major, contigouous, data matrix (with 3NPts elements).
        * @param NPts The number of points in the data matrix
        */
        void GetBoundedMeshRep(VerticesFaces& vf, const double *data, const uint32_t N) const;

        /**
        * \brief Return a vector of booleans indicating which points on the generated mesh 
        * are outside the maximum and minimum point cloud bounds.
        * @param StartN Generate N^2 vertices on the mesh
        * @param data A pointer to the colum major, contigouous, data matrix (with 3NPts elements).
        * @param N The number of points in the data matrix
        * @return A vector of bools specifying which vertices are inside the max/min bounds
        */
        std::vector<bool> GetPointCloudBounds(const uint32_t StartN, const double* data, const uint32_t N) const;

        /**
        * \brief Return vertices and faces of the quadric, assuming that it is a cylinder.
        * For a quadric cylinder one of the eigenvalues of the parameter A will be zero.
        * @param N Generate N^2 vertices on the mesh
        * @param data A pointer to the colum major, contigouous, data matrix (with 3NPts elements).
        * @param NPts The number of points in the data matrix
        * @return A cylinder mesh
        */
        VerticesFaces GetPrincipleCylinder(const uint32_t N, const double *data, const uint32_t NPts) const; 
    private:

        /**
        * \brief Create a unit circle with N^2 vertices
        */
        VerticesFaces GenerateUnitCircle(const uint32_t N ) const;

        /**
         * \brief Create a unit cylinder with N^2 vertices
         * @param Create N^2 vertices
         * @param The dimension in which the cylinder is unbounded.
         * @return A vertices / faces object representing the mesh
         */
        VerticesFaces GenerateUnitCylinder(const uint32_t N , const uint32_t principleDimension) const;

        /**
        * \brief Generate an ellipsoid from the parameters
        * @param N Generate N^2 vertices on the mesh
        * @param data A pointer to the colum major, contigouous, data matrix (with 3NPts elements).
        * @param NPts The number of points in the data matrix
        */
        VerticesFaces GenerateEllipsoid( const uint32_t N, const double *data, const uint32_t NPts ) const;

        /**
        * \brief Generate a cylinder from the parameters
        * @param N Generate N^2 vertices on the mesh
        * @param data A pointer to the colum major, contigouous, data matrix (with 3NPts elements).
        * @param NPts The number of points in the data matrix
        */
        VerticesFaces GeneratePrincipleCylinder( const uint32_t N, const double *data, const uint32_t NPts) const; ///< Create an ellipsoid

        Parameters m_params; ///< The quadric parameters
        
    };
}

#endif // QUADRIC_H

