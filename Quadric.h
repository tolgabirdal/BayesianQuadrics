#ifndef QUADRIC_H
#define QUADRIC_H

#include "stdint.h"
#include <vector>

#include "Matrix.h"

namespace Geometry
{
    class Quadric
    {
    public:
        class Parameters
        {
        public:
            double m_A[3][3];
            double m_b[3];
            double m_c;
        };

        struct FastDataStore // Some matrices which have been allocated for fast computation of the quadric
        {
            Eigen::Matrix<double,13,13> Q, eye;
            Eigen::Matrix<double,13,1> z, sol;
            double Z[13];
            double a;
        };

        Quadric();
        Quadric(const double* data, const uint32_t N, const double sigma);
        Quadric(const Parameters& params);

        Parameters GetParameters() const; ///< Returns the parameters
        double Fit(const double* data, const uint32_t N, const double sigma); ///< Fit to the data
        void FitFast( const double* data, const uint32_t N, const double sigma, FastDataStore* fds );

        VerticesFaces GetMeshRep(const uint32_t N, const double *data, const uint32_t NPts) const; ///< Get a mesh representation
        void GetBoundedMeshRep(VerticesFaces& vf, const double *data, const uint32_t N) const;
        std::vector<bool> GetPointCloudBounds(const uint32_t StartN, const double* data, const uint32_t N) const;
        VerticesFaces GetPrincipleCylinder(const uint32_t N, const double *data, const uint32_t NPts) const; ///< Return a  mesh representation for the principle cylinder
    private:
        VerticesFaces GenerateUnitCircle(const uint32_t N ) const; ///< Create a unit circle
        VerticesFaces GenerateUnitCylinder(const uint32_t N , const uint32_t principleDimension) const; ///< Create a unit cylinder

        VerticesFaces GenerateEllipsoid( const uint32_t N, const double *data, const uint32_t NPts ) const; ///< Create an ellipsoid
        VerticesFaces GeneratePrincipleCylinder( const uint32_t N, const double *data, const uint32_t NPts) const; ///< Create an ellipsoid

        Parameters m_params;
        
    };
}

#endif // QUADRIC_H

