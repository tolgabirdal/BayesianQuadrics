#include "Matrix.h"

Cu::Matrix ConvertToEigen(const double *in, const uint32_t M, const uint32_t N )
{
    Cu::Matrix out(M,N);
    out = Eigen::Map<Cu::Matrix>(const_cast<double *>(in), M, N);
    return out;
}

namespace Cu
{
void print_matrix_octave(const Cu::Matrix& A, const std::string& name)
{
    std::ofstream ostream;
    ostream.open(name);
    ostream << " # Created by Curve Fit" << std::endl;
    ostream << " # name: " << name << std::endl;
    ostream << " # type: matrix " << std::endl;
    ostream << " # rows : " << A.rows() << std::endl;
    ostream << " # columns : " << A.cols() << std::endl;
    for( uint32_t i=0; i<A.rows(); i++ )
    {
        for( uint32_t j=0; j<A.cols(); j++ )
        {
            ostream << " " << A(i,j);
        }
        ostream << std::endl;
    }
     ostream.close();
}

void print_matrix_octave(const std::vector<std::vector<double> > &A, const std::string& name)
{
    std::ofstream ostream;
    ostream.open(name);
    ostream << " # Created by Curve Fit" << std::endl;
    ostream << " # name: " << name << std::endl;
    ostream << " # type: matrix " << std::endl;
    ostream << " # rows : " << A[0].size() << std::endl;
    ostream << " # columns : " << A.size() << std::endl;
    for( uint32_t i=0; i<A.size(); i++ )
    {
        for( uint32_t j=0; j<A[0].size(); j++ )
        {
            ostream << " " << A[i][j];
        }
        ostream << std::endl;
    }
    ostream.close();
}

void print_matrix_octave(const double* data, uint32_t M, uint32_t N, const std::string& name)
{
    std::ofstream ostream;
    ostream.open(name);
    ostream << " # Created by Curve Fit" << std::endl;
    ostream << " # name: " << name << std::endl;
    ostream << " # type: matrix " << std::endl;
    ostream << " # rows : " << M << std::endl;
    ostream << " # columns : " << N << std::endl;
    for( uint32_t j=0; j<M; j++ )
    {
        for( uint32_t i=0; i<N; i++ )
        {
            ostream << " " << data[M*i+j];
        }
        ostream << std::endl;
    }
    ostream.close();
}

#ifdef DEBUG
    void print_matrix(const Cu::Matrix& P, const std::string& name)
    {
        std::cout << name << std::endl;
        for( int i=0; i<P.rows(); i++)
        {
            std::cout << " ";
            for( int j=0; j<P.cols(); j++)
            {
                std::cout << P(i,j) << " ";
            }
            std::cout << std::endl;
        }
    }
#else
    void print_matrix(const Cu::Matrix&, const std::string&)
    {
    }
#endif
}

namespace MatrixFuncs
{
    void HomogeneousNormalise( Cu::Matrix& in )
    {
        Cu::Matrix bottomRow = in.bottomRows(1);
        Cu::Matrix brep = bottomRow.replicate(in.rows(),1);
        in = in.cwiseQuotient(brep);
    }

    std::vector<uint32_t> find( const Cu::Matrix& in, const std::function<bool(float)>& func  )
    {
        std::vector<uint32_t> out(in.cols()*in.rows());
        uint32_t ind=0;
        for(uint32_t i=0; i<in.cols()*in.rows(); i++)
        {
            if(func(in.data()[i]))
            {
                out[ind++]=i;
            }
        }
        out.resize(ind);
        return out;
    }

    std::pair<Cu::Matrix, Cu::Matrix> EigenValues( const Cu::Matrix& A  )
    {
        Eigen::EigenSolver<Cu::Matrix> es(A);
        return std::make_pair(es.eigenvectors().real() / es.eigenvectors().real().determinant(), es.eigenvalues().real());
    }
}

namespace MeshFuncs
{
    template<typename T>
    inline void EraseFromList( std::vector<std::vector<T>>& list, const std::vector<bool>& indicators )
    {
        uint32_t revfit =  list.size()-1;
        uint32_t remit = indicators.size()-1;
        for(; true; --revfit)
        {
            if(indicators[remit])
            {
                list.erase( list.begin()+revfit );
            }
            if(revfit==0) break;
            --remit;
        }
    }

    void RemoveVerticesFaces( VerticesFaces& vf, std::vector<uint32_t>& VerticesToRemove )
    {
        std::vector<bool> removev(vf.first.size(),false);
        std::vector<bool> removef(vf.second.size(),false);
        std::vector<uint32_t> newFaceNames(vf.first.size(),0);

        //VerticesToRemove.sort();
        std::sort(VerticesToRemove.begin(), VerticesToRemove.end());

        for(auto vtr : VerticesToRemove)
        {
            removev[vtr] = true;
        }

        uint32_t i=0,index=0;
        for(auto& nfn : newFaceNames)
        {
            if(!removev[i++])
                nfn = index++;
        }

        auto bf = removef.begin();
        for(auto& f : vf.second)
        {
            for(auto& fi : f)
            {
                if(removev[fi])
                {
                    *bf = true;
                    break;
                }
            }
            ++bf;
        }

        // Remove faces
        EraseFromList<uint32_t>(vf.second,removef);

        for(auto& f1 : vf.second)
            for(auto& f2 : f1)
              f2 = newFaceNames[f2];

        // Remove vertices
        EraseFromList<double>(vf.first,removev);
    }

    void WritePly(const VerticesFaces& vf, const std::string filename)
    {
        std::ofstream f(filename, std::ofstream::out);
        f << "ply" << std::endl;
        f << "format ascii 1.0" << std::endl;
        f << "comment CurveFit generated" << std::endl;
        f << "element vertex " << vf.first.size() << std::endl;
        f << "property float x" << std::endl;
        f << "property float y" << std::endl;
        f << "property float z" << std::endl;
        f << "element face " << vf.second.size() << std::endl;
        f << "property list uchar int vertex_indices" << std::endl;
        f << "end_header" << std::endl;
        for(auto& v : vf.first)
        {
            f << v[0] << " " << v[1] << " " << v[2] << std::endl;
        }
        for(auto& fa : vf.second)
        {
            f << 3 << " " << fa[0] << " " << fa[1] << " " << fa[2] << std::endl;
        }
        f.close();
    }
}

namespace Funcs
{
    double normal_pdf(const double x, const double m, const double s)
    {
        static const double inv_sqrt_2pi = 0.398942280401433;
        double a = (x - m) / s;

        return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
    }
}
