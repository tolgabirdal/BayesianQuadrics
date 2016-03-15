#include "Quadric.h"


namespace Geometry
{
    Quadric::Quadric()
    {
        for(uint32_t i=0; i<3; ++i)
        {
            for(uint32_t j=0; j<3; ++j)
            {
                m_params.m_A[i][j]=0;
            }
            m_params.m_b[i]=0;
        }
        m_params.m_c=0;
    }

    Quadric::Quadric(const double* data, const uint32_t N, const double sigma)
    {
        Fit(data, N, sigma);
    }

    Quadric::Quadric(const Parameters& params) :
        m_params(params)
    {

    }

    Quadric::Parameters Quadric::GetParameters() const
    {
        return m_params;
    }

    void Quadric::FitFast( const double* data, const uint32_t N, const double sigma, FastDataStore* fds )
    {
        //Z.setZero(MN,N);
        fds->Q.setZero();
        fds->z.setZero();
        fds->sol.setZero();
        fds->eye.setZero();
        fds->z << 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0;

        uint32_t k, index, i, j;
        for( k=0; k<N; k++ )
        {
            fds->Z[0]  = data[3*k]*data[3*k];
            fds->Z[1]  = data[3*k]*data[1+3*k];
            fds->Z[2]  = data[3*k]*data[2+3*k];
            fds->Z[3]  = data[1+3*k]*data[3*k];
            fds->Z[4]  = data[1+3*k]*data[1+3*k];
            fds->Z[5]  = data[1+3*k]*data[2+3*k];
            fds->Z[6]  = data[2+3*k]*data[3*k];
            fds->Z[7]  = data[2+3*k]*data[1+3*k];
            fds->Z[8]  = data[2+3*k]*data[2+3*k];
            fds->Z[9]  = data[3*k];
            fds->Z[10] = data[1+3*k];
            fds->Z[11] = data[2+3*k];
            fds->Z[12] = 1;
            for(i=0;i<13;++i)
            for(j=0;j<13;j++)
                fds->Q(i,j) += fds->Z[i]*fds->Z[j];
        }

        fds->a=0;
        for(i=0;i<13;++i)
        for(j=0;j<13;j++)
        {
            fds->a += fds->sol(i)*fds->Q(i,j)*fds->sol(j);
            fds->Q(i,j)*=sigma*sigma;
            if(i==j) fds->Q(i,j)+=1;
        }

        fds->sol = static_cast<Eigen::Matrix<double,13,1>>(fds->Q.colPivHouseholderQr().solve(fds->z));

        index=0;
        for(i=0; i<3;i++)
        for(j=0;j<3;j++)
        {
            m_params.m_A[i][j] = std::real(fds->sol(index));
            index++;
        }
        for(i=0;i<3;i++)
        {
            m_params.m_b[i] = std::real(fds->sol(index));
            index++;
        }
        m_params.m_c = std::real(fds->sol(index));

        for(i=0;i<3;++i)
        {
            fds->a += m_params.m_A[i][i]*m_params.m_A[i][i] - 2*m_params.m_A[i][i];
        }
        fds->a += fds->a/(sigma*sigma);
    }

    double Quadric::Fit(const double* data, const uint32_t N, const double sigma)
    {
        const uint32_t M = 3;
        const uint32_t MN = M*M+M+1;

        Cu::Matrix Z(MN,N), Q(MN,MN), z(MN,1), sol(MN,1), eye(MN,MN);
        Z.setZero(MN,N);Q.setZero(MN,MN); z.setZero(MN,1);
        sol.setZero(MN,1); eye.setZero(MN,MN);

        uint32_t k, index, i, j;
        for( k=0; k<N; k++ )
        {
            index=0;
            for(i=0; i<M;i++)
            {
                for(j=0;j<M;j++)
                {
                    Z(index++,k) = data[i+M*k]*data[j+M*k];
                }
            }
            for( i=0; i<M; i++)
            {
                Z(index++,k) = data[i+M*k];
            }
            Z(index,k) = 1;
        }

        Q = Z*Z.transpose();

        for(k=0;k<3;k++)
        {
            z((M*k)+k)=1;
        }

        for(k=0;k<3*3;++k)
        {
            eye(k,k)=1;
        }

        sol = Cu::Matrix((sigma*sigma*Q+eye).colPivHouseholderQr().solve(z));

        index=0;
        for(i=0; i<M;i++)
        {
            for(j=0;j<M;j++)
            {
                m_params.m_A[i][j] = std::real(sol(index));
                index++;
            }
        }
        for(i=0;i<M;i++)
        {
            m_params.m_b[i] = std::real(sol(index));
            index++;
        }
        m_params.m_c = std::real(sol(index));

        Cu::Matrix r = sol.transpose()*Q*sol;
        double a=0;
        for(i=0;i<M;++i)
        {
            a+=m_params.m_A[i][i]*m_params.m_A[i][i] - 2*m_params.m_A[i][i];
        }
        return r(0) + a/(sigma*sigma);
    }

    inline std::vector<double> GeneratePoint(const double& theta, const double& phi)
    {
        std::vector<double> td(3,0.0);
        auto tdit = td.begin();
        *(tdit++) = std::sin(phi)*std::cos(theta);
        *(tdit++) = std::sin(phi)*std::sin(theta),
        *(tdit++) = std::cos(phi);
        return td;
    }

    inline std::vector<double> GenerateCPoint(const double& theta, const double& h, const uint32_t dimension)
    {
        std::vector<double> td(3,0.0);
        auto tdit = td.begin();
        switch(dimension)
        {
        case 0:
            *(tdit++) = h;
            *(tdit++) = std::cos(theta);
            *(tdit++) = std::sin(theta);
            break;
        case 1:
            *(tdit++) = std::sin(theta);
            *(tdit++) = h;
            *(tdit++) = std::cos(theta);
            break;
        case 2:
            *(tdit++) = std::cos(theta);
            *(tdit++) = std::sin(theta);
            *(tdit++) = h;
            break;
        }

        return td;
    }

    VerticesFaces Quadric::GenerateUnitCylinder( const uint32_t N, const uint32_t principleDimension ) const
    {
        std::vector<double> h(N+1);
        std::vector<double> theta(N+1);
        for(uint32_t i=0;i<=N;i++)
        {
          h[i] = 2*(static_cast<double>(i)/static_cast<double>(N))-1;
          theta[i] = ((2*M_PI) / static_cast<double>(N))*static_cast<double>(i);
        }

        std::vector<std::vector<double>> vertices(N*N*4, std::vector<double>(3,0.0));
        std::vector<std::vector<uint32_t>> faces(N*N*2, std::vector<uint32_t>(3,0));
        auto vit = vertices.begin();
        auto fit = faces.begin();

        uint32_t oi=0;
        uint32_t j,i;

        std::vector<uint32_t> f1(3), f2(3);
        for(j=0; j<N; j++)
        {
          for(i=0; i<N; i++)
          {
            *(vit++) = GenerateCPoint(theta[i], h[j],principleDimension);
            *(vit++) = GenerateCPoint(theta[i+1], h[j],principleDimension);
            *(vit++) = GenerateCPoint(theta[i], h[j+1],principleDimension);
            *(vit++) = GenerateCPoint(theta[i+1], h[j+1],principleDimension);

            f1[0]=oi;  f1[1]=oi+1;f1[2]=oi+2;
            f2[0]=oi+2;f2[1]=oi+1;f2[2]=oi+3;
            *(fit++) = f1;
            *(fit++) = f2;
            oi+=4;
          }
        }
        return std::make_pair( vertices, faces );
    }

    VerticesFaces Quadric::GenerateUnitCircle( const uint32_t N ) const
    {
        std::vector<double> phi(N+1);
        std::vector<double> theta(N+1);
        for(uint32_t i=0;i<=N;i++)
        {
          phi[i] = ((M_PI)/ static_cast<double>(N))*static_cast<double>(i);
          theta[i] = ((2*M_PI) / static_cast<double>(N))*static_cast<double>(i);
        }

        std::vector<std::vector<double>> vertices(N*N*4, std::vector<double>(3,0.0));
        std::vector<std::vector<uint32_t>> faces(N*N*2, std::vector<uint32_t>(3,0));

        uint32_t oi=0;
        uint32_t j,i;
        auto vit = vertices.begin();
        auto fit = faces.begin();
        std::vector<uint32_t> f1(3), f2(3);
        for(j=0; j<N; j++)
        {
          for(i=0; i<N; i++)
          {
            *(vit++) = GeneratePoint(theta[i], phi[j]);
            *(vit++) = GeneratePoint(theta[i+1], phi[j]);
            *(vit++) = GeneratePoint(theta[i], phi[j+1]);
            *(vit++) = GeneratePoint(theta[i+1], phi[j+1]);

            f1[0]=oi;  f1[1]=oi+1;f1[2]=oi+2;
            f2[0]=oi+2;f2[1]=oi+1;f2[2]=oi+3;
            *(fit++) = f1;
            *(fit++) = f2;
            oi+=4;
          }
        }
        return std::make_pair( vertices, faces );
    }

    inline void FixMinimumeigenValue( const double* data, std::pair<Cu::Matrix, Cu::Matrix>& US, Cu::Matrix& A, Cu::Matrix& b, Cu::Matrix& mu, Cu::Matrix& d, double c, uint32_t minrow, uint32_t NPts)
    {
        Cu::Matrix D = Eigen::Map<Cu::Matrix>(const_cast<double*>(data),3,NPts);
        //auto mu2 = Cu::Matrix(D.rowwise().mean());
        D.colwise() -= Cu::Vector(mu);
        Cu::Matrix reorient = US.first.transpose()*D;

        auto max = Cu::Matrix(reorient.rowwise().maxCoeff());
        auto min = Cu::Matrix(reorient.rowwise().minCoeff());
        US.second(minrow) = d(0)/std::max(std::pow(max(minrow),2), std::pow(min(minrow),2));

        for(uint32_t n,l=0;l<3;l++)
        for(n=0;n<3;++n)
            A(n,l)=US.first(n,l)*US.second(l);
        A = A*US.first.transpose();

        mu = -(A.fullPivHouseholderQr().solve(b))/2.0;
        d = mu.transpose() * A * mu;
        d(0) -= c;

        if(d(0)<0)
        {
            d = -d;
            mu=-mu;
            US.first = -US.first;
        }
    }

    inline Cu::Matrix FindC(const std::pair<Cu::Matrix, Cu::Matrix>& US, const Cu::Matrix& d)
    {
        Cu::Matrix C = ( US.first * Cu::Matrix(Cu::Vector(US.second).asDiagonal()).cwiseSqrt() ).transpose();
        C = std::sqrt(d(0))*C.inverse();
        return C;
    }

    inline void CreateVerticesFromC(std::vector<std::vector<double>>& vertex, const VerticesFaces& vf, const Cu::Matrix& C, const Cu::Matrix& mu)
    {
        uint32_t k,j,index=0;
        for(auto& vi : vf.first)
        {
            for(k=0;k<3;++k)
            {
                for(j=0;j<3;++j)
                {
                    vertex[index][k] += C(k,j)*vi[j];
                }
                vertex[index][k] += mu(k);
            }
            ++index;
        }
    }

    void Mud(const Cu::Matrix& A, const Cu::Matrix& b, const double c, Cu::Matrix& mu, Cu::Matrix& d)
    {
        mu = -(A.fullPivHouseholderQr().solve(b))/2.0;
        d = mu.transpose() * A * mu;
        d(0)-=c;
    }

    VerticesFaces Quadric::GeneratePrincipleCylinder( uint32_t N, const double *data, const uint32_t NPts ) const
    {
        const double* ma = reinterpret_cast<const double*>(m_params.m_A);
        const double* mb = reinterpret_cast<const double*>(m_params.m_b);
        Cu::Matrix A = Eigen::Map<Cu::Matrix>(const_cast<double *>(ma),3,3).transpose();
        Cu::Matrix b = Eigen::Map<Cu::Matrix>(const_cast<double *>(mb),3,1);

        auto US = MatrixFuncs::EigenValues(A);
        uint32_t minrow;
        double mins = Cu::Vector(US.second).minCoeff(&minrow);
        auto vf = GenerateUnitCylinder(N, minrow);

        Cu::Matrix mu,d;
        Mud(A,b,m_params.m_c,mu,d);

        if(d(0)<0)
        {
            d = -d;
            mu=-mu;
            US.first = -US.first;
        }

        Cu::Matrix I(3,3);
        I.setIdentity();

        if(mins < 0)
        {
            FixMinimumeigenValue(data, US, A, b, mu, d, m_params.m_c, minrow, NPts);
        }

        Cu::Matrix C = FindC(US, d);

        std::vector<std::vector<double>> vertex(vf.first.size(),std::vector<double>(3,0.0));
        CreateVerticesFromC(vertex, vf, C, mu);

        return std::make_pair(vertex, vf.second);
    }

    VerticesFaces Quadric::GenerateEllipsoid( uint32_t N, const double *data, const uint32_t NPts ) const
    {
        auto vf = GenerateUnitCircle(N);
        const double* ma = reinterpret_cast<const double*>(m_params.m_A);
        const double* mb = reinterpret_cast<const double*>(m_params.m_b);
        Cu::Matrix A = Eigen::Map<Cu::Matrix>(const_cast<double *>(ma),3,3).transpose();
        Cu::Matrix b = Eigen::Map<Cu::Matrix>(const_cast<double *>(mb),3,1);
        auto US = MatrixFuncs::EigenValues(A);
        US.first = -US.first; // Turns the quadric inside out (for consistency)
        uint32_t minrow;
        double mins = Cu::Vector(US.second).minCoeff(&minrow);

        Cu::Matrix I(3,3);
        I.setIdentity();

        Cu::Matrix mu,d;
        Mud(A,b,m_params.m_c,mu,d);

        if(d(0)<0)
        {
            d = -d;
            mu=-mu;
            US.first = -US.first;
        }

        if( mins < 0)
        { // @todo this eigenpair should be scaled properly.
            FixMinimumeigenValue(data, US, A, b, mu, d, m_params.m_c, minrow, NPts);
        }

        Cu::Matrix C = FindC(US, d);

        std::vector<std::vector<double>> vertex(vf.first.size(),std::vector<double>(3,0.0));
        CreateVerticesFromC(vertex, vf, C, mu);

        //Cu::print_matrix_octave(vertex, "vertexlist");
        return std::make_pair(vertex, vf.second);
    }

    VerticesFaces Quadric::GetPrincipleCylinder(const uint32_t N, const double *data, const uint32_t NPts) const
    {
        return GeneratePrincipleCylinder(N, data, NPts);
    }

    VerticesFaces Quadric::GetMeshRep(const uint32_t N, const double *data, const uint32_t NPts) const
    {
        return GenerateEllipsoid(N, data,  NPts);
    }

    std::vector<bool> Quadric::GetPointCloudBounds(const uint32_t StartN, const double* data, const uint32_t N) const
    {
        auto vf = GetMeshRep( StartN, data,  N );
        const double* ma = reinterpret_cast<const double*>(m_params.m_A);
        const double* mb = reinterpret_cast<const double*>(m_params.m_b);
        Cu::Matrix A = Eigen::Map<Cu::Matrix>(const_cast<double *>(ma),3,3).transpose();
        Cu::Matrix b = Eigen::Map<Cu::Matrix>(const_cast<double *>(mb),3,1);
        Cu::Matrix mu = -0.5*(A.inverse()*b);
        auto US = MatrixFuncs::EigenValues(A);

        Cu::Matrix D = Eigen::Map<Cu::Matrix>(const_cast<double*>(data),3,N);
        D.colwise() -= Cu::Vector(mu);
        Cu::Matrix reorient = US.first.transpose()*D;
        auto max = Cu::Matrix(reorient.rowwise().maxCoeff());
        auto min = Cu::Matrix(reorient.rowwise().minCoeff());

        uint32_t k,j,ind=0;
        double t[3];
        std::vector<bool> out(vf.first.size(), false);
        for(const auto& v : vf.first)
        {
            for(j=0;j<3;++j)
            {
                t[j]=0.0;
                for(k=0; k<3;++k)
                {
                    t[j]+= US.first.data()[3*j+k]*(v[k]-mu.data()[k]);
                }
                if(t[j]>max.data()[j] || t[j] < min.data()[j])
                {
                    out[ind]=true;
                    break;
                }
            }
            ++ind;
        }
        return out;
    }

    void Quadric::GetBoundedMeshRep(VerticesFaces& vf, const double *data, const uint32_t N) const
    {
        //auto vf = GetMeshRep( StartN );
        const double* ma = reinterpret_cast<const double*>(m_params.m_A);
        const double* mb = reinterpret_cast<const double*>(m_params.m_b);
        Cu::Matrix A = Eigen::Map<Cu::Matrix>(const_cast<double *>(ma),3,3).transpose();
        Cu::Matrix b = Eigen::Map<Cu::Matrix>(const_cast<double *>(mb),3,1);

        Cu::Matrix mu = -0.5*(A.inverse()*b);
        auto US = MatrixFuncs::EigenValues(A);

        Cu::Matrix D = Eigen::Map<Cu::Matrix>(const_cast<double*>(data),3,N);
        D.colwise() -= Cu::Vector(mu);
        Cu::Matrix reorient = US.first.transpose()*D;
        auto max = Cu::Matrix(reorient.rowwise().maxCoeff());
        auto min = Cu::Matrix(reorient.rowwise().minCoeff());

        std::vector<uint32_t> removev( vf.first.size());
        uint32_t k,j,ind=0,added=0;
        double t[3];
        for(const auto& v : vf.first)
        {
            for(j=0;j<3;++j)
            {
                t[j]=0.0;
                for(k=0; k<3;++k)
                {
                    t[j]+= US.first.data()[3*j+k]*(v[k]-mu.data()[k]);
                }
                if(t[j]>max.data()[j] || t[j] < min.data()[j])
                {
                    removev[added++] = ind;
                    break;
                }
            }
            ++ind;
        }
        removev.resize(added);

        MeshFuncs::RemoveVerticesFaces(vf, removev);

        //return vf;
    }
}
