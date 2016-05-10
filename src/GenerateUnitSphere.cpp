
#include <octave/oct.h>
#include <vector>
#include <math.h>

struct Td
{
  double x;
  double y;
  double z;
};

struct Ti
{
  uint32_t x;
  uint32_t y;
  uint32_t z;
};

Td GeneratePoint(const double& theta, const double& phi)
{
  Td pt
  {
    std::sin(phi)*std::cos(theta),
    std::sin(phi)*std::sin(theta),
    cos(phi)
  };
  return pt;
}

DEFUN_DLD(GenerateUnitSphere, args, , 
          "[v f] = GenerateUnitSphere( N )")
{
  const uint32_t N = args(0).uint_value();
  const uint32_t Ntot = N*N;
  
  std::vector<double> phi(N+1);
  std::vector<double> theta(N+1);
  for(uint32_t i=0;i<=N;i++)
  {
    phi[i] = ((M_PI)/ static_cast<double>(N))*static_cast<double>(i);
    theta[i] = ((2*M_PI) / static_cast<double>(N))*static_cast<double>(i);
  }
  
  std::vector<Td> vertices;
  std::vector<Ti> faces;
  
  uint32_t oi=0;
  uint32_t j,i;
  Td pt1,pt2,pt3,pt4;
  for(j=0; j<N; j++)
  {
    for(i=0; i<=N; i++)
    {
      pt1 = GeneratePoint(theta[j], phi[i]);
      vertices.push_back(pt1);
    }
  }
  
  uint32_t NT=N*(N+1);
  for(j=0; j<N; ++j)
  {
    for(i=0;i<N; ++i)
    {
      faces.push_back({ oi, (oi+1) %NT, (oi+N)%NT  });
      faces.push_back({ (oi+1)%NT, (oi+N+1)%NT , (oi+N)%NT });
      oi++;
    }
    oi++;
  }
  
  
  Matrix v(3,vertices.size());
  v.fill(0.0);
  Matrix f(3,faces.size());
  f.fill(0.0);
  i=0;
  for( auto vi : vertices )
  {
    v(0,i) = vi.x;
    v(1,i) = vi.y;
    v(2,i) = vi.z;
    i++;
  }
  
  i=0;
  for( auto fi : faces )
  {
    f(0,i) = fi.x;
    f(1,i) = fi.y;
    f(2,i) = fi.z;
    i++;
  }
  
  octave_value_list ovl;
  ovl(0) = v;
  ovl(1) = f;
  return ovl;
}
