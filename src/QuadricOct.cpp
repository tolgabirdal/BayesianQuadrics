
#include <octave/oct.h>
#include <vector>

#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <cassert>

#include "Matrix.h"
#include "Quadric.h"

void handler(int sig) 
{
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}

DEFUN_DLD(Quadric, args, , "[A b c] = Quadric( Data, sigma ) ")
{
  signal(SIGSEGV, handler);   // install our handler

  int nargin = args.length ();
  if( nargin != 2 )
  {
    print_usage();
  }

  const Matrix Data = args(0).matrix_value();
  const double sigma = args(1).scalar_value();
  
  if(Data.rows() != 3 )
  {
    std::cerr << "Quadric: There must be three rows to the data matrix" << std::endl;
    assert(false);
  }
  if(Data.cols() < 12)
  {
    std::cerr << "Quadric: There must be more than 12 observations" << std::endl;
    assert(false);
  }
  
 
  Geometry::Quadric::Quadric Q(Data.data(), Data.cols(), sigma);
  Geometry::Quadric::Parameters params = Q.GetParameters();
  
  Matrix A(3,3), b(3,1);
  for(uint32_t i,j=0;j<3;++j)
  {
    for(i=0;i<3;++i)
    {
      A(i,j) = params.m_A[i][j];
    }
    b(j) = params.m_b[j];
  }
  double c = params.m_c;
  
  octave_value_list ovl;
  ovl(0) = A;
  ovl(1) = b;
  ovl(2) = c;
  
  return ovl;
}
