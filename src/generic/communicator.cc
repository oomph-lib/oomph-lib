// MPI headers
#ifdef OOMPH_HAS_MPI
#include "mpi.h"
#endif

// Oomph-lib error handler
#include "communicator.h"
#include "matrices.h"

namespace oomph{

#ifdef OOMPH_HAS_MPI

//=============================================================================
/// \short A broadcast function for DenseMatrix<double>
//=============================================================================
 void OomphCommunicator::broadcast(const int& source, DenseMatrix<double>& x)
 {
   // Get number of entries on processor source (where the matrix exists)
   unsigned nrow,ncol;
   if (this->my_rank()==source)
    {
     nrow=x.nrow();
     if(nrow>0)
      {ncol=x.ncol();}
     else
      {ncol=0;}
    }

   // Broadcast to everybody how many entries to expect
   MPI_Bcast(&nrow,1,MPI_UNSIGNED_LONG,source,this->mpi_comm());
   MPI_Bcast(&ncol,1,MPI_UNSIGNED_LONG,source,this->mpi_comm());

   if(ncol!=0 && nrow!=0)
    {
     // convert into a C-style array
     double* x_bcast=new double[nrow*ncol];

     if(this->my_rank()==source)
      for (unsigned long i=0;i<nrow;i++)
       {
        for (unsigned long j=0;j<ncol;j++)
         {
          x_bcast[i*ncol+j]=x(i,j);
         }
       }
   
     // broadcast the array
     MPI_Bcast(x_bcast,ncol*nrow,MPI_DOUBLE,source,this->mpi_comm());

     // Now convert back into matrix (everywhere apart from source)
     if (this->my_rank()!=source)
      {
       x.resize(nrow,ncol);
       for (unsigned long i=0;i<nrow;i++)
        for (unsigned long j=0;j<ncol;j++)
         {
          x(i,j)=x_bcast[i*ncol+j];
         }
      }
     // delete C-style array
     delete[] x_bcast;
    }
   else
    {
     x.resize(0,0);
    }
  }

#endif

} // end of oomph namespace
