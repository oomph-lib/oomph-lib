//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
#include "oomph_utilities.h"
#include "Vector.h"
#include "matrices.h"
#include <algorithm>
#include <limits.h>

#ifdef OOMPH_HAS_MPI
#include "mpi.h"
#endif

namespace oomph
{

//==============================================
/// Namespace for error messages for broken
/// copy constructors and assignment operators
//==============================================
namespace BrokenCopy
{

 /// Issue error message and terminate execution
 void broken_assign(const std::string& class_name)
  {
   //Write the error message into a string
   std::string error_message = "Assignment operator for class\n\n";
   error_message += class_name;
   error_message += "\n\n";
   error_message += "is deliberately broken to avoid the accidental \n";
   error_message += "use of the inappropriate C++ default.\n";
   error_message += "If you really need an assignment operator\n"; 
   error_message += "for this class, write it yourself...\n";
   
   throw OomphLibError(error_message,"broken_assign()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 
 
 /// Issue error message and terminate execution
 void broken_copy(const std::string& class_name)
  {
   //Write the error message into a string
   std::string error_message = "Copy constructor for class\n\n";
   error_message += class_name;
   error_message += "\n\n";
   error_message += "is deliberately broken to avoid the accidental\n";
   error_message += "use of the inappropriate C++ default.\n";
   error_message += 
    "All function arguments should be passed by reference or\n"; 
   error_message += 
    "constant reference. If you really need a copy constructor\n";
   error_message += "for this class, write it yourself...\n";

   throw OomphLibError(error_message,"broken_copy()",
                       OOMPH_EXCEPTION_LOCATION);
  }
}


//======================================================================
/// Namespace for mathematical constants
//======================================================================
namespace MathematicalConstants
{
 /// Guess what...
 double Pi=4.0*std::atan(1.0);
}



//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////



//====================================================================
/// Namespace for global (cumulative) timings
//====================================================================
namespace CumulativeTimings
{

 /// (Re-)start i-th timer
 void start(const unsigned& i)
 {
  Start_time[i]=clock();
 }

 /// Halt i-th timer
 void halt(const unsigned& i)
 {
  Timing[i]+=clock()-Start_time[i];
 }

 /// Report time accumulated by i-th timer
 double cumulative_time(const unsigned& i)
 {
  return double(Timing[i])/CLOCKS_PER_SEC;
 }

 /// Reset i-th timer
 void reset(const unsigned& i)
 {
  Timing[i]=clock_t(0.0);
 }

 /// Reset all timers
 void reset()
 {
  unsigned n=Timing.size();
  for (unsigned i=0;i<n;i++)
   {
    Timing[i]=clock_t(0.0);
   }
 }

 /// Set number of timings that can be recorded in parallel
 void set_ntimers(const unsigned& ntimers)
 {
  Timing.resize(ntimers,clock_t(0.0));
  Start_time.resize(ntimers,clock_t(0.0));
 }

 /// Cumulative timings
 Vector<clock_t> Timing;

 /// Start times of active timers
 Vector<clock_t> Start_time;

}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////



//=================================================================
/// Evaluate the fitting function and its derivatives
/// w.r.t. fitting parameters (done by FD by default; can be
/// overloaded)
//=================================================================
double LevenbergMarquardtFittingFunctionObject::fitting_function(
 const double& x,
 Vector<double>& dfit_dparam)
{
 // Get reference value
 double fct=fitting_function(x);
 
 // FD step
 double eps_fd=1.0e-8;
 
 // Do FD loop
 unsigned n_param=Parameter.size();
 for (unsigned i=0;i<n_param;i++)
  {
   double backup=Parameter[i];
   Parameter[i]+=eps_fd;
   double fct_pls=fitting_function(x);
   dfit_dparam[i]=(fct_pls-fct)/eps_fd;
   Parameter[i]=backup;
  }
 return fct;
}




//=====================================================================
/// Fit the parameters to the pairs of (x,y) data specified, 
/// using max_iter Levenberg Marquardt iterations
//=====================================================================
void LevenbergMarquardtFitter::fit_it(
 const Vector<std::pair<double,double> >& fitting_data,
 const unsigned& max_iter, const bool& quiet)
{
 if (Fitting_function_object_pt==0)
  {
   throw OomphLibError("Fitting_function_object_pt==0",
                       "Problem::distribute()",
                       OOMPH_EXCEPTION_LOCATION);    
  }

 // Number of parameters:
 unsigned nparam=Fitting_function_object_pt->nparameter();

 // By default regard all parameters as fittable -- can generalise
 // this at some point.
 std::vector<bool> ia(nparam,true);
   
 // Chi squared
 double chisq=0.0;

 // Number of data pairs
 unsigned ndata=fitting_data.size();

 // Vector of standard deviations -- just set to one
 Vector<double> sig(ndata,1.0);

 // Move to vectors (as required by Num Rec interface
 Vector<double> x(ndata), y(ndata);
 for (unsigned i=0;i<ndata;i++)
  {
   x[i]=fitting_data[i].first;
   y[i]=fitting_data[i].second;
  }

 // "Workspace" for numerical recipes
 int ma=nparam; 
 DenseDoubleMatrix covar(ma,ma);
 DenseDoubleMatrix alpha(ma,ma);


 if (!quiet)
  {
   oomph_info << "Chi_squared" << " ";
   for (unsigned i=0;i<unsigned(ma);i++)
    { 
     oomph_info << " parameter " << i << " ";
    }       
   oomph_info << std::endl;
  }

 // Start iteration with almda negative for setup
 double alamda=-0.1;
 for (unsigned iter=0;iter<max_iter;iter++)
  {
   // This is where Num Rec code starts so now it gets really ugly
   static int mfit;
   static double ochisq;
   int j,k,l;
          
   static Vector<double> oneda(ma);
   static Vector<double> atry(ma), beta(ma), da(ma);

   // Initialisation
   if (alamda < 0.0)
    {
     mfit=0;
     for (j=0;j<ma;j++)
      {
       if (ia[j]) mfit++;
      }
     alamda=0.001;
     mrqcof(x,y,sig,Fitting_function_object_pt->parameter(),ia,
            alpha,beta,chisq);
     ochisq=chisq;
     for (j=0;j<ma;j++)
      {
       atry[j]=Fitting_function_object_pt->parameter(j); 
      }
    } 
     
   DenseDoubleMatrix temp(mfit,mfit);
   for (j=0;j<mfit;j++)
    {
     for (k=0;k<mfit;k++)
      {
       covar(j,k)=alpha(j,k);
      }
     covar(j,j)=alpha(j,j)*(1.0+alamda);
     for (k=0;k<mfit;k++)
      {
       temp(j,k)=covar(j,k);
      }
     oneda[j]=beta[j];
    }
     
   // Linear solver
   temp.solve(oneda); 
     
   for (j=0;j<mfit;j++)
    {
     for (k=0;k<mfit;k++)
      {
       covar(j,k)=temp(j,k);
      }
     da[j]=oneda[j];
    }
     
   // Converged
   if (alamda == 0.0)
    {
     return;
    }
     
          
   for (j=0,l=0;l<ma;l++)
    {
     if (ia[l]) atry[l]=Fitting_function_object_pt->parameter(l)+da[j++];
    }
   mrqcof(x,y,sig,atry,ia,covar,da,chisq);
   if (chisq < ochisq)
    {
     alamda *= 0.1;
     ochisq=chisq;
     for (j=0;j<mfit;j++)
      {
       for (k=0;k<mfit;k++) 
        {
         alpha(j,k)=covar(j,k);
        }
       beta[j]=da[j];
      }
     
     for (l=0;l<ma;l++)
      { 
       Fitting_function_object_pt->parameter(l)=atry[l];
      }
       
     if (!quiet)
      {
       // Output with fixed width
       std::cout.setf(std::ios_base::scientific,std::ios_base::floatfield);
       std::cout.width(15);
       std::cout << chisq << " ";
       for (l=0;l<ma;l++)
        { 
         std::cout << atry[l] << " ";
        }       
       std::cout << std::endl;
       // Reset
       std::cout.setf(std::_Ios_Fmtflags(0), std::ios_base::floatfield);
       std::cout.width(0);
      }
    }
   else
    {
     alamda *= 10.0;
     chisq=ochisq;
    }
     
  }

}

//==================================================================
/// Private helper function -- don't look into it...
//==================================================================
void LevenbergMarquardtFitter::mrqcof(Vector<double>& x, 
                                      Vector<double>& y, 
                                      Vector<double>& sig, 
                                      Vector<double>& a,
                                      std::vector<bool>& ia, 
                                      DenseDoubleMatrix& alpha, 
                                      Vector<double>& beta, 
                                      double& chisq)
{
 int i,j,k,l,m,mfit=0;
 double ymod,wt,sig2i,dy;
 
 int ndata=x.size();
 int ma=a.size();
 Vector<double> dyda(ma);
 for (j=0;j<ma;j++)
  {
   if (ia[j]) mfit++;
  }
 
 for (j=0;j<mfit;j++)
  {
   for (k=0;k<=j;k++)
    {
     alpha(j,k)=0.0;
    }
   beta[j]=0.0;
  }
 
 chisq=0.0;
 for (i=0;i<ndata;i++) 
  {
   Vector<double> backup=Fitting_function_object_pt->parameter();
   Fitting_function_object_pt->parameter()=a;
   ymod=Fitting_function_object_pt->fitting_function(x[i],dyda);
   Fitting_function_object_pt->parameter()=backup;
   sig2i=1.0/(sig[i]*sig[i]);
   dy=y[i]-ymod;
   for (j=0,l=0;l<ma;l++)
    {
     if (ia[l])
      {
       wt=dyda[l]*sig2i;
       for (k=0,m=0;m<l+1;m++)
        {
         if (ia[m]) alpha(j,k++) += wt*dyda[m];
        }
       beta[j++] += dy*wt;
      }
    }
   chisq += dy*dy*sig2i;
  }
 
 
 for (j=1;j<mfit;j++)
  {
   for (k=0;k<j;k++)
    {
     alpha(k,j)=alpha(j,k);
    }
  }
}




//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////



//======================================================================
/// Namespace for black-box FD Newton solver.
//======================================================================
namespace BlackBoxFDNewtonSolver
{


/// Function pointer for residual function: Parameters, unknowns, residuals
typedef void (*ResidualFctPt)(const Vector<double>&,
                              const Vector<double>&,
                              Vector<double>&);

 /// Max. # of Newton iterations
 unsigned Max_iter=20;

 /// \short Flag to indicate if progress of Newton iteration is to be 
 /// documented (defaults to false) 
 bool Doc_Progress=false;

 /// FD step
 double FD_step=1.0e-8;

 /// Tolerance
 double Tol=1.0e-8;

/// \short Black-box FD Newton solver:
/// Calling sequence for residual function is
/// \code residual_fct(parameters,unknowns,residuals) \endcode
/// where all arguments are double Vectors.
/// unknowns.size() = residuals.size() 
void black_box_fd_newton_solve(ResidualFctPt residual_fct, 
                               const Vector<double>& params, 
                               Vector<double>& unknowns)
{
 // Jacobian, current and advanced residual Vectors
 unsigned ndof=unknowns.size();
 DenseDoubleMatrix jacobian(ndof);
 Vector<double> residuals(ndof);
 Vector<double> residuals_pls(ndof);
 Vector<double> dx(ndof);

 // Newton iterations
 for (unsigned iloop=0;iloop<Max_iter;iloop++)
  {
   // Evaluate current residuals
   residual_fct(params,unknowns,residuals);

   // Check max. residuals
   double max_res = std::abs(*std::max_element(residuals.begin(),
                                               residuals.end(),
                                               AbsCmp<double>()));

   // Doc progress?
   if (Doc_Progress)
    {
     oomph_info << "\nNewton iteration iter=" << iloop 
                << "\ni residual[i] unknown[i] " << std::endl;
     for (unsigned i=0;i<ndof;i++)
      {
       oomph_info << i << " " << residuals[i] 
                  << " " << unknowns[i] << std::endl;
      }
    }
   
   // Converged?
   if (max_res<Tol) 
    {
     //oomph_info << "Converged " << std::endl;
     return;
    }

   // FD loop for Jacobian
   for (unsigned i=0;i<ndof;i++)
    {
     double backup=unknowns[i];
     unknowns[i]+=FD_step;
  
     // Evaluate advanced residuals
     residual_fct(params,unknowns,residuals_pls);

     // Do FD
     for (unsigned j=0;j<ndof;j++)
      {
       jacobian(j,i)=(residuals_pls[j]-residuals[j])/FD_step;
      }
     
     // Reset fd step
     unknowns[i]=backup;
    }
   
   // Solve (overwrites residuals)
   jacobian.solve(residuals);
   
   // Update:
   for (unsigned i=0;i<ndof;i++)
    {
     unknowns[i]-=residuals[i];
    }
  }
 
 // Failed to converge
 std::ostringstream error_stream;
 error_stream<< "Newton solver did not converge in " 
             << Max_iter << " steps " << std::endl;

 throw OomphLibError(error_stream.str(),
                     "black_box_fd_newton_solve()",
                     OOMPH_EXCEPTION_LOCATION);
}
}

//======================================================================
/// \short Set output directory (we try to open a file in it
/// to see if the directory exists -- if it doesn't we'll
/// issue a warning -- or, if directory_must_exist()==true,
/// die by throwing and OomphLibError
//======================================================================
void DocInfo::set_directory(const std::string& directory)
{
 // Try to open a file in output directory
 char filename[100];
 sprintf(filename,"%s/.dummy_check.dat",directory.c_str());
 std::ofstream some_file;
 some_file.open(filename);
 if (!some_file.is_open())
  {
   //Construct the error message
   std::string error_message = "Problem opening output file.\n";
   error_message += "I suspect you haven't created the output directory ";
   error_message += directory;
   error_message += "\n";
   
   //Issue a warning if the directory does not have to exist
   if(!Directory_must_exist)
    {
     //Create an Oomph Lib warning
     OomphLibWarning(error_message,"set_directory()",
                     OOMPH_EXCEPTION_LOCATION);
    }
   //Otherwise throw an erro
   else
    {
     error_message += "and the Directory_must_exist flag is true.\n";
     throw OomphLibError(error_message,"set_directory()",
                         OOMPH_EXCEPTION_LOCATION);
    }
  }
 //Write to the dummy file
 some_file << "Dummy file, opened to check if output directory " << std::endl;
 some_file << "exists. Can be deleted...." << std::endl;
 some_file.close();
 // Set directory
 Directory=directory;
}



//====================================================================
/// Namespace for command line arguments
//====================================================================
namespace CommandLineArgs
{

 /// Number of arguments + 1
 int Argc;

 /// Arguments themselves
 char** Argv;

 /// Set values
 void setup(int argc, char** argv)
 {
  Argc=argc;
  Argv=argv;
 }

 /// Doc the command line arguments
 void output()
 {
  oomph_info << "You are running the program: " 
            << CommandLineArgs::Argv[0] << std::endl;
  oomph_info << "with the following command line args: " << std::endl;
  for (int i=1;i<CommandLineArgs::Argc;i++)
   {
    oomph_info << CommandLineArgs::Argv[i] << std::endl;
   }
 }

}

#ifdef OOMPH_HAS_MPI
//========================================================================
/// Single (global) instantiation of the mpi output modifier
//========================================================================
MPIOutputModifier oomph_mpi_output;

//========================================================================
/// Precede the output by the processor ID but output everything
//========================================================================
 bool MPIOutputModifier::operator()(std::ostream &stream)
  {
   int my_rank = Communicator_pt->my_rank();
   
   if (!Output_from_single_processor)
    {
     stream << "Processor " << my_rank << ":   ";
     // Continue processing 
     return true;
    }
   else
    {
     if (unsigned(my_rank)==Output_rank)
      {
       stream << "Processor " << my_rank << ":   ";
       // Continue processing        
       return true;
      }
     else
      {
       return false;
      }
    }
  }

#endif

//=========================================================================
/// Basic namespace for MPI helper data and functions; in serial
/// just a basic version, containing default assignments for
/// My_rank and Nproc (simulating the run on a single processor
/// which is appropriate for the serial execution).
/// Namespace extended by copy from mpi_helpers.h in mpi/mpi_src/mpi_generic
//==========================================================================
namespace MPI_Helpers
{
 /// Processor rank
 int My_rank=0;
 
 /// Total number of processors
 int Nproc=1;
 
 /// Has MPI been initialised? (default: no)
 bool MPI_has_been_initialised=false;

#ifdef OOMPH_HAS_MPI

 /// the global communicator
 OomphCommunicator* Communicator_pt;

 /// initialize mpi
 void init(int argc, char **argv)
  {
    // call mpi int
    MPI_Init(&argc,&argv);
    MPI_has_been_initialised=true;

    // create the oomph-lib communicator using MPI_Comm_dup.
    // the communicator has the same group of processes but a new context
    MPI_Comm oomph_comm_world;
    MPI_Comm_dup(MPI_COMM_WORLD,&oomph_comm_world);

    // create the oomph-lib communicator
    // note: oomph_comm_world is deleted when the destructor of 
    // Communicator_pt is called
    Communicator_pt = new OomphCommunicator(oomph_comm_world,true);
    
    // Change MPI error handler so that error will return
    // rather than aborting
    MPI_Errhandler_set(oomph_comm_world, MPI_ERRORS_RETURN);
    
    // Use MPI output modifier: Each processor preceeds its output
    // by its rank
    oomph_mpi_output.communicator_pt() = Communicator_pt;
    oomph_info.output_modifier_pt() = &oomph_mpi_output;

    // LEGACY - until My_rank and Nproc are deleted
    // store My_rank and Nproc
    My_rank = Communicator_pt->my_rank();
    Nproc = Communicator_pt->nproc();
  }

 /// finalize mpi
 void finalize()
  {
    // delete the communicator
    delete Communicator_pt;

    // and call MPI_Finalize
    MPI_Finalize();
  }



 /// LEGACY
 /// Setup the namespace
 void setup()
 {
  std::ostringstream warning_stream;
  warning_stream  <<"WARNING:\n\n"
                  <<"The method MPI_Helpers::setup() has been replaced "
                  <<"with:\n\n       MPI_Helpers::init(...)\n\n"
                  <<"Notes:\n\n1. It is no longer necessary to call "
                  <<"MPI_Init(...)\n   at the begining of main(...) "
                  <<"because it is called\n   by MPI_Helpers::init(...)."
                  <<"\n\n2. The use of MPI_COMM_WORLD should be replaced "
                  <<"with:"
                  <<"\n\n       MPI_Helpers::Communicator_pt->mpi_comm()"
                  <<"\n\n   or (for example)\n\n"
                  <<"       problem_pt->communicator_pt()->mpi_comm()\n\n"
                  <<"3. Calls to MPI_Finalize() should be replaced with:"
                  <<"\n\n       MPI_Helpers::finalize()\n";
  ObsoleteCode::obsolete(warning_stream.str());

  // Check that MPI_Init has been called and throw an error if it's not
  int flag = 0;
  MPI_Initialized(&flag);
  if (!flag)
   {
    std::ostringstream error_message_stream;
    error_message_stream << "MPI_Init must be called before using "
                         << "MPI_Helpers::setup!!!\n";
    throw OomphLibError(error_message_stream.str(),
                        "MPI_Helpers::setup()",
                        OOMPH_EXCEPTION_LOCATION);
   }  

  // Set the bool to say MPI has been initialised
  // NB: make sure MPI_Init is called BEFORE MPI_Helpers::setup()
  MPI_has_been_initialised=true;

  Communicator_pt = new OomphCommunicator(MPI_COMM_WORLD,true);

  // Figure out number of processes
  MPI_Comm_size(MPI_COMM_WORLD,&Nproc);
  
  // Which processor am I?
  MPI_Comm_rank(MPI_COMM_WORLD,&My_rank);

  // Use MPI output modifier: Each processor preceeds its output
  // by its rank
  oomph_mpi_output.communicator_pt() = Communicator_pt;
  oomph_info.output_modifier_pt() = &oomph_mpi_output;
  
 }
#endif // OOMPH_HAS_MPI
}

//====================================================================
/// Namespace for flagging up obsolete parts of the code
//====================================================================
namespace ObsoleteCode
{
 
 /// Flag up obsolete parts of the code
 bool FlagObsoleteCode=true;

 /// Output warning message
 void obsolete()
 {
  if (FlagObsoleteCode)
   {
    std::string junk;
    oomph_info << "\n\n--------------------------------------------\n";
    oomph_info << "You are using obsolete code " << std::endl;
    oomph_info << "--------------------------------------------\n\n";
    oomph_info << "Enter: \"s\" to suppress further messages" << std::endl;
    oomph_info << 
     "       \"k\" to crash the code to allow a trace back in the debugger" 
         << std::endl;

    oomph_info << "       any other key to continue\n \n";
    oomph_info << 
     "                    [Note: Insert \n \n ";
    oomph_info << 
     "                            ObsoleteCode::FlagObsoleteCode=false;\n \n";
    oomph_info << 
     "                     into your code to suppress these messages \n";
    oomph_info << 
     "                     altogether.] \n";

    std::cin >> junk;
    if (junk=="s")
     {
      FlagObsoleteCode=false;
     }
    if (junk=="k")
     {
      throw OomphLibError("Killed","ObsoleteCode::obsolete()",
                          OOMPH_EXCEPTION_LOCATION);
     }
   }
 }


 ///Ouput a warning message with a string argument
 void obsolete(const std::string &message)
 {
  if(FlagObsoleteCode)
   {
    oomph_info << "\n\n------------------------------------" << std::endl;
    oomph_info << message << std::endl;
    oomph_info << "----------------------------------------" << std::endl;
    
    obsolete();
   }
 }

}




//====================================================================
/// Namespace for tecplot stuff
//====================================================================
namespace TecplotNames
{

  /// Tecplot colours 
  Vector<std::string> colour;


 /// Setup namespace
 void setup()
   {
    colour.resize(5);
    colour[0]="RED";
    colour[1]="GREEN";
    colour[2]="BLUE";
    colour[3]="CYAN";
    colour[4]="BLACK";
   }


}



#ifdef LEAK_CHECK

//====================================================================
/// Namespace for leak check: Keep a running count of all instantiated
/// objects -- add your own if you want to...
//====================================================================
namespace LeakCheckNames
{

 long QuadTree_build;
 long OcTree_build;
 long QuadTreeForest_build;
 long OcTreeForest_build;
 long RefineableQElement<2>_build;
 long RefineableQElement<3>_build;
 long MacroElement_build;
 long HangInfo_build;
 long Node_build;
 long GeomReference_build;
 long AlgebraicNode_build;

 void reset()
  {
   QuadTree_build=0;
   OcTree_build=0;
   QuadTreeForest_build=0;
   OcTreeForest_build=0;
   RefineableQElement<2>_build=0;
   RefineableQElement<3>_build=0;
   MacroElement_build=0;
   HangInfo_build=0;
   Node_build=0;
   GeomReference_build=0;
   AlgebraicNode_build=0;
  }

 void doc()
 {
  oomph_info << 
   "\n Leak check: # of builds - # of deletes for the following objects: \n\n";
   oomph_info << "LeakCheckNames::QuadTree_build " 
        << LeakCheckNames::QuadTree_build << std::endl;
   oomph_info << "LeakCheckNames::QuadTreeForest_build "  
        << LeakCheckNames::QuadTreeForest_build << std::endl;
   oomph_info << "LeakCheckNames::OcTree_build " 
        << LeakCheckNames::OcTree_build << std::endl;
   oomph_info << "LeakCheckNames::OcTreeForest_build "  
        << LeakCheckNames::OcTreeForest_build << std::endl;
   oomph_info << "LeakCheckNames::RefineableQElement<2>_build " 
        << LeakCheckNames::RefineableQElement<2>_build << std::endl;
   oomph_info << "LeakCheckNames::RefineableQElement<3>_build " 
        << LeakCheckNames::RefineableQElement<3>_build << std::endl;
   oomph_info << "LeakCheckNames::MacroElement_build " 
        <<  LeakCheckNames::MacroElement_build<< std::endl;
   oomph_info << "LeakCheckNames::HangInfo_build " 
        <<  LeakCheckNames::HangInfo_build<< std::endl;
   oomph_info << "LeakCheckNames::Node_build " 
        <<  LeakCheckNames::Node_build<< std::endl;
   oomph_info << "LeakCheckNames::GeomReference_build " 
        <<  LeakCheckNames::GeomReference_build<< std::endl;
   oomph_info << "LeakCheckNames::AlgebraicNode_build " 
        <<  LeakCheckNames::AlgebraicNode_build<< std::endl;
   oomph_info << std::endl;
 }



}


#endif



///////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


//====================================================================
/// Namespace for pause() command
//====================================================================
namespace PauseFlags
{
 
 /// Flag to enable pausing code -- pause the code by default
 bool PauseFlag=true;

}

//======================================================================
/// Pause and display message
//======================================================================
void pause(std::string message)
{
   std::string junk;
   if (PauseFlags::PauseFlag)
    {
     oomph_info << message 
          << "\n hit any key to continue [hit \"S\" "
          << "to suppress further interruptions]\n";
     std::cin >> junk;
     if (junk=="S")
      {
       PauseFlags::PauseFlag=false;
      }
    }
   else
    {
     oomph_info << "\n[Suppressed pause message] \n";
    }
 }





///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

//=============================================================================
/// Helper for recordning execution time.
//=============================================================================
namespace TimingHelpers
{

 /// returns the time in seconds after some point in past
 double timer()
  {
#ifdef OOMPH_HAS_MPI
   return MPI_Wtime();
#else
   time_t t = clock();
   return  double(t) / double(CLOCKS_PER_SEC);
#endif
  }
}//end of namespace TimingHelpers


// //====================================================================
// // Half-arsed attempt at docmenting memory usage
// //====================================================================
// void init_doc_memory_usage()
// {
//  int my_pid=int(getpid());
 
//  //oomph_info << "My pid " << my_pid << std::endl;
 
//  char check_mem_command[100];
//  sprintf(check_mem_command,
//    "rm -f memory_usage_%i.dat; date >>memory_usage_%i.dat ",my_pid,my_pid);
//  system(check_mem_command);
// }






// //====================================================================
// // Half-arsed attempt at docmenting memory usage
// //====================================================================
// void doc_memory_usage()
// {
//    int my_pid=int(getpid());

//    char check_mem_command[100];
//    sprintf(check_mem_command,
//     "top n1 -p %i | awk '{MEMORY=$12} END {print MEMORY}' >> memory_usage_%i.dat",my_pid,my_pid);

//    system(check_mem_command);
//  }

}
