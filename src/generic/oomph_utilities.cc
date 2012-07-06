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
#ifdef OOMPH_HAS_MPI
#include "mpi.h"
#endif

#include <algorithm>
#include <limits.h>
#include <cstring>


#include "oomph_utilities.h"
#include "Vector.h"
#include "matrices.h"

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
 int ma = static_cast<int>(nparam); 
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
       //Store the current output flags
       std::ios_base::fmtflags ff = std::cout.flags();
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
       std::cout.setf(ff, std::ios_base::floatfield);
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

 /// Use steplength control do make globally convergent (default false)
 bool Use_step_length_control=false; 

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
 Vector<double> gradient(ndof);
 Vector<double> newton_direction(ndof);
    
 double half_residual_squared=0.0;
 double max_step=0.0;
 
 // Newton iterations
 for (unsigned iloop=0;iloop<Max_iter;iloop++)
  {
   // Evaluate current residuals
   residual_fct(params,unknowns,residuals);
   
   // Get half of squared residual and find maximum step length
   // for step length control
   if (Use_step_length_control)
    {
     half_residual_squared=0.0;
     double sum=0.0;
     for (unsigned i=0;i<ndof;i++)
      {
       sum += sqrt(unknowns[i]);
       half_residual_squared+=residuals[i]*residuals[i];
      }
     half_residual_squared*=0.5;
     max_step=100.0*std::max(sqrt(sum),double(ndof));
    }
   
   
   // Check max. residuals
   double max_res = std::fabs(*std::max_element(residuals.begin(),
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
   
   
   // Get gradient
   if (Use_step_length_control)
    {
     for (unsigned i=0;i<ndof;i++)
      {
       double sum=0.0;
       for (unsigned j=0;j<ndof;j++)
        {
         sum += jacobian(j,i)*residuals[j];
        }
       gradient[i]=sum;
      }
    }
   
   
   // Solve 
   jacobian.solve(residuals,newton_direction);
   
   // Update
   if (Use_step_length_control)
    {     
     for (unsigned i=0;i<ndof;i++)
      {
       newton_direction[i]*=-1.0;
      }
     // Update with steplength control
     Vector<double> unknowns_old(unknowns);
     double half_residual_squared_old=half_residual_squared;
     bool check;
     line_search(unknowns_old,
                 half_residual_squared_old,
                 gradient,
                 residual_fct, 
                 params,
                 newton_direction,
                 unknowns,
                 half_residual_squared,
                 max_step,
                 check);
    }
   else
    {
     // Direct Newton update:
     for (unsigned i=0;i<ndof;i++)
      {
       unknowns[i]-=newton_direction[i];
      }
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









//=======================================================================
/// Line search helper for globally convergent Newton method
//=======================================================================
 void line_search(const Vector<double>& x_old, 
                  const double half_residual_squared_old, 
                  const Vector<double>& gradient, 
                  ResidualFctPt residual_fct, 
                  const Vector<double>& params,
                  Vector<double>& newton_dir,
                  Vector<double>& x, 
                  double& half_residual_squared,
                  const double stpmax,
                  bool& check)
 {
  const double ALF=1.0e-4;
  double TOLX=1.0e-16;
  double f2=0.0;
  double alam2=0.0;
  double tmplam;
  unsigned n=x_old.size();
  check=false;
  double sum=0.0;
  for (unsigned i=0;i<n;i++)
   {
    sum += newton_dir[i]*newton_dir[i];
   }
  sum=sqrt(sum);
  if (sum > stpmax)
   {
    for (unsigned i=0;i<n;i++)
     {
      newton_dir[i] *= stpmax/sum;
     }
   }
  double slope=0.0;
  for (unsigned i=0;i<n;i++)
   {
    slope += gradient[i]*newton_dir[i];
   }
  if (slope >= 0.0)
   {
    std::ostringstream error_stream;
    error_stream << "Roundoff problem in lnsrch: slope=" << slope << "\n";
    throw OomphLibError(error_stream.str(),
                        "BlackBoxFDNewtonSolver::line_search()",
                        OOMPH_EXCEPTION_LOCATION);
   }
  double test=0.0;
  for (unsigned i=0;i<n;i++)
   {
    double temp=std::fabs(newton_dir[i])/std::max(std::fabs(x_old[i]),1.0);  
    if (temp > test) test=temp;
   }
  double alamin=TOLX/test;
  double alam=1.0;
  for (;;) // hmmm
   {
    for (unsigned i=0;i<n;i++)
     {
      x[i]=x_old[i]+alam*newton_dir[i];
     }
    
   // Evaluate current residuals
    Vector<double> residuals(n);
   residual_fct(params,x,residuals);
    half_residual_squared=0.0;
    for (unsigned i=0;i<n;i++)
     {
      half_residual_squared+=residuals[i]*residuals[i];
     }
    half_residual_squared*=0.5;

    if (alam < alamin)
     {
      for (unsigned i=0;i<n;i++) x[i]=x_old[i];
      check=true;
      oomph_info << "Warning: Line search converged on x only!\n";
      return;
     } 
    else if (half_residual_squared <= half_residual_squared_old+ALF*alam*slope)
     {
      return;
     }
    else
     {
      if (alam == 1.0)
       {
        tmplam = -slope/(2.0*(half_residual_squared-
                              half_residual_squared_old-slope));
       }
      else
       {
        double rhs1=half_residual_squared-half_residual_squared_old-alam*slope;
        double rhs2=f2-half_residual_squared_old-alam2*slope;
        double a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
        double b=(-alam2*rhs1/(alam*alam)+alam*rhs2/
                  (alam2*alam2))/(alam-alam2);
        if (a == 0.0)
         {
          tmplam = -slope/(2.0*b);
         }
        else
         {
          double disc=b*b-3.0*a*slope;
          if (disc < 0.0)
           { 
            tmplam=0.5*alam;
           }       
          else if (b <= 0.0)
           {
            tmplam=(-b+sqrt(disc))/(3.0*a);
           }
          else
           {
            tmplam=-slope/(b+sqrt(disc));
           }
         }
        if (tmplam>0.5*alam)
         {
          tmplam=0.5*alam;
         }
       }
     }
    alam2=alam;
    f2 = half_residual_squared;
    alam=std::max(tmplam,0.1*alam);
   }
 }
}

//======================================================================
/// \short Set output directory (we try to open a file in it
/// to see if the directory exists -- if it doesn't we'll
/// issue a warning -- or, if directory_must_exist()==true,
/// die by throwing and OomphLibError
//======================================================================
void DocInfo::set_directory(const std::string& directory_)
{
 // Try to open a file in output directory
 std::ostringstream filename;
 filename << directory_ << "/.dummy_check.dat";
 std::ofstream some_file;
 some_file.open(filename.str().c_str());
 if (!some_file.is_open())
  {
   //Construct the error message
   std::string error_message = "Problem opening output file.\n";
   error_message += "I suspect you haven't created the output directory ";
   error_message += directory_;
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
 Directory=directory_;
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

 /// Map to indicate an input flag as having been set
 std::map<std::string,bool> Specified_command_line_flag;

 /// Map to associate an input flag with a double -- specified via pointer
 std::map<std::string,std::pair<bool,double*> > 
 Specified_command_line_double_pt;

 /// Map to associate an input flag with an int -- specified via pointer
 std::map<std::string,std::pair<bool,int*> > 
 Specified_command_line_int_pt;

 /// Map to associate an input flag with an unsigned -- specified via pointer
 std::map<std::string,std::pair<bool,unsigned*> > 
 Specified_command_line_unsigned_pt;

 /// Map to associate an input flag with a string -- specified via pointer
 std::map<std::string,std::pair<bool,std::string*> > 
 Specified_command_line_string_pt;

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
  std::stringstream str;
  for (int i=1;i<CommandLineArgs::Argc;i++)
   {
    str << CommandLineArgs::Argv[i] << " ";
   }
  oomph_info << str.str() << std::endl;
 }


 /// Specify possible argument-free command line flag
 void specify_command_line_flag(const std::string& command_line_flag)
 {
  Specified_command_line_flag[command_line_flag]=false;
 }

 /// \short Specify possible command line flag that specifies a double, accessed
 /// via pointer
 void specify_command_line_flag(const std::string& command_line_flag,
                                double* arg_pt)
 {
  Specified_command_line_double_pt[command_line_flag]=
   std::make_pair(false,arg_pt);
 }

 /// \short Specify possible command line flag that specifies an int, accessed
 /// via pointer
 void specify_command_line_flag(const std::string& command_line_flag,
                                int* arg_pt)
 {
  Specified_command_line_int_pt[command_line_flag]=
   std::make_pair(false,arg_pt);
 }

 /// \short Specify possible command line flag that specifies an unsigned, 
 /// accessed via pointer
 void specify_command_line_flag(const std::string& command_line_flag,
                                unsigned* arg_pt)
 {
  Specified_command_line_unsigned_pt[command_line_flag]=
   std::make_pair(false,arg_pt);
 }

 /// \short Specify possible command line flag that specifies a string, 
 /// accessed via pointer
 void specify_command_line_flag(const std::string& command_line_flag,
                                std::string* arg_pt)
 {
  Specified_command_line_string_pt[command_line_flag]=
   std::make_pair(false,arg_pt);
 }


 /// \short Check if command line flag has been set (value will have been
 /// assigned directly
 bool command_line_flag_has_been_set(const std::string& flag)
 {
  for (std::map<std::string,bool>::iterator it=
        Specified_command_line_flag.begin();
       it!=Specified_command_line_flag.end();it++)
   {
    if ((*it).first==flag) return (*it).second;
   }
  for (std::map<std::string,std::pair<bool,double*> >::iterator it=
        Specified_command_line_double_pt.begin();
       it!=Specified_command_line_double_pt.end();it++)
   {
    if ((*it).first==flag) return ((*it).second).first;
   }
  for (std::map<std::string,std::pair<bool,int*> >::iterator it=
        Specified_command_line_int_pt.begin();
       it!=Specified_command_line_int_pt.end();it++)
   {
    if ((*it).first==flag) return ((*it).second).first;
   }
  for (std::map<std::string,std::pair<bool,unsigned*> >::iterator it=
        Specified_command_line_unsigned_pt.begin();
       it!=Specified_command_line_unsigned_pt.end();it++)
   {
    if ((*it).first==flag) return ((*it).second).first;
   }
  for (std::map<std::string,std::pair<bool,std::string*> >::iterator it=
        Specified_command_line_string_pt.begin();
       it!=Specified_command_line_string_pt.end();it++)
   {
    if ((*it).first==flag) return ((*it).second).first;
   }
  return false;
 }

 /// Document specified command line flags
 void doc_specified_flags()
 {
  oomph_info << std::endl;
  oomph_info << "Specified (and recognised) command line flags:\n" ;
  oomph_info << "----------------------------------------------\n";
  for (std::map<std::string,bool>::iterator it=
        Specified_command_line_flag.begin();
       it!=Specified_command_line_flag.end();it++)
   {
    if ((*it).second) oomph_info << (*it).first << std::endl;
   }
  for (std::map<std::string,std::pair<bool,double*> >::iterator it=
        Specified_command_line_double_pt.begin();
       it!=Specified_command_line_double_pt.end();it++)
   {
    if (((*it).second).first) oomph_info << (*it).first << " " 
                                         << *(((*it).second).second) 
                                         << std::endl;
   }
  for (std::map<std::string,std::pair<bool,int*> >::iterator it=
        Specified_command_line_int_pt.begin();
       it!=Specified_command_line_int_pt.end();it++)
   {
    if (((*it).second).first) oomph_info << (*it).first << " " 
                                         << *(((*it).second).second) 
                                         << std::endl;
   }
  for (std::map<std::string,std::pair<bool,unsigned*> >::iterator it=
        Specified_command_line_unsigned_pt.begin();
       it!=Specified_command_line_unsigned_pt.end();it++)
   {
    if (((*it).second).first) oomph_info << (*it).first << " " 
                                         << *(((*it).second).second)
                                         << std::endl;
   }
  for (std::map<std::string,std::pair<bool,std::string*> >::iterator it=
        Specified_command_line_string_pt.begin();
       it!=Specified_command_line_string_pt.end();it++)
   {
    if (((*it).second).first) oomph_info << (*it).first << " " 
                                         << *(((*it).second).second)
                                         << std::endl;
   }
  oomph_info << std::endl;
 }


 /// Document available command line flags
 void doc_available_flags()
 {
  oomph_info << std::endl;
  oomph_info << "Available command line flags:\n" ;
  oomph_info << "-----------------------------\n";
  for (std::map<std::string,bool>::iterator it=
        Specified_command_line_flag.begin();
       it!=Specified_command_line_flag.end();it++)
   {
    oomph_info << (*it).first << std::endl;
   }
  for (std::map<std::string,std::pair<bool,double*> >::iterator it=
        Specified_command_line_double_pt.begin();
       it!=Specified_command_line_double_pt.end();it++)
   {
    oomph_info << (*it).first << " <double> " << std::endl;
   }
  for (std::map<std::string,std::pair<bool,int*> >::iterator it=
        Specified_command_line_int_pt.begin();
       it!=Specified_command_line_int_pt.end();it++)
   {
    oomph_info << (*it).first << " <int> " << std::endl;
   }
  for (std::map<std::string,std::pair<bool,unsigned*> >::iterator it=
        Specified_command_line_unsigned_pt.begin();
       it!=Specified_command_line_unsigned_pt.end();it++)
   {
    oomph_info << (*it).first << " <unsigned> " << std::endl;
   }
  for (std::map<std::string,std::pair<bool,std::string*> >::iterator it=
        Specified_command_line_string_pt.begin();
       it!=Specified_command_line_string_pt.end();it++)
   {
    oomph_info << (*it).first << " <string> " << std::endl;
   }
  oomph_info << std::endl;
 }



 /// Helper function to check if command line index is legal
 void check_arg_index(const int& argc,const int& arg_index)
 {
  if (arg_index>=argc)
   {    
    output();
    doc_available_flags();
    std::stringstream error_stream;
    error_stream << "Tried to read more command line arguments than\n"
                 << "specified. This tends to happen if a required argument\n"
                 << "to a command line flag was omitted, e.g. by running \n\n"
                 << "     ./a.out -some_double \n\n rather than\n\n"
                 << "     ./a.out -some_double 1.23 \n\n"
                 << "To aid the debugging I've output the available\n"
                 << "command line arguments above.\n";
    throw OomphLibError(
     error_stream.str(),
     "check_arg_index()",
     OOMPH_EXCEPTION_LOCATION);
   }
 }



 /// \short Parse command line, check for recognised flags and assign 
 /// associated values
 void parse_and_assign(int argc, char *argv[])
 {

  //Keep looping over command line arguments
  int arg_index = 1;
  while (arg_index < argc)
   {

    bool found_match=false;

    if ( (strcmp(argv[arg_index], "-help") == 0 ) ||
         (strcmp(argv[arg_index], "--help")== 0 ) )
     {
      oomph_info << "NOTE: You entered --help or -help on the commmand line\n";
      oomph_info << "so I'm going to tell you about this code's\n";
      oomph_info << "available command line flags and then return.\n";
      doc_available_flags();

#ifdef OOMPH_HAS_MPI
      int flag;
      MPI_Initialized(&flag);
      if (flag!=0) MPI_Helpers::finalize();
#endif
      oomph_info << "Shutting down...\n";
      exit(0);

     }


    //Check if the flag has been previously specified as a simple argument free
    //command line argument 
    for (std::map<std::string,bool>::iterator it=
          Specified_command_line_flag.begin();
         it!=Specified_command_line_flag.end();it++)
     {
      if ((*it).first==argv[arg_index])
       {
        Specified_command_line_flag[argv[arg_index]]=true;
        found_match=true;      
        break;
       }
     }
     
    if (!found_match)
     {      
      //Check if the flag has been previously specified as a 
      //command line argument that specifies (and is followed by) a double
      for (std::map<std::string,std::pair<bool,double*> >::iterator it=
            Specified_command_line_double_pt.begin();
           it!=Specified_command_line_double_pt.end();it++)
       {
        if ((*it).first==argv[arg_index])
         {
          //Next command line argument specifies the double. Set it via
          //the pointer
          arg_index++;
          check_arg_index(argc,arg_index);
          ((*it).second).first=true;
          *(((*it).second).second)=atof(argv[arg_index]);
          found_match=true;      
          break;
         }
       }
     }


    if (!found_match)
     {      
      //Check if the flag has been previously specified as a 
      //command line argument that specifies (and is followed by) an int
      for (std::map<std::string,std::pair<bool,int*> >::iterator it=
            Specified_command_line_int_pt.begin();
           it!=Specified_command_line_int_pt.end();it++)
       {
        if ((*it).first==argv[arg_index])
         {
          //Next command line argument specifies the int. Set it via
          //the pointer
          arg_index++;
          check_arg_index(argc,arg_index);
          ((*it).second).first=true;
          *(((*it).second).second)=atoi(argv[arg_index]);
          found_match=true;      
          break;
         }
       }
     }


    if (!found_match)
     {      
      //Check if the flag has been previously specified as a 
      //command line argument that specifies (and is followed by) an unsigned
      for (std::map<std::string,std::pair<bool,unsigned*> >::iterator it=
            Specified_command_line_unsigned_pt.begin();
           it!=Specified_command_line_unsigned_pt.end();it++)
       {
        if ((*it).first==argv[arg_index])
         {
          //Next command line argument specifies the unsigned. Set it via
          //the pointer
          arg_index++;
          check_arg_index(argc,arg_index);
          ((*it).second).first=true;
          *(((*it).second).second)=unsigned(atoi(argv[arg_index]));
          found_match=true;      
          break;
         }
       }
     }


    if (!found_match)
     {      
      //Check if the flag has been previously specified as a 
      //command line argument that specifies (and is followed by) a string
      for (std::map<std::string,std::pair<bool,std::string*> >::iterator it=
            Specified_command_line_string_pt.begin();
           it!=Specified_command_line_string_pt.end();it++)
       {
        if ((*it).first==argv[arg_index])
         {
          //Next command line argument specifies the string. Set it via
          //the pointer
          arg_index++;
          check_arg_index(argc,arg_index);
          ((*it).second).first=true;
          *(((*it).second).second)=argv[arg_index];
          found_match=true;      
          break;
         }
       }
     }


    // Oh dear, we still haven't found the argument in the list.
    // Maybe it was specified wrongly -- issue warning.
    if (!found_match)
     {      

      //Construct the error message
      std::string error_message = "Command line argument\n\n";
      error_message += argv[arg_index];
      error_message += "\n\nwas not recognised. This may be legal\n";
      error_message += "but seems sufficiently suspicious to flag up.\n";
      error_message += "If it should have been recognised, make sure you\n";
      error_message += "used the right number of \"-\" signs...\n";

      //Create an Oomph Lib warning
      OomphLibWarning(error_message,
                      "CommandLineArgs::parse_and_assign()",
                      OOMPH_EXCEPTION_LOCATION);
      
     }



    arg_index++;
   }
 }



 /// \short Parse previously specified command line, check for 
 /// recognised flags and assign associated values
 void parse_and_assign()
 {
  parse_and_assign(Argc,Argv);
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

//=============================================================================
/// Initialize mpi. If optional boolean flag is set to false, we use
/// MPI_COMM_WORLD itself as oomph-lib's communicator. Defaults to true.
//=============================================================================
 void MPI_Helpers::init(int argc, char **argv, 
                        const bool& make_duplicate_of_mpi_comm_world)
{
#ifdef OOMPH_HAS_MPI
 // call mpi int
 MPI_Init(&argc,&argv);


 // By default, create the oomph-lib communicator using MPI_Comm_dup so that
 // the communicator has the same group of processes but a new context
 MPI_Comm oomph_comm_world=MPI_COMM_WORLD;
 if (make_duplicate_of_mpi_comm_world)
  {
   MPI_Comm_dup(MPI_COMM_WORLD,&oomph_comm_world); 
  }

 if (MPI_COMM_WORLD!=oomph_comm_world)
  {
   oomph_info << "Oomph-lib communiator is a duplicate of MPI_COMM_WORLD\n";
  }
 else
  {
   oomph_info << "Oomph-lib communiator is MPI_COMM_WORLD\n";
  }
 
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
#else
 // create a serial communicator
 Communicator_pt = new OomphCommunicator;
#endif
 MPI_has_been_initialised=true;
}

//=============================================================================
/// finalize mpi
//=============================================================================
void MPI_Helpers::finalize()
{
 // delete the communicator
 delete Communicator_pt;
 
 // and call MPI_Finalize
#ifdef OOMPH_HAS_MPI
 MPI_Finalize();
#endif
}

//=============================================================================
/// access to the global oomph-lib communicator
//=============================================================================
OomphCommunicator* MPI_Helpers::communicator_pt()
{
#ifdef PARANOID
 if (!MPI_has_been_initialised)
  {
   std::ostringstream error_message_stream;
   error_message_stream 
    << "MPI has not been initialised.\n Call MPI_Helpers::init(...)";
   throw OomphLibError(error_message_stream.str(),
                       "MPI_Helpers::communicator_pt()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif
 return Communicator_pt;
}

bool MPI_Helpers::MPI_has_been_initialised = false;
OomphCommunicator* MPI_Helpers::Communicator_pt = 0;


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
   if(MPI_Helpers::mpi_has_been_initialised())
    {
     return MPI_Wtime();
    }
   else
#endif
    {
     time_t t = clock();
     return  double(t) / double(CLOCKS_PER_SEC);
    }
  }
}//end of namespace TimingHelpers


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////



//===============================================================
/// Namespace with helper functions to assess total memory usage
/// on the fly using system() -- details are very machine specific! This just
/// provides the overall machinery with default settings for
/// our own (linux machines). Uses the system command
/// to spawn a command that computes the total memory usage
/// on the machine where this is called.  [Disclaimer: works 
/// on my machine(s) -- no guarantees for any other platform; 
/// linux or not. MH]
//===============================================================
namespace MemoryUsage
{
 
  /// \short String containing system command that obtains memory usage
  /// of all processes.
  /// Default assigment for linux. [Disclaimer: works on my machine(s) --
  /// no guarantees for any other platform; linux or not. MH]
  std::string My_memory_usage_system_string="ps aux";

 /// \short Bool allowing quick bypassing of ALL operations related
 /// to memory usage monitoring -- this allows the code to remain
 /// "instrumented" without incurring the heavy penalties associated
 /// with the system calls and i/o. Default setting: false.
 bool Bypass_all_memory_usage_monitoring=false;
 
  /// \short String containing name of file in which we document
  /// my memory usage -- you may want to change this to allow different
  /// processors to write to separate files (especially in mpi 
  /// context). Note that file is appended to 
  /// so it ought to be emptied (either manually or by calling
  /// helper function empty_memory_usage_file()
  std::string My_memory_usage_filename="my_memory_usage.dat";

  /// \short Function to empty file that records my memory usage in
  /// file whose name is specified by My_memory_usage_filename.
  void empty_my_memory_usage_file()
  {
   // bail out straight away?
   if (Bypass_all_memory_usage_monitoring) return;

   // Open without appending and write header
   std::ofstream the_file;
   the_file.open(My_memory_usage_filename.c_str());
   the_file << "# My memory usage: \n";
   the_file.close();
  }


  
  /// \short Doc my memory usage, prepended by string (which allows
  /// identification from where the function is called, say) that records 
  /// memory usage in file whose name is specified by My_memory_usage_filename.
  /// Data is appended to that file; wipe it with empty_my_memory_usage_file().
  void doc_my_memory_usage(const std::string& prefix_string)
  {
   // bail out straight away?
   if (Bypass_all_memory_usage_monitoring) return;

   // Write prefix
   std::ofstream the_file;
   the_file.open(My_memory_usage_filename.c_str(),std::ios::app);
   the_file << prefix_string << " "; 
   the_file.close();
  
   // Sync all processors if in parallel
#ifdef OOMPH_HAS_MPI
   if (MPI_Helpers::mpi_has_been_initialised())
    {
     MPI_Barrier(MPI_Helpers::communicator_pt()->mpi_comm());
    }
#endif

   // Process memory usage command and write to file
   std::stringstream tmp;
   tmp << My_memory_usage_system_string 
       << " | awk '{if ($2==" << getpid() << "){print $4}}' >> " 
       << My_memory_usage_filename;
   int success=system(tmp.str().c_str());

   // Dummy command to shut up compiler warnings
   success+=1;
  }

 /// \short String containing system command that obtains total memory usage.
 /// Default assigment for linux. [Disclaimer: works on my machine(s) --
 /// no guarantees for any other platform; linux or not. MH]
 std::string Total_memory_usage_system_string=
  "ps aux | awk 'BEGIN{sum=0}{sum+=$4}END{print sum}'";
 
 /// \short String containing name of file in which we document total
 /// memory usage -- you may want to change this to allow different
 /// processors to write to separate files (especially in mpi 
 /// context). Note that file is appended to 
 /// so it ought to be emptied (either manually or by calling
 /// helper function empty_memory_usage_file()
 std::string Total_memory_usage_filename="memory_usage.dat";
 
 /// \short Function to empty file that records total memory usage in
 /// file whose name is specified by Total_memory_usage_filename.
 void empty_total_memory_usage_file()
 {
  // bail out straight away?
  if (Bypass_all_memory_usage_monitoring) return;

  // Open without appending and write header
  std::ofstream the_file;
  the_file.open(Total_memory_usage_filename.c_str());
  the_file << "# Total memory usage: \n";
  the_file.close();
 }
 
 /// \short Doc total memory usage, prepended by string (which allows
 /// identification from where the function is called, say) that records 
 /// memory usage in file whose name is specified by Total_memory_usage_filename.
 /// Data is appended to that file; wipe it with empty_memory_usage_file().
 void doc_total_memory_usage(const std::string& prefix_string)
 {
  // bail out straight away?
  if (Bypass_all_memory_usage_monitoring) return;

  // Write prefix
  std::ofstream the_file;
  the_file.open(Total_memory_usage_filename.c_str(),std::ios::app);
  the_file << prefix_string << " "; 
  the_file.close();
  
  // Sync all processors if in parallel
#ifdef OOMPH_HAS_MPI
  if (MPI_Helpers::mpi_has_been_initialised())
   {
    MPI_Barrier(MPI_Helpers::communicator_pt()->mpi_comm());
   }
#endif
  
  // Process memory usage command and write to file
  std::stringstream tmp;
  tmp << Total_memory_usage_system_string << " >> " 
      << Total_memory_usage_filename;
  int success=system(tmp.str().c_str());
  
  // Dummy command to shut up compiler warnings
  success+=1;
 }



 /// \short Function to empty file that records total and local memory usage
 /// in appropriate files
 void empty_memory_usage_files()
 {
  // bail out straight away?
  if (Bypass_all_memory_usage_monitoring) return;

  empty_my_memory_usage_file();
  empty_total_memory_usage_file();
  empty_top_file();
 }
 
 /// \short Doc total and local memory usage, prepended by string (which allows
 /// identification from where the function is called, say)
 void doc_memory_usage(const std::string& prefix_string)
 {
  // bail out straight away?
  if (Bypass_all_memory_usage_monitoring) return;

  doc_my_memory_usage(prefix_string);
  doc_total_memory_usage(prefix_string);
 }




 /// \short String containing system command that runs "top" (or equivalent)
 /// "indefinitely" and writes to file specified in Top_output_filename.
 /// Default assigment for linux. [Disclaimer: works on my machine(s) --
 /// no guarantees for any other platform; linux or not. MH]
 std::string Top_system_string="while true; do top -b -n 2 ; done  ";
 
 /// \short  String containing name of file in which we document "continuous"
 /// output from "top" (or equivalent)-- you may want to change this to 
 /// allow different processors to write to separate files (especially in mpi 
 /// context). Note that file is appended to 
 /// so it ought to be emptied (either manually or by calling
 /// helper function empty_top_file()
         std::string Top_output_filename="top_output.dat";
 
 /// \short Function to empty file that records continuous output from top in
 /// file whose name is specified by Top_output_filename
 void empty_top_file()
 {
  // bail out straight away?
  if (Bypass_all_memory_usage_monitoring) return;
  
  // Open without appending and write header
  std::ofstream the_file;
  the_file.open(Top_output_filename.c_str());
  the_file << "# Continuous output from top obtained with: \n";
  the_file << "# " << Top_system_string << "\n";
  the_file.close();
 }
 

 /// \short Start running top continously and output (append) into 
 /// file specified by Top_output_filename. Wipe that file  with 
 /// empty_top_file() if you wish. Note that this is again quite linux specific
 /// and unlikely to work on other operating systems.
 /// Insert optional comment into output file before starting.
 void run_continous_top(const std::string& comment)
  {
   // bail out straight away?
   if (Bypass_all_memory_usage_monitoring) return;
  
   // Sync all processors if in parallel
   std::string modifier="";

#ifdef OOMPH_HAS_MPI
   if (MPI_Helpers::mpi_has_been_initialised())
    {
     MPI_Barrier(MPI_Helpers::communicator_pt()->mpi_comm());     
     std::stringstream tmp;
     tmp << "_proc" << MPI_Helpers::communicator_pt()->my_rank();
     modifier=tmp.str();
    }
#endif

   // Process memory usage command and write to file
   std::stringstream tmp;

   // String stream seems unhappy about direct insertions of these
   std::string backslash="\\";
   std::string dollar="$";

   // Create files that spawn and kill continuous top
   tmp << "echo \"#/bin/bash\" > run_continuous_top" << modifier 
       << ".bash; "
       << "echo \" echo " << backslash << "\" kill " 
       << backslash << dollar 
       << backslash << dollar << " " << backslash << "\" > kill_continuous_top" 
       << modifier << ".bash; chmod a+x kill_continuous_top" 
       << modifier << ".bash; "  << Top_system_string 
       << " \" >> run_continuous_top" << modifier 
       << ".bash; chmod a+x run_continuous_top" << modifier << ".bash ";
   int success=system(tmp.str().c_str());

   // Add comment to annotate?
   if (comment!="")
    {
     insert_comment_to_continous_top(comment);
    }

   // Start spawning
   std::stringstream tmp2;
   tmp2 << "./run_continuous_top" << modifier << ".bash  >> " 
        << Top_output_filename << " & ";
   success=system(tmp2.str().c_str());

   // Dummy command to shut up compiler warnings
   success+=1;
  }



 /// \short Stop running top continously. Note that this is 
 /// again quite linux specific and unlikely to work on other operating 
 /// systems. Insert optional comment into output file before stopping.
 void stop_continous_top(const std::string& comment)
  {
   // bail out straight away?
   if (Bypass_all_memory_usage_monitoring) return;
  
   // Sync all processors if in parallel
   std::string modifier="";

#ifdef OOMPH_HAS_MPI
   if (MPI_Helpers::mpi_has_been_initialised())
    {
     MPI_Barrier(MPI_Helpers::communicator_pt()->mpi_comm());     
     std::stringstream tmp;
     tmp << "_proc" << MPI_Helpers::communicator_pt()->my_rank();
     modifier=tmp.str();
    }
#endif

   // Add comment to annotate?
   if (comment!="")
    {
     insert_comment_to_continous_top(comment);
    }
 
   // Kill
   std::stringstream tmp2;
   tmp2 << "./kill_continuous_top" << modifier << ".bash  >> " 
        << Top_output_filename << " & ";
   int success=system(tmp2.str().c_str());

   // Dummy command to shut up compiler warnings
   success+=1;
  }


 /// \short Insert comment into running continous top output
 void insert_comment_to_continous_top(const std::string& comment)
 {
  std::stringstream tmp;
  tmp << " echo \"" << comment << "\"  >> " 
      << Top_output_filename;
  int success=system(tmp.str().c_str());
  
  // Dummy command to shut up compiler warnings
  success+=1;
 }



} // end of namespace MemoryUsage


}
