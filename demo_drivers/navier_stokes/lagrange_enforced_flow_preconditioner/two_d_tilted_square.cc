//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision$
//LIC//
//LIC// $LastChangedDate$
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
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
//LIC//=====================================================================
//Driver for 2D tilted rectangle

#include <fenv.h>
#include <sstream>
#include <iomanip>
#include <ios>


// Generic includes
#include "generic.h"
#include "navier_stokes.h"

// The 2D mesh
#include "meshes/simple_rectangular_quadmesh.h"
#include "meshes/rectangular_quadmesh.h"


#include "meshes/simple_cubic_mesh.h"
#include "meshes/simple_cubic_tet_mesh.h"

using namespace std;

using namespace oomph;



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#ifdef OOMPH_HAS_HYPRE

//==========================================================================
/// helper method for the block diagonal F block preconditioner to allow 
/// hypre to be used for as a subsidiary block preconditioner
//==========================================================================
namespace Hypre_Subsidiary_Preconditioner_Helper
{
 Preconditioner* set_hypre_preconditioner()
 {
  return new HyprePreconditioner;
 }
}

#endif

 
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////
///
/// Namespace for all things related to the problem where the domain is a 
/// unit square.
///
////////////////////////////////////////////////////////////////////////////
namespace SquareLagrange
{
// Prob id, set by main method
//  const int* Prob_id_pt = 0;

//  std::string Prob_str = "";
//  std::string Ang_deg_str = "";
//  std::string Noel_str = "";


  ///////////////////////
  // Domain dimensions.//
  ///////////////////////
  //
  // This is a square domain: x,y \in [0,1]
  //

  // Min and max x value respectively.
  static const double X_min = 0.0;
  static const double X_max = 1.0;

  // Min and max y value respectively.
  static const double Y_min = 0.0;
  static const double Y_max = 1.0;

  // The length in the x and y direction respectively.
  static const double Lx = X_max - X_min;
  static const double Ly = Y_max - Y_min;

////////////////////////////////////////////////////////////////////////////

  // CL - set directly from the command line.
  // To set from CL - a CL value is set, this is changed depending on that
  // value.
  //
  // Problem parameter overview:
  //
  // // Solvers:
  //
  //
  // F_ns + L^T inv(W) L | B^T
  // --------------------------
  //                     | W
  //
  // W = 0 (SuperLU)
  // NS_solver = 0 (SuperLU) or 1 (LSC)
  // 
  // If NS_solver = 1, then we have:
  //
  // | F | B^T |   |
  // |----------   |
  // |   |-M_s |   |
  // |-------------|
  // |         | W |
  //
  // F_solver = 0 (SuperLU) or 1 (AMG)
  // P_solver = 0 (SuperLU) or 1 (AMG)
  // 

  // These are self explanatory:
  double Ang_deg = 30.0; //CL, Angle in degrees
  double Ang_rad = -1.0; // set by formula
  unsigned Noel = 4; //CL, Number of elements in 1D
  // the default is the norm of the momentum block.

  // Function to turn degrees into radians.
  inline double degtorad(const double& ang_deg)
  {
    return ang_deg * (MathematicalConstants::Pi / 180.0);
  }


  inline void setup_commandline_flags()
  {
    CommandLineArgs::specify_command_line_flag("--ang", &Ang_deg);

    CommandLineArgs::specify_command_line_flag("--noel", &Noel);
  }

//  inline void set_ang_str()
//  {
//    if(Prob_id_pt == 0)
//    {
//      std::ostringstream err_msg;
//      err_msg << "Oh dear, Prob_id_pt is null. Please set this in main().\n"
//        << "This should be stored in NSPP::Prob_id, and set by cmd via\n"
//        << "--prob_id \n"; 
//      throw OomphLibError(err_msg.str(),
//          OOMPH_CURRENT_FUNCTION,
//          OOMPH_EXCEPTION_LOCATION);
//    }
//
//    // If this is the vanilla problem, we set the angle as -1 and set the
//    // string as "A_". This would indicate that no angle is used.
//    // Furthermore, we ensure that no --ang is set.
//    if(Prob_str.compare("SqVa") == 0)
//    {
//      if(CommandLineArgs::command_line_flag_has_been_set("--ang"))
//      {
//        std::ostringstream err_msg;
//        err_msg << "prob_id is 88, doing vanilla LSC with no tilt.\n"
//          << "But you have set --ang, please do not set this."; 
//        throw OomphLibError(err_msg.str(),
//            OOMPH_CURRENT_FUNCTION,
//            OOMPH_EXCEPTION_LOCATION);
//      }
//      Ang = -1.0;
//      // Now we set the Ang_deg_str.
//      std::ostringstream strs;
//      strs << "A_";
//      Ang_deg_str = strs.str();
//    }
//    else
//    // This problem requires tilting, thus we set the Ang and Ang_deg_str.
//    {
//      // But first we check that --ang has been set.
//      // Check that Ang has been set.
//      if(!CommandLineArgs::command_line_flag_has_been_set("--ang"))
//      {
//        std::ostringstream err_msg;
//        err_msg << "Angle has not been set. Set (in degrees) with: \n"
//          << "--ang \n"; 
//        throw OomphLibError(err_msg.str(),
//            OOMPH_CURRENT_FUNCTION,
//            OOMPH_EXCEPTION_LOCATION);
//      }
//
//      // Now we need to convert Ang_deg into radians.
//      Ang = Ang_deg * (MathematicalConstants::Pi / 180.0);
//
//      // Now we set the Ang_deg_str.
//      std::ostringstream strs;
//      strs << "A" << Ang_deg;
//      Ang_deg_str = strs.str();
//    }
//  } // set_ang_str

//  inline void set_noel_str()
//  {
//    // Set Noel_str, used for book keeping.
//    if(CommandLineArgs::command_line_flag_has_been_set("--noel"))
//    {
//      std::ostringstream strs;
//      strs << "N" <<Noel;
//      Noel_str = strs.str();
//    }
//    else
//    {
//      std::ostringstream err_msg;
//      err_msg << "Please supply the number of elements in 1D using --noel.\n"
//        << std::endl;
//      throw OomphLibError(err_msg.str(),
//          OOMPH_CURRENT_FUNCTION,
//          OOMPH_EXCEPTION_LOCATION);
//    }
//  }

//  inline void generic_setup()
//  {
//    set_ang_str();
//    set_noel_str();
//  }

//  inline std::string prob_str()
//  {
//    set_prob_str();
//    return Prob_str;
//  }

//  inline std::string ang_deg_str()
//  {
//    set_ang_str();
//    return Ang_deg_str;
//  }

//  inline std::string noel_str()
//  {
//    set_noel_str();
//    return Noel_str;
//  }

//  inline std::string create_label()
//  {
//    std::string label = prob_str() + ang_deg_str() + noel_str();
//    return label; 
//  } // inlined function create_label

} // Namespace SquareLagrange

//=============================================================================
/// Namespace to hold variables common to all Navier Stokes problems.
/// Contains the
//=============================================================================
//namespace NavierStokesProblemParameters
//{
//
//  typedef std::map<int,std::string>::iterator int_string_map_it_type;
//
//
//  const static int Solver_type_DIRECT_SOLVE = 0;
//  const static int Solver_type_OOMPHLIB_GMRES = 1;
//  const static int Solver_type_TRILINOS_GMRES = 2;
//
//  // To fill
//  std::map<int,std::string> valid_solver_type_map;
//
//////////////////////////////////////
//
//  // STEADY wins, even if --dt, --time_start and --time_end is set,
//  // if Time_type == 0, then the time parameters will be ignored, steady
//  // state will be attempted. This is because no timing parameters is 
//  // required for steady state, thus we can ignore them.
//  const static int Time_type_STEADY = 0;
//
//  // Adaptive time stepping takes second precedence, this is because only
//  // two of the three time variables needs to be set 
//  // (--time_start and --time_end).
//  const static int Time_type_ADAPT = 1;
//
//  // This takes last precedence. We have done it this way so you can have
//  // all three time variables set and just change the Time_type to switch 
//  // between the different time stepping states.
//  const static int Time_type_FIXED = 2;
//
//  // This is what we set:
//  int Time_type = -1;
//
/////////////////////////////////
////  const static int MeshType_TETRAHEDRAL = 0;
////  const static int MeshType_HEXAHEDRAL = 1;
////
////  int Mesh_type = -1;
//
//  //////////////////////////////////////////
//
//  // From Commandline
//  int Solver_type = -1;
////  int Prob_id = -1;
//  bool Distribute_problem = false;
//  int Vis = -1;
//  double Rey = -1.0;
//  double Re_invFr = -1.0;
//
////  double Rey_start = -1.0;
////  double Rey_incre = -1.0;
////  double Rey_end = -1.0;
//
//  int Max_solver_iteration = -1;
//
//  double Delta_t = -1.0;
//
//
////  std::string Soln_dir_str = "";
////  std::string Itstime_dir_str = "";
//
//  // Inferred from commandline:
////  bool Doc_soln = false;
////  std::string Label_str = "";
//
//  // From code:
//  DocLinearSolverInfo* Doc_linear_solver_info_pt = 0;
//
//  // This is set by the main function.
//  // And --dt is set by the commandline. Then we time step from
//  // Time_start to Time_end in time step sizes of Delta_t
//  double Time_start = -1.0;
//  double Time_end = -1.0;
//
//  // Additional stuff required for quarter circle, but putting it
//  // here since it may be required elsewhere
//  
//  // Functional body force
////  void body_force(const double& time, const Vector<double>& x,
////                   Vector<double>& result)
////  {
////    result[0] = 0.0;
////    result[1] = -Re_invFr;
////  }
////
////  void zero_body_force(const double& time, const Vector<double>& x,
////                       Vector<double>& result)
////  {
////    result[0] = 0.0;
////    result[1] = 0.0;
////  }
//
////  Vector<double> Gravity(2);
//
////  inline void setup_commandline_flags()
////  {
////    CommandLineArgs::specify_command_line_flag("--dist_prob");
////
////    // A problem ID, there are eight different types of problems.
////    // Check the header file.
////    CommandLineArgs::specify_command_line_flag("--prob_id",&Prob_id);
////
////    // Flag to output the solution.
////    CommandLineArgs::specify_command_line_flag("--doc_soln", 
////        &Soln_dir_str);
////    CommandLineArgs::specify_command_line_flag("--visc", 
////        &Vis);
////    CommandLineArgs::specify_command_line_flag("--rey", &Rey);
////    CommandLineArgs::specify_command_line_flag("--rey_start", &Rey_start);
////    CommandLineArgs::specify_command_line_flag("--rey_incre", &Rey_incre);
////    CommandLineArgs::specify_command_line_flag("--rey_end", &Rey_end);
////    CommandLineArgs::specify_command_line_flag("--max_solver_iter", 
////                                               &Max_solver_iteration);
////    // Iteration count and times directory.
////    CommandLineArgs::specify_command_line_flag("--itstimedir", 
////        &Itstime_dir_str);
////
////    CommandLineArgs::specify_command_line_flag("--solver_type",
////        &Solver_type);
////
////    CommandLineArgs::specify_command_line_flag("--time_type",
////        &Time_type);
////
////    CommandLineArgs::specify_command_line_flag("--dt", &Delta_t);
////    CommandLineArgs::specify_command_line_flag("--time_start", &Time_start);
////    CommandLineArgs::specify_command_line_flag("--time_end", &Time_end);
////
////    CommandLineArgs::specify_command_line_flag("--mesh_type", &Mesh_type);
////  }
//
////  inline void generic_problem_setup(const unsigned& dim)
////  {
////    if(CommandLineArgs::command_line_flag_has_been_set("--dt"))
////    {
////      if(!CommandLineArgs::command_line_flag_has_been_set("--time_start"))
////      {
////      std::ostringstream err_msg;
////      err_msg << "Please set --time_start" << std::endl;
////      throw OomphLibError(err_msg.str(),
////          OOMPH_CURRENT_FUNCTION,
////          OOMPH_EXCEPTION_LOCATION);
////      }
////
////      if(!CommandLineArgs::command_line_flag_has_been_set("--time_end"))
////      {
////      std::ostringstream err_msg;
////      err_msg << "Please set --time_end" << std::endl;
////      throw OomphLibError(err_msg.str(),
////          OOMPH_CURRENT_FUNCTION,
////          OOMPH_EXCEPTION_LOCATION);
////      }
////    }
////
////    if(CommandLineArgs::command_line_flag_has_been_set("--mesh_type"))
////    {
////      // Check that the mesh valid
////      if(!((Mesh_type != 0) ||
////          (Mesh_type != 1)   )  )
////      {
////      std::ostringstream err_msg;
////      err_msg << "Unrecognised Mesh_type " << Mesh_type  << "\n" 
////              << "0 for triangle / tetrahedral \n"
////              << "1 for quads / hexahedral" << std::endl;
////      throw OomphLibError(err_msg.str(),
////          OOMPH_CURRENT_FUNCTION,
////          OOMPH_EXCEPTION_LOCATION);
////      }
////    }
////
////    // Do we have to distribute the problem?
////    if(CommandLineArgs::command_line_flag_has_been_set("--dist_prob"))
////    {
////      Distribute_problem = true;
////    }
////    else
////    {
////      Distribute_problem = false;
////    }
////
////    if(!CommandLineArgs::command_line_flag_has_been_set("--prob_id"))
////    {
////      std::ostringstream err_msg;
////      err_msg << "Please set --prob_id." << std::endl;
////      throw OomphLibError(err_msg.str(),
////          OOMPH_CURRENT_FUNCTION,
////          OOMPH_EXCEPTION_LOCATION);
////    }
////
////    // Document the solution? Default is false.
////    Doc_soln = false;
////    if(CommandLineArgs::command_line_flag_has_been_set("--doc_soln"))
////    {
////      // The argument immediately after --doc_soln is put into NSPP::Soln_dir_str.
////      // If this begins with "--", then no solution directory has been provided.
////      std::size_t found = Soln_dir_str.find("--");
////
////      // Check if they have set the solution directory.
////      if(found != std::string::npos)
////      {
////        std::ostringstream err_msg;
////        err_msg << "Please provide the doc_soln directory "
////          << "after the argument --doc_soln.\n" 
////          << "This must not start with \"--\"." << std::endl;
////
////        throw OomphLibError(err_msg.str(),
////            OOMPH_CURRENT_FUNCTION,
////            OOMPH_EXCEPTION_LOCATION);
////      }
////      else
////      {
////        Doc_soln = true;
////      }
////    }
////
////
////    // Set the viscuous term.
////    // Default: 0, Sim
////    if(CommandLineArgs::command_line_flag_has_been_set("--visc"))
////    {
////      if (Vis == 0)
////      {
////        if(dim == 2)
////        {
////          for (unsigned d = 0; d < 2; d++) 
////          {
////            NavierStokesEquations<2>::Gamma[d] = 0.0;
////          }
////        }
////        else
////        {
////          for (unsigned d = 0; d < 3; d++) 
////          {
////            NavierStokesEquations<3>::Gamma[d] = 0.0;
////          }
////        }
////      }
////      else if (Vis == 1)
////      {
////        if(dim == 2)
////        {
////          for (unsigned d = 0; d < 2; d++) 
////          {
////            NavierStokesEquations<2>::Gamma[d] = 1.0;
////          }
////        }
////        else
////        {
////          for (unsigned d = 0; d < 3; d++) 
////          {
////            NavierStokesEquations<3>::Gamma[d] = 1.0;
////          }
////        }
////      }
////      else
////      {
////        std::ostringstream err_msg;
////        err_msg << "Do not recognise viscuous term: " << Vis << ".\n"
////          << "Vis = 0 for simple form\n"
////          << "Vis = 1 for stress divergence form\n"
////          << std::endl;
////        throw OomphLibError(err_msg.str(),
////            OOMPH_CURRENT_FUNCTION,
////            OOMPH_EXCEPTION_LOCATION);
////      }
////    }
////    else
////    {
////      std::ostringstream err_msg;
////      err_msg << "Please set --visc to either 0 or 1.\n"
////        << std::endl;
////      throw OomphLibError(err_msg.str(),
////          OOMPH_CURRENT_FUNCTION,
////          OOMPH_EXCEPTION_LOCATION);
////    }
////
////    // Check if the Reynolds numbers have been set.
////    if(  CommandLineArgs::command_line_flag_has_been_set("--rey_start")
////        &&CommandLineArgs::command_line_flag_has_been_set("--rey_incre")
////        &&CommandLineArgs::command_line_flag_has_been_set("--rey_end")
////        &&CommandLineArgs::command_line_flag_has_been_set("--rey"))
////    {
////      std::ostringstream err_msg;
////      err_msg << "You have set all --rey* argument, please choose carefully!\n"
////        << std::endl;
////      throw OomphLibError(err_msg.str(),
////          OOMPH_CURRENT_FUNCTION,
////          OOMPH_EXCEPTION_LOCATION);
////    }
////    else if(  CommandLineArgs::command_line_flag_has_been_set("--rey_start")
////        &&CommandLineArgs::command_line_flag_has_been_set("--rey_incre")
////        &&CommandLineArgs::command_line_flag_has_been_set("--rey_end"))
////    {
////      oomph_info << "Looping Reynolds: \n"
////        << "Rey_start = " << Rey_start << std::endl; 
////      oomph_info << "Rey_incre = " << Rey_incre << std::endl; 
////      oomph_info << "Rey_end = " << Rey_end << std::endl; 
////    }
////    else if(!CommandLineArgs::command_line_flag_has_been_set("--rey"))
////    {
////      std::ostringstream err_msg;
////      err_msg << "No Reynolds numbers have been set.\n"
////        << "For a single Reynolds number, use --rey.\n"
////        << "For looping through Reynolds numbers, use:\n"
////        << "--rey_start --rey_incre --rey_end.\n"
////        << std::endl;
////      throw OomphLibError(err_msg.str(),
////          OOMPH_CURRENT_FUNCTION,
////          OOMPH_EXCEPTION_LOCATION);
////    }
////
////    // Store the iteration and timing results in a file?
////    // The its and time are always outputted in cout, but maybe we would like
////    // to output it to a file.
////    if(CommandLineArgs::command_line_flag_has_been_set("--itstimedir"))
////    {
////      // The argument immediately after --itstimedir is put into 
////      // NSPP::Itstime_dir_str.
////      // If this begins with "--", then no solution directory has been provided.
////      std::size_t found = Itstime_dir_str.find("--");
////
////      // Check if they have set the solution directory.
////      if(found != std::string::npos)
////      {
////        std::ostringstream err_msg;
////        err_msg << "Please provide the itstimedir directory "
////          << "after the argument --itstimedir.\n" 
////          << "This must not start with \"--\"." << std::endl;
////        throw OomphLibError(err_msg.str(),
////            OOMPH_CURRENT_FUNCTION,
////            OOMPH_EXCEPTION_LOCATION);
////      }
////    }
////
////    ////////////////////////////////////////////
////    //Set up the solver types
////    valid_solver_type_map.insert(
////        std::pair<int,std::string>(Solver_type_DIRECT_SOLVE, 
////                                   "Direct solve"));
////    valid_solver_type_map.insert(
////        std::pair<int,std::string>(Solver_type_OOMPHLIB_GMRES, 
////                                   "OOMPH-LIB's GMRES"));
////    valid_solver_type_map.insert(
////        std::pair<int,std::string>(Solver_type_TRILINOS_GMRES, 
////                                   "Trilinos Aztec00 GMRES"));
////
////    if(CommandLineArgs::command_line_flag_has_been_set("--solver_type"))
////    {
////#ifndef OOMPH_HAS_TRILINOS
////      if(Solver_type == Solver_type_TRILINOS_GMRES)
////      {
////        std::ostringstream err_msg;
////        err_msg << "Have set --solver_type to: " << Solver_type << "\n"
////          << "But OOMPH-LIB does not have trilinos!" << std::endl;
////        throw OomphLibError(err_msg.str(),
////            OOMPH_CURRENT_FUNCTION,
////            OOMPH_EXCEPTION_LOCATION);
////      }
////#endif
////
////
////      int_string_map_it_type solver_type_it;
////
////      // Check that the solver type is valid.
////      solver_type_it = valid_solver_type_map.find(Solver_type);
////
////      if(solver_type_it == valid_solver_type_map.end())
////      {
////        std::ostringstream err_msg;
////        err_msg << "Please provide a valid solver type "
////          << "after the argument --solver_type:\n"
////          << "Acceptable IDs are:\n";
////          // Loop through the solver types
////
////
////for(int_string_map_it_type iterator = valid_solver_type_map.begin(); 
////    iterator != valid_solver_type_map.end(); iterator++) 
////{
////  err_msg << iterator->first << " = " << iterator->second << "\n";
////}
////          err_msg << std::endl;
////
////        throw OomphLibError(err_msg.str(),
////            OOMPH_CURRENT_FUNCTION,
////            OOMPH_EXCEPTION_LOCATION);
////      }
////
////    }
////    else
////    {
////        std::ostringstream err_msg;
////        err_msg << "Please set --solver_type\n"
////          << "Acceptable IDs are:\n";
////          // Loop through the solver types
////
////for(int_string_map_it_type iterator = valid_solver_type_map.begin(); 
////    iterator != valid_solver_type_map.end(); iterator++) 
////{
////  err_msg << iterator->first << " = " << iterator->second << "\n";
////}
////          err_msg << std::endl;
////
////        throw OomphLibError(err_msg.str(),
////            OOMPH_CURRENT_FUNCTION,
////            OOMPH_EXCEPTION_LOCATION);
////    }
////
////    if(!CommandLineArgs::command_line_flag_has_been_set("--max_solver_iter"))
////    {
////      std::ostringstream err_msg;
////      err_msg << "Please set --max_solver_iter." << std::endl;
////      throw OomphLibError(err_msg.str(),
////          OOMPH_CURRENT_FUNCTION,
////          OOMPH_EXCEPTION_LOCATION);
////    }
////
////  } // NSPP::generic_problem_setup()
//
////  // NavierStokesProblemParameters::create_label();
////  inline std::string create_label()
////  {
////    std::string mesh_type_str="";
////    if(CommandLineArgs::command_line_flag_has_been_set("--mesh_type"))
////    {
////      if(Mesh_type == MeshType_TETRAHEDRAL)
////      {
////        mesh_type_str="Tet";
////      }
////      else if(Mesh_type == MeshType_HEXAHEDRAL)
////      {
////        mesh_type_str="Hex";
////      }
////      else
////      {
////      std::ostringstream err_msg;
////      err_msg << "Unrecognised Mesh_type: " << Mesh_type << std::endl;
////      throw OomphLibError(err_msg.str(),
////          OOMPH_CURRENT_FUNCTION,
////          OOMPH_EXCEPTION_LOCATION);
////
////      }
////    }
////
////    // Set the string for the Reynolds number.
////    std::string rey_str = "";
////    if(Rey >= 0)
////    {
////      std::ostringstream strs;
////      strs << "R" << Rey;
////      rey_str = strs.str();
////    }
////    else
////    {
////      std::ostringstream err_msg;
////      err_msg << "Something has gone wrong, the Reynolds number is negative\n"
////        << "Rey = " << Rey << "\n"
////        << "Please set it again using --rey.\n"
////        << std::endl;
////      throw OomphLibError(err_msg.str(),
////          OOMPH_CURRENT_FUNCTION,
////          OOMPH_EXCEPTION_LOCATION);
////    }
////
////    // Set the string for viscous term.
////    std::string vis_str = "";
////    if (Vis == 0)
////    {
////      vis_str = "Sim";
////    }
////    else if (Vis == 1)
////    {
////      vis_str = "Str";
////    } // else - setting viscuous term.
////    else
////    {
////      std::ostringstream err_msg;
////      err_msg << "Do not recognise viscuous term: " << Vis << ".\n"
////        << "Vis = 0 for simple form\n"
////        << "Vis = 1 for stress divergence form\n"
////        << std::endl;
////      throw OomphLibError(err_msg.str(),
////          OOMPH_CURRENT_FUNCTION,
////          OOMPH_EXCEPTION_LOCATION);
////    }
////
////    std::string label = mesh_type_str + vis_str + rey_str;
////
////    return label;
////  } // create_label()
//
////  inline void doc_iter_times(Problem* problem_pt,
////                             DocLinearSolverInfo* doc_linear_solver_info_pt)
////  {
////     unsigned iters = 0;
////     double preconditioner_setup_time = 0.0;
////     double solver_time = 0.0;
////
////     // Get the iteration counts and preconditioner setup time
////#ifdef PARANOID
////     IterativeLinearSolver* iterative_solver_pt
////       = dynamic_cast<IterativeLinearSolver*>
////       (problem_pt->linear_solver_pt());
////     if(iterative_solver_pt == 0)
////     {
////       std::ostringstream error_message;
////       error_message << "Cannot cast the solver pointer." << std::endl;
////
////       throw OomphLibError(error_message.str(),
////           OOMPH_CURRENT_FUNCTION,
////           OOMPH_EXCEPTION_LOCATION);
////     }
////     else
////     {
////       iters = iterative_solver_pt->iterations();
////       preconditioner_setup_time 
////         = iterative_solver_pt->preconditioner_pt()->setup_time();
////     }
////#else
////     iters = static_cast<IterativeLinearSolver*>
////       (problem_pt->linear_solver_pt())->iterations();
////     preconditioner_setup_time = static_cast<IterativeLinearSolver*>
////       (problem_pt->linear_solver_pt())->preconditioner_pt()->setup_time();
////#endif
////
////     // Set the solver time.
////     if(Solver_type == Solver_type_TRILINOS_GMRES)
////     {
////       TrilinosAztecOOSolver* trilinos_solver_pt 
////         = dynamic_cast<TrilinosAztecOOSolver*>(problem_pt->linear_solver_pt());
////       solver_time = trilinos_solver_pt->linear_solver_solution_time();
////     }
////     else
////     {
////       solver_time 
////         = problem_pt->linear_solver_pt()->linear_solver_solution_time();
////     }
////
////     doc_linear_solver_info_pt->add_iteration_and_time
////       (iters,preconditioner_setup_time,solver_time);
////  }
//
//  
//} // NavierStokesProblemParameters

//namespace PreconditionerHelpers
//{
//  //////////////////////////////////////////////////////////////////////
//
//
//  //////////////////////////////////////////////////////////////////////
//
//  // Set my problem constructor.
////  Vector<Mesh*> Mesh_pt;
//
//  // Set by main()
////  DocLinearSolverInfo* Doc_linear_solver_info_pt = 0;
//
//  // Set by main()
////  Problem* Problem_pt = 0;
//
//  // Set by main(), from NSPP
////  int* Vis_pt = 0;
//
////  std::string* Label_str_pt = 0;
//
//  ///////////////////////////////////////////////////
//  // Set directly by CL
////  int W_solver = -1;
////  int NS_solver = -1;
////  int F_solver = -1;
////  int P_solver = -1;
//
////  double F_amg_strength = -1.0;
////  double F_amg_damping = -1.0;
////  int F_amg_coarsening = -1;
////  int F_amg_simple_smoother = -1;
////  int F_amg_complex_smoother = -1;
////  int F_amg_iterations = -1;
////  int F_amg_smoother_iterations = -1;
//
////  double P_amg_strength = -1.0;
////  double P_amg_damping = -1.0;
////  int P_amg_coarsening = -1;
////  int P_amg_simple_smoother = -1;
////  int P_amg_complex_smoother = -1;
////  int P_amg_iterations = -1;
////  int P_amg_smoother_iterations = -1;
//
////  std::string Doc_prec_dir_str = "";
////  double Scaling_sigma = 0.0;
//
//  // Set indirectly by CL
//
//  // Set to true if --doc_prec is provided with a directory.
////  bool Doc_prec = false;
//
//  // Set to true if --print_hypre is set.
////  bool Print_hypre = false;
//
//  // Set to false if --sigma value is set.
////  bool Use_axnorm = true;
//
//  // Set to true if --bdw is set.
////  bool Use_block_diagonal_w = false;
//
//  // Set to true if --lsc_only is set.
////  bool Lsc_only = false;
//
//
////  Preconditioner* Lgr_preconditioner_pt = 0;
////  Preconditioner* NS_preconditioner_pt = 0;
////  Preconditioner* F_preconditioner_pt = 0;
////  Preconditioner* P_preconditioner_pt = 0;
//
//
//  inline void setup_commandline_flags()
//  {
//    // Flag to output the preconditioner, used for debugging.
//    // string
//    CommandLineArgs::specify_command_line_flag(
//        "--doc_prec",&Doc_prec_dir_str);
//
//    // No parameter after.
//    CommandLineArgs::specify_command_line_flag(
//        "--lsc_only");
//
//    // double
//    CommandLineArgs::specify_command_line_flag(
//        "--sigma",&Scaling_sigma);
//
//    // int
//    CommandLineArgs::specify_command_line_flag(
//        "--w_solver",&W_solver);
//
//    // Nothing set
//    CommandLineArgs::specify_command_line_flag("--bdw");
//
//    // int
//    CommandLineArgs::specify_command_line_flag(
//        "--ns_solver",&NS_solver);
//
//    // int
//    CommandLineArgs::specify_command_line_flag(
//        "--p_solver",&P_solver);
//
//    // int
//    CommandLineArgs::specify_command_line_flag(
//        "--f_solver",&F_solver);
//
//    // NS_F block AMG parameters
//    // double
//    CommandLineArgs::specify_command_line_flag(
//        "--f_amg_str",&f_amg_strength);
//    // double
//    CommandLineArgs::specify_command_line_flag(
//        "--f_amg_damp",&f_amg_damping);
//
//    // int
//    CommandLineArgs::specify_command_line_flag("--f_amg_coarse",
//        &f_amg_coarsening);
//
//    // int
//    CommandLineArgs::specify_command_line_flag(
//        "--f_amg_sim_smoo",&f_amg_simple_smoother);
//
//    // int
//    CommandLineArgs::specify_command_line_flag(
//        "--f_amg_com_smoo",&f_amg_complex_smoother);
//
//    // int
//    CommandLineArgs::specify_command_line_flag(
//        "--f_amg_iter",&f_amg_iterations);
//
//    //int
//    CommandLineArgs::specify_command_line_flag("--f_amg_smiter",
//        &f_amg_smoother_iterations);
//
//    // NS_P block AMG parameters
//    // double
//    CommandLineArgs::specify_command_line_flag(
//        "--p_amg_str",&p_amg_strength);
//    // double
//    CommandLineArgs::specify_command_line_flag(
//        "--p_amg_damp",&p_amg_damping);
//    // int
//    CommandLineArgs::specify_command_line_flag(
//        "--p_amg_coarse",&p_amg_coarsening);
//    // int
//    CommandLineArgs::specify_command_line_flag(
//        "--p_amg_sim_smoo",&p_amg_simple_smoother);
//    // int
//    CommandLineArgs::specify_command_line_flag(
//        "--p_amg_com_smoo",&p_amg_complex_smoother);
//
//    // int
//    CommandLineArgs::specify_command_line_flag(
//        "--p_amg_iter",&p_amg_iterations);
//    // int
//    CommandLineArgs::specify_command_line_flag("--p_amg_smiter",
//        &p_amg_smoother_iterations);
//
//    CommandLineArgs::specify_command_line_flag("--print_hypre");
//  }
//
//  inline void generic_setup()
//  {
//    if(CommandLineArgs::command_line_flag_has_been_set("--lsc_only"))
//    {
//      Lsc_only = true;
//    }
//    else
//    {
//      Lsc_only = false;
//    }
//
//    // Document the preconditioner? Default is false.
//    if(CommandLineArgs::command_line_flag_has_been_set("--doc_prec"))
//    {
//      // The argument immediately after --doc_prec is put into SL::Doc_prec_dir.
//      // If this begins with "--", then no prec directory has been provided.
//      std::size_t found = Doc_prec_dir_str.find("--");
//
//      // Check if they have set the doc_prec directory.
//      if(found != std::string::npos)
//      {
//        std::ostringstream err_msg;
//        err_msg 
//          << "Please provide the doc_prec directory "
//          << "after the argument --doc_prec.\n" 
//          << "This must not start with \"--\"." << std::endl;
//
//        throw OomphLibError(err_msg.str(),
//            OOMPH_CURRENT_FUNCTION,
//            OOMPH_EXCEPTION_LOCATION);
//      }
//      else
//      {
//        Doc_prec = true;
//      }
//    }
//    else
//    {
//      Doc_prec = false;
//    }
//
//
//    // If we are using SuperLU for the Navier-Stokes block (NS_solver = 0), 
//    // the F_solver and P_solver should not be set.
//    if((NS_solver == 0) && 
//        (CommandLineArgs::command_line_flag_has_been_set("--p_solver") ||
//         CommandLineArgs::command_line_flag_has_been_set("--f_solver")))
//    {
//      std::ostringstream err_msg;
//      err_msg << "NS_solver = 0, using SuperLU for the Navier-Stokes block.\n"
//        << "But you have either --f_solver or --p_solver.\n"
//        << "These should NOT be set, unless you want to use LSC.\n"
//        << "In which case, set --ns_solver 1.\n"
//        << std::endl;
//      throw OomphLibError(err_msg.str(),
//          OOMPH_CURRENT_FUNCTION,
//          OOMPH_EXCEPTION_LOCATION);
//    }
//
//    if(CommandLineArgs::command_line_flag_has_been_set("--sigma"))
//    {
//      Use_axnorm = false;
//    }
//    else
//    {
//      Use_axnorm = true;
//    }
//
//    if(CommandLineArgs::command_line_flag_has_been_set("--bdw"))
//    {
//      Use_block_diagonal_w = true;
//    }
//    else
//    {
//      Use_block_diagonal_w = false;
//    }
//
//    if(CommandLineArgs::command_line_flag_has_been_set("--print_hypre"))
//    {
//      Print_hypre = true;
//    }
//    else
//    {
//      Print_hypre = false;
//    }
//  } // LPH::generic_setup()
//
//
//  inline Preconditioner* get_lsc_preconditioner()
//  {
//    if(Lsc_only && (Mesh_pt.size() != 1))
//    {
//      std::ostringstream err_msg;
//      err_msg << "You have chosen to use only the LSC preconditioner\n"
//        << "Thus the Vector Mesh_pt must contain exactly one mesh,\n"
//        << "the Navier Stokes bulk mesh.\n"
//        << std::endl;
//      throw OomphLibError(err_msg.str(),
//          OOMPH_CURRENT_FUNCTION,
//          OOMPH_EXCEPTION_LOCATION);
//    }
//
//
//    // If ns_solver is 1, this means we want to use LSC.
//    // So the F_solver and P_solver must be set.
//    if((F_solver == -1) || (P_solver == -1))
//    {
//      std::ostringstream err_msg;
//      err_msg << "Getting LSC preconditioner for NS block.\n"
//        << "But --f_solver and --p_solver have not been set.\n"
//        << "0 - Exact (SuperLU)\n"
//        << "xy - please check the code for more details.\n"
//        << std::endl;
//      throw OomphLibError(err_msg.str(),
//          OOMPH_CURRENT_FUNCTION,
//          OOMPH_EXCEPTION_LOCATION);
//    }
//
//    // Check that the problem pointer is set (not null).
//    // LSC requires a problem pointer.
//    if(Problem_pt == 0)
//    {
//      std::ostringstream err_msg;
//      err_msg << "Please set the Problem_pt variable.\n"
//        << std::endl;
//      throw OomphLibError(err_msg.str(),
//          OOMPH_CURRENT_FUNCTION,
//          OOMPH_EXCEPTION_LOCATION);
//    }
//
//    // Create the NS LSC preconditioner.
//    NavierStokesSchurComplementPreconditioner* ns_preconditioner_pt =
//      new NavierStokesSchurComplementPreconditioner(Problem_pt);
//
//
//    // Give LSC the bulk mesh (Navier-Stokes mesh).
//    ns_preconditioner_pt->set_navier_stokes_mesh(Mesh_pt[0]);
//
//
//
//    //// Setting the F solver within the NS block
//    /////////////////////////////////////////////
//
//    // Preconditioner for the F block:
//    Preconditioner* f_preconditioner_pt = 0;
//
//    // f_solver == 0 is default, so do nothing.
//    //
//    // AMG depends on the Reynolds number so we check that the Reynolds number
//    // is set if LSC is switched on - even if we do not want to use AMG...
//    // for consistency.
//    if(Vis_pt == 0)
//    {
//      std::ostringstream err_msg;
//      err_msg << "Please set your Vis_pt variable.\n"
//        << "E.g. in main method do: LPH::Vis_pt = &NSPP::Vis;"
//        << std::endl;
//      throw OomphLibError(err_msg.str(),
//          OOMPH_CURRENT_FUNCTION,
//          OOMPH_EXCEPTION_LOCATION);
//    }
//
//    if((*Vis_pt) == -1)
//    {
//      std::ostringstream err_msg;
//      err_msg << "Please set your viscuous term in NSPP header.\n"
//        << std::endl;
//      throw OomphLibError(err_msg.str(),
//          OOMPH_CURRENT_FUNCTION,
//          OOMPH_EXCEPTION_LOCATION);
//    }
//
//
//    ///////////////////////////////////////// FFFFFFFFFFFFFFFFFF
//
//
//    if(F_solver == 11)
//    {
//#ifdef OOMPH_HAS_HYPRE
//      // LSC takes type "Preconditioner".
//      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//        set_hypre_for_2D_poison_problem();
//#endif
//    }
//    else if(F_solver == 12)
//    {
//#ifdef OOMPH_HAS_HYPRE
//      // LSC takes type "Preconditioner".
//      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//        set_hypre_for_navier_stokes_momentum_block();
//#endif
//    }
//    else if(F_solver == 13)
//    {
//#ifdef OOMPH_HAS_HYPRE
//      // LSC takes type "Preconditioner".
//      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//        set_hypre_for_CLJPGSStrn075();
//#endif
//    }
//    else if(F_solver == 14)
//    {
//#ifdef OOMPH_HAS_HYPRE
//      // LSC takes type "Preconditioner".
//      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//        set_hypre_for_RSGSStrn075();
//#endif
//    }
//    else if(F_solver == 15)
//    {
//#ifdef OOMPH_HAS_HYPRE
//      // LSC takes type "Preconditioner".
//      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//        set_hypre_for_CLJPPilutStrn075();
//#endif
//    }
//    else if(F_solver == 16)
//    {
//#ifdef OOMPH_HAS_HYPRE
//      // LSC takes type "Preconditioner".
//      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//        set_hypre_for_RSPilutStrn075();
//#endif
//    }
//    else if(F_solver == 17)
//    {
//#ifdef OOMPH_HAS_HYPRE
//      // LSC takes type "Preconditioner".
//      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//        set_hypre_for_augmented_momentum_block();
//#endif
//    }
//    else if(F_solver == 81)
//    {
//#ifdef OOMPH_HAS_HYPRE
//      // LSC takes type "Preconditioner".
//      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//        set_hypre_for_CLJPGSStrn0668();
//#endif
//    }
//    else if(F_solver == 82)
//    {
//#ifdef OOMPH_HAS_HYPRE
//      // LSC takes type "Preconditioner".
//      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//        set_hypre_for_CLJPJStrn0668();
//#endif
//    }
//    else if(F_solver == 83)
//    {
//#ifdef OOMPH_HAS_HYPRE
//      // LSC takes type "Preconditioner".
//      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//        set_hypre_for_CLJPPilutStrn0668();
//#endif
//    }
//    else if(F_solver == 84)
//    {
//#ifdef OOMPH_HAS_HYPRE
//      // LSC takes type "Preconditioner".
//      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//        set_hypre_for_RSGSStrn0668();
//#endif
//    }
//    else if(F_solver == 85)
//    {
//#ifdef OOMPH_HAS_HYPRE
//      // LSC takes type "Preconditioner".
//      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//        set_hypre_for_RSJStrn0668();
//#endif
//    }
//    else if(F_solver == 86)
//    {
//#ifdef OOMPH_HAS_HYPRE
//      // LSC takes type "Preconditioner".
//      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//        set_hypre_for_RSPilutStrn0668();
//#endif
//    }
//    else if(F_solver == 2)
//    {
//      //f_preconditioner_pt = new RayBlockDiagonalPreconditioner<CRDoubleMatrix>;
//      f_preconditioner_pt = new BlockDiagonalPreconditioner<CRDoubleMatrix>;
//    }
//    else if(F_solver == 3)
//    {
//      f_preconditioner_pt = new BlockDiagonalPreconditioner<CRDoubleMatrix>;
//#ifdef OOMPH_HAS_HYPRE
//      dynamic_cast<BlockDiagonalPreconditioner<CRDoubleMatrix>* >
//        (f_preconditioner_pt)->set_subsidiary_preconditioner_function
//        (Hypre_Subsidiary_Preconditioner_Helper::set_hypre_for_2D_poison_problem);
//#endif
//    }
//    else if (F_solver == 69)
//    {
//#ifdef OOMPH_HAS_HYPRE
//      // This is what set_defaults_for_2D_poisson_problem() does:
//      f_amg_iterations = 1;
//      //f_amg_simple_smoother = 1; // GS - commented out since it is reset below
//      //f_amg_strength = 0.25; // commented out since it is reset below, depending on the viscous term.
//      //f_amg_coarsening = 0; // CLJP - commented out since it is reset below.
//      // END OF 2D poisson stuff.
//
//      // Set this to -1 to make sure complex isn't used
//      f_amg_complex_smoother = -1;
//
//      // The default smoother iterations is two.
//      f_amg_smoother_iterations = 2;
//
//      // Now set my own stuff.
//      f_amg_simple_smoother = 1; // GS
//      f_amg_coarsening = 1; // RSs.
//
//      // There is no damping with GS, otherwise we set the parameter:
//      f_amg_damping = -1.0;
//
//      // Different amg strength for simple/stress divergence for viscous term.
//      // This is only set to the below defaults if the strength parameter is
//      // not already set.
//      if(f_amg_strength < 0)
//      {
//        const int vis = *Vis_pt;
//        if(vis == 0)
//        {
//          // Simple form
//          f_amg_strength = 0.25;
//        }
//        else if (vis == 1)
//        {
//          // Stress divergence form
//          f_amg_strength = 0.668;
//        }
//        else
//        {
//          std::ostringstream err_msg;
//          err_msg << "Do not recognise viscuous term: " << vis << std::endl;
//
//          throw OomphLibError(err_msg.str(),
//              OOMPH_CURRENT_FUNCTION,
//              OOMPH_EXCEPTION_LOCATION);
//        }
//      }
//
//      // Setup the preconditioner.
//      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//        get_custom_hypre_preconditioner(
//            f_amg_iterations, f_amg_smoother_iterations, 
//            f_amg_simple_smoother, f_amg_complex_smoother,
//            f_amg_damping, f_amg_strength,
//            f_amg_coarsening);
//
//      if(Print_hypre)
//      {
//        Hypre_Subsidiary_Preconditioner_Helper::print_hypre_settings(
//            f_preconditioner_pt);
//      }
//#endif
//    }
//    else if (F_solver == 96)
//    {
//#ifdef OOMPH_HAS_HYPRE
//      // In this method we pass everything to the global AMG parameters and
//      // let that handle the work.
//
//      // AMG coarsening:
//      // Set: RayGlobalAMGParam::amg_coarsening = 
//      //
//      // 0 - CLJP
//      // 1 - Classical RS
//      // 3 - modified RS
//      // 6 - Falgout - default on documentation
//      // 8 - PMIS
//      // 10 - HMIS
//      // 11 - One pass on RS coarsening on each processor, not recommended.
//
//      // AMG smoother:
//      // Set: RayGlobalAMGParam::amg_smoother = 
//      // 0 - Jacobi (Need to set damping as well)
//      // 1 - Gauss-Seidel, sequential, very slow in parallel
//      // 2 - GS - interior parallel, serial on boundary.
//      // 3 - hybrid GS or SOR, forward solve
//      // 4 - hybrid GS or SOR, backwards solve.
//      // 6 - hybrid symmetric GS or SSOR.
//      // REMEMBER TO SET SIMPLE SMOOTHING TO TRUE. This should be handled
//      // by your hypre preconditioner creation function.
//      //
//      // Complex smoothing:
//      // 6 - Schwarz
//      // 7 - Pilut
//      // 8 - ParaSails
//      // 9 - Euclid
//      // REMEMBER TO SET SIMPLE SMOOTHING TO FALSE, this should be handled
//      // by the hypre subsidiary helper function.
//
//      // First check that the user has selected a coarsening strategy
//      // using --f_amg_coarse.
//      // If the user has, then f_amg_coarsening should no longer be -1.
//      bool coarsening_ok 
//        = RayGen::check_amg_coarsener(f_amg_coarsening);
//      if((f_amg_coarsening < 0) || (!coarsening_ok))
//      {
//        std::ostringstream err_msg;
//        err_msg 
//          << "Please set a coarsening strategy with --f_amg_coarse.\n"
//          << "Current coarsening strategy (f_amg_coarsening) is: " 
//          << f_amg_coarsening << "\n\n"
//
//          << "You have either not set it or it is not valid.\n\n"
//
//          << "Valid IDs are:\n"
//          << "0 - CLJP\n"
//          << "1 - Classical RS\n"
//          << "3 - modified RS\n"
//          << "6 - Falgout (default in documentation)\n"
//          << "8 - PMIS\n"
//          << "10 - HMIS\n"
//          << "11 - One pass on RS coarsening on each processor, "
//          << "not recommended.\n"
//          << std::endl;
//
//        throw OomphLibError(err_msg.str(),
//            OOMPH_CURRENT_FUNCTION,
//            OOMPH_EXCEPTION_LOCATION);
//      }
//
//      // First check if the user has selected both simple and complex
//      // smoothing, this is not allowed, either one or the other.
//      if((f_amg_simple_smoother >= 0) && 
//          (f_amg_complex_smoother >= 0))
//      {
//        std::ostringstream err_msg;
//        err_msg 
//          << "Both simple and complex smoothers are set.\n"
//          << "Please choose one or the other.\n\n" 
//          << "f_amg_simple_smoother is " << f_amg_simple_smoother << "\n"
//          << "f_amg_complex_smoother is " << f_amg_complex_smoother << "\n\n"
//
//          << "Simple smoother IDs, set with --f_amg_sim_smoo\n"
//          << "0 - Jacobi (Need to set damping as well)\n"
//          << "1 - Gauss-Seidel, sequential, very slow in parallel\n"
//          << "2 - GS - interior parallel, serial on boundary\n"
//          << "3 - hybrid GS or SOR, forward solve\n"
//          << "4 - hybrid GS or SOR, backwards solve\n"
//          << "6 - hybrid symmetric GS or SSOR.\n\n"
//
//          << "Complex smoother IDs, set with --f_amg_com_smoo\n"
//          << "6 - Schwarz (default in documentation)\n"
//          << "7 - Pilut\n"
//          << "8 - ParaSails\n"
//          << "9 - Euclid\n" 
//          << std::endl;
//
//        throw OomphLibError(err_msg.str(),
//            OOMPH_CURRENT_FUNCTION,
//            OOMPH_EXCEPTION_LOCATION);
//      }
//      else if((f_amg_simple_smoother < 0) && (f_amg_complex_smoother < 0))
//      {
//        std::ostringstream err_msg;
//        err_msg 
//          << "Please select a smoother for the f block.\n"
//          << "Use --f_amg_sim_smoo or --f_amg_com_smoo flag.\n\n"
//
//          << "Simple smoother IDs, set with --f_amg_sim_smoo\n"
//          << "0 - Jacobi (Need to set damping as well)\n"
//          << "1 - Gauss-Seidel, sequential, very slow in parallel\n"
//          << "2 - GS - interior parallel, serial on boundary\n"
//          << "3 - hybrid GS or SOR, forward solve\n"
//          << "4 - hybrid GS or SOR, backwards solve\n"
//          << "6 - hybrid symmetric GS or SSOR.\n\n"
//
//          << "Complex smoother IDs, set with --f_amg_com_smoo\n"
//          << "6 - Schwarz (default in documentation)\n"
//          << "7 - Pilut\n"
//          << "8 - ParaSails\n"
//          << "9 - Euclid\n" 
//          << std::endl;
//
//        throw OomphLibError(err_msg.str(),
//            OOMPH_CURRENT_FUNCTION,
//            OOMPH_EXCEPTION_LOCATION);
//      }
//
//
//      // Only one of the two smoothers have been set. We see if these are
//      // valid smoothing IDs.
//      if(f_amg_simple_smoother >= 0)
//      {
//        // check if simple smoother is okay.
//        bool sim_smoother_ok
//          = RayGen::check_amg_sim_smoother(f_amg_simple_smoother);
//
//        if(!sim_smoother_ok)
//        {
//          std::ostringstream err_msg;
//          err_msg 
//            << "Please provide a valid simple smoother \n"
//            << "using --f_amg_sim_smoo. You have provided: " 
//            << f_amg_simple_smoother <<"\n\n"
//            << "Valid IDs are:\n"
//            << "0 - Jacobi (Need to set damping as well)\n"
//            << "1 - Gauss-Seidel, sequential, very slow in parallel\n"
//            << "2 - GS - interior parallel, serial on boundary\n"
//            << "3 - hybrid GS or SOR, forward solve\n"
//            << "4 - hybrid GS or SOR, backwards solve\n"
//            << "6 - hybrid symmetric GS or SSOR.\n"
//            << std::endl;
//
//          throw OomphLibError(err_msg.str(),
//              OOMPH_CURRENT_FUNCTION,
//              OOMPH_EXCEPTION_LOCATION);
//        }
//      }
//      else if(f_amg_complex_smoother >=0)
//      {
//        // check if complex smoother is valid.
//        bool com_smoother_ok
//          = RayGen::check_amg_com_smoother(f_amg_complex_smoother);
//        if(!com_smoother_ok)
//        {
//          std::ostringstream err_msg;
//          err_msg 
//            << "Please provide a valid complex smoother \n"
//            << "using --f_amg_com_smoo. You have provided: " 
//            << f_amg_complex_smoother <<"\n\n"
//            << "Valid IDs are:\n"
//            << "6 - Schwarz (default in documentation)\n"
//            << "7 - Pilut\n"
//            << "8 - ParaSails\n"
//            << "9 - Euclid\n"
//            << std::endl;
//
//          throw OomphLibError(err_msg.str(),
//              OOMPH_CURRENT_FUNCTION,
//              OOMPH_EXCEPTION_LOCATION);
//        }
//      }
//      else
//      {
//        // Should never get here, something has gone wrong.
//        std::ostringstream err_msg;
//        err_msg 
//          << "Something when wrong with your selection of smoother.\n"
//          << "Choose to use either a simple smoother or complex smoother\n"
//          << "with the flags --f_amg_sim_smoo and --f_amg_com_smoo.\n\n"
//
//          << "Current f_amg_simple_smoother is " 
//          << f_amg_simple_smoother << "\n"
//          << "Current f_amg_complex_smoother is " 
//          << f_amg_complex_smoother << "\n\n"
//
//          << "Simple smoother IDs, set with --f_amg_sim_smoo\n"
//          << "0 - Jacobi (Need to set damping as well)\n"
//          << "1 - Gauss-Seidel, sequential, very slow in parallel\n"
//          << "2 - GS - interior parallel, serial on boundary\n"
//          << "3 - hybrid GS or SOR, forward solve\n"
//          << "4 - hybrid GS or SOR, backwards solve\n"
//          << "6 - hybrid symmetric GS or SSOR.\n\n"
//
//          << "Complex smoother IDs, set with --f_amg_com_smoo\n"
//          << "6 - Schwarz (default in documentation)\n"
//          << "7 - Pilut\n"
//          << "8 - ParaSails\n"
//          << "9 - Euclid\n" 
//          << std::endl;
//
//        throw OomphLibError(err_msg.str(),
//            OOMPH_CURRENT_FUNCTION,
//            OOMPH_EXCEPTION_LOCATION);
//      }
//
//
//      // Damping is required for Jacobi or hybrid SOR smoothers.
//      // First check if amg_simple_smoother is one of those.
//      bool damping_required
//        = RayGen::check_amg_sim_smoother_damping_required(f_amg_simple_smoother);
//
//      if(damping_required && (f_amg_damping < 0))
//      {
//        std::ostringstream err_msg;
//        err_msg 
//          << "Have selected simple smoother: " << f_amg_simple_smoother << "\n"
//          << "Damping parameter is required for this smoother.\n"
//          << "Please set it via --f_amg_damp.\n"
//          << std::endl;
//
//        throw OomphLibError(err_msg.str(),
//            OOMPH_CURRENT_FUNCTION,
//            OOMPH_EXCEPTION_LOCATION);
//      }
//      // Note that we haven't checked the reverse, i.e. if a smoother
//      // which does not required damping is set, and the damping parameter
//      // is set, we allow this... if implemented properly, the hypre
//      // code should just ignore it. It will also be interesting to find
//      // smoothers which the damping parameter affects the result.
//
//
//      // Now check if the strength parameter is set. 
//      // But don't throw an error, just a warning.
//      if(f_amg_strength < 0)
//      {
//        std::ostringstream warning_msg;
//        warning_msg 
//          << "You have not set the strength parameter for --f_amg_str\n"
//          << "The default settings will be used.\n"
//          << "0.25 for simple form of the viscous term or \n"
//          << "0.668 for stress divergence form of the viscous term.\n"
//          << std::endl;
//
//        throw OomphLibWarning(warning_msg.str(),
//            OOMPH_CURRENT_FUNCTION,
//            OOMPH_EXCEPTION_LOCATION);
//
//        const int vis = (*Vis_pt);
//        if(vis == 0)
//        {
//          // Simple form
//          f_amg_strength = 0.25;
//        }
//        else if(vis == 1)
//        {
//          // Stress divergence form
//          f_amg_strength = 0.668;
//        }
//        else
//        {
//          std::ostringstream err_msg;
//          err_msg 
//            << "Do not recognise viscuous term: " << vis << "\n"
//            << "Please point it to NSPP::Vis\n" 
//            << std::endl;
//
//          throw OomphLibError(err_msg.str(),
//              OOMPH_CURRENT_FUNCTION,
//              OOMPH_EXCEPTION_LOCATION);
//        }
//      }
//
//      // Now check the f_amg_iterations
//      if(f_amg_iterations < 0)
//      {
//        std::ostringstream err_msg;
//        err_msg 
//          << "Have not set f_amg_iterations via --f_amg_iter\n"
//          << std::endl;
//
//        throw OomphLibError(err_msg.str(),
//            OOMPH_CURRENT_FUNCTION,
//            OOMPH_EXCEPTION_LOCATION);
//      }
//
//      if(f_amg_smoother_iterations < 0)
//      {
//        std::ostringstream err_msg;
//        err_msg 
//          << "Have not set f_amg_smoother_iterations via --f_amg_smiter\n"
//          << std::endl;
//
//        throw OomphLibError(err_msg.str(),
//            OOMPH_CURRENT_FUNCTION,
//            OOMPH_EXCEPTION_LOCATION);
//      }
//
//      // Setup the preconditioner.
//      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//        get_custom_hypre_preconditioner(
//            f_amg_iterations, f_amg_smoother_iterations, 
//            f_amg_simple_smoother, f_amg_complex_smoother,
//            f_amg_damping, f_amg_strength,
//            f_amg_coarsening);
//
//      if(Print_hypre)
//      {
//        Hypre_Subsidiary_Preconditioner_Helper::print_hypre_settings(
//            f_preconditioner_pt);
//      }
//#endif
//    }
//    ////////////////////////////////////////////////////////////////////////
//    // For the below, we have 8090, 8091, 8092 and 8093 with the following 
//    // configuration:
//    // For the LSC F block:
//    // 8090 - block diagonal with SuperLU
//    // 8091 - upper block triangular with Super LU
//    // 8092 - lower block triangular with Super LU
//    // 8093 - Full SuperLU f block.
//    else if(F_solver == 8090)
//    {
//#ifdef OOMPH_HAS_HYPRE
//      // Create a block diagonal preconditioner.
//      f_preconditioner_pt =
//        new BlockDiagonalPreconditioner<CRDoubleMatrix>;
//#endif
//    }
//    else if(F_solver == 8091)
//    {
//#ifdef OOMPH_HAS_HYPRE
//      // Create a block triangular preconditioner.
//      f_preconditioner_pt =
//        new BlockTriangularPreconditioner<CRDoubleMatrix>;
//
//      // Use upper triangular preconditioner.
//      dynamic_cast<BlockTriangularPreconditioner<CRDoubleMatrix>* >
//      (f_preconditioner_pt)->upper_triangular();
//#endif
//    }
//    else if(F_solver == 8092)
//    {
//#ifdef OOMPH_HAS_HYPRE
//      // Create a triangular preconditioner.
//      f_preconditioner_pt =
//        new BlockTriangularPreconditioner<CRDoubleMatrix>;
//
//      // Use lower triangular preconditioner.
//      dynamic_cast<BlockTriangularPreconditioner<CRDoubleMatrix>* >
//      (f_preconditioner_pt)->lower_triangular();
//#endif
//    }
//    else if(F_solver == 8093)
//    {
//#ifdef OOMPH_HAS_HYPRE
//    // This is using super LU for the full F block.
//    // Since SuperLU is the default behaviour, we simply set the pointer 
//    // to 0.
//    f_preconditioner_pt = 0;
//#endif
//    }
//    ////////////////////////////////////////////////////////////////////////
//    // For the below, we have 9090, 9091, 9092 and 9093 with the following 
//    // configuration:
//    // For the LSC F block:
//    // 9090 - block diagonal with Hypre
//    // 9091 - upper block triangular with Hypre
//    // 9092 - lower block triangular with Hypre
//    // 9093 - full AMG.
//    else if(F_solver == 9090)
//    {
//#ifdef OOMPH_HAS_HYPRE
//      // Create a block diagonal preconditioner.
//      f_preconditioner_pt =
//        new BlockDiagonalPreconditioner<CRDoubleMatrix>;
//
//      // Now, since f_precondtioner_pt is a Preconditioner*, it needs to be
//      // caste to a block diagonal one if we want to call functions from
//      // that class.
//      dynamic_cast<BlockDiagonalPreconditioner<CRDoubleMatrix>* >
//      (f_preconditioner_pt)->set_subsidiary_preconditioner_function
//      (Hypre_Subsidiary_Preconditioner_Helper::set_hypre_JhalfStrnSimOneVTwoTwoRS);
//
//      // All done. The below is prints the Hypre settings.
//
//      // Check the Hypre values used, we encapsulate this so we can easily 
//      // take it out later.
//      {
//        // Create a new preconditioner with the above function we set.
//        Preconditioner* check_prec_pt = 
//          Hypre_Subsidiary_Preconditioner_Helper::
//          set_hypre_JhalfStrnSimOneVTwoTwoRS();
//
//        // Now print it out to see the settings!
//        Hypre_Subsidiary_Preconditioner_Helper::
//          print_hypre_settings(check_prec_pt);
//      }
//#endif
//    }
//    else if(F_solver == 9091)
//    {
//#ifdef OOMPH_HAS_HYPRE
//      // Create a block triangular preconditioner.
//      f_preconditioner_pt =
//        new BlockTriangularPreconditioner<CRDoubleMatrix>;
//
//      // Use upper triangular preconditioner.
//      dynamic_cast<BlockTriangularPreconditioner<CRDoubleMatrix>* >
//      (f_preconditioner_pt)->upper_triangular();
//
//      // Set the Hypre preconditioner.
//      dynamic_cast<BlockTriangularPreconditioner<CRDoubleMatrix>* >
//      (f_preconditioner_pt)->set_subsidiary_preconditioner_function
//      (Hypre_Subsidiary_Preconditioner_Helper::set_hypre_JhalfStrnSimOneVTwoTwoRS);
//
//      // All done. The below is prints the Hypre settings.
//
//      // Check the Hypre values used, we encapsulate this so we can easily 
//      // take it out later.
//      {
//        // Create a new preconditioner with the above function we set.
//        Preconditioner* check_prec_pt = 
//          Hypre_Subsidiary_Preconditioner_Helper::
//          set_hypre_JhalfStrnSimOneVTwoTwoRS();
//
//        // Now print it out to see the settings!
//        Hypre_Subsidiary_Preconditioner_Helper::
//          print_hypre_settings(check_prec_pt);
//      }
//
//#endif
//    }
//    else if(F_solver == 9092)
//    {
//#ifdef OOMPH_HAS_HYPRE
//      // Create a triangular preconditioner.
//      f_preconditioner_pt =
//        new BlockTriangularPreconditioner<CRDoubleMatrix>;
//
//      // Use lower triangular preconditioner.
//      dynamic_cast<BlockTriangularPreconditioner<CRDoubleMatrix>* >
//      (f_preconditioner_pt)->lower_triangular();
//
//      // Set the hypre preconditioner.
//      dynamic_cast<BlockTriangularPreconditioner<CRDoubleMatrix>* >
//      (f_preconditioner_pt)->set_subsidiary_preconditioner_function
//      (Hypre_Subsidiary_Preconditioner_Helper::set_hypre_JhalfStrnSimOneVTwoTwoRS);
//
//      // All done. The below is prints the Hypre settings.
//
//      // Check the Hypre values used, we encapsulate this so we can easily 
//      // take it out later.
//      {
//        // Create a new preconditioner with the above function we set.
//        Preconditioner* check_prec_pt = 
//          Hypre_Subsidiary_Preconditioner_Helper::
//          set_hypre_JhalfStrnSimOneVTwoTwoRS();
//
//        // Now print it out to see the settings!
//        Hypre_Subsidiary_Preconditioner_Helper::
//          print_hypre_settings(check_prec_pt);
//      }
//#endif
//    }
//    else if(F_solver == 9093)
//    {
//#ifdef OOMPH_HAS_HYPRE
//    // Create a new hypre preconditioner
//    f_preconditioner_pt =
//      Hypre_Subsidiary_Preconditioner_Helper::
//      set_hypre_JhalfStrnSimOneVTwoTwoRS();
//    
//    // Print it to check the settings.
//    Hypre_Subsidiary_Preconditioner_Helper::
//      print_hypre_settings(f_preconditioner_pt);
//#endif
//    }
//    ////////////////////////////////////////////////////////////////////////
//    // For the below, we have 9090, 9091, 9092 and 9093 with the following 
//    // configuration:
//    // For the LSC F block:
//    // 9090 - block diagonal with Hypre
//    // 9091 - upper block triangular with Hypre
//    // 9092 - lower block triangular with Hypre
//    // 9093 - full AMG.
//    else if(F_solver == 9190)
//    {
//#ifdef OOMPH_HAS_HYPRE
//      // Create a block diagonal preconditioner.
//      f_preconditioner_pt =
//        new BlockDiagonalPreconditioner<CRDoubleMatrix>;
//
//      // Now, since f_precondtioner_pt is a Preconditioner*, it needs to be
//      // caste to a block diagonal one if we want to call functions from
//      // that class.
//      dynamic_cast<BlockDiagonalPreconditioner<CRDoubleMatrix>* >
//      (f_preconditioner_pt)->set_subsidiary_preconditioner_function
//      (Hypre_Subsidiary_Preconditioner_Helper::set_hypre_JhalfStrnStrOneVTwoTwoRS);
//
//      // All done. The below is prints the Hypre settings.
//
//      // Check the Hypre values used, we encapsulate this so we can easily 
//      // take it out later.
//      {
//        // Create a new preconditioner with the above function we set.
//        Preconditioner* check_prec_pt = 
//          Hypre_Subsidiary_Preconditioner_Helper::
//          set_hypre_JhalfStrnStrOneVTwoTwoRS();
//
//        // Now print it out to see the settings!
//        Hypre_Subsidiary_Preconditioner_Helper::
//          print_hypre_settings(check_prec_pt);
//      }
//#endif
//    }
//    else if(F_solver == 9191)
//    {
//#ifdef OOMPH_HAS_HYPRE
//      // Create a block triangular preconditioner.
//      f_preconditioner_pt =
//        new BlockTriangularPreconditioner<CRDoubleMatrix>;
//
//      // Use upper triangular preconditioner.
//      dynamic_cast<BlockTriangularPreconditioner<CRDoubleMatrix>* >
//      (f_preconditioner_pt)->upper_triangular();
//
//      // Set the Hypre preconditioner.
//      dynamic_cast<BlockTriangularPreconditioner<CRDoubleMatrix>* >
//      (f_preconditioner_pt)->set_subsidiary_preconditioner_function
//      (Hypre_Subsidiary_Preconditioner_Helper::set_hypre_JhalfStrnStrOneVTwoTwoRS);
//
//      // All done. The below is prints the Hypre settings.
//
//      // Check the Hypre values used, we encapsulate this so we can easily 
//      // take it out later.
//      {
//        // Create a new preconditioner with the above function we set.
//        Preconditioner* check_prec_pt = 
//          Hypre_Subsidiary_Preconditioner_Helper::
//          set_hypre_JhalfStrnStrOneVTwoTwoRS();
//
//        // Now print it out to see the settings!
//        Hypre_Subsidiary_Preconditioner_Helper::
//          print_hypre_settings(check_prec_pt);
//      }
//
//#endif
//    }
//    else if(F_solver == 9192)
//    {
//#ifdef OOMPH_HAS_HYPRE
//      // Create a triangular preconditioner.
//      f_preconditioner_pt =
//        new BlockTriangularPreconditioner<CRDoubleMatrix>;
//
//      // Use lower triangular preconditioner.
//      dynamic_cast<BlockTriangularPreconditioner<CRDoubleMatrix>* >
//      (f_preconditioner_pt)->lower_triangular();
//
//      // Set the hypre preconditioner.
//      dynamic_cast<BlockTriangularPreconditioner<CRDoubleMatrix>* >
//      (f_preconditioner_pt)->set_subsidiary_preconditioner_function
//      (Hypre_Subsidiary_Preconditioner_Helper::set_hypre_JhalfStrnStrOneVTwoTwoRS);
//
//      // All done. The below is prints the Hypre settings.
//
//      // Check the Hypre values used, we encapsulate this so we can easily 
//      // take it out later.
//      {
//        // Create a new preconditioner with the above function we set.
//        Preconditioner* check_prec_pt = 
//          Hypre_Subsidiary_Preconditioner_Helper::
//          set_hypre_JhalfStrnStrOneVTwoTwoRS();
//
//        // Now print it out to see the settings!
//        Hypre_Subsidiary_Preconditioner_Helper::
//          print_hypre_settings(check_prec_pt);
//      }
//#endif
//    }
//    else if(F_solver == 9193)
//    {
//#ifdef OOMPH_HAS_HYPRE
//    // Create a new hypre preconditioner
//    f_preconditioner_pt =
//      Hypre_Subsidiary_Preconditioner_Helper::
//      set_hypre_JhalfStrnStrOneVTwoTwoRS();
//    
//    // Print it to check the settings.
//    Hypre_Subsidiary_Preconditioner_Helper::
//      print_hypre_settings(f_preconditioner_pt);
//#endif
//    }
//
//    // Now set the F preconditioner.
//    F_preconditioner_pt = f_preconditioner_pt;
//
//    // Set the preconditioner in the LSC preconditioner.
//    ns_preconditioner_pt->set_f_preconditioner(f_preconditioner_pt);
//
//    ///////////////////////////////////////// FFFFFFFFFFFFFFFFFF
//
//
//    // P block solve.
//    ///////////////////////////////////////////
//
//    // Pointer to the preconditioner.
//    Preconditioner * p_preconditioner_pt = 0;
//
//    //SL::P_solver == 0 is default, so do nothing.
//    if(P_solver == 1) 
//    {
//#ifdef OOMPH_HAS_HYPRE
//
//      p_preconditioner_pt = new HyprePreconditioner;
//
//      // Cast it to a Hypre preconditioner so we can set AMG settings.
//      HyprePreconditioner* hypre_preconditioner_pt =
//        static_cast<HyprePreconditioner*>(p_preconditioner_pt);
//
//      Hypre_default_settings::
//        set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);
//
//      if(Print_hypre)
//      {
//        Hypre_Subsidiary_Preconditioner_Helper::print_hypre_settings(
//            p_preconditioner_pt);
//      }
//
//      // Set it as the p preconditioner for LSC
//      //     ns_preconditioner_pt->set_p_preconditioner(p_preconditioner_pt);
//#endif
//    }
//    else if(P_solver == 13)
//    {
// #ifdef OOMPH_HAS_HYPRE
//
//      p_preconditioner_pt = new HyprePreconditioner;
//
//      // Cast it to a Hypre preconditioner so we can set AMG settings.
//      HyprePreconditioner* hypre_preconditioner_pt =
//        static_cast<HyprePreconditioner*>(p_preconditioner_pt);
//
//      Hypre_default_settings::
//        set_defaults_for_3D_poisson_problem(hypre_preconditioner_pt);
//
//      if(Print_hypre)
//      {
//        Hypre_Subsidiary_Preconditioner_Helper::print_hypre_settings(
//            p_preconditioner_pt);
//      }
//
//      // Set it as the p preconditioner for LSC
//      //     ns_preconditioner_pt->set_p_preconditioner(p_preconditioner_pt);
//#endif
//     
//    }
//    else if(P_solver == 96)
//    {
//#ifdef OOMPH_HAS_HYPRE
//
//      p_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//        get_custom_hypre_preconditioner(
//            p_amg_iterations, p_amg_smoother_iterations, 
//            p_amg_simple_smoother, p_amg_complex_smoother,
//            p_amg_damping, p_amg_strength,
//            p_amg_coarsening);
//
//      if(Print_hypre)
//      {
//        Hypre_Subsidiary_Preconditioner_Helper::print_hypre_settings(
//            p_preconditioner_pt);
//      }
//#endif
//    }
//    else if(P_solver == 2)
//    {
//#ifdef OOMPH_HAS_HYPRE
//      p_preconditioner_pt = new HyprePreconditioner;
//
//      HyprePreconditioner* hypre_preconditioner_pt =
//        static_cast<HyprePreconditioner*>(p_preconditioner_pt);
//
//      hypre_preconditioner_pt->hypre_method() = HyprePreconditioner::BoomerAMG;
//
//      // Setup v-cycles
//      hypre_preconditioner_pt->set_amg_iterations(2);
//      hypre_preconditioner_pt->amg_smoother_iterations() = 2;
//
//      // Setup smoother
//      // simple: 0 - DJ, 1 - GS
//      // compelx: Pilut - 7
//      hypre_preconditioner_pt->amg_using_simple_smoothing();
//      hypre_preconditioner_pt->amg_simple_smoother() = 0;
//      // only applicable for DJ
//      hypre_preconditioner_pt->amg_damping() = 0.8;
//
//      // Setup coarsening
//      // 0 - CLJP
//      // 1 - RS
//      hypre_preconditioner_pt->amg_coarsening() = 1;
//#endif
//    }
//    P_preconditioner_pt = P_preconditioner_pt;
//
//    ns_preconditioner_pt->set_p_preconditioner(p_preconditioner_pt);
//
//
//    return ns_preconditioner_pt;
//
//
//
//  } // EoFunc get_lsc_preconditioner(...)
//
//  inline Preconditioner* get_lgr_preconditioner()
//  {
//
//    LagrangeEnforcedflowPreconditioner* prec_pt
//      = new LagrangeEnforcedflowPreconditioner;
//
////    SimpleAugmentationPreconditioner* prec_pt
////      = new SimpleAugmentationPreconditioner;
//   // Set the mesh
//    if(Mesh_pt.size() < 2)
//    {
//      std::ostringstream err_msg;
//      err_msg << "There must be at least two meshes.\n"
//        << "Since the mesh size is 1, did you mean to use lsc only?\n"
//        << "If so, set --lsc_only and --f_solver --p_solver\n"
//        << std::endl;
//      throw OomphLibError(err_msg.str(),
//          OOMPH_CURRENT_FUNCTION,
//          OOMPH_EXCEPTION_LOCATION);
//    }
//
//    prec_pt->set_meshes(Mesh_pt);
//
//    // Set W solver.
//    if(W_solver == -1)
//    {
//      std::ostringstream err_msg;
//      err_msg << "There W_solver has not been set.\n"
//        << std::endl;
//      throw OomphLibError(err_msg.str(),
//          OOMPH_CURRENT_FUNCTION,
//          OOMPH_EXCEPTION_LOCATION);
//    }
//    else if(W_solver == 0)
//    {
//      // Using SuperLU, this is the default, do nothing.
//    }
//    else if (W_solver ==1)
//    {
//      prec_pt->set_lagrange_multiplier_subsidiary_preconditioner
//        (Lagrange_Enforced_Flow_Preconditioner_Subsidiary_Operator_Helper
//         ::get_lagrange_multiplier_preconditioner);
//    }
//    else
//    {
//      std::ostringstream err_msg;
//      err_msg << "There is no other W solver set.\n"
//        << std::endl;
//      throw OomphLibError(err_msg.str(),
//          OOMPH_CURRENT_FUNCTION,
//          OOMPH_EXCEPTION_LOCATION);
//    }
//
//    // The preconditioner for the fluid block:
//    if(NS_solver == -1)
//    {
//      std::ostringstream err_msg;
//      err_msg << "The NS solver has not been set.\n"
//        << std::endl;
//      throw OomphLibError(err_msg.str(),
//          OOMPH_CURRENT_FUNCTION,
//          OOMPH_EXCEPTION_LOCATION);
//    }
//    else if(NS_solver == 0) // Exact solve.
//    {
//      // This is the default, do nothing.
//      // But the param's F_solver and P_solver should not have been set,
//      // i.e. it should stay as -1.
//
//      if((F_solver != -1) || (P_solver != -1))
//      {
//        std::ostringstream err_msg;
//        err_msg << "Doing exact NS solve. (NS_solver is 0)\n"
//          << "but you have set F_solver and P_solver as well."
//          << "Please leave these as -1.\n"
//          << std::endl;
//        throw OomphLibError(err_msg.str(),
//            OOMPH_CURRENT_FUNCTION,
//            OOMPH_EXCEPTION_LOCATION);
//      }
//    }
//    else if(NS_solver == 1) // LSC
//    {
//
//
//      NS_preconditioner_pt = get_lsc_preconditioner();
//      // Set the NS preconditioner as LSC.
//      prec_pt->set_navier_stokes_lsc_preconditioner(NS_preconditioner_pt);
//
//
//     } // if for using LSC as NS prec.
//    else
//    {
//      pause("There is no solver for NS.");
//    }
//
//    if(!Use_axnorm)
//    {
//      prec_pt->scaling_sigma() = Scaling_sigma;
//    }
//
//    // Set the doc info for book keeping purposes. First we check that it is
//    // actually set.
//    if(Doc_linear_solver_info_pt == 0)
//    {
//      std::ostringstream err_msg;
//      err_msg << "Please set Doc_linear_solver_info_pt\n"
//        << std::endl;
//      throw OomphLibError(err_msg.str(),
//          OOMPH_CURRENT_FUNCTION,
//          OOMPH_EXCEPTION_LOCATION);
//    }
//
//    prec_pt->set_doc_linear_solver_info_pt(Doc_linear_solver_info_pt);
//
//    if(Use_block_diagonal_w)
//    {
//      prec_pt->use_block_diagonal_w_block();
//    }
//    else
//    {
//      prec_pt->use_diagonal_w_block();
//    }
//
//    if(Doc_prec)
//    {
//      prec_pt->enable_doc_prec();
//    }
//
//    //     Set the label, use to output information from the preconditioner, such
//    //     as the block matrices and the rhs vector
//    if(Label_str_pt == 0)
//    {
//      std::ostringstream err_msg;
//      err_msg << "Please set Label_str_pt, this should point\n"
//        << "to the Label_str in NSPP namespace." << std::endl;
//      throw OomphLibError(err_msg.str(),
//          OOMPH_CURRENT_FUNCTION,
//          OOMPH_EXCEPTION_LOCATION);
//    }
//
//    prec_pt->set_label_pt(Label_str_pt);
//    prec_pt->set_doc_prec_directory_pt(&Doc_prec_dir_str);
//
//    Lgr_preconditioner_pt = prec_pt;
//    return prec_pt;
//  }
//
//  inline Preconditioner* get_preconditioner()
//  {
//    if(Lsc_only)
//    {
//      return get_lsc_preconditioner();
//    }
//    else
//    {
//      return get_lgr_preconditioner();
//    }
//  } // LPH::setup_preconditioner
//
//  inline void clean_up_memory()
//  {
//    if(Lgr_preconditioner_pt != 0)
//    {
//      delete Lgr_preconditioner_pt;
//    }
//
//    if(NS_preconditioner_pt != 0)
//    {
//      delete NS_preconditioner_pt;
//    }
//
//    if(F_preconditioner_pt != 0)
//    {
//      delete F_preconditioner_pt;
//    }
//    if(P_preconditioner_pt != 0)
//    {
//      delete P_preconditioner_pt;
//    }
//
//   
//  }
//
//  inline std::string create_lsc_label()
//  {
//    std::string f_str = "";
//    std::string p_str = "";
//    // Now we continue with setting the string for the solvers.
//    // Only set the f_str if NS_solver > 0
//    if(NS_solver == 1 || Lsc_only)
//    {
//      switch(F_solver)
//      {
//        case 0:
//          f_str = "Fe";
//          break;
//        case 69:
//          f_str = "Fa";
//          break;
//        case 96:
//          f_str = "Fray";
//          break;
//        case 11:
//          f_str = "Fh2dp";
//          break;
//        case 12:
//          f_str = "Fhns";
//          break;
//        case 13:
//          f_str = "CLJPGSStrn075";
//          break;
//        case 14:
//          f_str = "FRSGSStrn075";
//          break;
//        case 15:
//          f_str = "FCLJPPilutStrn075";
//          break;
//        case 16:
//          f_str = "FRSPilutStrn075";
//          break;
//        case 17:
//          f_str = "Fray_old"; // I have no short hand for this...
//          break;
//        case 81:
//          f_str = "CLJPGSStrn0668";
//          break;
//        case 82:
//          f_str = "CLJPJStrn0668";
//          break;
//        case 83:
//          f_str = "CLJPPilutStrn0668";
//          break;
//        case 84:
//          f_str = "RSGSStrn0668";
//          break;
//        case 85:
//          f_str = "RSJStrn0668";
//          break;
//        case 86:
//          f_str = "RSPilutStrn0668";
//          break;
//        case 2:
//          f_str = "Fde";
//          break;
//        case 3:
//          f_str = "Fda";
//          break;
//        default:
//          {
//            std::ostringstream err_msg;
//            err_msg << "There is an unrecognised F_solver, recognised F_solver:\n"
//              << "Look at rayheader.h\n"
//              << std::endl;
//            throw OomphLibError(err_msg.str(),
//                OOMPH_CURRENT_FUNCTION,
//                OOMPH_EXCEPTION_LOCATION);
//          }
//      }  // switch for f_solver
//
//      switch(P_solver)
//      {
//        case 0:
//          p_str = "Pe";
//          break;
//        case 1:
//          p_str = "Pa";
//          break;
//        case 13:
//          p_str = "Pa3d";
//          break;
//        case 96:
//          p_str = "Pray";
//          break;
//        default:
//          {
//            std::ostringstream err_msg;
//            err_msg << "There is an unrecognised P_solver, recognised P_solver:\n"
//              << "Look at rayheader.h\n"
//              << std::endl;
//            throw OomphLibError(err_msg.str(),
//                OOMPH_CURRENT_FUNCTION,
//                OOMPH_EXCEPTION_LOCATION);
//          }
//      } // switch for p_solver
//    } // if ns_solve > 0
//
//
//    std::string prec_str = f_str + p_str;
//    return prec_str;
//  }
//
//  inline std::string create_lgr_label()
//  {
//    std::string w_str = "";
//
//    std::string ns_str = "";
//
//    std::string sigma_str = "";
//
//    // Set the string for W_solver.
//    switch(W_solver)
//    {
//      case 0:
//        w_str = "We";
//        break;
//      case 1:
//        w_str = "Wc";
//        break;
//      default:
//        {
//          std::ostringstream err_msg;
//          err_msg << "There is an unrecognised W_solver,\n"
//            << "recognised W_solver:\n"
//            << "0 = (We) SuperLU solve\n"
//            << std::endl;
//          throw OomphLibError(err_msg.str(),
//              OOMPH_CURRENT_FUNCTION,
//              OOMPH_EXCEPTION_LOCATION);
//        }
//    } // switch
//
//    if(Use_block_diagonal_w)
//    {
//      w_str += "bd";
//    }
//    else
//    {
//      w_str += "d";
//    }
//
//    // Set the string for NS_solver
//    switch(NS_solver)
//    {
//      case 0:
//        {
//          ns_str = "Ne";
//        }
//        break;
//      case 1:
//        ns_str = "Nl";
//        break;
//      default:
//        {
//          std::ostringstream err_msg;
//          err_msg << "There is an unrecognised NS_solver.\n"
//            << "Recognised NS_solver:\n"
//            << "0 = (Ne) SuperLU\n"
//            << "1 = (Nl) LSC preconditioner\n"
//            << std::endl;
//          throw OomphLibError(err_msg.str(),
//              OOMPH_CURRENT_FUNCTION,
//              OOMPH_EXCEPTION_LOCATION);
//        }
//    } // switch NS_solver
//
//    if(CommandLineArgs::command_line_flag_has_been_set("--sigma"))
//    {
//      std::ostringstream strs;
//      strs << "S" << Scaling_sigma;
//      sigma_str = strs.str();
//    }
//
//    std::string prec_str = w_str + ns_str + create_lsc_label() + sigma_str;
//    return prec_str;
//  }
//
//  inline std::string create_label()
//  {
//    std::string prec_str;
//    if(Lsc_only)
//    {
//      prec_str = create_lsc_label();
//    }
//    else
//    {
//      prec_str = create_lgr_label();
//    }
//
//    return prec_str;
//  } // LPH::create_label()
//
//} // end of namespace LagrangianPreconditionerHelpers


//==start_of_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Variables
{

 // Problem dimension
 static const unsigned Dim = 2;

 // Min and max x value respectively.
 static const double X_min = 0.0;
 static const double X_max = 1.0;

 // Min and max y value respectively.
 static const double Y_min = 0.0;
 static const double Y_max = 1.0;

 // The domain length in the x and y direction respectively.
 static const double Lx = X_max - X_min;
 static const double Ly = Y_max - Y_min;

 /// Reynolds number
 double Re = 100.0;
 double Ang_deg = 30.0;
 double Ang_rad = -1.0;
 unsigned Noel = 6;
 bool Use_lsc = false;

 
 inline double degtorad(const double& ang_deg)
 {
   return ang_deg * (MathematicalConstants::Pi / 180.0);
 }



 /// Storage for number of iterations during Newton steps 
// Vector<unsigned> Iterations;

 /// Storage for linear solver times during Newton steps 
// Vector<double> Linear_solver_time;

 /// Traction at the outflow boundary
// void prescribed_traction(const double& t,
//                          const Vector<double>& x,
//                          const Vector<double> &n,
//                          Vector<double>& traction)
// {
//  traction.resize(3);
//  traction[0]=1.0;
//  traction[1]=0.0;
//  traction[2]=0.0;
// } 

} // namespace Global_Variables



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


namespace oomph
{
//========================================================================
/// \short A Sloping Mesh  class.
///
/// Derived from RectangularQuadMesh:
/// the same mesh rotated with an angle phi (in radians)
//========================================================================
 template<class ELEMENT>
 class SlopingQuadMesh : public RectangularQuadMesh<ELEMENT>
 {
 public:

  /// Constructor.
  SlopingQuadMesh(const unsigned& nx, const unsigned& ny,
                  const double& lx,  const double& ly, const double& phi ) :
   RectangularQuadMesh<ELEMENT>(nx,ny,lx,ly)
   {
    // Find out how many nodes there are
    const unsigned n_node=this->nnode();

    // Loop over all nodes
    for (unsigned n=0;n<n_node;n++)
     {
      // Pointer to node:
      Node* nod_pt=this->node_pt(n);

      // Get the x/y coordinates
      const double x=nod_pt->x(0);
      const double y=nod_pt->x(1);

      // Set new nodal coordinates
      nod_pt->x(0)=x*cos(phi)-y*sin(phi);
      nod_pt->x(1)=x*sin(phi)+y*cos(phi);
     }
   }
 };
} // end of namespace oomph

//===start_of_problem_class=============================================
//======================================================================

template<class ELEMENT>
class TiltedCavityProblem : public Problem
{
public:

 /// \short Constructor: Pass number of elements in x and y directions and
 /// lengths
 TiltedCavityProblem();

 /// Update before solve is empty
 void actions_before_newton_solve()
 {
//   namespace NSPP = NavierStokesProblemParameters;
//   if(NSPP::Solver_type != NSPP::Solver_type_DIRECT_SOLVE)
//   {
    // Initialise counters for each newton solve.
//    Doc_linear_solver_info_pt->setup_new_time_step();
//   }
 }

 /// \short Update after solve is empty
 void actions_after_newton_solve()
 {
 }

 void actions_after_newton_step()
 {
//   namespace NSPP = NavierStokesProblemParameters;
//   if(NSPP::Solver_type != NSPP::Solver_type_DIRECT_SOLVE)
//   {
//     NSPP::doc_iter_times(this,Doc_linear_solver_info_pt);
//   }
 }

 void actions_before_distribute()
 {
//   namespace NSPP = NavierStokesProblemParameters;
//   namespace SL = SquareLagrange;

//   if(NSPP::Distribute_problem)
//   {
//     if(NSPP::Prob_id == SL::PID_SQ_PO)
//     {
//       GenericProblemSetup::delete_flux_elements(Surface_mesh_P_pt);

//      rebuild_global_mesh();
//     }
//     else
//     {
//   std::ostringstream err_msg;
//   err_msg << "Please set up the distributed bit for problem id: "
//           << NSPP::Prob_id << ".\n"
//           << std::endl;

//   throw OomphLibError(err_msg.str(),
//       OOMPH_CURRENT_FUNCTION,
//       OOMPH_EXCEPTION_LOCATION);
//     }
//   }
 }

 void actions_after_distribute()
 {
//   namespace NSPP = NavierStokesProblemParameters;
//   namespace SL = SquareLagrange;

//   if(NSPP::Distribute_problem)
//   {
//   if(NSPP::Prob_id == SL::PID_SQ_PO)
//   {
//     create_parall_outflow_lagrange_elements(1,
//         Bulk_mesh_pt,Surface_mesh_P_pt);
//     rebuild_global_mesh();
//   }
//   else
//   {
//   std::ostringstream err_msg;
//   err_msg << "Please set up the distributed bit for problem id: "
//           << NSPP::Prob_id << ".\n"
//           << std::endl;

//   throw OomphLibError(err_msg.str(),
//       OOMPH_CURRENT_FUNCTION,
//       OOMPH_EXCEPTION_LOCATION);
//   }
//   }
 }

 /// Doc the solution
 void doc_solution();

 /// \short Create lagrange elements on boundary b of the Mesh pointed
 /// to by bulk_mesh_pt and add them to the Mesh object pointed to by
 /// surface_mesh_pt
 void create_parall_outflow_lagrange_elements(const unsigned &b,
                                              Mesh* const &bulk_mesh_pt,
                                              Mesh* const &surface_mesh_pt);

private:


 /// Pointer to the "bulk" mesh
 Mesh* Bulk_mesh_pt;

 /// Pointer to the "surface" mesh
 Mesh* Surface_mesh_P_pt;

 // Preconditioner
 Preconditioner* Prec_pt;

 // Preconditioner for the Navier-Stokes block
 Preconditioner* Navier_stokes_prec_pt;

 // Solver
 IterativeLinearSolver* Solver_pt;

 DocLinearSolverInfo* Doc_linear_solver_info_pt;

 // Enumeration for the sides
 //         2
 //    -------------
 //    |           |
 //  3 |           | 1
 //    |           |
 //    |           |
 //    -------------
 //         0
 //
 unsigned Bottom_bound; // 0
 unsigned Right_bound;  // 1
 unsigned Top_bound;    // 2
 unsigned Left_bound;   // 3

};


//===start_of_constructor=============================================
/// Problem constructor
//====================================================================
template<class ELEMENT> // rrrback - changed here.
TiltedCavityProblem<ELEMENT>::TiltedCavityProblem()
{
  // Alias the namespace for convenience
  namespace GV = Global_Variables;

  // Assign the enumerations for the boundaries.
  Bottom_bound = 0;
  Right_bound = 1;
  Top_bound = 2;
  Left_bound = 3;
  
//  Doc_linear_solver_info_pt = NSPP::Doc_linear_solver_info_pt;

  // First we set the tilted cavity mesh.
  Bulk_mesh_pt = new SlopingQuadMesh<ELEMENT>(GV::Noel,GV::Noel,
                                              GV::Lx,GV::Ly,
                                              GV::Ang_rad); 

  // Set the boundary conditions, recall that the boundaries are:
  //
  //             2 non slip
  //         ----------
  //         |        |
  // 3 Inflow|        |1 P.O. (parallel outflow)
  //         |        |
  //         ----------
  //             0 non slip

  const unsigned po_bound = Right_bound; // parallel outflow

  // Create a "surface mesh" that will contain only
  // ImposeParallelOutflowElements in boundary 1
  // The constructor just creates the mesh without
  // giving it any elements, nodes, etc.
  Surface_mesh_P_pt = new Mesh;

  // Create ImposeParallelOutflowElement from all elements that are
  // adjacent to the Neumann boundary.
  create_parall_outflow_lagrange_elements(po_bound,
                                          Bulk_mesh_pt,
                                          Surface_mesh_P_pt);

  // Add the two meshes to the problem.
  add_sub_mesh(Bulk_mesh_pt);
  add_sub_mesh(Surface_mesh_P_pt);

  // combine all sub-meshes into a single mesh.
  build_global_mesh();

  // All nodes are free by default
  // just pin the ones that have Dirichlet conditions
  const unsigned num_bound = mesh_pt()->nboundary();
  for(unsigned ibound=0;ibound<num_bound;ibound++)
  { 
    if(ibound != po_bound)
    {
      unsigned num_nod=mesh_pt()->nboundary_node(ibound);
      for (unsigned inod=0;inod<num_nod;inod++)
      {
        // Get node
        Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
 
        nod_pt->pin(0);
        nod_pt->pin(1);
        
        nod_pt->set_value(0,0);
        nod_pt->set_value(1,0);
      }
    }
  }

  // Which is the inflow boundary?
  const unsigned in_bound = Left_bound;

  // The number of nodes on a boundary
  unsigned num_nod;

  num_nod = mesh_pt()->nboundary_node(in_bound);
  for(unsigned inod=0;inod<num_nod;inod++)
  {
    Node* nod_pt=mesh_pt()->boundary_node_pt(in_bound,inod);

    // Pin both velocity components
    nod_pt->pin(0);
    nod_pt->pin(1);
    
    // Get the x and y cartesian coordinates
    double x0=nod_pt->x(0);
    double x1=nod_pt->x(1);

    // Tilt x1 by -SL::Ang, this will give us the original coordinate.
    double x1_old = x0*sin(-GV::Ang_rad) + x1*cos(-GV::Ang_rad);

    // Now calculate the parabolic inflow at this point.
    double u0_old = (x1_old - GV::Y_min)*(GV::Y_max - x1_old);
   
    // Now apply the rotation to u0_old, using rotation matrices.
    // with x = u0_old and y = 0, i.e. R*[u;0] since we have the
    // velocity in the x direction only. There is no velocity
    // in the y direction.
    double u0=u0_old*cos(GV::Ang_rad);
    double u1=u0_old*sin(GV::Ang_rad);

    nod_pt->set_value(0,u0);
    nod_pt->set_value(1,u1);
  }

 
  //Complete the problem setup to make the elements fully functional

  //Loop over the elements
  unsigned n_el = Bulk_mesh_pt->nelement();
  for(unsigned e=0;e<n_el;e++)
  {
    //Cast to a fluid element
    ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

    //Set the Reynolds number, etc
    el_pt->re_pt() = &GV::Re;
  } // for(unsigned e=0;e<n_el;e++)

  //Assign equation numbers
  oomph_info << "\n equation numbers : "
             << assign_eqn_numbers() << std::endl;
 

  // Only do this bit if we do NOT have a direct solver.
  // Create the vector of mesh pointers!
  Vector<Mesh*> mesh_pt;
  mesh_pt.resize(2,0);
  mesh_pt[0] = Bulk_mesh_pt;
  mesh_pt[1] = Surface_mesh_P_pt;

  // Create the preconditioner.
  LagrangeEnforcedflowPreconditioner* prec_pt
    = new LagrangeEnforcedflowPreconditioner;

  prec_pt->set_meshes(mesh_pt);
  
  NavierStokesSchurComplementPreconditioner* lsc_prec_pt = 0;
  if(GV::Use_lsc)
  {
    // Create the NS LSC preconditioner.
    lsc_prec_pt = new NavierStokesSchurComplementPreconditioner(this);
    lsc_prec_pt->set_navier_stokes_mesh(Bulk_mesh_pt);
  }

  prec_pt->set_navier_stokes_lsc_preconditioner(lsc_prec_pt);

  Navier_stokes_prec_pt = lsc_prec_pt;

  Prec_pt = prec_pt;


  IterativeLinearSolver* solver_pt = new GMRES<CRDoubleMatrix>;
  // We use RHS preconditioning. Note that by default,
  // left hand preconditioning is used.
  static_cast<GMRES<CRDoubleMatrix>*>(solver_pt)->set_preconditioner_RHS();

  solver_pt->tolerance() = 1.0e-6;
  solver_pt->max_iter() = 100;
  solver_pt->preconditioner_pt() = Prec_pt;
  this->linear_solver_pt() = solver_pt;
  this->newton_solver_tolerance() = 1.0e-6;
}

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::doc_solution()
{

//  namespace NSPP = NavierStokesProblemParameters;
  
//  std::ofstream some_file;
//  std::stringstream filename;
//  filename << NSPP::Soln_dir_str<<"/"<<NSPP::Label_str<<".dat";

  // Number of plot points
//  unsigned npts=5;

  // Output solution
//  some_file.open(filename.str().c_str());
//  Bulk_mesh_pt->output(some_file,npts);
//  some_file.close();
}

//============start_of_create_parall_outflow_lagrange_elements===========
/// Create ImposeParallelOutflowElement on the b-th boundary of the
/// Mesh object pointed to by bulk_mesh_pt and add the elements to the
/// Mesh object pointeed to by surface_mesh_pt.
//=======================================================================
template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::
create_parall_outflow_lagrange_elements(const unsigned &b,
                                        Mesh* const &bulk_mesh_pt,
                                        Mesh* const &surface_mesh_pt)
{
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = bulk_mesh_pt->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    bulk_mesh_pt->boundary_element_pt(b,e));

   // What is the index of the face of element e along boundary b?
   int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);

   // Build the corresponding impose_impenetrability_element
   ImposeParallelOutflowElement<ELEMENT>* flux_element_pt = new
    ImposeParallelOutflowElement<ELEMENT>(bulk_elem_pt,
                                          face_index,0);


   // Add the prescribed-flux element to the surface mesh
   surface_mesh_pt->add_element_pt(flux_element_pt);

   // Loop over the nodes
   unsigned nnod=flux_element_pt->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt = flux_element_pt->node_pt(j);

     // Is the node also on boundary 0 or 2?
     if ((nod_pt->is_on_boundary(0))||(nod_pt->is_on_boundary(2)))
      {
       // How many nodal values were used by the "bulk" element
       // that originally created this node?
       unsigned n_bulk_value=flux_element_pt->nbulk_value(j);

       // The remaining ones are Lagrange multipliers and we pin them.
       unsigned nval=nod_pt->nvalue();
       for (unsigned j=n_bulk_value;j<nval;j++)
        {
         nod_pt->pin(j);
        }
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//==start_of_main======================================================
/// Driver for Lagrange enforced flow preconditioner
//=====================================================================
int main(int argc, char **argv)
{
#ifdef OOMPH_HAS_MPI
 MPI_Helpers::init(argc,argv);
#endif


 // Alias the namespace for convenience.
 namespace GV = Global_Variables;
// namespace LPH = LagrangianPreconditionerHelpers;
// namespace SL = SquareLagrange;

 // Set up doc info - used to store information on solver and 
 // iteration time.
// DocLinearSolverInfo doc_linear_solver_info;
// NSPP::Doc_linear_solver_info_pt = &doc_linear_solver_info;



 //Label for output
// DocInfo doc_info;
 
// //Set output directory
// doc_info.set_directory("RESLT");
 
 //Doc number of gmres iterations
// ofstream out_file;

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

// CommandLineArgs::specify_command_line_flag("--dist_prob");

 // Flag to output the solution.
// CommandLineArgs::specify_command_line_flag("--doc_soln", 
//        &Soln_dir_str);

 // int, ini to -1
// CommandLineArgs::specify_command_line_flag("--visc", &NSPP::Vis);
 
 // Default is 100.0
 CommandLineArgs::specify_command_line_flag("--re", &GV::Re);
 
 // int init to -1
// CommandLineArgs::specify_command_line_flag("--max_solver_iter", 
//                                            &NSPP::Max_solver_iteration);

//    // Iteration count and times directory.
//    CommandLineArgs::specify_command_line_flag("--itstimedir", 
//        &Itstime_dir_str);

 // int init to -1
// CommandLineArgs::specify_command_line_flag("--solver_type",
//                                            &NSPP::Solver_type);

 // int init to -1, takes values 0, 1, 2
// CommandLineArgs::specify_command_line_flag("--time_type",
//                                            &NSPP::Time_type);

 // double init to -1
// CommandLineArgs::specify_command_line_flag("--dt", &NSPP::Delta_t);

 // double init to -1
// CommandLineArgs::specify_command_line_flag("--time_start", 
//                                            &NSPP::Time_start);

 // double init to -1
// CommandLineArgs::specify_command_line_flag("--time_end", 
//                                            &NSPP::Time_end);

 ///////////////////////////////////////////////////////////////
 // Default to 30 degrees
 CommandLineArgs::specify_command_line_flag("--ang", &GV::Ang_deg);

 // Default to 6
 CommandLineArgs::specify_command_line_flag("--noel", &GV::Noel);


 ///////////////////////////////////////////////////////////////
    // Flag to output the preconditioner, used for debugging.
    // string
//    CommandLineArgs::specify_command_line_flag(
//        "--doc_prec",&Doc_prec_dir_str);

    // No parameter after.
//    CommandLineArgs::specify_command_line_flag(
//        "--lsc_only");

    // double
//    CommandLineArgs::specify_command_line_flag(
//        "--sigma",&Scaling_sigma);

    // int
//    CommandLineArgs::specify_command_line_flag(
//        "--w_solver",&W_solver);

    // Nothing set
//    CommandLineArgs::specify_command_line_flag("--bdw");

   CommandLineArgs::specify_command_line_flag("--use_lsc");

    // int
//    CommandLineArgs::specify_command_line_flag(
//        "--ns_solver",&NS_solver);

    // int
//    CommandLineArgs::specify_command_line_flag(
//        "--p_solver",&P_solver);

    // int
//    CommandLineArgs::specify_command_line_flag(
//        "--f_solver",&F_solver);

    // NS_F block AMG parameters
    // double
//    CommandLineArgs::specify_command_line_flag(
//        "--f_amg_str",&F_amg_strength);
    // double
//    CommandLineArgs::specify_command_line_flag(
//        "--f_amg_damp",&F_amg_damping);

    // int
//    CommandLineArgs::specify_command_line_flag("--f_amg_coarse",
//        &F_amg_coarsening);

    // int
//    CommandLineArgs::specify_command_line_flag(
//        "--f_amg_sim_smoo",&F_amg_simple_smoother);

    // int
//    CommandLineArgs::specify_command_line_flag(
//        "--f_amg_com_smoo",&F_amg_complex_smoother);

    // int
//    CommandLineArgs::specify_command_line_flag(
//        "--f_amg_iter",&F_amg_iterations);

    //int
//    CommandLineArgs::specify_command_line_flag("--f_amg_smiter",
//        &F_amg_smoother_iterations);

    // NS_P block AMG parameters
    // double
//    CommandLineArgs::specify_command_line_flag(
//        "--p_amg_str",&P_amg_strength);
    // double
//    CommandLineArgs::specify_command_line_flag(
//        "--p_amg_damp",&P_amg_damping);
    // int
//    CommandLineArgs::specify_command_line_flag(
//        "--p_amg_coarse",&P_amg_coarsening);
    // int
//    CommandLineArgs::specify_command_line_flag(
//        "--p_amg_sim_smoo",&P_amg_simple_smoother);
    // int
//    CommandLineArgs::specify_command_line_flag(
//        "--p_amg_com_smoo",&P_amg_complex_smoother);

    // int
//    CommandLineArgs::specify_command_line_flag(
//        "--p_amg_iter",&P_amg_iterations);
    // int
//    CommandLineArgs::specify_command_line_flag("--p_amg_smiter",
//        &P_amg_smoother_iterations);

//    CommandLineArgs::specify_command_line_flag("--print_hypre");
 


 //    CommandLineArgs::specify_command_line_flag("--mesh_type", &Mesh_type);

 // Parse the above flags.
 CommandLineArgs::parse_and_assign();
 CommandLineArgs::doc_specified_flags();

 GV::Ang_rad = GV::degtorad(GV::Ang_deg);
 if(CommandLineArgs::command_line_flag_has_been_set("--use_lsc"))
 {
   GV::Use_lsc = true;
 }
 else
 {
   GV::Use_lsc = false;
 }

 
 TiltedCavityProblem< QTaylorHoodElement<GV::Dim> > problem;

    // Solve the problem
    problem.newton_solve();

 

#ifdef OOMPH_HAS_MPI
   MPI_Helpers::finalize();
#endif
   
   
} // end_of_main


