//
// Created by iqraa on 13-9-24.
//
//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented,
//LIC// multi-physics finite-element library, available
//LIC// at http://www.oomph-lib.org.
//LIC//
//LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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

//Oomph-lib includes
#include "generic.h"
#include "solid.h"
#include "constitutive.h"

// The mesh form gmsh
#include "meshes/gmsh_brick_mesh.h"

using namespace std;
using namespace oomph;



/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////


//=======start_namespace==========================================
/// Global variables
//================================================================
namespace Input
{
    /// Elastic modulus
    double E=1.0e8;

    /// Poisson's ratio
    double Nu=0.3;


    /// Pointer to constitutive law
    ConstitutiveLaw* Constitutive_law_pt= nullptr;
    
    /// Density
    double Density = 2500;
    
    /// Gravity
    double Gravity=9.81;

    /// Gravity as body force
    void body_force(const double& time,
                    const Vector<double> &xi,
                    Vector<double> &b)
    {
        b[0]=0.0;
        b[1]=0.0;
        b[2]=-Gravity*Density;
    }
}


template<class ELEMENT>
class CantileverProblem : public Problem
{

public:

    /// Constructor:
    CantileverProblem(const string& mesh_file)
    {
        /// Verbose to print all Gmsh info/entities/elements.
        bool Verbose = false;

        Problem::mesh_pt() = new SolidGmshBrickMesh<ELEMENT>(mesh_file, Verbose);


        // Complete build of elements
        unsigned n_element=mesh_pt()->nelement();
        for(unsigned i=0;i<n_element;i++)
        {
            // Cast to a solid element
            auto *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

            // Set the constitutive law
            el_pt->constitutive_law_pt() = Input::Constitutive_law_pt;

            // Set the body force
            el_pt->body_force_fct_pt() = Input::body_force;
            
            // Set density
            el_pt->lambda_sq_pt() = &Input::Density;

        } // done build of elements


        // Pin the left boundary (boundary 0) in all directions
        unsigned b=0;
        unsigned n_side = mesh_pt()->nboundary_node(b);

        //Loop over the nodes
        for(unsigned i=0;i<n_side;i++)
        {
            mesh_pt()->boundary_node_pt(b,i)->pin_position(0);
            mesh_pt()->boundary_node_pt(b,i)->pin_position(1);
            mesh_pt()->boundary_node_pt(b,i)->pin_position(2);
        }

        // Pin the redundant solid pressures (if any)
        PVDEquationsBase<3>::pin_redundant_nodal_solid_pressures(mesh_pt()->element_pt());

        //Assign equation numbers
        assign_eqn_numbers();

        // Prepare output directory
        Doc_info.set_directory("RESLT");

    }

    /// Update function (empty)
    void actions_after_newton_solve() {}

    /// Update function (empty)
    void actions_before_newton_solve() {}

    /// Actions before adapt. Empty
    void actions_before_adapt(){}

    /// Actions after adapt
    void actions_after_adapt()
    {
        // Pin the redundant solid pressures (if any)
        PVDEquationsBase<3>::pin_redundant_nodal_solid_pressures(mesh_pt()->element_pt());
    }

    /// Doc the solution
    void doc_solution()
    {

        ofstream some_file;
        char filename[100];

        // Number of plot points
        unsigned n_plot = 2;

        // Output shape of deformed body
        sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),Doc_info.number());
        some_file.open(filename);
        mesh_pt()->output(some_file,n_plot);
        some_file.close();

        // Increment label for output files
        Doc_info.number()++;
    };


    /// Access function for the mesh
    SolidGmshBrickMesh<ELEMENT>* mesh_pt()
    {
        return dynamic_cast<SolidGmshBrickMesh<ELEMENT>*>(Problem::mesh_pt());
    }


    void solve()
    {
        // Set output directory
        char dirname[100];
        sprintf(dirname,"RESLT");

        // Prepare output
        Doc_info.set_directory(dirname);

        // Doc solution
        doc_solution();

        //Parameter incrementation
        unsigned nstep=2;

        double g_increment=25.0e-2;
        for(unsigned i=0;i<nstep;i++)
        {
            // Increment load
            Input::Gravity+=g_increment;

            // Solve it
            newton_solve();

            // Doc solution
            doc_solution();
        }
    };

private:

    /// DocInfo object for output
    DocInfo Doc_info;
    
    /// Bulk mesh pointer
    SolidGmshBrickMesh<ELEMENT>* Solid_mesh_pt;
};




//=======start_of_main==================================================
/// Driver for 3D cantilever beam loaded by gravity
//======================================================================
int main(int argc, char* argv[])
{
    // Store command line arguments
    CommandLineArgs::setup(argc,argv);

    // Check number of command line arguments: Need exactly three.
    if (argc!=2)
    {
        std::string error_message =
            "Wrong number of command line arguments.\n";
        error_message +=
            "Must specify the following file name  \n";
        error_message += 
            "cantilever.msh\n";

   throw OomphLibError(error_message,
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }

    // Convert arguments to strings that specify the input file names
    string Mesh_file(argv[1]);


    // Generalised Hookean constitutive equations
    Input::Constitutive_law_pt = new GeneralisedHookean(&Input::Nu, &Input::E);

    //Set up the problem with pure displacement based elements
    CantileverProblem<QPVDElement<3,2> > problem(Mesh_file);
    problem.solve();


    delete Input::Constitutive_law_pt;
    Input::Constitutive_law_pt= nullptr;


} //end of main






