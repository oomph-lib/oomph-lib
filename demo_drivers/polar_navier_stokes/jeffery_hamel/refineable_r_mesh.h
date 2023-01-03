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
//Header file for my Refineable_r_mesh class

namespace oomph
{

//======================================================================
///  My log r spaced mesh
//======================================================================
template <class ELEMENT>
class Refineable_r_mesh : public RectangularQuadMesh<ELEMENT>,
public RefineableQuadMesh<ELEMENT>
{
protected:

 //Storage for pointers to fluid elements
 Vector<GeneralisedElement*> Fluid_elt_pt;

public:

 /// Return  pointer to fluid element e
 ELEMENT* fluid_elt_pt(unsigned e) 
   {return Fluid_elt_pt[e];}
 /// Return length of fluid element vector
 unsigned fluid_elt_length() {return Fluid_elt_pt.size();}
 /// Return pointer to vector of all Fluid elements
 Vector<GeneralisedElement*> fluid_elt_vector() {return Fluid_elt_pt;}

 /// Constructor, which "builds" the mesh. The arguments are the number
 /// of elements in each direction.
 Refineable_r_mesh(const unsigned int &nx,const unsigned int &ny) :
 RectangularQuadMesh<ELEMENT>(nx,ny,Global_Physical_Variables::R_l,
                              Global_Physical_Variables::R_r,
                              -1.0,1.0)
  {  
   //Call the generic mesh construction routine
   //(Defined in RectangularQuadMesh<ELEMENT>)
   //this->build_mesh();

   //Function to log space the mesh
   if(Global_Physical_Variables::log_mesh) stretch_mesh();

   // Nodal positions etc. were created in constructor for
   // RectangularMesh<...>. Only need to setup quadtree forest
   this->setup_quadtree_forest();
  }

private:

  void stretch_mesh()
  {
    unsigned long n_node = this -> nnode();
    //Get number of elements in x
    unsigned Nx = this->Nx;

    //Look up the size of the mesh
    double R_min = this->Xmin;
    double R_max = this->Xmax;

    //Work out a best guess for the last 7% of elements in x
    //Always rounds down due to truncation
    int end=static_cast<int>(0.07*Nx);
    //int end=static_cast<int>(0.2*Nx);
    if(end==0)
     {
      std::cout << std::endl << "Warning, no elements in outlet region" << std::endl;
      std::cout << "Adjusting to one element" << std::endl;
      end=1;
     }

    //Work out how where that element boundary will fall:
    double cut=R_min+(R_max-R_min)*((static_cast<double>(Nx-end))/Nx);

    //Now spread those final 7% of elements out over the final 10% of the domain
    double a=cut;
    double b=R_min+0.9*(R_max-R_min);
       
    if(Global_Physical_Variables::new_outlet_region)
     {
      //Output what I'm doing
      std::cout << std::endl << "Number of end elements in r: " << end << std::endl;
      std::cout << "Start of outlet region: " << cut << std::endl << std::endl;
     }

    for(unsigned long n=0;n<n_node;n++)
      {
       //Work out old nodal position:
       double r_pos = this -> Node_pt[ n ] -> x( 0 );

       if(Global_Physical_Variables::new_outlet_region)
	{
	 if(r_pos<a)
	  {
	   //Log space until old_r = a
	   this -> Node_pt[ n ] -> x( 0 ) = log_spacing( r_pos,a,b,R_min);
	  }
	 else
	  {
	   //Linear spacing from old_r = a to 1
	   this -> Node_pt[ n ] -> x( 0 ) = linear_spacing(r_pos,a,b,R_max);
	  }	    
	}
       else
	{
	 //Log space across whole domain
	 this -> Node_pt[ n ] -> x( 0 ) = log_spacing( r_pos,R_max,R_max,R_min);
	}
      }
	
  }//End of stretch mesh

  double log_spacing( const double& r,const double& a,const double& b,double R_min)
   {
    return R_min*exp(((r-R_min)/(a-R_min))*log(b/R_min));
   }

  double linear_spacing( const double& r,const double& a,const double& b,double R_max)
   {
    if((R_max-a)<1.e-14) 
    {std::cout << "Warning - no elements in outlet region.  Attempting to divide by zero" << std::endl;}
    return R_max+(r - R_max)*((R_max-b)/(R_max-a));
   }

}; //End of Refineable_r_mesh class

} //End of namespace oomph
