// Biot Savart

// Generic oomph-lib header
#include "generic.h"

using namespace std;

using namespace oomph;


//=========================================================================
///  1D line in 3D space 
//=========================================================================
class MyLine : public GeomObject
{

public:
 
 /// Constructor
 MyLine(const double& ampl_bend_x, const double start_bend_x,
        const double& l_bend_x) :
  GeomObject(1,3), Ampl_bend_x(ampl_bend_x), 
  Start_bend_x(start_bend_x), L_bend_x(l_bend_x)
  {}
 
 
 /// Broken copy constructor
 MyLine(const MyLine& dummy) 
  { 
   BrokenCopy::broken_copy("MyLine");
  } 
 
 /// Broken assignment operator
 void operator=(const MyLine&) 
  {
   BrokenCopy::broken_assign("MyLine");
  }


 /// Destructor:  Empty
 ~MyLine() {}
 
 /// \short Position Vector at Lagrangian coordinate zeta 
 void position(const Vector<double>& zeta, Vector<double>& r) const
  {
   // Position Vector
   r[0] = 0.0;
   r[1] = 0.0;
   r[2] = zeta[0];

   // Bend
   if (zeta[0]>Start_bend_x)
    {
     double excess=(zeta[0]-Start_bend_x)/L_bend_x;
     if (excess>0.5) excess=0.5;
     r[0]+=Ampl_bend_x*0.5*(1.0-cos(2.0*MathematicalConstants::Pi*excess));
    }
  }



 /// \short Derivative of position Vector w.r.t. to coordinates: 
 /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i). 
 /// Evaluated at current time.
 virtual void dposition(const Vector<double>& zeta, 
                        DenseMatrix<double> &drdzeta) const
  {
   // Tangent vector
   drdzeta(0,0)=0.0;
   drdzeta(0,1)=0.0;
   drdzeta(0,2)=1.0;

   // Bend:
   if (zeta[0]>Start_bend_x)
    {
     double excess=(zeta[0]-Start_bend_x)/L_bend_x;
     if (excess>0.5) excess=0.5;
     drdzeta(0,0)+=
      Ampl_bend_x*sin(2.0*MathematicalConstants::Pi*excess)*
      MathematicalConstants::Pi/L_bend_x;
    }
  }

 /// Amplitude for bending in x-direction
 double Ampl_bend_x;

 /// Start point for bending in x-direction
 double Start_bend_x;

 /// Wavelength for bending in x-direction
 double L_bend_x;

};





//==start_of_main======================================================
/// Potential at position r, due to unit line charge distributed
/// along GeomObject with specified start and end points
//=====================================================================
double potential(const Vector<double> r, GeomObject* line_pt, 
                 const double& s_min, const double& s_max)
               
{

 // Arclength and positon on GeomObject
 Vector<double> s(1);
 Vector<double> r_line(3);
 DenseMatrix<double> dr_line(3,1);

 // Number of steps for trapezoidal rule
 unsigned nstep=1000;
 double ds=(s_max-s_min)/double(nstep-1);
 
 // First point
 s[0]=s_min;
 line_pt->position(s,r_line);
 line_pt->dposition(s,dr_line);

 // Distance and differential
 double dist=0.0;
 double drds=0.0;
 for (unsigned i=0;i<3;i++)
  {
   dist+= (r_line[i]-r[i])*(r_line[i]-r[i]);
   drds+=dr_line(0,i)*dr_line(0,i);
  }
 

 // Initialise potential
 double pot=0.5*drds/dist*ds;

 // Move along line and increment
 for (unsigned j=0;j<nstep-1;j++)
  {
   // Next point
   s[0]+=ds;
   line_pt->position(s,r_line);
   line_pt->dposition(s,dr_line);

   // Distance and differential
   dist=0.0;
   drds=0.0;
   for (unsigned i=0;i<3;i++)
    {
     dist+= (r_line[i]-r[i])*(r_line[i]-r[i]);
     drds+=dr_line(0,i)*dr_line(0,i);
    }
   
   // Add to potential
   pot+=drds/dist*ds;
   
  }
 

 // Final point
 s[0]=s_max;
 line_pt->position(s,r_line);
 line_pt->dposition(s,dr_line);

 // Distance and differential
 dist=0.0;
 drds=0.0;
 for (unsigned i=0;i<3;i++)
  {
   dist+= (r_line[i]-r[i])*(r_line[i]-r[i]);
   drds+=dr_line(0,i)*dr_line(0,i);
  }

 // Add to potential
 pot+=0.5*drds/dist*ds;

 return pot;

}

//==start_of_main======================================================
/// 
//=====================================================================
int main()
{

 // Charged line
 double ampl_bend1_x= 0.5;
 double ampl_bend2_x=-0.5;
 double start_bend_x=0.0;
 double l_bend_x=2.0;
 GeomObject* line1_pt=new MyLine(ampl_bend1_x,start_bend_x,l_bend_x);
 GeomObject* line2_pt=new MyLine(ampl_bend2_x,start_bend_x,l_bend_x);
 double s_min=-10.5; 
 double s_max= 10.5;
 
 /// Loop over regularly spaced 3D plot points
 unsigned nplot=20;
 double limit=1.0;

 double x_min=-limit;
 double x_max= limit;
 
 double y_min=-limit;
 double y_max= limit;

 double z_min=-limit;
 double z_max= 2.0*limit;


 Vector<double> x(3);

 ofstream outfile("potential.dat");
 outfile << "ZONE I=" << nplot << ",J=" << nplot << ",K=" << nplot << std::endl;
 
 for (unsigned ix=0;ix<nplot;ix++)
  {
   x[0]=x_min+(x_max-x_min)*double(ix)/double(nplot-1);
   for (unsigned iy=0;iy<nplot;iy++)
    {
     x[1]=y_min+(y_max-y_min)*double(iy)/double(nplot-1);
     for (unsigned iz=0;iz<nplot;iz++)
      {
       x[2]=z_min+(z_max-z_min)*double(iz)/double(nplot-1);
              
       double pot=potential(x, line1_pt, s_min,s_max);
       pot+=potential(x, line2_pt, s_min,s_max);

       outfile << x[0] << " " << x[1] << " " << x[2] << " " << pot << std::endl;
      }
    }
  }

 outfile.close();

} // end_of_main











