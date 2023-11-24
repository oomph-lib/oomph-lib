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

 private:

  /// dimension
  unsigned el_dim;

  /// Storage for elemental error estimate -- used for post-processing
  double Error;

  /// Number of previous history values to be used for  extrapolation
  /// of strain rate
  unsigned Nprev_for_extrapolation_of_strain_rate;

  /// Boolean to indicate if we're using a fixed point iteration for
  /// the strain rate that forms the basis for the invariant
  bool Use_fixed_point_for_strain_rate;

  /// Boolean to indicate whether we use Aitken extrapolation
  /// during the fixed point iterations
  bool Use_aitken_extrapolation;



 public:

  /// Constructor initialise error and set default for number of 
  /// previous history values to be used for extrapolation of strain rate
  MyTaylorHoodElement()
   {
    el_dim=2;
    Error=0.0;
    Nprev_for_extrapolation_of_strain_rate=4;
    Use_fixed_point_for_strain_rate=false;
    Use_aitken_extrapolation=false;

    // Make space for fixed point iteration on strain rate
    unsigned fp_n_val = 3;
    Fixed_point_iteration_guess_for_strain_rate.resize(fp_n_val);
    unsigned n_intpt = integral_pt()->nweight();
    for (unsigned val=0;val<fp_n_val;val++)
     {
      Fixed_point_iteration_guess_for_strain_rate[val].resize(n_intpt);
      for (unsigned ipt=0;ipt<n_intpt;ipt++)
       {
        Fixed_point_iteration_guess_for_strain_rate[val][ipt].resize(el_dim,
                                                                     el_dim,
                                                                     0.0);
       }
     }

    Aitken_index=0;
   }
  
  
  /// Number of previous history values to be used for  extrapolation
  /// of strain rate
  unsigned& nprev_for_extrapolation_of_strain_rate()
   {
    return Nprev_for_extrapolation_of_strain_rate;
   }
  

  /// Enable use of fixed point iteration (sets
  /// current best guess based on extrapolation)
  void enable_fixed_point_iteration_for_strain_rate()
   {
    Aitken_index=0;
    update_latest_fixed_point_iteration_guess_for_strain_rate();
    Use_fixed_point_for_strain_rate=true;
   }

  /// Disable use of fixed point iteration
  void disable_fixed_point_iteration_for_strain_rate()
   {
    Use_fixed_point_for_strain_rate=false;
    Use_aitken_extrapolation=false;
   }

  /// Enable use of Aitken extrapolation
  void enable_aitken_extrapolation()
   {
    Use_aitken_extrapolation=true;
   }

  /// Disable use of Aitken extrapolation
  void disable_aitken_extrapolation()
   {
    Use_aitken_extrapolation=false;
   }
  
  /// Return latest guess (obtained via fixed point iteration)
  /// for strain rate at integration point ipt
  void latest_fixed_point_iteration_guess_for_strain_rate(
   const unsigned& ipt, DenseMatrix<double>& strainrate) const
   {
    for (unsigned i=0;i<el_dim;i++)
     {
      for (unsigned j=0;j<el_dim;j++)
       {
        strainrate(i,j)=
         Fixed_point_iteration_guess_for_strain_rate[Aitken_index-1][ipt](i,j);
       }
     }
   }
 
  
  /// Update latest guess (obtained via fixed point iteration)
  /// for strain rate from current actual strain rate
  void update_latest_fixed_point_iteration_guess_for_strain_rate()
   {
    Vector<double> s(el_dim);
    DenseMatrix<double> strain_rate(el_dim,el_dim,0.0);
    unsigned n_intpt = integral_pt()->nweight();
    for (unsigned ipt=0;ipt<n_intpt;ipt++)
     {
      for(unsigned i=0;i<el_dim;i++) s[i] = integral_pt()->knot(ipt,i);
      this->strain_rate(s,strain_rate);
      for (unsigned i=0;i<el_dim;i++)
       {
        for (unsigned j=0;j<el_dim;j++)
         {
          Fixed_point_iteration_guess_for_strain_rate[Aitken_index][ipt](i,j)
           =strain_rate(i,j);
         }
       }
     }

    if(Aitken_index == 2)
     {
      Aitken_index = 0;
      for (unsigned ipt=0;ipt<n_intpt;ipt++)
       {
        for (unsigned i=0;i<el_dim;i++)
         {
          for (unsigned j=0;j<el_dim;j++)
           {
            double v0 = 
             Fixed_point_iteration_guess_for_strain_rate[0][ipt](i,j);
            double v1 =
             Fixed_point_iteration_guess_for_strain_rate[1][ipt](i,j);
            double v2 =
             Fixed_point_iteration_guess_for_strain_rate[2][ipt](i,j);

            double new_value=v2;

            if(Use_aitken_extrapolation)
             {
              double max_diff=std::max(std::fabs(v1-v0),std::fabs(v2-v1));
              
              if(max_diff > 1.0e-16)
               {
                new_value=v2-std::pow(v2-v1,2.0)/(v2-2.0*v1+v0);
               }
             }

            Fixed_point_iteration_guess_for_strain_rate[Aitken_index][ipt](i,j)=
             new_value;

           }
         }
       }
     }

    Aitken_index++;
   }


  /// Get strain-rate tensor: \f$ e_{ij} \f$  where 
  /// \f$ i,j = r,z,\theta \f$ (in that order). Extrapolated
  /// from history values evaluated at integration point ipt. Overloaded
  /// version from base class.
  void extrapolated_strain_rate(const unsigned& ipt,
                                DenseMatrix<double>& strainrate) const
   {
    if (Use_fixed_point_for_strain_rate)
     {
      latest_fixed_point_iteration_guess_for_strain_rate(ipt,strainrate);
     }
    else
     {
      Vector<double> s(el_dim);
      for(unsigned i=0;i<el_dim;i++) s[i] = integral_pt()->knot(ipt,i);
      extrapolated_strain_rate(s,strainrate);
     }
   }


  /// Get strain-rate tensor: \f$ e_{ij} \f$  where 
  /// \f$ i,j = r,z,\theta \f$ (in that order). Extrapolated
  /// from history values evaluated at local coordinate s. Overloaded
  /// version from base class.
  void extrapolated_strain_rate(const Vector<double>& s, 
                                DenseMatrix<double>& strainrate) const
   {

#ifdef PARANOID
    if ((strainrate.ncol()!=2)||(strainrate.nrow()!=2))
     {
      std::ostringstream error_message;
      error_message  << "The strain rate has incorrect dimensions " 
                     << strainrate.ncol() << " x " 
                     << strainrate.nrow() << " Not 2" << std::endl;
      
      throw OomphLibError(error_message.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }
#endif
    
    
    // Get required previous strain rates
    Vector<DenseMatrix<double> > 
     previous_strain_rate(Nprev_for_extrapolation_of_strain_rate);
    for (unsigned it=0;it<Nprev_for_extrapolation_of_strain_rate;it++)
     {
      previous_strain_rate[it].resize(el_dim,el_dim);
      strain_rate(it+1,s,previous_strain_rate[it]);
     }
    
    // Get timestepper from first node
    TimeStepper* time_stepper_pt=node_pt(0)->time_stepper_pt();

// hierher #ifdef PARANOID
    if (time_stepper_pt->nprev_values()<
        Nprev_for_extrapolation_of_strain_rate)
     {
      oomph_info << "Won't work: " << time_stepper_pt->nprev_values() 
                 << " < " 
                 << Nprev_for_extrapolation_of_strain_rate << std::endl;
      abort();
     }
// hierher #endif
    
    
    // Which extrapolation are we doing?
    switch(Nprev_for_extrapolation_of_strain_rate)
     {
      
      // Zero-th order extrapolation; use one previous value
     case 1:
     {
      strainrate=previous_strain_rate[0];
     }
     break;
     
     
     // First order extrapolation -- two history values
    case 2:
    {
     // Current and previous timesteps
     double dt=time_stepper_pt->time_pt()->dt(0);
     double dt_minus_1=time_stepper_pt->time_pt()->dt(1);
     
     // Extrapolate
     for (unsigned i=0;i<el_dim;i++)
      {
       for (unsigned j=0;j<el_dim;j++)
        {
         double u_minus_one=previous_strain_rate[0](i,j);
         double u_minus_two=previous_strain_rate[1](i,j);
         
         // Rate of changed based on previous two solutions
         double slope=(u_minus_one-u_minus_two)/dt_minus_1;
         
         // Extrapolated value from previous computed one to current one
         strainrate(i,j)=u_minus_one+slope*dt;
        }
      }
    }
    break;


     // Four history values
    case 4:
    {
     // Extrapolate
     for (unsigned i=0;i<el_dim;i++)
      {
       for (unsigned j=0;j<el_dim;j++)
        {
         double u_minus_one  =previous_strain_rate[0](i,j);
         double u_minus_two  =previous_strain_rate[1](i,j);
         double u_minus_three=previous_strain_rate[2](i,j);
         double u_minus_four =previous_strain_rate[3](i,j);
         
         // Current and previous timesteps
         double dt=time_stepper_pt->time_pt()->dt(0);
         double dt_minus_1=time_stepper_pt->time_pt()->dt(1);
         double dt_minus_2=time_stepper_pt->time_pt()->dt(2);
         double dt_minus_3=time_stepper_pt->time_pt()->dt(3);
         

         double MapleGenVar1=0.0;
         double MapleGenVar2=0.0;
         double MapleGenVar3=0.0;
         double MapleGenVar4=0.0;
         double MapleGenVar5=0.0;
         double MapleGenVar6=0.0;
         double MapleGenVar7=0.0;
         double t0=0.0;

           MapleGenVar2 = -1.0;
      MapleGenVar7 = dt*dt*dt*dt_minus_1*dt_minus_1*dt_minus_2*u_minus_four-dt*
dt*dt*dt_minus_1*dt_minus_1*dt_minus_2*u_minus_three-dt*dt*dt*dt_minus_1*
dt_minus_1*dt_minus_3*u_minus_three+dt*dt*dt*dt_minus_1*dt_minus_1*dt_minus_3*
u_minus_two+dt*dt*dt*dt_minus_1*dt_minus_2*dt_minus_2*u_minus_four-dt*dt*dt*
dt_minus_1*dt_minus_2*dt_minus_2*u_minus_three-2.0*dt*dt*dt*dt_minus_1*
dt_minus_2*dt_minus_3*u_minus_three+2.0*dt*dt*dt*dt_minus_1*dt_minus_2*
dt_minus_3*u_minus_two-dt*dt*dt*dt_minus_1*dt_minus_3*dt_minus_3*u_minus_three+
dt*dt*dt*dt_minus_1*dt_minus_3*dt_minus_3*u_minus_two;
      MapleGenVar6 = -dt*dt*dt*dt_minus_2*dt_minus_2*dt_minus_3*u_minus_one+dt*
dt*dt*dt_minus_2*dt_minus_2*dt_minus_3*u_minus_two-dt*dt*dt*dt_minus_2*
dt_minus_3*dt_minus_3*u_minus_one+dt*dt*dt*dt_minus_2*dt_minus_3*dt_minus_3*
u_minus_two+2.0*dt*dt*dt_minus_1*dt_minus_1*dt_minus_1*dt_minus_2*u_minus_four
-2.0*dt*dt*dt_minus_1*dt_minus_1*dt_minus_1*dt_minus_2*u_minus_three-2.0*dt*dt*
dt_minus_1*dt_minus_1*dt_minus_1*dt_minus_3*u_minus_three+2.0*dt*dt*dt_minus_1*
dt_minus_1*dt_minus_1*dt_minus_3*u_minus_two+3.0*dt*dt*dt_minus_1*dt_minus_1*
dt_minus_2*dt_minus_2*u_minus_four-3.0*dt*dt*dt_minus_1*dt_minus_1*dt_minus_2*
dt_minus_2*u_minus_three+MapleGenVar7;
      MapleGenVar7 = -6.0*dt*dt*dt_minus_1*dt_minus_1*dt_minus_2*dt_minus_3*
u_minus_three+6.0*dt*dt*dt_minus_1*dt_minus_1*dt_minus_2*dt_minus_3*u_minus_two
-3.0*dt*dt*dt_minus_1*dt_minus_1*dt_minus_3*dt_minus_3*u_minus_three+3.0*dt*dt*
dt_minus_1*dt_minus_1*dt_minus_3*dt_minus_3*u_minus_two+dt*dt*dt_minus_1*
dt_minus_2*dt_minus_2*dt_minus_2*u_minus_four-dt*dt*dt_minus_1*dt_minus_2*
dt_minus_2*dt_minus_2*u_minus_three-3.0*dt*dt*dt_minus_1*dt_minus_2*dt_minus_2*
dt_minus_3*u_minus_one-3.0*dt*dt*dt_minus_1*dt_minus_2*dt_minus_2*dt_minus_3*
u_minus_three+6.0*dt*dt*dt_minus_1*dt_minus_2*dt_minus_2*dt_minus_3*u_minus_two
-3.0*dt*dt*dt_minus_1*dt_minus_2*dt_minus_3*dt_minus_3*u_minus_one;
      MapleGenVar5 = -3.0*dt*dt*dt_minus_1*dt_minus_2*dt_minus_3*dt_minus_3*
u_minus_three+6.0*dt*dt*dt_minus_1*dt_minus_2*dt_minus_3*dt_minus_3*u_minus_two
-dt*dt*dt_minus_1*dt_minus_3*dt_minus_3*dt_minus_3*u_minus_three+dt*dt*
dt_minus_1*dt_minus_3*dt_minus_3*dt_minus_3*u_minus_two-2.0*dt*dt*dt_minus_2*
dt_minus_2*dt_minus_2*dt_minus_3*u_minus_one+2.0*dt*dt*dt_minus_2*dt_minus_2*
dt_minus_2*dt_minus_3*u_minus_two-3.0*dt*dt*dt_minus_2*dt_minus_2*dt_minus_3*
dt_minus_3*u_minus_one+3.0*dt*dt*dt_minus_2*dt_minus_2*dt_minus_3*dt_minus_3*
u_minus_two-dt*dt*dt_minus_2*dt_minus_3*dt_minus_3*dt_minus_3*u_minus_one+dt*dt
*dt_minus_2*dt_minus_3*dt_minus_3*dt_minus_3*u_minus_two+MapleGenVar6+
MapleGenVar7;
      MapleGenVar7 = dt*dt_minus_1*dt_minus_1*dt_minus_1*dt_minus_1*dt_minus_2*
u_minus_four-dt*dt_minus_1*dt_minus_1*dt_minus_1*dt_minus_1*dt_minus_2*
u_minus_three-dt*dt_minus_1*dt_minus_1*dt_minus_1*dt_minus_1*dt_minus_3*
u_minus_three+dt*dt_minus_1*dt_minus_1*dt_minus_1*dt_minus_1*dt_minus_3*
u_minus_two+2.0*dt*dt_minus_1*dt_minus_1*dt_minus_1*dt_minus_2*dt_minus_2*
u_minus_four-2.0*dt*dt_minus_1*dt_minus_1*dt_minus_1*dt_minus_2*dt_minus_2*
u_minus_three-4.0*dt*dt_minus_1*dt_minus_1*dt_minus_1*dt_minus_2*dt_minus_3*
u_minus_three+4.0*dt*dt_minus_1*dt_minus_1*dt_minus_1*dt_minus_2*dt_minus_3*
u_minus_two-2.0*dt*dt_minus_1*dt_minus_1*dt_minus_1*dt_minus_3*dt_minus_3*
u_minus_three+2.0*dt*dt_minus_1*dt_minus_1*dt_minus_1*dt_minus_3*dt_minus_3*
u_minus_two;
      MapleGenVar6 = dt*dt_minus_1*dt_minus_1*dt_minus_2*dt_minus_2*dt_minus_2*
u_minus_four-dt*dt_minus_1*dt_minus_1*dt_minus_2*dt_minus_2*dt_minus_2*
u_minus_three-3.0*dt*dt_minus_1*dt_minus_1*dt_minus_2*dt_minus_2*dt_minus_3*
u_minus_one-3.0*dt*dt_minus_1*dt_minus_1*dt_minus_2*dt_minus_2*dt_minus_3*
u_minus_three+6.0*dt*dt_minus_1*dt_minus_1*dt_minus_2*dt_minus_2*dt_minus_3*
u_minus_two-3.0*dt*dt_minus_1*dt_minus_1*dt_minus_2*dt_minus_3*dt_minus_3*
u_minus_one-3.0*dt*dt_minus_1*dt_minus_1*dt_minus_2*dt_minus_3*dt_minus_3*
u_minus_three+6.0*dt*dt_minus_1*dt_minus_1*dt_minus_2*dt_minus_3*dt_minus_3*
u_minus_two-dt*dt_minus_1*dt_minus_1*dt_minus_3*dt_minus_3*dt_minus_3*
u_minus_three+dt*dt_minus_1*dt_minus_1*dt_minus_3*dt_minus_3*dt_minus_3*
u_minus_two+MapleGenVar7;
      MapleGenVar7 = -4.0*dt*dt_minus_1*dt_minus_2*dt_minus_2*dt_minus_2*
dt_minus_3*u_minus_one+4.0*dt*dt_minus_1*dt_minus_2*dt_minus_2*dt_minus_2*
dt_minus_3*u_minus_two-6.0*dt*dt_minus_1*dt_minus_2*dt_minus_2*dt_minus_3*
dt_minus_3*u_minus_one+6.0*dt*dt_minus_1*dt_minus_2*dt_minus_2*dt_minus_3*
dt_minus_3*u_minus_two-2.0*dt*dt_minus_1*dt_minus_2*dt_minus_3*dt_minus_3*
dt_minus_3*u_minus_one+2.0*dt*dt_minus_1*dt_minus_2*dt_minus_3*dt_minus_3*
dt_minus_3*u_minus_two-dt*dt_minus_2*dt_minus_2*dt_minus_2*dt_minus_2*
dt_minus_3*u_minus_one+dt*dt_minus_2*dt_minus_2*dt_minus_2*dt_minus_2*
dt_minus_3*u_minus_two-2.0*dt*dt_minus_2*dt_minus_2*dt_minus_2*dt_minus_3*
dt_minus_3*u_minus_one+2.0*dt*dt_minus_2*dt_minus_2*dt_minus_2*dt_minus_3*
dt_minus_3*u_minus_two-dt*dt_minus_2*dt_minus_2*dt_minus_3*dt_minus_3*
dt_minus_3*u_minus_one;
      MapleGenVar4 = dt*dt_minus_2*dt_minus_2*dt_minus_3*dt_minus_3*dt_minus_3*
u_minus_two-dt_minus_1*dt_minus_1*dt_minus_1*dt_minus_2*dt_minus_2*dt_minus_3*
u_minus_one-dt_minus_1*dt_minus_1*dt_minus_1*dt_minus_2*dt_minus_3*dt_minus_3*
u_minus_one-2.0*dt_minus_1*dt_minus_1*dt_minus_2*dt_minus_2*dt_minus_2*
dt_minus_3*u_minus_one-3.0*dt_minus_1*dt_minus_1*dt_minus_2*dt_minus_2*
dt_minus_3*dt_minus_3*u_minus_one-dt_minus_1*dt_minus_1*dt_minus_2*dt_minus_3*
dt_minus_3*dt_minus_3*u_minus_one-dt_minus_1*dt_minus_2*dt_minus_2*dt_minus_2*
dt_minus_2*dt_minus_3*u_minus_one-2.0*dt_minus_1*dt_minus_2*dt_minus_2*
dt_minus_2*dt_minus_3*dt_minus_3*u_minus_one-dt_minus_1*dt_minus_2*dt_minus_2*
dt_minus_3*dt_minus_3*dt_minus_3*u_minus_one+MapleGenVar5+MapleGenVar6+
MapleGenVar7;
      MapleGenVar5 = 1/dt_minus_1;
      MapleGenVar3 = MapleGenVar4*MapleGenVar5;
      MapleGenVar1 = MapleGenVar2*MapleGenVar3;
      MapleGenVar2 = 1/dt_minus_2/dt_minus_3/(dt_minus_1*dt_minus_1*dt_minus_2+
dt_minus_1*dt_minus_1*dt_minus_3+2.0*dt_minus_1*dt_minus_2*dt_minus_2+3.0*
dt_minus_1*dt_minus_2*dt_minus_3+dt_minus_1*dt_minus_3*dt_minus_3+dt_minus_2*
dt_minus_2*dt_minus_2+2.0*dt_minus_2*dt_minus_2*dt_minus_3+dt_minus_2*
dt_minus_3*dt_minus_3);
      t0 = MapleGenVar1*MapleGenVar2;

      // Extrapolated value from previous computed ones to current one
      strainrate(i,j)=t0;
      
        }
      }
    }
    break;
    
     default:
     {
      oomph_info << "Never get here\n";
      abort();
     }
     break;
     
     }
   }
  
 

  /// Set error value for post-processing
  void set_error(const double& error){Error=error;}
  
  /// Return variable identifier
  std::string variable_identifier()
   {
    std::string txt="VARIABLES=";
    txt+="\"x\","; // 1
    txt+="\"y\","; // 2
    txt+="\"u\","; // 3
    txt+="\"v\","; // 4
    txt+="\"p\","; // 5
    txt+="\"du/dt\","; // 6
    txt+="\"dv/dt\","; // 7
    //txt+="\"u_m\","; // 8
    //txt+="\"v_m\","; // 9
    //txt+="\"x_h1\","; // 10
    //txt+="\"y_h1\","; // 11
    //txt+="\"x_h2\","; // 12
    //txt+="\"y_h2\","; // 13
    //txt+="\"u_h1\","; // 14
    //txt+="\"v_h1\","; // 15
    //txt+="\"u_h2\","; // 16
    //txt+="\"v_h2\","; // 17
    //txt+="\"strain_rate(0,0)\","; // 18
    //txt+="\"strain_rate(1,1)\","; // 19
    //txt+="\"strain_rate(0,1)\","; // 20
    //txt+="\"stress(0,0)\","; // 21
    //txt+="\"stress(1,1)\","; // 22
    //txt+="\"stress(0,1)\","; // 23
    txt+="\"invariant_of_strain_rate\","; // 24
    txt+="\"invariant_of_strain_rate_extrapol\","; // 25
    txt+="\"invariant_of_stress\","; // 26
    txt+="\"viscosity\","; // 27
    txt+="\"yield\","; // 28
    txt+="\"error\","; // 29  
    txt+="\"size\",";  // 30
    //txt+="\"strain_rate_m1(0,0)\","; //31
    //txt+="\"strain_rate_m1(1,1)\","; //32
    //txt+="\"strain_rate_m1(0,1)\","; //33
    //txt+="\"strain_rate_m2(0,0)\","; //v34
    //txt+="\"strain_rate_m2(1,1)\","; //v35
    //txt+="\"strain_rate_m2(0,1)\","; //v36
    //txt+="\"strain_rate_extrapol(0,0)\","; //v37
    //txt+="\"strain_rate_extrapol(1,1)\","; //v38
    //txt+="\"strain_rate_extrapol(0,1)\","; //v39
 
    txt+="\n";
    return txt;
   }

  
  /// Overload output function
  void output(std::ostream &outfile, 
              const unsigned &nplot)
   {
    
    // Vector of local coordinates
    Vector<double> s(el_dim);
    
    // Acceleration
    Vector<double> dudt(el_dim);
    
    // Mesh velocity
    Vector<double> mesh_veloc(el_dim,0.0);
   
    // Tecplot header info
    outfile << tecplot_zone_string(nplot);
   
    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Get continuous time from timestepper of first node
    //double time=node_pt(0)->time_stepper_pt()->time_pt()->time();
   
    //Set up memory for the shape functions
    Shape psif(n_node);
    DShape dpsifdx(n_node,el_dim);
   
    // Loop over plot points
    unsigned num_plot_points=nplot_points(nplot);
    for (unsigned iplot=0;iplot<num_plot_points;iplot++)
     {
     
      // Get local coordinates of plot point
      get_s_plot(iplot,nplot,s);
     
      //Call the derivatives of the shape and test functions
      dshape_eulerian(s,psif,dpsifdx);
     
      //Allocate storage
      Vector<double> mesh_veloc(el_dim,0.0);
      Vector<double> int_x(el_dim,0.0);
      Vector<double> dudt(el_dim,0.0);
      Vector<double> dudt_ALE(el_dim,0.0);
      DenseMatrix<double> interpolated_dudx(el_dim,el_dim,0.0);
     
      // //Initialise everything to zero
      // for(unsigned i=0;i<el_dim;i++)
      //  {
      //   mesh_veloc[i]=0.0;
      //   dudt[i]=0.0;
      //   dudt_ALE[i]=0.0;
      //   for(unsigned j=0;j<el_dim;j++)
      //    {
      //     interpolated_dudx(i,j) = 0.0;
      //    }
        
      //   int_x[i] = interpolated_x(s,i);
      //  }
     
      //Calculate velocities and derivatives

      //Loop over directions
      for(unsigned i=0;i<el_dim;i++)
       {
        //Get the index at which velocity i is stored
        unsigned u_nodal_index = u_index_nst(i);
        // Loop over nodes
        for(unsigned l=0;l<n_node;l++) 
         {
          dudt[i]+=du_dt_nst(l,u_nodal_index)*psif[l];
          mesh_veloc[i]+=dnodal_position_dt(l,i)*psif[l];
         }
       }
     
      //Loop over directions
      for(unsigned i=0;i<el_dim;i++)
       {
        //Get the index at which velocity i is stored
        unsigned u_nodal_index = u_index_nst(i);
        // Loop over nodes
        for(unsigned l=0;l<n_node;l++) 
         {
          //Loop over derivative directions for velocity gradients
          for(unsigned j=0;j<el_dim;j++)
           {                               
            interpolated_dudx(i,j) += nodal_value(l,u_nodal_index)*
             dpsifdx(l,j);
           }
         }
       }
     
      // Get dudt in ALE form (incl mesh veloc)
      for(unsigned i=0;i<el_dim;i++)
       {
        dudt_ALE[i]=dudt[i];
        for (unsigned k=0;k<el_dim;k++)
         {
          dudt_ALE[i]-=mesh_veloc[k]*interpolated_dudx(i,k);
         }
       }
     

      // Actual rate of strain
      DenseMatrix<double> rate_of_strain(el_dim,el_dim,0.0);
      this->strain_rate(s,rate_of_strain);

      // Extrapolated  (or recycle  actual one if there's no extrapol)
      DenseMatrix<double> rate_of_strain_extrapol(rate_of_strain);
      if (Use_extrapolated_strainrate_to_compute_second_invariant)
       {
        this->extrapolated_strain_rate(s,rate_of_strain_extrapol);
       }


      // Get associated invariants
      double second_invariant_strain=
       SecondInvariantHelper::second_invariant(rate_of_strain);
      double second_invariant_strain_extrapol=
       SecondInvariantHelper::second_invariant(rate_of_strain_extrapol);

      // Get viscosity hierher see above; compute this through the
      // equation/element (not directly via the constitutive equation)
      double viscosity=this->Constitutive_eqn_pt
       ->viscosity(second_invariant_strain_extrapol);
      
      // Stress
      DenseMatrix<double> stress(el_dim,el_dim,0.0);

      // Calculate the stress
      for(unsigned i=0;i<el_dim;i++)
       {
        for(unsigned j=0;j<el_dim;j++)
         {
          stress(i,j)=viscosity*rate_of_strain(i,j);
         }
       }

      double second_invariant_stress=
       SecondInvariantHelper::second_invariant(stress);    

      // Flag indicating if material has yielded (1.0) or not (0.0)
      double yield=0.0;

      // identify whether material is yielded (1) or not (0)
      if (sqrt(std::fabs(second_invariant_stress)) 
          >= Problem_Parameter::Yield_stress)
       {
        yield=1.0;
       }
      
      // Coordinates
      for(unsigned i=0;i<el_dim;i++) 
       {
        outfile << interpolated_x(s,i) << " ";
       }
     
      // Velocities
      for(unsigned i=0;i<el_dim;i++) 
       {
        outfile << interpolated_u_nst(s,i) << " ";
       }

      // Pressure
      outfile << interpolated_p_nst(s)  << " ";

      // Accelerations
      for(unsigned i=0;i<el_dim;i++) 
       {
        outfile << dudt_ALE[i] << " ";
       }
     
      // Mesh velocity
      //for(unsigned i=0;i<el_dim;i++) 
      // {
      //  outfile << mesh_veloc[i] << " ";
      // }
     
      // History values of coordinates
      //unsigned n_prev=node_pt(0)->position_time_stepper_pt()->ntstorage();
      //for (unsigned t=1;t<n_prev-2;t++)
      // {
      //  for(unsigned i=0;i<el_dim;i++) 
      //   {
      //    outfile << interpolated_x(t,s,i) << " ";
      //   }
      // }
     
      // History values of velocities
      //n_prev=node_pt(0)->time_stepper_pt()->ntstorage();
      //for (unsigned t=1;t<n_prev-2;t++)
      // {
      //  for(unsigned i=0;i<el_dim;i++) 
      //   {
      //    outfile << interpolated_u_nst(t,s,i) << " ";
      //   }
      // }

      //outfile << rate_of_strain(0,0) << " ";

      //outfile << rate_of_strain(1,1) << " ";

      //outfile << rate_of_strain(0,1) << " ";

      //outfile << stress(0,0) << " ";

      //outfile << stress(1,1) << " ";

      //outfile << stress(0,1) << " ";

      outfile << second_invariant_strain << " ";
      outfile << second_invariant_strain_extrapol << " ";
      outfile << second_invariant_stress << " ";

      outfile << viscosity << " ";

      outfile << yield << " ";

      outfile << Error << " ";

      outfile << size() << " ";

      // Rate of strain at previous timestep
      //this->strain_rate(1,s,rate_of_strain);
      //outfile << rate_of_strain(0,0) << " ";
      //outfile << rate_of_strain(1,1) << " ";
      //outfile << rate_of_strain(0,1) << " ";


      // Rate of strain at second previous timestep
      //this->strain_rate(2,s,rate_of_strain);
      //outfile << rate_of_strain(0,0) << " ";
      //outfile << rate_of_strain(1,1) << " ";
      //outfile << rate_of_strain(0,1) << " ";

      // Extrapolated rate of strain 
      //outfile << rate_of_strain_extrapol(0,0) << " ";
      //outfile << rate_of_strain_extrapol(1,1) << " ";
      //outfile << rate_of_strain_extrapol(0,1) << " ";

      outfile << std::endl;        

     }
    
    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile,nplot); 
    }




  /// Get 'flux' for Z2 error recovery
  void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
   {
#ifdef PARANOID
    unsigned num_entries=el_dim+(el_dim*(el_dim-1))/2;
    if (flux.size() < num_entries)
     {
      std::ostringstream error_message;
      error_message << "The flux vector has the wrong number of entries, " 
                    << flux.size() << ", whereas it should be at least " 
                    << num_entries << std::endl;
      throw OomphLibError(error_message.str(),
                          "RefineableNavierStokesEquations::get_Z2_flux()",
                          OOMPH_EXCEPTION_LOCATION);
     }
#endif
    
    // Get strain rate matrix
    DenseMatrix<double> strainrate(el_dim);
    DenseMatrix<double> stress(el_dim,el_dim,0.0);

    // get strainrate
    this->strain_rate(s,strainrate);
    double second_invariant_strain=
     SecondInvariantHelper::second_invariant(strainrate);

    // get viscosity
    double viscosity=this->Constitutive_eqn_pt
     ->viscosity(second_invariant_strain);

    // Calculate the stress
    for(unsigned i=0;i<el_dim;i++)
     {
      for(unsigned j=0;j<el_dim;j++)
       {
        stress(i,j)=viscosity*strainrate(i,j);
       }
     }

    // get the second invariant of the stress tensor
    //double second_invariant_stress=
    // SecondInvariantHelper::second_invariant(stress);

    // Pack into flux Vector
    unsigned icount=0;
    
    // identify whether material is yielded (0) or not (1)
    double yield=0.0;
    //if(sqrt(std::fabs(second_invariant_stress)) 
    //   >= Problem_Parameter::Yield_stress)
    // obacht
    if(std::fabs(second_invariant_strain) > 
       Problem_Parameter::Critical_strain_rate)
     {
      yield=1.0;
     }

    
    // Add bias to create lots of refinement near yield surface
    // by creating a jump in the z2 flux (Note: only makes sense if 
    // refinement happens at every timestep, otherwise really fine 
    // region gets immediately left behind by rapidly moving
    // yield surface
    bool add_yield_bias=false;
    if (!add_yield_bias) yield=0.0;

    // Start with diagonal terms
    for(unsigned i=0;i<el_dim;i++)
     {
      // Force more refinement at yield surface by adding yield value
      // hierher: This does the right thing but triangle can't handle the refinement
      // properly. Mesh seems ok though. MH will investigate triangle mesh adaptation
      // separately.
      flux[icount]=strainrate(i,i)+yield; 
      icount++;
     }
    
    //Off diagonals row by row
    for(unsigned i=0;i<el_dim;i++)
     {
      for(unsigned j=i+1;j<el_dim;j++)
       {
        // force more refinement at yield surface by adding yield value
        flux[icount]=strainrate(i,j)+yield;
        icount++;
       }
     }
   }
  


  /// Get square of L2 norms of (i) strain invariant, (ii) its 
  /// extrapolated value, (iii) difference between the two. Returns area
  /// as a check
  double square_of_norm_of_strain_invariant(double& norm_squared,
                                            double& extrapolated_norm_squared,
                                            double& error_norm_squared)
   {
    // Initialise
    norm_squared=0.0;
    extrapolated_norm_squared=0.0;
    error_norm_squared=0.0;        
    double area=0.0;

    //Number of integration points
    unsigned n_intpt = integral_pt()->nweight();
    
    //Set the Vector to hold local coordinates
    Vector<double> s(el_dim);

    //Loop over the integration points
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
     {
      //Assign values of s
      for(unsigned i=0;i<el_dim;i++) s[i] = integral_pt()->knot(ipt,i);
      
      //Get the integral weight
      double w = integral_pt()->weight(ipt);
      
      // Get Jacobain of mapping between local and Eulerian coords
      double J = this->J_eulerian(s);
      
      //Premultiply the weights and the Jacobian
      double W = w*J;
      
      // Invariant of actual strain rate tensor
      DenseMatrix<double> strainrate(el_dim,el_dim,0.0);
      this->strain_rate(s,strainrate);
      double invariant=SecondInvariantHelper::second_invariant(strainrate);

      // Invariant of extrapolated strain rate tensor
      // (or recycle  actual one if there's no extrapol)
      DenseMatrix<double> extrapolated_strainrate(strainrate);
      if (Use_extrapolated_strainrate_to_compute_second_invariant)
       {
        this->extrapolated_strain_rate(s,extrapolated_strainrate);
       }
      double extrapolated_invariant=
       SecondInvariantHelper::second_invariant(extrapolated_strainrate);

      //Assemble norms
      norm_squared+=
       invariant*invariant*W;
      extrapolated_norm_squared+=
       extrapolated_invariant*extrapolated_invariant*W;
      error_norm_squared+=
       (extrapolated_invariant-invariant)*
       (extrapolated_invariant-invariant)*W;
      area+=W;
     }

    return area;

   }

  /// Get square of L2 norms of (i) viscosity, (ii) its 
  /// extrapolated value, (iii) difference between the two. Returns area
  /// as a check
  double square_of_norm_of_viscosity(double& norm_squared,
                                     double& extrapolated_norm_squared,
                                     double& error_norm_squared)
   {
    // Initialise
    norm_squared=0.0;
    extrapolated_norm_squared=0.0;
    error_norm_squared=0.0;        
    double area=0.0;

    //Number of integration points
    unsigned n_intpt = integral_pt()->nweight();
    
    //Set the Vector to hold local coordinates
    Vector<double> s(el_dim);

    //Loop over the integration points
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
     {
      //Assign values of s
      for(unsigned i=0;i<el_dim;i++) s[i] = integral_pt()->knot(ipt,i);
      
      //Get the integral weight
      double w = integral_pt()->weight(ipt);
      
      // Get Jacobain of mapping between local and Eulerian coords
      double J = this->J_eulerian(s);
      
      //Premultiply the weights and the Jacobian
      double W = w*J;
      
      // Invariant of actual strain rate tensor
      DenseMatrix<double> strainrate(el_dim,el_dim,0.0);
      this->strain_rate(s,strainrate);
      double invariant=SecondInvariantHelper::second_invariant(strainrate);

      double viscosity=this->Constitutive_eqn_pt->viscosity(invariant);

      // Invariant of extrapolated strain rate tensor
      // (or recycle  actual one if there's no extrapol)
      DenseMatrix<double> extrapolated_strainrate(strainrate);
      if (Use_extrapolated_strainrate_to_compute_second_invariant)
       {
        this->extrapolated_strain_rate(s,extrapolated_strainrate);
       }
      double extrapolated_invariant=
       SecondInvariantHelper::second_invariant(extrapolated_strainrate);

      double extrapolated_viscosity=
       this->Constitutive_eqn_pt->viscosity(extrapolated_invariant);

      //Assemble norms
      norm_squared+=
       viscosity*viscosity*W;
      extrapolated_norm_squared+=
       extrapolated_viscosity*extrapolated_viscosity*W;
      error_norm_squared+=
       (extrapolated_viscosity-viscosity)*
       (extrapolated_viscosity-viscosity)*W;
      area+=W;
     }

    return area;

   }

  /// Get square of L2 norms of (i) current strainrate, (ii) its 
  /// latest guess from fixed point iteration, (iii) difference 
  /// between the two. Returns area as a check
  double square_of_norm_of_fixed_point(double& norm_squared,
                                       double& latest_guess_norm_squared,
                                       double& error_norm_squared)
   {
    // Initialise
    norm_squared=0.0;
    latest_guess_norm_squared=0.0;
    error_norm_squared=0.0;        
    double area=0.0;

    //Number of integration points
    unsigned n_intpt = integral_pt()->nweight();
    
    //Set the Vector to hold local coordinates
    Vector<double> s(2);

    //Loop over the integration points
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
     {
      //Assign values of s
      for(unsigned i=0;i<2;i++) s[i] = integral_pt()->knot(ipt,i);
      
      //Get the integral weight
      double w = integral_pt()->weight(ipt);
      
      // Get Jacobain of mapping between local and Eulerian coords
      double J = this->J_eulerian(s);
      
      //Premultiply the weights and the Jacobian
      double W = w*J;
      
      // Invariant of actual strain rate tensor
      DenseMatrix<double> current_strainrate(el_dim,el_dim,0.0);
      this->strain_rate(s,current_strainrate);
      double invariant=
       SecondInvariantHelper::second_invariant(current_strainrate);

      // Invariant of extrapolated strain rate tensor
      // (or recycle  actual one if there's no extrapol)
      DenseMatrix<double> latest_guess_strainrate(el_dim,el_dim,0.0);
      latest_fixed_point_iteration_guess_for_strain_rate(
       ipt, latest_guess_strainrate);
      
      double latest_guess_invariant=
       SecondInvariantHelper::second_invariant(latest_guess_strainrate);

      //Assemble norms
      norm_squared+=
       invariant*invariant*W;
      latest_guess_norm_squared+=
       latest_guess_invariant*latest_guess_invariant*W;
      error_norm_squared+=
       (latest_guess_invariant-invariant)*
       (latest_guess_invariant-invariant)*W;
      area+=W;
     }

    return area;

   }

  /// Get square of L2 norm of velocity
  double square_of_l2_norm()
   {

    // Initalise
    double sum=0.0;
    
    //Find out how many nodes there are
    unsigned n_node = nnode();
    
    //Find the indices at which the local velocities are stored
    unsigned u_nodal_index[el_dim+1];
    for(unsigned i=0;i<el_dim;i++) {u_nodal_index[i] = u_index_nst(i);}
    
    //Set up memory for the velocity shape fcts
    Shape psif(n_node);
    DShape dpsidx(n_node,el_dim);
    
    //Number of integration points
    unsigned n_intpt = integral_pt()->nweight();
    
    //Set the Vector to hold local coordinates
    Vector<double> s(el_dim);
    
    //Loop over the integration points
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
     {
      //Assign values of s
      for(unsigned i=0;i<el_dim;i++) s[i] = integral_pt()->knot(ipt,i);
      
      //Get the integral weight
      double w = integral_pt()->weight(ipt);
      
      // Call the derivatives of the veloc shape functions
      // (Derivs not needed but they are free)
      double J = this->dshape_eulerian_at_knot(ipt,psif,dpsidx);
      
      //Premultiply the weights and the Jacobian
      double W = w*J;
      
      //Calculate velocities
      Vector<double> interpolated_u(el_dim,0.0);
      
      // Loop over nodes
      for(unsigned l=0;l<n_node;l++)
       {
        //Loop over directions
        for(unsigned i=0;i<el_dim;i++)
         {
          //Get the nodal value
          double u_value = raw_nodal_value(l,u_nodal_index[i]);
          interpolated_u[i] += u_value*psif[l];
         }
       }

      //Assemble square of L2 norm
      for(unsigned i=0;i<el_dim;i++)
       {
        sum+=interpolated_u[i]*interpolated_u[i]*W;
       }
     }

    return sum;

   }


  /* /// Output solution in data vector at local cordinates s: */
  /* /// r,z,u_r,u_z,u_phi,p */
  /* void point_output_data(const Vector<double> &s, Vector<double>& data) */
  /*  { */

  /*   // get the time */
  /*   double time=this->node_pt(0)->time_stepper_pt()->time(); */
    
  /*   /// get the viscosity ratio */
  /*   // hierher not used double visc_ratio=this->viscosity_ratio(); */

  /*   DenseMatrix<double> strainrate(el_dim,el_dim,0.0); */
  /*   this->strain_rate(s,strainrate); */

  /*   DenseMatrix<double> extrapolated_strainrate(el_dim,el_dim,0.0); */
    
  /*   double second_invariant=0.0; */
    
  /*   if(!Use_extrapolated_strainrate_to_compute_second_invariant) */
  /*    { */
  /*     second_invariant=SecondInvariantHelper::second_invariant(strainrate); */
  /*    } */
  /*   else */
  /*    { */
  /*     this->extrapolated_strain_rate(s,extrapolated_strainrate); */
  /*     second_invariant=SecondInvariantHelper::second_invariant( */
  /*      extrapolated_strainrate); */
  /*    } */
    
  /*   // Get the viscosity according to the constitutive equation */
  /*   double viscosity=Constitutive_eqn_pt->viscosity(second_invariant); */

  /*   // Output the components of the position */
  /*   for(unsigned i=0;i<el_dim;i++) */
  /*    { */
  /*     data.push_back(interpolated_x(s,i)); // 1 2 */
  /*    } */
    
  /*   // Output the components of the FE representation of u at s */
  /*   for(unsigned i=0;i<el_dim;i++) */
  /*    { */
  /*     data.push_back(interpolated_u_nst(s,i)); // 3 4 */
  /*    } */
    
  /*   double x=interpolated_x(s,0); */

  /*   data.push_back(Problem_Parameter::exact_soln(time,x)); //5 */

  /*   // Output FE representation of p at s */
  /*   data.push_back(interpolated_p_nst(s)); // 6 */

  /*   data.push_back(strainrate(0,0)); // 7 */
    
  /*   data.push_back(strainrate(1,1)); // 8 */
    
  /*   data.push_back(strainrate(0,1)); // 9 */
   
  /*   data.push_back(second_invariant); // 10 */
    
  /*   data.push_back(viscosity); // 11 */
    
  /*  } */


 private:

  /// Current best guess for strain rate tensor (fixed point iteration)
  Vector<Vector<DenseMatrix<double> > > Fixed_point_iteration_guess_for_strain_rate;
 /// unsigned storing the number of fixed point iterations after the last
 /// Aitken extrapolation
 unsigned Aitken_index;


// hierher moved additional member functions out I don't need them.
// If required they can be re-instated (but only after merging them
// properly -- there's so much code duplication!
//#include "leftover_bits.h"
