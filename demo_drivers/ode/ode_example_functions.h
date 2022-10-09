//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
#ifndef OOMPH_ODE_EXAMPLE_FUNCTIONS_H
#define OOMPH_ODE_EXAMPLE_FUNCTIONS_H

#include "generic.h"
#include "ode.h"

namespace oomph
{
 using namespace MathematicalConstants;
 using namespace StringConversion;


 namespace deriv_functions
 {

  inline Vector<double> cos(const double& time, const Vector<double>& x)
  {
   Vector<double> values(1);
   values[0] = std::cos(time);
   return values;
  }
  inline Vector<double> dcos(const double& t, const Vector<double>& x,
                             const Vector<double>& u)
  {
   Vector<double> deriv(1, 0.0);
   deriv[0] = -1*std::sin(t);
   return deriv;
  }


  inline Vector<double> sin(const double& time, const Vector<double>& x)
  {
   Vector<double> values(1);
   values[0] = std::sin(time);
   return values;
  }
  inline Vector<double> dsin(const double& t, const Vector<double>& x,
                             const Vector<double>& u)
  {
   Vector<double> deriv(1, 0.0);
   deriv[0] = std::cos(t);
   return deriv;
  }


  inline Vector<double> exp(const double& time, const Vector<double>& x)
  {
   Vector<double> values(1);
   values[0] = std::exp(time);
   return values;
  }
  inline Vector<double> dexp(const double& t, const Vector<double>& x,
                             const Vector<double>& u)
  {
   Vector<double> deriv(1, 0.0);
   deriv[0] = u[0];
   return deriv;
  }


  // A polynomial of degree 2
  inline Vector<double> poly2(const double& time, const Vector<double>& x)
  {
   double b0 = 0.5, b1 = 0, b2 = 1;

   Vector<double> values(1);
   values[0] =  b2*time*time + b1*time +b0;
   return values;
  }
  inline Vector<double> dpoly2(const double& t, const Vector<double>& x,
                               const Vector<double>& u)
  {
   double b1 = 0, b2 = 1;

   Vector<double> deriv(1, 0.0);
   deriv[0] = 2*t*b2 + b1;
   return deriv;
  }


  // A polynomial of degree 3
  inline Vector<double> poly3(const double& time, const Vector<double>& x)
  {
   double a0 = 0.5, a1 = 0, a2 = 0, a3 = 1;

   Vector<double> values(1);
   values[0] = a3*time*time*time + a2*time*time + a1*time +a0;
   return values;
  }
  inline Vector<double> dpoly3(const double& t, const Vector<double>& x,
                               const Vector<double>& u)
  {
   double a1 = 0, a2 = 0, a3 = 1;

   Vector<double> deriv(1, 0.0);
   deriv[0] = 3*t*t*a3 + 2*t*a2 + a1;
   return deriv;
  }

  // stiff ode, example from Iserles pg. 54
  inline Vector<double> stiff_test(const double& time, const Vector<double>& x)
  {
   Vector<double> x1(2), x2(2);
   x1[0] = 0; x1[1] = 0;
   x2[0] = 1.0; x2[1] = 999.0/10;

   Vector<double> values(2);
   values[0] = std::exp(-100*time)*x1[0] + std::exp(-0.1*time)*x2[0];
   values[1] = std::exp(-100*time)*x1[1] + std::exp(-0.1*time)*x2[1];
   return values;
  }
  inline Vector<double> dstiff_test(const double& t, const Vector<double>& x,
                                    const Vector<double>& u)
  {
   Vector<double> deriv(2, 0.0);
   deriv[0] = -100*u[0] + u[1];
   deriv[1] = 0*u[0] - 0.1*u[1];
   return deriv;
  }

  /// A damped oscillation solution
  class DampedOscillation : public SolutionFunctorBase
  {
  public:
   /// Constructor
   DampedOscillation()
   {
    Beta = 0.5;
    Omega = 2;
   }

   /// Virtual destructor
   virtual ~DampedOscillation() {}

   /// Function call
   Vector<double> operator()(const double& t, const Vector<double>& x) const
   {
    Vector<double> values(1);
    values[0] = std::exp(-Beta*t) * std::sin(Omega*t);
    return values;
   }

   /// Derivative call
   Vector<double> derivative(const double& t, const Vector<double>& x,
                             const Vector<double>& u) const
   {
    Vector<double> deriv(1, 0.0);
    deriv[0] = -Beta * std::exp(-Beta*t) * std::sin(Omega*t)
     + Omega * std::exp(-Beta*t) * std::cos(Omega*t);
    return deriv;
   }

   double Beta;
   double Omega;
  };

  /// Another stiff solution: Atkinson equation (8.1) pg 128
  class SimpleStiffTest : public SolutionFunctorBase
  {
  public:
   /// Constructor
   SimpleStiffTest()
   {
    Lambda = 100;
    Y_intial = 1;
   }

   /// Virtual destructor
   virtual ~SimpleStiffTest() {}

   /// Function call
   Vector<double> operator()(const double& t, const Vector<double>& x) const
   {
    Vector<double> values(1);
    values[0] = std::exp(-Lambda * t) * Y_intial;
    return values;
   }

   /// Derivative call
   Vector<double> derivative(const double& t, const Vector<double>& x,
                             const Vector<double>& u) const
   {
    Vector<double> deriv(1, 0.0);
    deriv[0] = -Lambda * u[0];
    return deriv;
   }

   double Lambda;
   double Y_intial;
  };

  /// Another stiff solution: Atkinson pg. 158, also example 8.2 pg 129.
  class OrderReductionTest : public SolutionFunctorBase
  {
  public:
   /// Constructor
   OrderReductionTest()
   {
    Lambda = -100;
   }

   /// Virtual destructor
   virtual ~OrderReductionTest() {}

   /// Function call
   Vector<double> operator()(const double& t, const Vector<double>& x) const
   {
    Vector<double> values(1);
    values[0] = std::sin(t);
    return values;
   }

   /// Derivative call
   Vector<double> derivative(const double& t, const Vector<double>& x,
                             const Vector<double>& u) const
   {
    Vector<double> deriv(1, 0.0);
    deriv[0] = Lambda*u[0] - Lambda*std::sin(t) + std::cos(t);
    return deriv;
   }

   double Lambda;
  };

 }

 namespace ODEFactories
 {

  SolutionFunctorBase* exact_solutions_factory(const std::string& exact_name)
  {
   if(exact_name == "damped_oscillation")
    {
     return new deriv_functions::DampedOscillation;
    }
   else if(exact_name == "simple_stiff")
    {
     return new deriv_functions::SimpleStiffTest;
    }
   else if(exact_name == "order_reduction")
    {
     return new deriv_functions::OrderReductionTest;
    }
   else if(exact_name == "strong_order_reduction")
    {
     return new deriv_functions::OrderReductionTest;
    }

   SolutionFunctorBase::TimeSpaceToDoubleVectFctPt fpt;
   SolutionFunctorBase::TimeSpaceValueToDoubleVectFctPt dfpt;

   if(exact_name == "sin")
    {
     fpt = &deriv_functions::sin;
     dfpt = &deriv_functions::dsin;
    }
   else if(exact_name == "cos")
    {
     fpt = &deriv_functions::cos;
     dfpt = &deriv_functions::dcos;
    }
   else if(exact_name == "exp")
    {
     fpt = &deriv_functions::exp;
     dfpt = &deriv_functions::dexp;
    }
   else if(exact_name == "poly3")
    {
     fpt = &deriv_functions::poly3;
     dfpt = &deriv_functions::dpoly3;
    }
   else if(exact_name == "stiff_test")
    {
     fpt = &deriv_functions::stiff_test;
     dfpt = &deriv_functions::dstiff_test;
    }
   else if(exact_name == "poly2")
    {
     fpt = &deriv_functions::poly2;
     dfpt = &deriv_functions::dpoly2;
    }
   else
    {
     throw OomphLibError("Unrecognised exact solution " + exact_name,
                         OOMPH_EXCEPTION_LOCATION,
                         OOMPH_CURRENT_FUNCTION);
    }

   return new SolutionFunctor(fpt, dfpt);
  }

 }


} // End of oomph namespace

#endif
