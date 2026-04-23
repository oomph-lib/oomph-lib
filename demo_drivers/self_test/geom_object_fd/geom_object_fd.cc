//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2026 Matthias Heil and Andrew Hazel
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
// Driver for adaptive 2D rectangular driven cavity. Solved with black
// box adaptation, using Taylor Hood and Crouzeix Raviart elements.

// Generic oomph-lib header
#include "generic.h"

// Geometric object with boundary header
#include "generic/geom_obj_with_boundary.h"

// Geometric objects header
#include "generic/geom_objects.h"

using namespace std;

using namespace oomph;


//========================================================================
/// Test GeomObject for the code of Finite difference method 
///
/// \f[ {\bf r} = 
/// (zeta_1+zeta_2+zeta_1^2+zeta_2^2+zeta_1 zeta_2+sin(zeta_1+zeta_2), 
/// zeta_1^2 zeta_2+zeta_1 zeta_2^2+sin(zeta_1 zeta_2)+cos(zeta_1+zeta_2),
/// zeta_1 zeta_2+zeta_1^2 zeta_2^2+sin(zeta_1)+cos(zeta_2)
/// +sin(zeta_1+zeta_2))^T \f]
///
//========================================================================
class TestFDMethodShape : public GeomObject
{
public:
  /// Constructor: 
  TestFDMethodShape(): GeomObject(2, 3){}

  /// Broken copy constructor
  TestFDMethodShape(const TestFDMethodShape &dummy) = delete;

  /// Broken assignment operator
  void operator=(const TestFDMethodShape &) = delete;

  /// Position vector
  void position(const Vector<double> &zeta, Vector<double> &r) const
  {
    r[0] = zeta[0] + zeta[1] + zeta[0] * zeta[0] + zeta[1] * zeta[1] 
           + zeta[0] * zeta[1] + sin(zeta[0] + zeta[1]);
    r[1] = zeta[0] * zeta[0] * zeta[1] + zeta[0] * zeta[1] * zeta[1] 
           + sin(zeta[0] * zeta[1]) + cos(zeta[0] + zeta[1]);
    r[2] = zeta[0] * zeta[1] + zeta[0] * zeta[0] * zeta[1] * zeta[1]
    + sin(zeta[0]) + cos(zeta[1]) + sin(zeta[0] + zeta[1]);
  }

  /// Position vector (dummy unsteady version returns steady version)
  void position(const unsigned &t,
                const Vector<double> &zeta,
                Vector<double> &r) const
  {
    position(zeta, r);
  }

  /// How many items of Data does the shape of the object depend on?
  unsigned ngeom_data() const
  {
    return 0;
  }

  /// Analytical 1st derivs w.r.t. zeta.
  void dposition(const Vector<double> &zeta,
                 DenseMatrix<double> &drdzeta) const
  {
    // dr_dzeta_1
    drdzeta(0, 0) = 0.1e1 + 0.2e1 * zeta[0] + zeta[1] + cos(zeta[0] + zeta[1]);
    drdzeta(0, 1) = 0.2e1 * zeta[0] * zeta[1] + zeta[1] * zeta[1] + 
                    zeta[1] * cos(zeta[0] * zeta[1]) - sin(zeta[0] + zeta[1]);
    drdzeta(0, 2) = zeta[1] + 0.2e1 * zeta[0] * zeta[1] * zeta[1] + cos(zeta[0]) 
                    + cos(zeta[0] + zeta[1]);

    // dr_dzeta_2
    drdzeta(1, 0) = 0.1e1 + 0.2e1 * zeta[1] + zeta[0] + cos(zeta[0] + zeta[1]);
    drdzeta(1, 1) = zeta[0] * zeta[0] + 0.2e1 * zeta[0] * zeta[1] + 
                    zeta[0] * cos(zeta[0] * zeta[1]) - sin(zeta[0] + zeta[1]);
    drdzeta(1, 2) = zeta[0] + 0.2e1 * zeta[0] * zeta[0] * zeta[1] - sin(zeta[1]) 
                    + cos(zeta[0] + zeta[1]);
  }

  /// Analytical 2nd derivs w.r.t. zeta.
  void d2position(const Vector<double> &zeta,
                  RankThreeTensor<double> &ddrdzeta) const
  {
    // d2r_dzeta_1zeta_1
    ddrdzeta(0, 0, 0) = 0.2e1 - sin(zeta[0] + zeta[1]);
    ddrdzeta(0, 0, 1) = 0.2e1 * zeta[1] - zeta[1] * zeta[1] * 
                        sin(zeta[0] * zeta[1]) - cos(zeta[0] + zeta[1]);
    ddrdzeta(0, 0, 2) = 0.2e1 * zeta[1] * zeta[1] - sin(zeta[0]) 
                        - sin(zeta[0] + zeta[1]);

    // d2r_dzeta_2zeta_2
    ddrdzeta(1, 1, 0) = 0.2e1 - sin(zeta[0] + zeta[1]);
    ddrdzeta(1, 1, 1) = 0.2e1 * zeta[0] - zeta[0] * zeta[0] * 
                        sin(zeta[0] * zeta[1]) - cos(zeta[0] + zeta[1]);
    ddrdzeta(1, 1, 2) = 0.2e1 * zeta[0] * zeta[0] - cos(zeta[1]) 
                        - sin(zeta[0] + zeta[1]);

    // d2r_dzeta_1zeta_2
    ddrdzeta(0, 1, 0) = 0.1e1 - sin(zeta[0] + zeta[1]);
    ddrdzeta(0, 1, 1) = 0.2e1 * zeta[0] + 0.2e1 * zeta[1] + 
                        cos(zeta[0] * zeta[1]) - zeta[1] * zeta[0] * 
                        sin(zeta[0] * zeta[1]) - cos(zeta[0] + zeta[1]);
    ddrdzeta(0, 1, 2) = 0.1e1 + 0.4e1 * zeta[0] * zeta[1] 
                        - sin(zeta[0] + zeta[1]);

    // d2r_dzeta_2zeta_1
    ddrdzeta(1, 0, 0) = 0.1e1 - sin(zeta[0] + zeta[1]);
    ddrdzeta(1, 0, 1) = 0.2e1 * zeta[0] + 0.2e1 * zeta[1] + 
                        cos(zeta[0] * zeta[1]) - zeta[1] * zeta[0] * 
                        sin(zeta[0] * zeta[1]) - cos(zeta[0] + zeta[1]);
    ddrdzeta(1, 0, 2) = 0.1e1 + 0.4e1 * zeta[0] * zeta[1] 
                        - sin(zeta[0] + zeta[1]);
  }

  /// Position Vector and analytical 1st and 2nd derivs w.r.t. zeta.
  void d2position(const Vector<double> &zeta,
                  Vector<double> &r,
                  DenseMatrix<double> &drdzeta,
                  RankThreeTensor<double> &ddrdzeta) const
  {
    // Get the position
    position(zeta, r);

    // Do the first derivative
    dposition(zeta, drdzeta);

    // Do the second derivative
    d2position(zeta, ddrdzeta);
  }

};


//==start_of_main======================================================
/// Driver for validation of FD method code
//=====================================================================
int main()
{

 // Set output directory
 DocInfo doc_info;
 doc_info.set_directory("RESLT");

 // Create the geometric object
 TestFDMethodShape* geom_object_pt = new TestFDMethodShape();

 // Initialize two variables
 Vector<double> zeta(2, 0.0);

 // Assign values for zeta
 zeta[0] = 0.1;
 zeta[1] = 0.3;

 // ---------------------------------
 // 1st derivative (dposition) check
 // ---------------------------------

 // Compute analytical 1st derivative
 DenseMatrix<double> drdzeta_exact(2, 3);
 geom_object_pt->dposition(zeta, drdzeta_exact);

 // Compute FD 1st derivative
 DenseMatrix<double> drdzeta_FD(2, 3);
 geom_object_pt->GeomObject::dposition(zeta, drdzeta_FD);

 // Get the number of rows
 unsigned nrow = drdzeta_exact.nrow();

 // Get the number of columns
 unsigned ncol = drdzeta_exact.ncol();

 // Open a file to record the analytical and FD results
 std::ofstream file1_analytical_and_FD("RESLT/first_derivs.dat");

 // Output the analytical and FD results
 for (unsigned j = 0; j < nrow; j++)
 {
   for (unsigned k = 0; k < ncol; k++)
   {
     file1_analytical_and_FD << drdzeta_exact(j, k) << " " << 
     drdzeta_FD(j, k) << std::endl;
   }
 }
 // Close the file
 file1_analytical_and_FD.close();

 // ------------------------------------------------------------------
 // This part of the code is used to assess the optimal value of
 // the FD step
 // ------------------------------------------------------------------
 // Open a file for plotting the difference between
 // the analytical and FD results
 std::ofstream file1("RESLT/first_derivs_difference.dat");

 // Loop over a range of finite difference step sizes
 unsigned npts = 1000;
 for (unsigned i = 0; i < npts; i++)
 {
   // Define logarithmically spaced step size
   double exponent = -16.0 + 16.0 * double(i) / double(npts - 1);
   double dzeta = std::pow(10.0, exponent);

   // Output step size
   file1 << dzeta << " ";

   // Set FD step for 1st derivative
   geom_object_pt->first_derivative_fd_step() = dzeta;

   // Compute FD 1st derivative
   DenseMatrix<double> drdzeta_FD(2, 3);
   geom_object_pt->GeomObject::dposition(zeta, drdzeta_FD);

   // Output the analytical and FD results
   for (unsigned j = 0; j < nrow; j++)
   {
     for (unsigned k = 0; k < ncol; k++)
     {
       file1 << drdzeta_FD(j, k) << " " << drdzeta_exact(j, k) << " ";
     }
   }
   file1 << std::endl;
 }
 // Close the file
 file1.close();
 // ------------------------------------------------------------------

 // -------------------------------
 // 2nd derivative (d2position)
 // -------------------------------

 // Compute analytical 2nd derivative
 RankThreeTensor<double> ddrdzeta_exact(2, 2, 3);
 geom_object_pt->d2position(zeta, ddrdzeta_exact);

 // Compute FD 2nd derivative
 RankThreeTensor<double> ddrdzeta_FD(2, 2, 3);
 geom_object_pt->GeomObject::d2position(zeta, ddrdzeta_FD);

 // Get the number of index1
 unsigned nindex1 = ddrdzeta_exact.nindex1();

 // Get the number of index2
 unsigned nindex2 = ddrdzeta_exact.nindex2();

 // Get the number of index3
 unsigned nindex3 = ddrdzeta_exact.nindex3();

 // Open a file to record the analytical and FD results
 std::ofstream file2_analytical_and_FD("RESLT/second_derivs.dat");

 // Output the analytical and FD results
 for (unsigned j = 0; j < nindex1; j++)
 {
   for (unsigned k = 0; k < nindex2; k++)
   {
     for (unsigned h = 0; h < nindex3; h++)
     {
       file2_analytical_and_FD << ddrdzeta_exact(j, k, h) << " " <<
       ddrdzeta_FD(j, k, h) << std::endl;
     }
   }
 }
 // Close the file
 file2_analytical_and_FD.close();

 // ------------------------------------------------------------------
 // This part of the code is used to assess the optimal value of
 // the FD step
 // ------------------------------------------------------------------
 // Open a file
 std::ofstream file2("RESLT/second_derivs_difference.dat");

 // Loop over a range of finite difference step sizes
 npts = 1000;
 for (unsigned i = 0; i < npts; i++)
 {
   // Define logarithmically spaced step size
   double exponent = -8.0 + 8.0 * double(i) / double(npts - 1);
   double dzeta = std::pow(10.0, exponent);

   // Output step size
   file2 << dzeta << " ";

   // Set FD step for 2nd derivative
   geom_object_pt->second_derivative_fd_step() = dzeta;

   // Compute FD 2nd derivative
   RankThreeTensor<double> ddrdzeta_FD(2, 2, 3);
   geom_object_pt->GeomObject::d2position(zeta, ddrdzeta_FD);

   // Output the analytical and FD results
   for (unsigned j = 0; j < nindex1; j++)
   {
     for (unsigned k = 0; k < nindex2; k++)
     {
       for (unsigned h = 0; h < nindex3; h++)
       {
         file2 << ddrdzeta_FD(j, k, h) << " " << ddrdzeta_exact(j, k, h) << " ";
       }
     }
   }
   file2 << std::endl;
 }
 // Close the file
 file2.close();
 // ------------------------------------------------------------------


} // end_of_main











