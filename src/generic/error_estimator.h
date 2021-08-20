// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================
#ifndef OOMPH_ERROR_ESTIMATOR_NAMESPACE_HEADER
#define OOMPH_ERROR_ESTIMATOR_NAMESPACE_HEADER

#include "mesh.h"
#include "quadtree.h"
#include "nodes.h"
#include "algebraic_elements.h"

namespace oomph
{
  //========================================================================
  /// Base class for spatial error estimators
  //========================================================================
  class ErrorEstimator
  {
  public:
    /// Default empty constructor
    ErrorEstimator() {}

    /// Broken copy constructor
    ErrorEstimator(const ErrorEstimator&)
    {
      BrokenCopy::broken_copy("ErrorEstimator");
    }

    /// Broken assignment operator
    void operator=(const ErrorEstimator&)
    {
      BrokenCopy::broken_assign("ErrorEstimator");
    }

    /// Empty virtual destructor
    virtual ~ErrorEstimator() {}

    /// \short Compute the elemental error-measures for a given mesh
    /// and store them in a vector.
    void get_element_errors(Mesh*& mesh_pt, Vector<double>& elemental_error)
    {
      // Create dummy doc info object and switch off output
      DocInfo doc_info;
      doc_info.disable_doc();
      // Forward call to version with doc.
      get_element_errors(mesh_pt, elemental_error, doc_info);
    }


    /// \short Compute the elemental error measures for a given mesh
    /// and store them in a vector. Doc errors etc.
    virtual void get_element_errors(Mesh*& mesh_pt,
                                    Vector<double>& elemental_error,
                                    DocInfo& doc_info) = 0;
  };


  //========================================================================
  /// \short Base class for finite elements that can compute the
  /// quantities that are required for the Z2 error estimator.
  //========================================================================
  class ElementWithZ2ErrorEstimator : public virtual FiniteElement
  {
  public:
    /// Default empty constructor
    ElementWithZ2ErrorEstimator() {}

    /// Broken copy constructor
    ElementWithZ2ErrorEstimator(const ElementWithZ2ErrorEstimator&)
    {
      BrokenCopy::broken_copy("ElementWithZ2ErrorEstimator");
    }

    /// Broken assignment operator
    void operator=(const ElementWithZ2ErrorEstimator&)
    {
      BrokenCopy::broken_assign("ElementWithZ2ErrorEstimator");
    }

    /// \short Number of 'flux' terms for Z2 error estimation
    virtual unsigned num_Z2_flux_terms() = 0;

    /// \short A stuitable error estimator for
    /// a multi-physics elements may require one Z2 error estimate for each
    /// field (e.g. velocity and temperature in a fluid convection problem).
    /// It is assumed that these error estimates will each use
    /// selected flux terms. The number of compound fluxes returns the number
    /// of such combinations of the flux terms. Default value is one and all
    /// flux terms are combined with equal weight.
    virtual unsigned ncompound_fluxes()
    {
      return 1;
    }

    /// \short Z2 'flux' terms for Z2 error estimation
    virtual void get_Z2_flux(const Vector<double>& s, Vector<double>& flux) = 0;

    /// \short Plot the error when compared against a given exact flux.
    /// Also calculates the norm of the error and that of the exact flux.
    virtual void compute_exact_Z2_error(
      std::ostream& outfile,
      FiniteElement::SteadyExactSolutionFctPt exact_flux_pt,
      double& error,
      double& norm)
    {
      std::string error_message =
        "compute_exact_Z2_error undefined for this element \n";
      outfile << error_message;

      throw OomphLibError(
        error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    /// \short Return the compound flux index of each flux component
    /// The default (do nothing behaviour) will mean that all indices
    /// remain at the default value zero.
    virtual void get_Z2_compound_flux_indices(Vector<unsigned>& flux_index) {}

    /// \short Number of vertex nodes in the element
    virtual unsigned nvertex_node() const = 0;

    /// \short Pointer to the j-th vertex node in the element. Needed for
    /// efficient patch assmbly
    virtual Node* vertex_node_pt(const unsigned& j) const = 0;

    /// \short Order of recovery shape functions
    virtual unsigned nrecovery_order() = 0;

    /// \short Return the geometric jacobian (should be overloaded in
    /// cylindrical and spherical geometries).
    /// Default value one is suitable for Cartesian coordinates
    virtual double geometric_jacobian(const Vector<double>& x)
    {
      return 1.0;
    }
  };


  //========================================================================
  /// Z2-error-estimator:  Elements that can
  /// be used with Z2 error estimation should be derived from
  /// the base class ElementWithZ2ErrorEstimator and implement its
  /// pure virtual member functions to provide the following functionality:
  /// - pointer-based access to the vertex nodes in the element
  ///   (this is required to facilitate setup of element patches).
  /// - the computation of a suitable "flux vector" which represents
  ///   a quantity whose FE representation is discontinuous across
  ///   element boundaries but would become continuous under infinite
  ///   mesh refinement.
  ///
  /// As an example consider the finite element solution of the Laplace problem,
  /// \f$ \partial^2 u/\partial x_i^2 = 0 \f$. If we approximate the
  /// unknown \f$ u \f$ on a finite element mesh with \f$ N \f$ nodes as
  /// \f[
  /// u^{[FE]}(x_k) = \sum_{j=1}^{N} U_j \ \psi_j(x_k),
  /// \f]
  /// where the \f$ \psi_j(x_k) \f$ are the (global) \f$ C^0 \f$ basis
  /// functions, the finite-element representation of the flux,
  /// \f$ f_i = \partial u/\partial x_i \f$,
  /// \f[
  /// f_i^{[FE]} = \sum_{j=1}^{N} U_j \ \frac{\partial \psi_j}{\partial x_i}
  /// \f]
  /// is discontinuous between elements but the magnitude of the jump
  /// decreases under mesh refinement.  We denote the number
  /// of flux terms by \f$N_{flux}\f$, so for a 2D (3D) Laplace problem,
  /// \f$N_{flux}=2 \ (3).\f$
  ///
  /// The idea behind Z2 error estimation is to compute an
  /// approximation for the flux that is more accurate than the value
  /// \f$ f_i^{[FE]} \f$  obtained directly from the finite element
  /// solution. We refer to the improved approximation for the flux
  /// as the "recovered flux" and denote it by \f$ f_i^{[rec]} \f$. Since
  /// \f$ f_i^{[rec]} \f$ is more accurate than \f$ f_i^{[FE]} \f$, the
  /// difference between the two provides a measure of the error.
  /// In Z2 error estimation, the "recovered flux" is determined by
  /// projecting the discontinuous, FE-based flux \f$ f_i^{[FE]} \f$
  /// onto a set of continuous basis functions, known as the "recovery shape
  /// functions". In principle, one could use the
  /// finite element shape functions \f$ \psi_j(x_k) \f$ as the
  /// recovery shape functions but then the determination of the
  /// recovered flux would involve the solution of a linear system
  /// that is as big as the original problem. To avoid this, the projection
  /// is performed over small patches of elements within which
  /// low-order polynomials (defined in terms of the global, Eulerian
  /// coordinates) are used as the recovery shape functions.
  ///
  /// Z2 error estimation is known to be "optimal" for many self-adjoint
  /// problems. See, e.g., Zienkiewicz & Zhu's original papers
  /// "The superconvergent patch recovery and a posteriori error estimates."
  /// International Journal for Numerical Methods in Engineering  \b 33 (1992)
  /// Part I: 1331-1364; Part II: 1365-1382.
  /// In non-self adjoint problems, the error estimator only
  /// analyses the "self-adjoint part" of the differential operator.
  /// However, in many cases, this still produces good error indicators
  /// since the error estimator detects under-resolved, sharp gradients
  /// in the solution.
  ///
  /// Z2 error estimation works as follows:
  /// -# We combine the elements in the mesh into a number of "patches",
  ///    which consist of all elements that share a common vertex node.
  ///    Most elements will therefore be members of multiple patches.
  /// -# Within each patch \f$p\f$, we expand the "recovered flux" as
  ///    \f[
  ///    f^{[rec](p)}_i(x_k) = \sum_{j=1}^{N_{rec}}
  ///    F^{(p)}_{ij} \ \psi^{[rec]}_j(x_k) \mbox{ \ \ \ for
  ///    $i=1,...,N_{flux}$,} \f] where the functions \f$ \psi^{[rec]}_j(x_k)\f$
  ///    are the recovery shape functions, which are functions of the global,
  ///    Eulerian coordinates. Typically, these are chosen to be low-order
  ///    polynomials. For instance, in 2D, a bilinear representation of \f$
  ///    f^{(p)}_i(x_0,x_1) \f$ involves the \f$N_{rec}=3\f$ recovery shape
  ///    functions \f$ \psi^{[rec]}_0(x_0,x_1)=1, \ \psi^{[rec]}_1(x_0,x_1)=x_0
  ///    \f$ and \f$ \psi^{[rec]}_2(x_0,x_1)=x_1\f$.
  ///
  ///    We determine the coefficients \f$ F^{(p)}_{ij} \f$ by enforcing
  ///    \f$ f^{(p)}_i(x_k) =  f^{[FE]}_i(x_k)\f$ in its weak form:
  ///    \f[
  ///    \int_{\mbox{Patch $p$}} \left(
  ///    f^{[FE]}_i(x_k) - \sum_{j=1}^{N_{rec}}
  ///    F^{(p)}_{ij} \ \psi^{[rec]}_j(x_k) \right) \psi^{[rec]}_l(x_k)\ dv = 0
  ///    \mbox{ \ \ \ \ for $l=1,...,N_{rec}$ and $i=1,...,N_{flux}$}.
  ///    \f]
  ///    Once the \f$ F^{(p)}_{ij} \f$  are determined in a given patch,
  ///    we determine the values of the recovered flux at
  ///    all nodes that are part of the patch. We denote the
  ///    value of the recovered flux at node \f$ k \f$ by
  ///    \f$ {\cal F}^{(p)}_{ik}\f$.
  ///
  ///    We repeat this procedure for every patch. For nodes that are part of
  ///    multiple patches, the procedure
  ///    will provide multiple, slightly different nodal values for the
  ///    recovered flux. We average these values via \f[
  ///    {\cal F}_{ij} = \frac{1}{N_p(j)}
  ///    \sum_{\mbox{Node $j \in $ patch $p$}}
  ///    {\cal F}^{(p)}_{ij},
  ///    \f]
  ///    where \f$N_p(j)\f$ denotes the number of patches that node \f$ j\f$ is
  ///    a member of. This allows us to obtain a globally-continuous,
  ///    finite-element based representation of the recovered flux as \f[
  ///    f_i^{[rec]} = \sum_{j=1}^{N} {\cal F}_{ij}\ \psi_j,
  ///    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ (1)
  ///    \f]
  ///    where the \f$ \psi_j \f$ are the (global) finite element
  ///    shape functions.
  /// -# Since the recovered flux in (1), is based on nodal values, we can
  ///    evaluate it locally within each of the \f$ N_e\f$ elements in the mesh
  ///    to obtain a normalised elemental error estimate via
  ///    \f[
  ///    E_{e} =  \sqrt{ \frac{
  ///             \int_{\mbox{$e$}} \left(
  ///              f_i^{[rec]}  - f_i^{[FE]} \right)^2 dv}
  ///              {\sum_{e'=1}^{N_e}
  ///               \int_{\mbox{$e'$}} \left(
  ///                f_i^{[rec]} \right)^2 dv} } \mbox{\ \ \ for $e=1,...,N_e$.}
  ///    \f]
  ///    In this (default) form, mesh refinement, based on this error estimate
  ///    will lead to an equidistribution of the error across all elements.
  ///    Usually, this is the desired behaviour. However, there are
  ///    cases in which the solution evolves towards a state in which
  ///    the flux tends to zero while the solution itself becomes so simple
  ///    that it can be represented exactly on any finite element mesh
  ///    (e.g. in spin-up problems in which the fluid motion approaches
  ///    a rigid body rotation). In that case the mesh adaptation tries
  ///    to equidistribute any small roundoff errors in the solution,
  ///    causing strong, spatially uniform mesh refinement throughout
  ///    the domain, even though the solution is already fully converged.
  ///    For such cases, it is possible to replace the denominator in the
  ///    above expression by a prescribed reference flux, which may be
  ///    specified via the access function
  ///    \code
  ///    reference_flux_norm()
  ///    \endcode
  ///
  ///
  //========================================================================
  class Z2ErrorEstimator : public virtual ErrorEstimator
  {
  public:
    /// \short Function pointer to combined error estimator function
    typedef double (*CombinedErrorEstimateFctPt)(const Vector<double>& errors);

    /// Constructor: Set order of recovery shape functions
    Z2ErrorEstimator(const unsigned& recovery_order)
      : Recovery_order(recovery_order),
        Recovery_order_from_first_element(false),
        Reference_flux_norm(0.0),
        Combined_error_fct_pt(0)
    {
    }


    /// \short Constructor: Leave order of recovery shape functions open
    /// for now -- they will be read out from first element of the mesh
    /// when the error estimator is applied
    Z2ErrorEstimator()
      : Recovery_order(0),
        Recovery_order_from_first_element(true),
        Reference_flux_norm(0.0),
        Combined_error_fct_pt(0)
    {
    }

    /// Broken copy constructor
    Z2ErrorEstimator(const Z2ErrorEstimator&)
    {
      BrokenCopy::broken_copy("Z2ErrorEstimator");
    }

    /// Broken assignment operator
    void operator=(const Z2ErrorEstimator&)
    {
      BrokenCopy::broken_assign("Z2ErrorEstimator");
    }

    /// Empty virtual destructor
    virtual ~Z2ErrorEstimator() {}

    /// \short Compute the elemental error measures for a given mesh
    /// and store them in a vector.
    void get_element_errors(Mesh*& mesh_pt, Vector<double>& elemental_error)
    {
      // Create dummy doc info object and switch off output
      DocInfo doc_info;
      doc_info.disable_doc();
      // Forward call to version with doc.
      get_element_errors(mesh_pt, elemental_error, doc_info);
    }

    /// \short Compute the elemental error measures for a given mesh
    /// and store them in a vector.
    /// If doc_info.enable_doc(), doc FE and recovered fluxes in
    /// - flux_fe*.dat
    /// - flux_rec*.dat
    void get_element_errors(Mesh*& mesh_pt,
                            Vector<double>& elemental_error,
                            DocInfo& doc_info);

    /// Access function for order of recovery polynomials
    unsigned& recovery_order()
    {
      return Recovery_order;
    }

    /// Access function for order of recovery polynomials (const version)
    unsigned recovery_order() const
    {
      return Recovery_order;
    }

    /// Access function: Pointer to combined error estimate function
    CombinedErrorEstimateFctPt& combined_error_fct_pt()
    {
      return Combined_error_fct_pt;
    }

    ///\short  Access function: Pointer to combined error estimate function.
    /// Const version
    CombinedErrorEstimateFctPt combined_error_fct_pt() const
    {
      return Combined_error_fct_pt;
    }

    /// \short Setup patches: For each vertex node pointed to by nod_pt,
    /// adjacent_elements_pt[nod_pt] contains the pointer to the vector that
    /// contains the pointers to the elements that the node is part of.
    void setup_patches(Mesh*& mesh_pt,
                       std::map<Node*, Vector<ElementWithZ2ErrorEstimator*>*>&
                         adjacent_elements_pt,
                       Vector<Node*>& vertex_node_pt);

    /// Access function for prescribed reference flux norm
    double& reference_flux_norm()
    {
      return Reference_flux_norm;
    }

    /// Access function for prescribed reference flux norm (const. version)
    double reference_flux_norm() const
    {
      return Reference_flux_norm;
    }

    /// Return a combined error estimate from all compound errors
    double get_combined_error_estimate(const Vector<double>& compound_error);

  private:
    /// \short Given the vector of elements that make up a patch,
    /// the number of recovery and flux terms, and the spatial
    /// dimension of the problem, compute
    /// the matrix of recovered flux coefficients and return
    /// a pointer to it.
    void get_recovered_flux_in_patch(
      const Vector<ElementWithZ2ErrorEstimator*>& patch_el_pt,
      const unsigned& num_recovery_terms,
      const unsigned& num_flux_terms,
      const unsigned& dim,
      DenseMatrix<double>*& recovered_flux_coefficient_pt);


    /// \short Return number of coefficients for expansion of recovered fluxes
    /// for given spatial dimension of elements.
    /// (We use complete polynomials of the specified given order.)
    unsigned nrecovery_terms(const unsigned& dim);


    /// \short Recovery shape functions as functions of the global, Eulerian
    /// coordinate x of dimension dim.
    /// The recovery shape functions are  complete polynomials of
    /// the order specified by Recovery_order.
    void shape_rec(const Vector<double>& x,
                   const unsigned& dim,
                   Vector<double>& psi_r);

    /// \short Integation scheme associated with the recovery shape functions
    /// must be of sufficiently high order to integrate the mass matrix
    /// associated with the recovery shape functions
    Integral* integral_rec(const unsigned& dim, const bool& is_q_mesh);

    /// Order of recovery polynomials
    unsigned Recovery_order;

    /// Bool to indicate if recovery order is to be read out from
    /// first element in mesh or set globally
    bool Recovery_order_from_first_element;

    /// Doc flux and recovered flux
    void doc_flux(Mesh* mesh_pt,
                  const unsigned& num_flux_terms,
                  MapMatrixMixed<Node*, int, double>& rec_flux_map,
                  const Vector<double>& elemental_error,
                  DocInfo& doc_info);

    /// Prescribed reference flux norm
    double Reference_flux_norm;

    /// Function pointer to combined error estimator function
    CombinedErrorEstimateFctPt Combined_error_fct_pt;
  };


  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////


  //========================================================================
  /// Dummy error estimator, allows manual specification of refinement
  /// pattern by forcing refinement in regions defined by elements in
  /// a reference mesh.
  //========================================================================
  class DummyErrorEstimator : public virtual ErrorEstimator
  {
  public:
    /// \short Constructor. Provide mesh and number of the elements that define
    /// the regions within which elements are to be refined subsequently.
    /// Also specify the node number of a central node
    /// within elements -- it's used to determine if an element is
    /// in the region where refinement is supposed to take place.
    /// Optional boolean flag (defaulting to false) indicates that
    /// refinement decision is based on Lagrangian coordinates -- only
    /// applicable to solid meshes.
    DummyErrorEstimator(Mesh* mesh_pt,
                        const Vector<unsigned>& elements_to_refine,
                        const unsigned& central_node_number,
                        const bool& use_lagrangian_coordinates = false)
      : Use_lagrangian_coordinates(use_lagrangian_coordinates),
        Central_node_number(central_node_number)
    {
#ifdef PARANOID
#ifdef OOMPH_HAS_MPI
      if (mesh_pt->is_mesh_distributed())
      {
        throw OomphLibError(
          "Can't use this error estimator on distributed meshes!",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
#endif
#endif

#ifdef PARANOID
      if (mesh_pt->nelement() == 0)
      {
        throw OomphLibError(
          "Can't build error estimator if there are no elements in mesh\n",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
#endif

      unsigned dim = mesh_pt->finite_element_pt(0)->node_pt(0)->ndim();
      if (use_lagrangian_coordinates)
      {
        SolidNode* solid_nod_pt =
          dynamic_cast<SolidNode*>(mesh_pt->finite_element_pt(0)->node_pt(0));
        if (solid_nod_pt != 0)
        {
          dim = solid_nod_pt->nlagrangian();
        }
      }
      unsigned nregion = elements_to_refine.size();
      Region_low_bound.resize(nregion);
      Region_upp_bound.resize(nregion);
      for (unsigned e = 0; e < nregion; e++)
      {
        Region_low_bound[e].resize(dim, 1.0e20);
        Region_upp_bound[e].resize(dim, -1.0e20);
        FiniteElement* el_pt =
          mesh_pt->finite_element_pt(elements_to_refine[e]);
        unsigned nnod = el_pt->nnode();
        for (unsigned j = 0; j < nnod; j++)
        {
          Node* nod_pt = el_pt->node_pt(j);
          for (unsigned i = 0; i < dim; i++)
          {
            double x = nod_pt->x(i);
            if (use_lagrangian_coordinates)
            {
              SolidNode* solid_nod_pt = dynamic_cast<SolidNode*>(nod_pt);
              if (solid_nod_pt != 0)
              {
                x = solid_nod_pt->xi(i);
              }
            }
            if (x < Region_low_bound[e][i])
            {
              Region_low_bound[e][i] = x;
            }
            if (x > Region_upp_bound[e][i])
            {
              Region_upp_bound[e][i] = x;
            }
          }
        }
      }
    }


    /// \short Constructor. Provide vectors to "lower left" and "upper right"
    /// vertices of refinement region
    /// Also specify the node number of a central node
    /// within elements -- it's used to determine if an element is
    /// in the region where refinement is supposed to take place.
    /// Optional boolean flag (defaulting to false) indicates that
    /// refinement decision is based on Lagrangian coordinates -- only
    /// applicable to solid meshes.
    DummyErrorEstimator(Mesh* mesh_pt,
                        const Vector<double>& lower_left,
                        const Vector<double>& upper_right,
                        const unsigned& central_node_number,
                        const bool& use_lagrangian_coordinates = false)
      : Use_lagrangian_coordinates(use_lagrangian_coordinates),
        Central_node_number(central_node_number)
    {
#ifdef PARANOID
      if (mesh_pt->nelement() == 0)
      {
        throw OomphLibError(
          "Can't build error estimator if there are no elements in mesh\n",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
#endif

      unsigned nregion = 1;
      Region_low_bound.resize(nregion);
      Region_upp_bound.resize(nregion);
      Region_low_bound[0] = lower_left;
      Region_upp_bound[0] = upper_right;
    }


    /// Broken copy constructor
    DummyErrorEstimator(const DummyErrorEstimator&)
    {
      BrokenCopy::broken_copy("DummyErrorEstimator");
    }


    /// Broken assignment operator
    void operator=(const DummyErrorEstimator&)
    {
      BrokenCopy::broken_assign("DummyErrorEstimator");
    }


    /// Empty virtual destructor
    virtual ~DummyErrorEstimator() {}


    /// \short Compute the elemental error measures for a given mesh
    /// and store them in a vector. Doc errors etc.
    virtual void get_element_errors(Mesh*& mesh_pt,
                                    Vector<double>& elemental_error,
                                    DocInfo& doc_info)
    {
#ifdef PARANOID
      if (doc_info.is_doc_enabled())
      {
        std::ostringstream warning_stream;
        warning_stream
          << "No output defined in DummyErrorEstimator::get_element_errors()\n"
          << "Ignoring doc_info flag.\n";
        OomphLibWarning(warning_stream.str(),
                        "DummyErrorEstimator::get_element_errors()",
                        OOMPH_EXCEPTION_LOCATION);
      }
#endif
      unsigned nregion = Region_low_bound.size();
      unsigned nelem = mesh_pt->nelement();
      for (unsigned e = 0; e < nelem; e++)
      {
        elemental_error[e] = 0.0;

        // Check if element is in the regions to be refined
        // (based on coords of its central node)
        Node* nod_pt =
          mesh_pt->finite_element_pt(e)->node_pt(Central_node_number);
        for (unsigned r = 0; r < nregion; r++)
        {
          bool is_inside = true;
          unsigned dim = Region_low_bound[r].size();
          for (unsigned i = 0; i < dim; i++)
          {
            double x = nod_pt->x(i);
            if (Use_lagrangian_coordinates)
            {
              SolidNode* solid_nod_pt = dynamic_cast<SolidNode*>(nod_pt);
              if (solid_nod_pt != 0)
              {
                x = solid_nod_pt->xi(i);
              }
            }
            if (x < Region_low_bound[r][i])
            {
              is_inside = false;
              break;
            }
            if (x > Region_upp_bound[r][i])
            {
              is_inside = false;
              break;
            }
          }
          if (is_inside)
          {
            elemental_error[e] = 1.0;
            break;
          }
        }
      }
    }

  private:
    /// \short Use Lagrangian coordinates to decide which element is to be
    /// refined
    bool Use_lagrangian_coordinates;

    /// \short Number of local node that is used to identify if an element
    /// is located in the refinement region
    unsigned Central_node_number;

    /// Upper bounds for the coordinates of the refinement regions
    Vector<Vector<double>> Region_upp_bound;

    /// Lower bounds for the coordinates of the refinement regions
    Vector<Vector<double>> Region_low_bound;
  };


} // namespace oomph

#endif
