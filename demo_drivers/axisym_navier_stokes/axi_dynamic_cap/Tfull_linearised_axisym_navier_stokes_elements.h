#ifndef TFULL_LINEARISED_AXISYM_NAVIER_STOKES_HEADER
#define TFULL_LINEARISED_AXISYM_NAVIER_STOKES_HEADER

//#include "generic.h"
#include "full_linearised_axisym_navier_stokes_elements.h"

namespace oomph
{
  class TaylorHoodFullLinearisedAxisymmetricNavierStokesEquations
    : public virtual FullLinearisedAxisymmetricNavierStokesEquations,
      public virtual TElement<2, 3>
  {
  private:
    static const unsigned Initial_Nvalue[];

  protected:
    /// Static array of ints to hold conversion from pressure
    /// node numbers to actual node numbers
    static const unsigned Pconv[];

    inline double dshape_and_dtest_eulerian_lin_axi_nst(const Vector<double>& s,
                                                        Shape& psi,
                                                        DShape& dpsidx,
                                                        Shape& test,
                                                        DShape& dtestdx) const
    {
      // Call the geometrical shape functions and derivatives
      double J = this->dshape_eulerian(s, psi, dpsidx);
      // Test functions are the shape functions
      test = psi;
      dtestdx = dpsidx;
      // Return the jacobian
      return J;
    }

    inline double dshape_and_dtest_eulerian_at_knot_lin_axi_nst(
      const unsigned& ipt,
      Shape& psi,
      DShape& dpsidx,
      Shape& test,
      DShape& dtestdx) const
    {
      // Call the geometrical shape functions and derivatives
      double J = this->dshape_eulerian_at_knot(ipt, psi, dpsidx);
      // Test functions are the shape functions
      test = psi;
      dtestdx = dpsidx;
      // Return the jacobian
      return J;
    }

    inline void pshape_lin_axi_nst(const Vector<double>& s, Shape& psi) const
    {
      psi[0] = s[0];
      psi[1] = s[1];
      psi[2] = 1.0 - s[0] - s[1];
    }

    inline void pshape_lin_axi_nst(const Vector<double>& s,
                                   Shape& psi,
                                   Shape& test) const
    {
      // Call the pressure shape functions
      this->pshape_lin_axi_nst(s, psi);
      // Test functions are shape functions
      test = psi;
    }


  public:
    inline virtual unsigned required_nvalue(const unsigned& n) const
    {
      return Initial_Nvalue[n];
    }

    virtual unsigned npres_lin_axi_nst() const
    {
      return 3;
    }

    /// Return number of velocity values at each node in the element
    unsigned n_u_lin_axi_nst() const
    {
      return 6;
    }

    double p_lin_axi_nst(const unsigned& n_p, const unsigned& i) const
    {
      return this->nodal_value(Pconv[n_p], 6 + i);
    }

    virtual void pin_pressure(const unsigned& i)
    {
      for (unsigned n_p = 0; n_p < 3; n_p++)
      {
        this->node_pt(Pconv[n_p])->pin(6 + i);
      }
    }

    /// Return the local equation numbers for the pressure values.
    virtual int p_local_eqn(const unsigned& n, const unsigned& i)
    {
      return this->nodal_local_eqn(Pconv[n], 6 + i);
    }

    virtual inline unsigned xhat_index_lin_axi_nst(const unsigned& n,
                                                   const unsigned& i) const
    {
      double index = 6;
      if (n < 3)
      {
        index += 2;
      }

      return i + index;
    }

    void pin()
    {
      for (unsigned i = 0; i < 6; i++)
      {
        for (unsigned j = 0; j < 6; j++)
        {
          this->node_pt(i)->pin(j);
        }
      }
      for (unsigned i = 0; i < 3; i++)
      {
        for (unsigned j = 0; j < 2; j++)
        {
          this->node_pt(i)->pin(6 + j);
        }
      }
    }

    virtual void output(std::ostream& outfile)
    {
      FullLinearisedAxisymmetricNavierStokesEquations::output(outfile);
    }

    virtual void output(std::ostream& outfile, const unsigned& n_plot)
    {
      FullLinearisedAxisymmetricNavierStokesEquations::output(outfile, n_plot);
    }

    virtual void output(FILE* file_pt)
    {
      FullLinearisedAxisymmetricNavierStokesEquations::output(file_pt);
    }

    virtual void output(FILE* file_pt, const unsigned& n_plot)
    {
      FullLinearisedAxisymmetricNavierStokesEquations::output(file_pt, n_plot);
    }
  };

  template<>
  class FaceGeometry<TaylorHoodFullLinearisedAxisymmetricNavierStokesEquations>
    : public TElement<1, 3>
  {
  public:
    /// Constructor: Call the constructor for the
    /// appropriate lower-dimensional QElement
    FaceGeometry() : TElement<1, 3>() {}
  };

} // namespace oomph

#endif
