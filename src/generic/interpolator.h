#ifndef OOMPH_INTERPOLATOR_H
#define OOMPH_INTERPOLATOR_H

namespace oomph
{

 // Forward decl.
 class FiniteElement;
 class Time;

 namespace InterpolatorHelpers
 {

  /// Determine if an element contains any hanging nodes
  inline bool has_hanging_nodes(const FiniteElement* ele_pt)
  {
   for(unsigned nd=0, nnode=ele_pt->nnode(); nd<nnode; nd++)
    {
     Node* nd_pt = ele_pt->node_pt(nd);
     if(nd_pt->is_hanging())
      {
       return true; 
      }
    }
   return false;
  }


  /// If paranoid check that timesteppers in all nodes in the element are
  /// the same.
  inline void check_timesteppers_same(const FiniteElement* ele_pt)
  {
#ifdef PARANOID
   for(unsigned nd=0, nnode=ele_pt->nnode(); nd<nnode; nd++)
    {
     Node* nd_pt = ele_pt->node_pt(nd);

     // Check ts_pts are all the same
     if(nd_pt->time_stepper_pt() != ele_pt->node_pt(0)->time_stepper_pt())
      {
       std::ostringstream error_msg;
       error_msg << "Time steppers should all be the same within an element!";
       throw OomphLibError(error_msg.str(), OOMPH_CURRENT_FUNCTION,
                           OOMPH_EXCEPTION_LOCATION);
      }

     //??ds external data? internal data? etc?
    }
#endif
  }


  /// If paranoid check that all nodes in the element have the same
  /// nvalues.
  inline void check_nvalues_in_element_same(const FiniteElement* ele_pt)
  {
#ifdef PARANOID
   unsigned Nvalue = ele_pt->node_pt(0)->nvalue();
   for(unsigned nd=0, nnode=ele_pt->nnode(); nd<nnode; nd++)
    {
     Node* nd_pt = ele_pt->node_pt(nd);

     // Check all number of values is the same for all nodes
     if(nd_pt->nvalue() != Nvalue)
      {
       std::ostringstream error_msg;
       error_msg << "Number of values must be the same at all nodes to use this interpolator.";
       throw OomphLibError(error_msg.str(), OOMPH_CURRENT_FUNCTION,
                           OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif
  }


  /// ??ds this is probably nasty, switch to using only pointers so we have
  /// a sensible null value?
  static double NotYetCalculatedValue;


  /// Check to see if variables are initialised
  inline bool uninitialised(double var) {return var == NotYetCalculatedValue;}


  /// \short Check to see if a vector is initialised. Recusive checking
  /// would be faster but takes a fairly large amount of time.
  template<typename T>
  inline bool uninitialised(Vector<T> var)
  {
   bool uninitialised_entry = false;
   for(unsigned j=0; j<var.size(); j++)
    {
     if(uninitialised(var[j]))
      {
       uninitialised_entry = true;
       break;
      }
    }
   return var.empty() || uninitialised_entry;
   // return var.empty();
  }


  /// Check to see if a pointer/array is uninitialised
  inline bool uninitialised(const double* var)
  {
   if(var == 0) return true;
   else return uninitialised(var[0]);
  }

 }

 using namespace InterpolatorHelpers;


 /* Implementation notes:

  * Use memoization to ensure we don't duplicate calculations.

  * Lots of code is repeated because c++ function pointers suck

  * Use Vector<Vector<double> > rather than matrix so that we can easily
  (and efficiently) return vectors of a single row (e.g. all derivatives of
  one value).

  * No need to template by DIM or Nnode: it has no effect on speed (tested).


  TODO:

  * We really have two classes here: a memoiser and an interpolator.  Split
  them. Memoising interface over an interpolator is probably best, then we
  can swap the interpolator out without breaking any memoising.

  * Pass in parameter for time deriv order to interpolation (instead of
  separate interpolate_dxxxxdt functions)?

  LIMITATIONS:

  * Can be slow due to creation/destruction of vectors. This could be fixed
  by adding the option to clear the memory or (as I'm now using for
  micromagnetics problems --David) by using raw C-arrays instead.

  * Only simple shape functions are used. Fancier things should be easy
  enough to implement if needed though.

  * Number of values and time stepper pt must be the same at all nodes
  within the element (so can't use it for Navier-Stokes). This can probably
  be allowed if needed, possibly at a small speed penalty.

 */

 // ============================================================
 /// Interpolator to easily calculate values, positions, derivatives
 /// etc. within elements.
 // ============================================================
 class GeneralInterpolator
 {

 public:

  /// Default constructor. Choose element to interpolate within and
  /// position in element to interpolate at.
  GeneralInterpolator(const FiniteElement* this_element,
                      const Vector<double> &s)
   :
   // Cache some loop end conditions
   Nnode(this_element->nnode()),
   Dim(this_element->dim()),
   Nprev_value_zeroth_derivative(this_element->node_pt(0)->
                                 time_stepper_pt()->
                                 nprev_values_for_value_at_evaluation_time()),
   Nprev_value_zeroth_pos_derivative(this_element->node_pt(0)->
                                     position_time_stepper_pt()->
                                     nprev_values_for_value_at_evaluation_time()),
   Nprev_value_derivative(this_element->node_pt(0)->
                          time_stepper_pt()->ntstorage()),
   Nprev_value_pos_derivative(this_element->node_pt(0)->
                              position_time_stepper_pt()->ntstorage()),
   Nvalue(this_element->node_pt(0)->nvalue()),

  // Initialise pointers
   This_element(this_element),
   Ts_pt(this_element->node_pt(0)->time_stepper_pt()),
   Ts_weights_pt(Ts_pt->weights_pt()),
   Position_ts_weights_pt(this_element->node_pt(0)->position_time_stepper_pt()
                          ->weights_pt()),

  // Initialise storage for shape functions
   Psi(Nnode),
   Test(Nnode),
   Dpsidx(Nnode, Dim),
   Dtestdx(Nnode, Dim),
   S(s),
   
  // Use negative time to signify that it has not been calculated yet
   Intp_time(NotYetCalculatedValue),

  // Initialise storage for value derivatives (because this is a two
  // dimensional array we need to do some initialisation now).
   Dvaluesdx(Nvalue),
   D2valuesdxdt(Nvalue)
   { 
    // If paranoid check that all nodes have the same nvalues.
    InterpolatorHelpers::check_nvalues_in_element_same(this_element);

    build();
   }

  /// More complex constructor. Choose element to interpolate within,
  /// position in element to interpolate at and a time stepper to do time
  /// interpolation. WARNING: if you use this construtor you need to make
  /// sure that the previous values stored (in nodes) have the same meaning
  /// in both the given time stepper and the nodes' time steppers.
  GeneralInterpolator(const FiniteElement* this_element,
                      const Vector<double> &s,
                      const TimeStepper* ts_pt)
   :
   // Cache some loop end conditions
   Nnode(this_element->nnode()),
   Dim(this_element->dim()),
   Nprev_value_zeroth_derivative(ts_pt->
                                 nprev_values_for_value_at_evaluation_time()),
   Nprev_value_zeroth_pos_derivative(this_element->node_pt(0)->
                                     position_time_stepper_pt()->
                                     nprev_values_for_value_at_evaluation_time()),
   Nprev_value_derivative(ts_pt->ntstorage()),
   Nprev_value_pos_derivative(this_element->node_pt(0)->
                              position_time_stepper_pt()->ntstorage()),
   Nvalue(this_element->node_pt(0)->nvalue()),

  // Initialise pointers
   This_element(this_element),
   Ts_pt(ts_pt),
   Ts_weights_pt(ts_pt->weights_pt()),
   Position_ts_weights_pt(this_element->node_pt(0)->position_time_stepper_pt()
                          ->weights_pt()),

  // Initialise storage for shape functions
   Psi(Nnode),
   Test(Nnode),
   Dpsidx(Nnode, Dim),
   Dtestdx(Nnode, Dim),
   S(s),
   
  // Use negative time to signify that it has not been calculated yet
   Intp_time(NotYetCalculatedValue),

  // Initialise storage for value derivatives (because this is a two
  // dimensional array we need to do some initialisation now).
   Dvaluesdx(Nvalue),
   D2valuesdxdt(Nvalue)
   {
    build();
   }


  virtual void build()
  {
   // Set up shape + test functions
   J = This_element->dshape_eulerian(S, Psi, Dpsidx);
   Test = Psi;
   Dtestdx = Dpsidx;

   // Find out if any nodes in this element are hanging
   Has_hanging_nodes = InterpolatorHelpers::has_hanging_nodes(This_element);

   // If paranoid check that all nodes have the same nvalues.
   InterpolatorHelpers::check_nvalues_in_element_same(This_element);
  }

  /// Destructor
  virtual ~GeneralInterpolator() {}

  const TimeStepper* ts_pt() const {return Ts_pt;}

  double time()
  {
   if(uninitialised(Intp_time)) {Intp_time = interpolate_time();}
   return Intp_time;
  }

  double x(const unsigned &i)
  {return x()[i];}

  const Vector<double> &x()
  {
   if(uninitialised(X)) X = interpolate_x();
   return X;
  }

  const Vector<double> &s()
   {
    return S;
   }

  double dxdt(const unsigned &i)
  {return dxdt()[i];}

  const Vector<double> &dxdt()
  {
   if(uninitialised(Dxdt)) Dxdt = interpolate_dxdt();
   return Dxdt;
  }

  double value(const unsigned &i_val)
  {return value()[i_val];}

  const Vector<double> &value()
  {
   if(uninitialised(Values)) Values = interpolate_values(0, Nvalue);
   return Values;
  }

  double dvaluedt(const unsigned &i_val)
  {return dvaluedt()[i_val];}

  const Vector<double> &dvaluedt()
  {
   if(uninitialised(Dvaluesdt)) Dvaluesdt = interpolate_dvaluesdt(0, Nvalue);
   return Dvaluesdt;
  }

  double dvaluedx(const unsigned &i_val, const unsigned &direction)
  {return dvaluedx(i_val)[direction];}


   const Vector<double>& dvaluedx(const unsigned &i_val)
   {
     if(uninitialised(Dvaluesdx[i_val]))
       {
         Dvaluesdx[i_val] = interpolate_dvaluesdx(i_val);
       }
     return Dvaluesdx[i_val];
   }

  const Vector<double>& d2valuedxdt(const unsigned &i_val)
   {
     if(uninitialised(D2valuesdxdt[i_val]))
       {
         D2valuesdxdt[i_val] = interpolate_d2valuesdxdt(i_val);
       }
     return D2valuesdxdt[i_val];
   }

  // Access functions for Jacobian and shape/test functions
  double j() const {return J;}
  double psi(const unsigned &i) const {return Psi(i);}
  double test(const unsigned &i) const {return Test(i);}
  double dpsidx(const unsigned &i, const unsigned &i_direc) const
  {return Dpsidx(i, i_direc);}
  double dtestdx(const unsigned &i, const unsigned &i_direc) const
  {return Dtestdx(i, i_direc);}


  /// Interpolate evaluation position
  virtual Vector<double> interpolate_x()
  {
   if(!Has_hanging_nodes) {return interpolate_x_raw();}
   else {return interpolate_x_with_hanging_nodes();}
  }

  /// Interpolate derivative of position
  virtual Vector<double> interpolate_dxdt() const
  {
   if(!Has_hanging_nodes) {return interpolate_dxdt_raw();}
   else {return interpolate_dxdt_with_hanging_nodes();}

  }

  /// Interpolate evaluation time
  virtual double interpolate_time()
  {
   double time = 0.0;
   Time* time_pt = This_element->node_pt(0)->time_stepper_pt()->time_pt();
   for(unsigned i_tm=0; i_tm<Nprev_value_zeroth_derivative; i_tm++)
    {
     // the number of pointers chained together here seems a bit silly... ds
     time += time_pt->time(i_tm) * (*Ts_weights_pt)(0,i_tm);
    }
   return time;
  }

  // Interpolate values (by default get all, optionally just get some range.
  virtual Vector<double> interpolate_values(const unsigned &start,
                                            const unsigned &end)
  {
   if(!Has_hanging_nodes) {return interpolate_values_raw(start, end);}
   else {return interpolate_values_with_hanging_nodes(start, end);}
  }

  virtual Vector<double> interpolate_dvaluesdt(const unsigned &start,
                                               const unsigned &end)
  {
   if(!Has_hanging_nodes) {return interpolate_dvaluesdt_raw(start, end);}
   else {return interpolate_dvaluesdt_with_hanging_nodes(start, end);}
  }

  virtual Vector<double> interpolate_dvaluesdx(const unsigned &i_val)
  {
   if(!Has_hanging_nodes) {return interpolate_dvaluesdx_raw(i_val);}
   else {return interpolate_dvaluesdx_with_hanging_nodes(i_val);}
  }

  virtual Vector<double> interpolate_d2valuesdxdt(const unsigned &i_val)
  {
   if(!Has_hanging_nodes) {return interpolate_d2valuesdxdt_raw(i_val);}
   else {return interpolate_d2valuesdxdt_with_hanging_nodes(i_val);}
  }

 protected:

  // Loop end conditions etc.
  const unsigned Nnode;
  const unsigned Dim;
  const unsigned Nprev_value_zeroth_derivative;
  const unsigned Nprev_value_zeroth_pos_derivative;
  const unsigned Nprev_value_derivative;
  const unsigned Nprev_value_pos_derivative;
  const unsigned Nvalue;

  bool Has_hanging_nodes;

  // Cached pointers
  const FiniteElement* This_element;
  const TimeStepper* Ts_pt;
  const DenseMatrix<double>* Ts_weights_pt;
  const DenseMatrix<double>* Position_ts_weights_pt;

  // Jacobian + shape/test functions
  double J;
  Shape Psi;
  Shape Test;
  DShape Dpsidx;
  DShape Dtestdx;
  const Vector<double> S;

  // Interpolated value storage (note we can't name the time variable
  // "Time" because that is used for the time class).
  double Intp_time;

  Vector<double> X;
  Vector<double> Dxdt;

  Vector<double> Values;
  Vector<double> Dvaluesdt;
  Vector<Vector<double> > Dvaluesdx;
  Vector<Vector<double> > D2valuesdxdt;


  // Checks to see if private variables are initialised
  bool uninitialised(double var) {return var == NotYetCalculatedValue;}

  /// \short Check to see if a vector is initialised. Recusive checking
  /// is more robust but takes a fairly large amount of time.
  template<typename T> bool uninitialised(Vector<T> var)
  {
    bool uninitialised_entry = false;
    for(unsigned j=0; j<var.size(); j++)
     {
      if(uninitialised(var[j]))
       {
        uninitialised_entry = true;
        break;
       }
     }
    return var.empty() || uninitialised_entry;
  }


  // Position interpolation
  // ============================================================

  Vector<double> interpolate_x_with_hanging_nodes() const
  {
   Vector<double> output(Dim, 0.0);
   for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
    {
     for(unsigned i_tm=0; i_tm<Nprev_value_zeroth_pos_derivative; i_tm++)
      {
       for(unsigned j=0; j<Dim; j++)
        {
         output[j] += This_element->nodal_position(i_tm, i_nd, j)
          * Psi(i_nd) * (*Position_ts_weights_pt)(0, i_tm);
        }
      }
    }
   return output;
  }

  /// \short Use if no nodes are hanging (raw access functions are faster).
  Vector<double> interpolate_x_raw() const
  {
   Vector<double> output(Dim, 0.0);
   for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
    {
     for(unsigned i_tm=0; i_tm<Nprev_value_zeroth_pos_derivative; i_tm++)
      {
       for(unsigned j=0; j<Dim; j++)
        {
         output[j] += This_element->raw_nodal_position(i_tm, i_nd, j)
          * Psi(i_nd) * (*Position_ts_weights_pt)(0, i_tm);
        }
      }
    }
   return output;
  }

  // Position derivatives interpolation
  // ============================================================

  Vector<double> interpolate_dxdt_raw() const
  {
   Vector<double> output(Dim);
   for(unsigned l=0;l<Nnode;l++)
    {
     for(unsigned j=0;j<Dim;j++)
      {
       for(unsigned t=0;t<Nprev_value_pos_derivative;t++)
        {
         output[j] += (*Position_ts_weights_pt)(1,t)
          * This_element->raw_nodal_position(t, l, j) * Psi(l);
        }
      }
    }
   return output;
  }

  Vector<double> interpolate_dxdt_with_hanging_nodes() const
  {
   Vector<double> output(Dim);
   for(unsigned l=0;l<Nnode;l++)
    {
     for(unsigned j=0;j<Dim;j++)
      {
       for(unsigned t=0;t<Nprev_value_pos_derivative;t++)
        {
         output[j] += (*Position_ts_weights_pt)(1,t)
          * This_element->nodal_position(t, l, j) * Psi(l);
        }
      }
    }
   return output;
  }


  // Interpolate values
  // ============================================================

  //??ds copy this to all function's doc

  /// \short Interpolate [something] [with/without] considering hanging nodes (it
  /// is faster to use "raw" access functions which ignore hanging nodes
  /// where possible). Interpolates a range of values [start, end)
  /// (i.e. including start but not end).

  Vector<double> interpolate_values_raw(const unsigned &start,
                                        const unsigned &end) const
  {
   Vector<double> output(end - start);
   for(unsigned j=start; j<end; j++)
    {
     for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
      {
       for(unsigned i_tm=0; i_tm<Nprev_value_zeroth_derivative; i_tm++)
        {
         output[j-start] += This_element->raw_nodal_value(i_tm, i_nd, j)
          * Psi(i_nd) * (*Ts_weights_pt)(0, i_tm);
        }
      }
    }
   
   return output;
  }

  Vector<double> interpolate_values_with_hanging_nodes(const unsigned &start,
                                                       const unsigned &end) const
  {
   Vector<double> output(end - start);
   for(unsigned j=start; j<end; j++)
    {
     for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
      {
       for(unsigned i_tm=0; i_tm<Nprev_value_zeroth_derivative; i_tm++)
        {
         output[j-start] += This_element->nodal_value(i_tm, i_nd, j)
          * Psi(i_nd) * (*Ts_weights_pt)(0, i_tm);
        }
      }
    }
   return output;
  }


  // Interpolate derivatives of values w.r.t. time
  // ============================================================

  Vector<double> interpolate_dvaluesdt_raw(const unsigned &start,
                                           const unsigned &end) const
  {
   Vector<double> output(end - start);
   for(unsigned j=start; j<end; j++)
    {
     for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
      {
       for(unsigned i_tm=0; i_tm<Nprev_value_derivative; i_tm++)
        {
         output[j-start] += This_element->raw_nodal_value(i_tm, i_nd, j)
          * Psi(i_nd) * (*Ts_weights_pt)(1, i_tm);
        }
      }
    }
   return output;
  }

  Vector<double> interpolate_dvaluesdt_with_hanging_nodes(const unsigned &start,
                                                          const unsigned &end) const
  {
   Vector<double> output(end - start);
   for(unsigned j=start; j<end; j++)
    {
     for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
      {
       for(unsigned i_tm=0; i_tm<Nprev_value_derivative; i_tm++)
        {
         output[j-start] += This_element->nodal_value(i_tm, i_nd, j)
          * Psi(i_nd) * (*Ts_weights_pt)(1, i_tm);
        }
      }
    }
   return output;
  }


  // Interpolate derivatives of values w.r.t. position
  // ============================================================


  Vector<double> interpolate_dvaluesdx_raw(const unsigned &i_value) const
  {
    Vector<double> output(Dim, 0.0);
    for(unsigned i_direc=0; i_direc<Dim; i_direc++)
      {
        for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
          {
            for(unsigned i_tm=0; i_tm<Nprev_value_zeroth_derivative; i_tm++)
              {
                output[i_direc] += This_element->raw_nodal_value(i_tm, i_nd, i_value)
                  * Dpsidx(i_nd, i_direc) * (*Ts_weights_pt)(0, i_tm);
              }
          }
      }
    return output;
  }

  Vector<double> interpolate_dvaluesdx_with_hanging_nodes(const unsigned &i_value) const
  {
    Vector<double> output(Dim, 0.0);

    for(unsigned i_direc=0; i_direc<Dim; i_direc++)
      {
        for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
          {
            for(unsigned i_tm=0; i_tm<Nprev_value_zeroth_derivative; i_tm++)
              {
                output[i_direc] += This_element->nodal_value(i_tm, i_nd, i_value)
                  * Dpsidx(i_nd, i_direc) * (*Ts_weights_pt)(0, i_tm);
              }
          }
      }
    return output;
  }

  // Interpolate derivatives of values w.r.t. position and time
  // ============================================================


  Vector<double> interpolate_d2valuesdxdt_raw(const unsigned &i_value) const
  {
    Vector<double> output(Dim, 0.0);
    for(unsigned i_direc=0; i_direc<Dim; i_direc++)
      {
        for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
          {
            for(unsigned i_tm=0; i_tm<Nprev_value_derivative; i_tm++)
              {
                output[i_direc] += This_element->raw_nodal_value(i_tm, i_nd, i_value)
                  * Dpsidx(i_nd, i_direc) * (*Ts_weights_pt)(1, i_tm);
              }
          }
      }
    return output;
  }

  Vector<double> interpolate_d2valuesdxdt_with_hanging_nodes(const unsigned &i_value) const
  {
    Vector<double> output(Dim, 0.0);

    for(unsigned i_direc=0; i_direc<Dim; i_direc++)
      {
        for(unsigned i_nd=0; i_nd<Nnode; i_nd++)
          {
            for(unsigned i_tm=0; i_tm<Nprev_value_derivative; i_tm++)
              {
                output[i_direc] += This_element->nodal_value(i_tm, i_nd, i_value)
                  * Dpsidx(i_nd, i_direc) * (*Ts_weights_pt)(1, i_tm);
              }
          }
      }
    return output;
  }


 };


/// Modified version of the itnerpolator for face elements: 1) can't get x
/// derivatives very easily here so throw errors if we try. 2) Don't try to
/// get dpsi etc (because we can't), just get shape and Jacobian separately.
class FaceElementInterpolator : public GeneralInterpolator
{

public:

 /// Constructor, jsut use underlying interpolator
 FaceElementInterpolator(const FaceElement* this_element,
                                const Vector<double> &s)
  : GeneralInterpolator(this_element, s) {}

 /// More complex constructor. Choose element to interpolate within,
 /// position in element to interpolate at and a time stepper to do time
 /// interpolation. WARNING: if you use this construtor you need to make
 /// sure that the previous values stored (in nodes) have the same meaning
 /// in both the given time stepper and the nodes' time steppers.
 FaceElementInterpolator(const FaceElement* this_element,
                         const Vector<double> &s,
                         const TimeStepper* ts_pt)
  : GeneralInterpolator(this_element, s, ts_pt) {}

 virtual void build()
 {
  // Set up shape + test functions
  J = This_element->J_eulerian(s());
  This_element->shape(s(), Psi);
  Test = Psi;

  // Find out if any nodes in this element are hanging
  Has_hanging_nodes = InterpolatorHelpers::has_hanging_nodes(This_element);

  // If paranoid check that all nodes have the same nvalues.
  InterpolatorHelpers::check_nvalues_in_element_same(This_element);
 }

 double dpsidx(const unsigned &i, const unsigned &i_direc) const
 {
  broken_function_error();
  return 0;
 }

 double dtestdx(const unsigned &i, const unsigned &i_direc) const
 {
  broken_function_error();
  return 0;
 }

private:

 void broken_function_error() const
  {
   std::string err = "Cannot calculate derivatives w.r.t. x in face elements";
   throw OomphLibError(err, OOMPH_EXCEPTION_LOCATION,
                       OOMPH_CURRENT_FUNCTION);
  }

 // /// Broken constructors
 // // FaceElementInterpolator() {}
 // FaceElementInterpolator(FaceElementInterpolator& d) {}
 // void operator=(FaceElementInterpolator& d) {}

};


} // End of oomph namespace

#endif
