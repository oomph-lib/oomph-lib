#ifndef MY_ERROR_ESTIMATOR_HEADER
#define MY_ERROR_ESTIMATOR_HEADER

#include "generic.h"

namespace oomph
{
  //========================================================================
  /// Constant error estimator, allows manual specification of refinement
  /// pattern by forcing refinement in regions defined by elements in
  /// a reference mesh.
  //========================================================================
  class ConstantErrorEstimator : public virtual ErrorEstimator
  {
  public:
    /// Constructor. Provide mesh and number of the elements that define
    /// the regions within which elements are to be refined subsequently.
    /// Also specify the node number of a central node
    /// within elements -- it's used to determine if an element is
    /// in the region where refinement is supposed to take place.
    /// Optional boolean flag (defaulting to false) indicates that
    /// refinement decision is based on Lagrangian coordinates -- only
    /// applicable to solid meshes.
    ConstantErrorEstimator(const double& error_value = 0.0)
      : Error_value(error_value)
    {
    }

    /// Broken copy constructor
    ConstantErrorEstimator(const ConstantErrorEstimator&) = delete;

    /// Broken assignment operator
    void operator=(const ConstantErrorEstimator&) = delete;

    /// Empty virtual destructor
    virtual ~ConstantErrorEstimator() {}

    /// Compute the elemental error measures for a given mesh
    /// and store them in a vector. Doc errors etc.
    virtual void get_element_errors(Mesh*& mesh_pt,
                                    Vector<double>& elemental_error,
                                    DocInfo& doc_info)
    {
#ifdef PARANOID
      if (doc_info.is_doc_enabled())
      {
        std::ostringstream warning_stream;
        warning_stream << "No output defined in "
                          "ConstantErrorEstimator::get_element_errors()\n"
                       << "Ignoring doc_info flag.\n";
        OomphLibWarning(warning_stream.str(),
                        "ConstantErrorEstimator::get_element_errors()",
                        OOMPH_EXCEPTION_LOCATION);
      }
#endif
      unsigned nelem = mesh_pt->nelement();
      for (unsigned e = 0; e < nelem; e++)
      {
        elemental_error[e] = Error_value;
      }
    }

  private:
    double Error_value;
  };


  //========================================================================
  /// My error estimator, allows manual specification of refinement
  /// pattern by forcing greater refinement for lower elements.
  //========================================================================
  class MyErrorEstimator : public virtual ErrorEstimator
  {
  public:
    /// Constructor. Provide mesh and number of the elements that define
    /// the regions within which elements are to be refined subsequently.
    /// Also specify the node number of a central node
    /// within elements -- it's used to determine if an element is
    /// in the region where refinement is supposed to take place.
    /// Optional boolean flag (defaulting to false) indicates that
    /// refinement decision is based on Lagrangian coordinates -- only
    /// applicable to solid meshes.
    MyErrorEstimator(const double& refinement_gradient)
      : Refinement_gradient(refinement_gradient)
    {
    }

    /// Broken copy constructor
    MyErrorEstimator(const MyErrorEstimator&) = delete;

    /// Broken assignment operator
    void operator=(const MyErrorEstimator&) = delete;

    /// Empty virtual destructor
    virtual ~MyErrorEstimator() {}


    /// Compute the elemental error measures for a given mesh
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
          << "No output defined in MyErrorEstimator::get_element_errors()\n"
          << "Ignoring doc_info flag.\n";
        OomphLibWarning(warning_stream.str(),
                        "MyErrorEstimator::get_element_errors()",
                        OOMPH_EXCEPTION_LOCATION);
      }
#endif
      unsigned nelem = mesh_pt->nelement();
      for (unsigned e = 0; e < nelem; e++)
      {
        elemental_error[e] = 0.0;

        // Check if element is in the regions to be refined
        // (based on coords of its central node)

        Vector<double> s(2);
        s[0] = 0.5;
        s[1] = 0.5;
        unsigned direction = 1;
        double y = mesh_pt->finite_element_pt(e)->interpolated_x(s, direction);
        double area = mesh_pt->finite_element_pt(e)->size();
        elemental_error[e] = Refinement_gradient * (2.0 - y) * area;
      }
    }

  private:
    double Refinement_gradient;
  };

  //========================================================================
  /// Contact line error estimator, allows manual specification of refinement
  /// pattern by forcing greater refinement for elements closer to the contact
  /// line.
  //========================================================================
  class ContactlineErrorEstimator : public virtual ErrorEstimator
  {
  public:
    /// Constructor. Provide mesh and number of the elements that define
    /// the regions within which elements are to be refined subsequently.
    /// Also specify the node number of a central node
    /// within elements -- it's used to determine if an element is
    /// in the region where refinement is supposed to take place.
    /// Optional boolean flag (defaulting to false) indicates that
    /// refinement decision is based on Lagrangian coordinates -- only
    /// applicable to solid meshes.
    ContactlineErrorEstimator(Node* contact_line_node_pt,
                              const double& min_element_length,
                              const double& element_length_ratio)
      : Min_element_length(min_element_length),
        Element_length_ratio(element_length_ratio)
    {
      Contact_line_solid_node_pt = contact_line_node_pt;
    }

    /// Broken copy constructor
    ContactlineErrorEstimator(const ContactlineErrorEstimator&) = delete;

    /// Broken assignment operator
    void operator=(const ContactlineErrorEstimator&) = delete;

    /// Empty virtual destructor
    virtual ~ContactlineErrorEstimator() {}


    /// Compute the elemental error measures for a given mesh
    /// and store them in a vector. Doc errors etc.
    virtual void get_element_errors(Mesh*& mesh_pt,
                                    Vector<double>& elemental_error,
                                    DocInfo& doc_info)
    {
#ifdef PARANOID
      if (doc_info.is_doc_enabled())
      {
        std::ostringstream warning_stream;
        warning_stream << "No output defined in "
                          "ContactlineErrorEstimator::get_element_errors()\n"
                       << "Ignoring doc_info flag.\n";
        OomphLibWarning(warning_stream.str(),
                        "ContactlineErrorEstimator::get_element_errors()",
                        OOMPH_EXCEPTION_LOCATION);
      }
#endif
      unsigned nelem = mesh_pt->nelement();
      for (unsigned e = 0; e < nelem; e++)
      {
        elemental_error[e] = 0.0;

        // Get coords of the centre of the element
        Vector<double> s(2);
        s[0] = 0.5;
        s[1] = 0.5;
        double x_center = mesh_pt->finite_element_pt(e)->interpolated_x(s, 0);
        double y_center = mesh_pt->finite_element_pt(e)->interpolated_x(s, 1);

        // Get coords of the contact line
        double x_contact_line =
          Contact_line_solid_node_pt->x(0);
        double y_contact_line =
          Contact_line_solid_node_pt->x(1);

        // Get the area
        double area = mesh_pt->finite_element_pt(e)->size();

        // Get distance from contact line
        double L = pow(pow(x_center - x_contact_line, 2.0) +
                         pow(y_center - y_contact_line, 2.0),
                       0.5);

        // Compute target area
        double target_area =
          0.5 * pow((Min_element_length - L * (1 - Element_length_ratio)) /
                      Element_length_ratio,
                    2.0);

        // Error = Max permitted error * area / target area
        // If the area is too large, the error will be large and the region
        // refined
        // If the area is too small, the error will be small and the region
        // unrefined
        elemental_error[e] = area / target_area;
      }
    }

  private:
    double Min_element_length;
    double Element_length_ratio;
    Node* Contact_line_solid_node_pt;
  };

  //========================================================================
  /// Corner error estimator, allows manual specification of refinement
  /// pattern by forcing greater refinement for elements closer to the
  /// corners.
  //========================================================================
  class CornerErrorEstimator : public virtual ErrorEstimator
  {
  public:
    /// Constructor. Provide mesh and number of the elements that define
    /// the regions within which elements are to be refined subsequently.
    /// Also specify the node number of a central node
    /// within elements -- it's used to determine if an element is
    /// in the region where refinement is supposed to take place.
    /// Optional boolean flag (defaulting to false) indicates that
    /// refinement decision is based on Lagrangian coordinates -- only
    /// applicable to solid meshes.
    CornerErrorEstimator(Node*& contact_line_node_pt,
                         Node*& inner_corner_node_pt,
                         double* const& min_element_length_pt,
                         double* const& inner_min_element_length_pt,
                         const double& element_length_ratio)
      : Min_element_length_pt(min_element_length_pt),
        Inner_min_element_length_pt(inner_min_element_length_pt),
        Element_length_ratio(element_length_ratio)
    {
      Contact_line_solid_node_pt = contact_line_node_pt;
      Inner_corner_node_pt = inner_corner_node_pt;
    }

    /// Broken copy constructor
    CornerErrorEstimator(const CornerErrorEstimator&) = delete;

    /// Broken assignment operator
    void operator=(const CornerErrorEstimator&) = delete;

    /// Empty virtual destructor
    virtual ~CornerErrorEstimator() {}


    /// Compute the elemental error measures for a given mesh
    /// and store them in a vector. Doc errors etc.
    virtual void get_element_errors(Mesh*& mesh_pt,
                                    Vector<double>& elemental_error,
                                    DocInfo& doc_info)
    {
#ifdef PARANOID
      if (doc_info.is_doc_enabled())
      {
        std::ostringstream warning_stream;
        warning_stream << "No output defined in "
                          "CornerErrorEstimator::get_element_errors()\n"
                       << "Ignoring doc_info flag.\n";
        OomphLibWarning(warning_stream.str(),
                        "CornerErrorEstimator::get_element_errors()",
                        OOMPH_EXCEPTION_LOCATION);
      }
#endif
      unsigned nelem = mesh_pt->nelement();
      for (unsigned e = 0; e < nelem; e++)
      {
        elemental_error[e] = 0.0;

        // Get coords of the centre of the element
        Vector<double> s(2);
        s[0] = 0.5;
        s[1] = 0.5;
        double x_center = mesh_pt->finite_element_pt(e)->interpolated_x(s, 0);
        double y_center = mesh_pt->finite_element_pt(e)->interpolated_x(s, 1);

        // Get coords of the contact line
        double x_contact_line =
          Contact_line_solid_node_pt->x(0);
        double y_contact_line =
          Contact_line_solid_node_pt->x(1);

        double x_inner = Inner_corner_node_pt->x(0);
        double y_inner = Inner_corner_node_pt->x(1);

        // Get the area
        double area = mesh_pt->finite_element_pt(e)->size();

        // Get distance from contact line
        double L1 = pow(pow(x_center - x_contact_line, 2.0) +
                          pow(y_center - y_contact_line, 2.0),
                        0.5);

        double L2 =
          pow(pow(x_center - x_inner, 2.0) + pow(y_center - y_inner, 2.0), 0.5);

        // Compute target area
        double target_area =
          0.5 *
          pow(std::min(*Min_element_length_pt - L1 * (1 - Element_length_ratio),
                       *Inner_min_element_length_pt -
                         L2 * (1 - Element_length_ratio)) /
                Element_length_ratio,
              2.0);

        // Error = Max permitted error * area / target area
        // If the area is too large, the error will be large and the region
        // refined
        // If the area is too small, the error will be small and the region
        // unrefined
        elemental_error[e] = area / target_area;
      }
    }

  private:
    double* Min_element_length_pt;
    double* Inner_min_element_length_pt;
    double Element_length_ratio;
    Node* Contact_line_solid_node_pt;
    Node* Inner_corner_node_pt;
  };

}; // namespace oomph
#endif
