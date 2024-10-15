#ifndef MY_EIGENPROBLEM_HEADER
#define MY_EIGENPROBLEM_HEADER

#include "generic.h"
#include "complex_less.h"

namespace oomph
{
  // Provides the functionality to solve the problem's eigenproblem.
  class MyEigenproblem : public Problem
  {
  private:
    // Private variables
    DocInfo Doc_info;
    // Provide storage for the latest eigensolve. Not full implemented yet.
    Vector<DoubleVector> eigenvector_real_storage;
    Vector<DoubleVector> eigenvector_imag_storage;

  public:
    MyEigenproblem() : Problem(), Doc_info("RESLT") {}
    DocInfo& doc_info()
    {
      return Doc_info;
    }

    // Solve the eigenproblem for the n least stable eigensolution
    // Note that the eigensolutions are not documented
    Vector<std::complex<double>> solve_n_most_unstable_eigensolutions(
      const unsigned& n_eval)
    {
      oomph_info << "solve_n_most_unstable_eigensolutions" << std::endl;
      Vector<std::complex<double>> eigenvalues;


#ifdef OOMPH_HAS_TRILINOS
      EigenSolver* old_eigen_solver_pt = this->eigen_solver_pt();
      ANASAZI new_eigen_solver;
      this->eigen_solver_pt() = &new_eigen_solver;
#endif


      this->solve_eigenproblem(n_eval,
                               eigenvalues,
                               eigenvector_real_storage,
                               eigenvector_imag_storage);

      Vector<std::complex<double>> sorted_eigenvalues = eigenvalues;
      std::sort(sorted_eigenvalues.rbegin(),
                sorted_eigenvalues.rend(),
                ComplexLess<double>());

      // Index array
      Vector<unsigned> eig_index(eigenvalues.size());
      // Populate with ascending integers, 0, 1, 2, 3, 4, ...
      std::iota(eig_index.begin(), eig_index.end(), 0);

      ComplexLess<double> cl;

      std::sort(eig_index.rbegin(),
                eig_index.rend(),
                [&eigenvalues, &cl](unsigned i1, unsigned i2)
                { return cl(eigenvalues[i1], eigenvalues[i2]); });

      for (unsigned i = 0; i < n_eval; i++)
      {
        // Print to screen
        printf("%i, %.17g + %.17g i\n",
               eig_index[i],
               sorted_eigenvalues[i].real(),
               sorted_eigenvalues[i].imag());
      }
      doc_eigenvalues(n_eval, sorted_eigenvalues);

#ifdef OOMPH_HAS_TRILINOS
      this->eigen_solver_pt() = old_eigen_solver_pt;
#endif

      return sorted_eigenvalues;
    }

    // Solve the eigenproblem for the n least stable eigensolution
    // Document them as we go as we do not store the solutions in the problem
    // dofs.
    Vector<std::complex<double>> solve_and_document_n_most_unstable_eigensolutions(
      const unsigned& n_eval)
    {
      oomph_info << "solve_and_document_n_most_unstable_eigensolutions"
                 << std::endl;

      Vector<std::complex<double>> eigenvalues;
      Vector<DoubleVector> eigenvector_real;
      Vector<DoubleVector> eigenvector_imag;

#ifdef OOMPH_HAS_TRILINOS
      EigenSolver* old_eigen_solver_pt = this->eigen_solver_pt();
      ANASAZI new_eigen_solver;
      this->eigen_solver_pt() = &new_eigen_solver;
#endif

      this->solve_eigenproblem(
        n_eval, eigenvalues, eigenvector_real, eigenvector_imag);

      Vector<std::complex<double>> sorted_eigenvalues = eigenvalues;
      std::sort(sorted_eigenvalues.rbegin(),
                sorted_eigenvalues.rend(),
                ComplexLess<double>());

      // Index array
      Vector<unsigned> eig_index(eigenvalues.size());
      // Populate with ascending integers, 0, 1, 2, 3, 4, ...
      std::iota(eig_index.begin(), eig_index.end(), 0);

      ComplexLess<double> cl;

      std::sort(eig_index.rbegin(),
                eig_index.rend(),
                [&eigenvalues, &cl](unsigned i1, unsigned i2)
                { return cl(eigenvalues[i1], eigenvalues[i2]); });

      Vector<Vector<double>> eigenvector_magnitude;
      Vector<Vector<double>> eigenvector_argument;
      Vector<Vector<double>> eigenvector_plot;
      Vector<DoubleVector> eigenvector_plot2;
      eigenvector_magnitude.resize(n_eval);
      eigenvector_argument.resize(n_eval);
      eigenvector_plot.resize(n_eval);
      eigenvector_plot2.resize(n_eval);
      const unsigned eigenvector_length = this->ndof();
      for (unsigned i = 0; i < n_eval; i++)
      {
        eigenvector_magnitude[i].resize(eigenvector_length);
        eigenvector_argument[i].resize(eigenvector_length);
        eigenvector_plot[i].resize(eigenvector_length);
        for (unsigned j = 0; j < eigenvector_length; j++)
        {
          eigenvector_magnitude[i][j] =
            pow(eigenvector_real[i][j] * eigenvector_real[i][j] +
                  eigenvector_imag[i][j] * eigenvector_imag[i][j],
                0.5);
          eigenvector_argument[i][j] =
            atan2(eigenvector_imag[i][j], eigenvector_real[i][j]);
          eigenvector_plot[i][j] =
            eigenvector_magnitude[i][j] * cos(eigenvector_argument[i][j]);
        }
        eigenvector_plot2[i] = eigenvector_real[i];
        for (unsigned j = 0; j < eigenvector_length; j++)
        {
          eigenvector_plot2[i][j] = eigenvector_plot[i][j];
        }
      }

      // Print to screen
      for (unsigned i = 0; i < n_eval; i++)
      {
        printf("%i, %.17g + %.17g i\n",
               eig_index[i],
               sorted_eigenvalues[i].real(),
               sorted_eigenvalues[i].imag());
      }


      // Open an output file for the sorted eigenvalues
      doc_eigenvalues(n_eval, sorted_eigenvalues);

      DoubleVector backup_dofs;
      this->get_dofs(backup_dofs);

      for (unsigned i = 0; i < n_eval; i++)
      {
        const unsigned n_amplidute = 1;
        for (unsigned j = 0; j < n_amplidute; j++)
        {
          this->add_eigenvector_to_dofs(std::pow(ndof(), 0.5) * 0.2,
                                        eigenvector_plot2[eig_index[i]]);
          // Output result
          this->doc_solution();

          this->set_dofs(backup_dofs);
        }
      }

      for (unsigned i = 0; i < n_eval; i++)
      {
        std::string filename =
          Doc_info.directory() + "/eigenvector" + to_string(i) + ".dat";
        std::ofstream output_stream(filename);
        for (unsigned j = 0; j < ndof(); j++)
        {
          output_stream << eigenvector_real[eig_index[i]][j] << " ";
          output_stream << eigenvector_imag[eig_index[i]][j] << std::endl;
        }
        output_stream.close();
      }

#ifdef OOMPH_HAS_TRILINOS
      this->eigen_solver_pt() = old_eigen_solver_pt;
#endif

      return sorted_eigenvalues;
    }

    // Perturb the current state by an eigensolution
    void perturb_with_eigensolution(const double& amplitude)
    {
      // Find eigensolution
      Vector<std::complex<double>> eigenvalues;
      Vector<DoubleVector> eigenvector_real;
      Vector<DoubleVector> eigenvector_imag;

#ifdef OOMPH_HAS_TRILINOS
      EigenSolver* old_eigen_solver_pt = this->eigen_solver_pt();
      ANASAZI new_eigen_solver;
      this->eigen_solver_pt() = &new_eigen_solver;
#endif

      const int n_eval = 8;
      this->solve_eigenproblem(
        n_eval, eigenvalues, eigenvector_real, eigenvector_imag);

      // Find the most unstable eigenvalue
      Vector<unsigned> eig_index(eigenvalues.size());
      // Populate with ascending integers, 0, 1, 2, 3, 4, ...
      std::iota(eig_index.begin(), eig_index.end(), 0);

      ComplexLess<double> cl;
      Vector<unsigned>::iterator max1 =
        std::max_element(eig_index.begin(),
                         eig_index.end(),
                         [&eigenvalues, &cl](unsigned i1, unsigned i2)
                         { return cl(eigenvalues[i1], eigenvalues[i2]); });

      std::cout << "max: " << *max1 << std::endl;
      std::cout << "Eigenvalue: " << eigenvalues[*max1] << std::endl;

      // Add eigensolution to the current state
      this->add_eigenvector_to_dofs(amplitude, eigenvector_real[*max1]);

#ifdef OOMPH_HAS_TRILINOS
      this->eigen_solver_pt() = old_eigen_solver_pt;
#endif
    }

    // Sanity check to ensure that anasazi is solving our problem.
    void check_first_eigensolution_solves_the_problem()
    {
      // Find eigensolution
      Vector<std::complex<double>> eigenvalues;
      Vector<DoubleVector> eigenvector_real;
      Vector<DoubleVector> eigenvector_imag;

#ifdef OOMPH_HAS_TRILINOS
      EigenSolver* old_eigen_solver_pt = this->eigen_solver_pt();
      ANASAZI new_eigen_solver;
      this->eigen_solver_pt() = &new_eigen_solver;
#endif

      const int n_eval = 8;
      this->solve_eigenproblem(
        n_eval, eigenvalues, eigenvector_real, eigenvector_imag);

      Vector<std::complex<double>> sorted_eigenvalues = eigenvalues;
      std::sort(sorted_eigenvalues.rbegin(),
                sorted_eigenvalues.rend(),
                ComplexLess<double>());

      // Index array
      Vector<unsigned> eig_index(eigenvalues.size());
      // Populate with ascending integers, 0, 1, 2, 3, 4, ...
      std::iota(eig_index.begin(), eig_index.end(), 0);

      ComplexLess<double> cl;

      std::sort(eig_index.rbegin(),
                eig_index.rend(),
                [&eigenvalues, &cl](unsigned i1, unsigned i2)
                { return cl(eigenvalues[i1], eigenvalues[i2]); });

      check_eigensolution(sorted_eigenvalues[0],
                          eigenvector_real[eig_index[0]]);

#ifdef OOMPH_HAS_TRILINOS
      this->eigen_solver_pt() = old_eigen_solver_pt;
#endif
    }

    // Broken doc solution
    virtual void doc_solution() = 0;

  private:
    // Helper function to document the eigenvalues
    void doc_eigenvalues(const unsigned& n_doc,
                         Vector<std::complex<double>> const& sorted_eigenvalues)
    {
      std::string filename;
      filename = Doc_info.directory() + "/eigenvalues" +
                 to_string(Doc_info.number()) + ".dat";

      // Open an output file for the sorted eigenvalues
      std::ofstream output_stream(filename);
      output_stream.precision(16);

      for (unsigned i = 0; i < n_doc; i++)
      {
        // Send to file
        output_stream << sorted_eigenvalues[i].real() << " "
                      << sorted_eigenvalues[i].imag() << std::endl;
      }
      output_stream.close();
    }

    // Helper function to check the eigensolution
    void check_eigensolution(const std::complex<double>& evalue,
                             const DoubleVector& evector)
    {
      DoubleVector backup_dofs;
      get_dofs(backup_dofs);

      CRDoubleMatrix mass_matrix;
      CRDoubleMatrix jacobian;
      get_eigenproblem_matrices(mass_matrix, jacobian);
      // Add forcing term to residuals
      const unsigned N = ndof();
      DoubleVector residuals;
      get_dofs(residuals);
      for (unsigned i = 0; i < N; i++)
      {
        residuals[i] = 0;
        for (unsigned j = 0; j < N; j++)
        {
          residuals[i] += jacobian(i, j) * evector[j] -
                          evalue.real() * mass_matrix(i, j) * evector[j];
        }
      }
      set_dofs(residuals);
      // doc_solution();

      // Output residuals
      std::ofstream output_stream("eigenmode-residuals.dat");
      residuals.output(output_stream, 16);
      output_stream.close();

      // Restore problem state
      this->set_dofs(backup_dofs);
    }
  };

} // namespace oomph
#endif
