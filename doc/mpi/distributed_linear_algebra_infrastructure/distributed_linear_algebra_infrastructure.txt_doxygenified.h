/**

\mainpage Distributed Linear Algebra Infrastructure

In this document we provide an overview of \c oomph-lib's distributed linear
algebra framework, discussing its design and key functionality. 
The aim of this framework is to facilitate the parallel (distributed) 
execution of linear algebra type operations with little (or no) user 
intervention. This requires all linear algebra data and computations 
to be distributed over all available processes.


We begin by defining the \c OomphCommunicator, a class that is fundamental  
to distributed computing in \c oomph-lib. Next we discuss the class 
\c LinearAlgebraDistribution which specifies the  distribution of  the data 
and computations over the available processes. In
Sections \ref double_vector, \ref cr_double_matrix and \ref
distributed_linear_algebra_object we discuss \c oomph-lib's distributed linear
algebra objects including the key containers (matrices and vectors) 
and operators (solvers and preconditioners). Finally, we demonstrate how the 
distributed linear algebra framework is typically used in practice.


The primary aim of this document is to provide an overview of the design 
and functionality of distributed linear algebra capabilities in 
\c oomph-lib, and hence 
we do not discuss every method of every class; we refer the reader to the class
documentation for a complete specification.


\section oomph_communicator OomphCommunicator

\c Oomph-lib employs MPI for distributed memory
parallelisation. Fundamental to MPI is the communicator  
( \c MPI_Comm )
which determines which processes are involved in a parallel
computation. Although \c oomph-lib is implemented in C++, the C MPI
bindings are utilised. \c Oomph-lib provides the class \c
OomphCommunicator as an object-oriented wrapper for a \c
MPI_Comm.

Calling \c MPI_Helpers::init(argc,argv) is equivalent to calling \c
MPI_Init(argc,argv). It initialises MPI
( i.e. calls \c MPI_Init(argc,argv) ) and creates \c oomph-lib's global
communicator:

\dontinclude generic_mpi_test.cc
\skipline start_of_main
\until init(

 The newly created communicator is available via a \c MPI_Helpers
class static method:

\skipline Get the global 
\until communicator_pt

 and is equivalent to the \c MPI_COMM_WORLD
communicator in that it represents all the processes that \c oomph-lib
knows about. By default, this communicator contains exactly the same
set of processes as \c MPI_COMM_WORLD. 

The \c OomphCommunicator provides a number of access functions to 
communicator data including the rank of the process and the number of 
the processes:

\until std::endl;


\section linear_algebra_distribution LinearAlgebraDistribution


Distributed memory parallelisation requires data and computations to 
be distributed
over the available processes in some way. In this document we are
solely interested in distributed linear algebra. We choose to
distribute the linear algebra objects row-wise, and constrain the
distribution such that each process is associated with a single
contiguous set of rows. The class \c LinearAlgebraDistribution allows
the specification of such a distribution. 


The distribution is defined by two integers defining the first
global row and the number of local rows associated with a
process. This data is sufficient to allow a mapping between a global
row number and a local row number on a particular process.

To construct a \c LinearAlgebraDistribution in which 100 rows are
uniformly distributed across the set of processes specified by \c
comm_pt we write:

\until std::endl;

In this example, if run on four processes, the first 25 rows are associated 
with process 0, the next 25 rows are on process 1 and so on. In
general, in a uniform distribution of \f$n_r\f$ global rows over
\f$n_p\f$ processes the first row on process
\f$p=0,1,\ldots,(n_p-1)\f$ is \f$\lfloor
\frac{pn_r}{n_p}\rfloor\f$. It is also possible the specify
alternative user defined distributions; see the class documentation for details.


An optional third ( \c bool ) argument (default: true) in the
constructor indicates that we require a distributed linear algebra 
distribution. However, on some occasions we may want
to replicate all rows of a linear algebra object on all processes. This is
achieved by simply making the third argument \c false (non-distributed):

\until << std::endl;

This example illustrates two other features of \c
LinearAlgebraDistribution. Firstly, the default constructor creates
an empty distribution, and secondly for every (non-default)
constructor there is an equivalent \c build(...) method to 
"re-construct" the object.


The state of the object is accessible through a range of methods.

\until bool built

The \c built() method indicates if the object specifies a
distribution, or is empty.


\section double_vector DoubleVector

The simplest distributed linear algebra object is \c DoubleVector, a
distributed vector of \c doubles developed specifically for linear
algebra (It differs from a \c Vector<double> which simply provides a
container for \c doubles ). For example, the following command constructs a \c DoubleVector with a uniform distribution (specified by the distributed 
\c LinearAlgebraDistribution defined in the previous section) 
and unit elements:

\until my_vector

 
To access the vector elements the \c operator[] is implemented. 
For example to increment every element by one:

\until }

It is the \c oomph-lib convention that the data in \c DoubleVector (and all 
other distributed linear algebra object) is accessed using local indices. The 
following loop documents the local row number, the global row number, and the
value of the elements on each process:

\until }

To change the distribution of a \c DoubleVector while retaining the
data, we provide the \c redistribute(...) method. For example to
change \c my_vector from uniformly distributed to locally replicated:

\until redistribute

Just like the \c LinearAlgebraDistribution, we provide \c build() methods 
that mirror the behaviour of all non-default constructors. For example
to revert \c my_vector to a uniform distribution with unit elements:

\until distributed_distribution

It is important to differentiate between \c build(...) and \c
redistribute(...); calling \c build(...) deletes the existing data,
effectively re-constructing the object, whereas \c redistribute(...)
retains the vector's data.

Like the \c LinearAlgebraDistribution, a \c DoubleVector need not
contain any data. To generate an object in this state, we could 
instantiate an object using the default constructor or call the 
\c clear() method:


\until clear

Again the \c built() method returns the state of the object 
and indicates if it contains any data.


\section cr_double_matrix CRDoubleMatrix

\c CRDoubleMatrix is the only distributed matrix in \c oomph-lib. It
employs sparse compressed row storage to store \c double coefficients.

A \c CRDoubleMatrix has three fundamental states: 

- A \c CRDoubleMatrix can have no distribution or coefficients in
which case \c my_matrix->distribution_built() and
\c my_matrix->built() are both \c false.
- A (built) distribution but
no coefficients in which case \c my_matrix->distribution_built() is \c true but \c
 my_matrix->built() is still \c false.
- A (built) distribution and
 coefficients in which case \c my_matrix->distribution_built() and \c
 my_matrix->built() are both \c true.  

For example, to construct an empty matrix we call:

\until my_matrix

To specify the distribution as defined by the 
\c LinearAlgebraDistribution \c distributed_distribution we write:

\until my_matrix.build

The distribution has now been specified but the coefficients have not.
Like the \c DoubleVector, rows are indexed locally and hence the
coefficients rows must be indexed locally. For example, 
to populate \c my_matrix
as a square identity matrix, we write:

\until my_matrix.build

 We note that the column indices are global because only the rows are
 distributed. The assembly of a \c CRDoubleMatrix is now
 complete. 


We constructed the matrix in two stages by first specifying 
the distribution and then specifying the coefficients. However it is possible 
to perform this operation in just one step, by using the appropriate 
constructor or \c build(...) method, for example:

\until row_start);

\section distributed_linear_algebra_object DistributedLinearAlgebraObject


In this section we introduce the class \c
DistributedLinearAlgebraObject, 
a base class for all distributed linear algebra objects. 
This class encapsulates a \c LinearAlgebraDistribution, provides (protected)
access to derived classes to update ( \c build_distribution(...) ) and
clear ( \c clear() ) the stored distribution. Secondly, it provides
methods to simplify access to commonly used \c
LinearAlgebraDistribution data. For example, because a 
\c CRDoubleMatrix is a \c DistributedLinearAlgebraObject,

\until first_row();

can be replaced with

\until first_row();


 
\c DistributedLinearAlgebraObjects can be divided into two types:
containers and operators. We have already reviewed the containers \c
DoubleVector and \c CRDoubleMatrix. A wide range of operator classes
have been implemented in \c oomph-lib to operate on these
containers. In particular, all \c LinearSolvers, \c
IterativeLinearSolvers and \c Preconditioners (discussed in the
<a href="../../../linear_solvers/html/index.html">Linear Solvers
Tutorial</A>) are \c DistributedLinearAlgebraObjects. We finish this section 
by reviewing the key linear algebra operators:

- \c SuperLUSolver is a \c LinearSolver wrapper to both the 
<a href="http://crd.lbl.gov/~xiaoye/SuperLU/#superlu"><code>SuperLU</code></a>
  direct solver and the 
<a href="http://crd.lbl.gov/~xiaoye/SuperLU/#superlu_dist">
   <code>SuperLU Dist</code></a> distributed direct
  solver. By default, whenever possible this class will automatically
  perform distributed solves. 

- \c TrilinosAztecOOSolver is an \c IterativeLinearSolver wrapper to
  the  <a href="http://trilinos.sandia.gov/packages/aztecoo">
   <code>Trilinos AztecOO</code></a> package implementation of
  distributed Krylov methods including CG, GMRES and BiCGStab.


- \c TrilinosMLPreconditioner is a wrapper to the
  distributed <a href="http://trilinos.sandia.gov/packages/ml">
  Trilinos ML AMG preconditioners</a>.

- \c TrilinosIFPACKPreconditioner is a wrapper to the
  distributed  <a href="http://trilinos.sandia.gov/packages/ifpack">
  Trilinos IFPACK preconditioners</a>.

- \c HyprePreconditioner is a wrapper to the
  distributed <a href="https://computation.llnl.gov/casc/linear_solvers/sls_hypre.html">
  <code>Hypre Scalable Linear Solvers</code></a>
  package, of particular interest is the
  classical AMG implementation \c BoomerAMG.


- \c MatrixVectorProduct is a wrapper to the
  <a href="http://trilinos.sandia.gov/packages/epetra/">
  <code>Trilinos Epetra</code></a> distributed matrix-vector product 
  implementation. 


\section distributed_linear_algebra_in_practice Distributed Linear Algebra In Practice

Having discussed \c oomph-lib's linear algebra infrastructure,
we finally remark that \c oomph-lib is implemented such that
linear algebra in \c oomph-lib is automatically distributed if executed under MPI on multiple
processes. Specifically, a user should not need to specify either a \c
LinearAlgebraDistribution or a \c OomphCommunicator, unless they wish
to customise some aspect of the parallelisation.


All functionality is designed such that if a user does not specify a
\c LinearAlgebraDistribution, then as much data and computation as
possible will be uniformly distributed over all available processes.


As an example, we consider the \c Problem method \c get_jacobian(...). If the user does not specify a return distribution for the
Jacobian and residuals, then \c oomph-lib will uniformly distribute both containers.

\until std::endl;

On the other hand, a user can specify a return distribution by
 setting the distribution of the matrix
and vector prior to calling \c get_jacobian(...). 

\until std::endl;


We finally remark that because all linear algebra operations are 
automatically distributed, to parallelise \c oomph-lib's Newton solve 
phase, the user need only run their executable under MPI on multiple processes.



\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:\n\n
<CENTER>
<A HREF="../../../../self_test/mpi/generic_mpi/">
self_test/mpi/generic_mpi
</A>
</CENTER>\n
- The driver code is: \n\n
<CENTER>
<A HREF="../../../../self_test/mpi/generic_mpi/generic_mpi_test.cc">
self_test/mpi/generic_mpi/generic_mpi_test.cc
</A>
</CENTER>

<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/

