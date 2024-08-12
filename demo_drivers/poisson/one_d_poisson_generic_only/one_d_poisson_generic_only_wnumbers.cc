1  //LIC// ====================================================================
2  //LIC// This file forms part of oomph-lib, the object-oriented, 
3  //LIC// multi-physics finite-element library, available 
4  //LIC// at http://www.oomph-lib.org.
5  //LIC// 
6  //LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
7  //LIC// 
8  //LIC// This library is free software; you can redistribute it and/or
9  //LIC// modify it under the terms of the GNU Lesser General Public
10 //LIC// License as published by the Free Software Foundation; either
11 //LIC// version 2.1 of the License, or (at your option) any later version.
12 //LIC// 
13 //LIC// This library is distributed in the hope that it will be useful,
14 //LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
15 //LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
16 //LIC// Lesser General Public License for more details.
17 //LIC// 
18 //LIC// You should have received a copy of the GNU Lesser General Public
19 //LIC// License along with this library; if not, write to the Free Software
20 //LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
21 //LIC// 02110-1301  USA.
22 //LIC// 
23 //LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
24 //LIC// 
25 //LIC//====================================================================
26 //Demonstration code to solve a one-dimensional Poisson equation
27 //using the OOMPH-LIB base classes only
28 
29 //Include std C++ IO library
30 #include <iostream>
31 
32 //Include the OOMPH-LIB base classes
33 #include "generic.h"
34 
35 using namespace std;
36 
37 using namespace oomph;
38 
39 //------------------------GEOMETRIC ELEMENT----------------------------------
40 
41 //===========================================================================
42 /// A "Geometric" Element consists of Nodes only and contains all the 
43 /// information necessary to find the value of one, or more, global 
44 /// coordinates at the nodes and at any point within the element. It should 
45 /// also specify a spatial integration scheme that can be used to 
46 /// calculate the size of the element. 
47 ///
48 /// A linear, one dimensional geometric element consists of two Nodes and 
49 /// uses linear interpolation between the Nodes. The integration scheme is 
50 /// a two-point Gaussian quadrature.
51 //===========================================================================
52 class TwoNodeGeometricElement : public FiniteElement
53 {
54 private:
55  /// Integration scheme that will be used to integrate over the element.
56  /// Simple Gaussian quadrature in one dimension, with two Gauss points.
57  static Gauss<1,2> Default_spatial_integration_scheme;
58 
59 public:
60 
61  /// Element constructor; resize Node_pt to the number of nodes (two)
62  /// and set the integration scheme. 
63  TwoNodeGeometricElement()
64   {
65    //Linear interpolation requires two Nodes per element.
66    //In fact, calling this function merely provides storage for 
67    //the pointers to the Nodes and initialises the pointers to NULL. 
68    //The Nodes themselves are created during the mesh generation 
69    //process by the functions FiniteElement::construct_node(...) 
70    //which stores the pointers to the newly created Nodes 
71    //in the element's own internal storage.
72    this->set_n_node(2);
73 
74    //The element is one-dimensional 
75    this->set_dimension(1);
76 
77    //Set the pointer to the spatial integration scheme
78    set_integration_scheme(&Default_spatial_integration_scheme);
79   }
80 
81  /// Define the shape functions for interpolation across the element.
82  void shape(const Vector<double> &s, Shape &psi) const
83   {
84    //There are two shape functions (one per node)
85    //In terms of the local coordinate, s[0]:
86    //Node 0 is at s[0] = -1.0
87    //Node 1 is at s[0] =  1.0
88    
89    //psi[0] takes the value one at node 0 and zero at node 1
90    psi[0] = 0.5*(1.0 - s[0]);
91 
92    //psi[1] takes the value one at node 1 and zero at node 0
93    psi[1] = 0.5*(1.0 + s[0]);
94   }
95 
96  /// Define the derivatives of the shape functions w.r.t. the local coordinate
97  void dshape_local(const Vector<double> &s, Shape &psi, DShape &dpsids) const
98   {
99    //Call the shape functions
100   shape(s,psi);
101
102   //The derivative of psi[0] wrt s[0] is -0.5
103   dpsids(0,0) = -0.5;
104
105   //The derivative of psi[0] wrt s[0] is 0.5
106   dpsids(1,0) = 0.5;
107  }
108 
109 /// The element consists of two nodes along its one (and only) edge
110 unsigned nnode_1d() const {return 2;}
111
112};
113
114//-----------------------SPECIFIC FINITE ELEMENT---------------------------
115
116//=========================================================================
117/// A specific Element needs to implement the functions that compute 
118/// the elemental contributions to the global residual vector and the 
119/// global Jacobian matrix formed when using Newton's method. 
120/// By using calls to the interpolation scheme defined in an
121/// an underlying "Geometric" Element, the equations can be written 
122/// generically and Elements that solve the same equations, 
123/// but using different orders of interpolation may be easily 
124/// constructed, using a different "geometric" base class.
125///
126/// A linear, Poisson element solves a one-dimensional Poisson equation
127/// using the integration and interpolation scheme defined in the
128/// above "Geometric" element.
129//=========================================================================
130class TwoNodePoissonElement : public TwoNodeGeometricElement
131 {
132
133 public:
134
135  /// Constructor: Initialise the sign of the source function to +1
136  TwoNodePoissonElement()
137   {
138    // Initialise sign of source function
139    Sign=1;
140   }
141
142  /// Define a specific source function. For greater generality, and
143  /// inclusion in a library, this could be defined as a source function 
144  /// pointer, that would then be set externally.
145  double f(const double &x) const 
146   {
147    return double(Sign)*30.0*sin(sqrt(30.0)*x);
148   }
149
150  /// Access function to the sign in the source function
151  int& sign() {return Sign;}
152
153  /// Define an access function to the first data value stored 
154  /// at each node. In a more general "Equation" element, such abstraction
155  /// is essential, because different Elements will store the same variables
156  /// in different locations.
157  double u(const unsigned &n) {return node_pt(n)->value(0);}
158
159  /// For the Poisson equation, only one value is stored at each node
160  unsigned required_nvalue(const unsigned &n) const {return 1;}
161
162  /// Calculate the elemental contributions to the global 
163  /// residual vector for the weak form of the Poisson equation
164  void get_residuals(Vector<double> &residuals)
165   {
166    //Find the number of degrees of freedom (unpinned values) in the element
167    unsigned n_dof = ndof();
168
169    //Initialise all the residuals to zero
170    for(unsigned i=0;i<n_dof;i++) {residuals[i] = 0.0;}
171    
172    //Find the number of nodes in the element
173    unsigned n_node = nnode();
174
175    //Allocate memory for shape functions and their derivatives:
176    // There's one shape function for each node:
177    Shape psi(n_node);
178
179    // Each of the n_node shape functions has one derivative with 
180    // respect to the single local coordinate:
181    DShape dpsidx(n_node,1);
182
183    //Storage for the single local coordinate
184    Vector<double> s(1);
185    
186    //Find the number of integration points in the underlying 
187    //geometric element's integration scheme 
188    unsigned n_intpt = integral_pt()->nweight();
189
190    //Loop over the integration points
191    for(unsigned ipt=0;ipt<n_intpt;ipt++)
192     {
193      //Set the value of the local coordinate to be the integration 
194      //scheme's knot point
195      s[0] = integral_pt()->knot(ipt,0);
196
197      //Find the weight of the integration scheme at this knot point
198      double w = integral_pt()->weight(ipt);
199
200      //Find the shape functions and their derivatives at the knot point. 
201      //This function is implemented in FiniteElement.
202      //It also returns the Jacobian of the mapping from local to 
203      //global coordinates.
204      double J = dshape_eulerian(s,psi,dpsidx);
205
206      //Premultiply the weight and the Jacobian
207      double W = w*J;
208      
209      //Allocate storage for the value of the field variable u,
210      //its derivative and the global position at the knot point.
211      //Initialise them all to zero.
212      double interpolated_x=0.0, interpolated_u=0.0, interpolated_dudx=0.0;
213
214      //Calculate the interpolated values by  looping over the shape 
215      //functions and summing the appropriate contributions
216      for(unsigned n=0;n<n_node;n++) 
217       {
218        interpolated_x += nodal_position(n,0)*psi[n];
219        interpolated_u += u(n)*psi[n];
220        interpolated_dudx += u(n)*dpsidx(n,0);
221       }
222      
223      // Evaluate the source function
224      double source=f(interpolated_x);
225
226      //ASSEMBLE THE RESIDUALS
227      
228      //Loop over the test functions (same as the shape functions
229      //since we're implementing an isoparametric element)
230      for(unsigned l=0;l<n_node;l++)
231       {
232        //Get the local equation number
233        //The variable is the first (only) value stored at the nodes
234        int local_eqn_number = nodal_local_eqn(l,0);
235
236        //If the equation is not a Dirichlet boundary condition
237        if(local_eqn_number >= 0)
238         {
239          //Add body force/source term here 
240          residuals[local_eqn_number] += source*psi[l]*W;
241
242          //Add the Poisson bit itself
243          residuals[local_eqn_number] += interpolated_dudx*dpsidx(l,0)*W;
244         }
245       }
246
247     } //End of loop over the integration points
248
249   } //End of function
250
251
252  /// Calculate the elemental contribution to the global residual
253  /// vector and to the Jacobian matrix dR_{i}/du_{j} used in the Newton method
254  void get_jacobian(Vector<double> &residuals, DenseMatrix<double> &jacobian)
255   {
256    //First, calculate the residuals
257    get_residuals(residuals);
258
259    //Find the number of degrees of freedom (unpinned values) in the element
260    unsigned n_dof = ndof();
261
262    //Initialise all entries of the Jacobian matrix to zero
263    for(unsigned i=0;i<n_dof;i++) 
264     {
265      for(unsigned j=0;j<n_dof;j++) {jacobian(i,j) = 0.0;}
266     }
267    
268    //Find the number of nodes in the element
269    unsigned n_node = nnode();
270    //Allocate memory for shape functions and their derivatives
271    Shape psi(n_node);
272    DShape dpsidx(n_node,1);
273
274    //Storage for the local coordinate
275    Vector<double> s(1);
276    
277    //Find the number of integration points in the underlying
278    //geometric element's integration scheme 
279    unsigned n_intpt = integral_pt()->nweight();
280
281    //Loop over the integration points
282    for(unsigned ipt=0;ipt<n_intpt;ipt++)
283     {
284      //Set the value of the local coordinate to be the integration 
285      //scheme's knot point
286      s[0] = integral_pt()->knot(ipt,0);
287
288      //Find the weight of the integration scheme at this knot point
289      double w = integral_pt()->weight(ipt);
290
291      //Find the shape functions and their derivatives at the knot point. 
292      //This function is implemented in FiniteElement.
293      //It also returns the Jacobian of the mapping from local to 
294      //global coordinates.
295      double J = dshape_eulerian(s,psi,dpsidx);
296
297      //Premultiply the weight and the Jacobian
298      double W = w*J;
299            
300      //ASSEMBLE THE JACOBIAN TERMS
301      
302      //Loop over the test (shape) functions
303      for(unsigned l=0;l<n_node;l++)
304       {
305        //Get the local equation number
306        //The variable is the first (only) value stored at the nodes
307        int local_eqn_number = nodal_local_eqn(l,0);
308
309        //If the equation is not a Dirichlet boundary condition
310        if(local_eqn_number >= 0)
311         {
312          //Loop over the degrees of freedom
313          for(unsigned l2=0;l2<n_node;l2++)
314           {
315            //Get the local degree of freedom number
316            //The variable is the first (only) value stored at the nodes
317            int local_dof_number = nodal_local_eqn(l2,0);
318
319            //If the degree of freedom is not pinned
320            if(local_dof_number >= 0)
321             {
322              //Add the contribution to the Jacobian
323              jacobian(local_eqn_number,local_dof_number) += 
324               dpsidx(l2,0)*dpsidx(l,0)*W;
325             }
326           }
327         }
328       }
329     } //End of loop over the integration points
330
331   } //End of function
332  
333  //Define an output function for the element 
334  void output(ostream &output) 
335   {
336    //Read out the number of nodes in the element   
337    unsigned n_node = nnode();
338
339    //Loop over the nodes and print out the global coordinate 
340    //and value of the field variable, u, at each node
341    for(unsigned n=0;n<n_node;n++)
342     {
343      output << nodal_position(n,0) << " " << u(n) << std::endl;
344     }
345   } //End of function
346
347
348
349  /// Self test function: The sign in the source function should
350  /// only have the values +/- 1. Following the general oomph-lib convention,
351  /// the self_test() returns 0 for success, and 1 for failure:
352  unsigned self_test()
353   {
354    // Initialise success flag
355    unsigned success=0;
356
357    // Run the generic FiniteElement self test
358    success=FiniteElement::self_test();
359
360    // Do additional test for this function
361    if ((Sign!=1)&&(Sign!=-1))
362     {
363      cout << "Sign of source function should be +/- 1," << std::endl;
364      cout << "but it is:" << Sign << std::endl;
365      success=1;
366     }
367
368    // Return success flag
369    return success;
370
371   } // End of self test
372
373
374 private:
375
376  /// The sign of the source function
377  int Sign;
378
379}; //End of the class
380
381/// Define the static spatial integration scheme
382Gauss<1,2> TwoNodeGeometricElement::Default_spatial_integration_scheme;
383
384//----------------------ONE DIMENSIONAL (LINE) MESH---------------------------
385
386//============================================================================
387/// A simple one dimensional mesh: uniformly spaced nodes in the domain x=0,1
388//============================================================================
389template<class ELEMENT>
390class OneDimMesh : public Mesh
391{
392
393public:
394
395 /// Mesh Constructor. The argument is the desired number of elements
396 OneDimMesh(const unsigned &n_element)
397 {
398  //Resize the vector of pointers to elements: there are n_element elements
399  Element_pt.resize(n_element); 
400
401  //Construct the first element (Note the use of the template parameter)
402  Element_pt[0] = new ELEMENT;
403
404  //Construct the first node and add it to the Mesh::Node_pt vector
405  //Note: The FiniteElement::construct_boundary_node(j) function
406  //builds the element's j-th local node, and provides the functionality
407  //that allows it to be located on a Mesh boundary -- essentially this
408  //involves allocating additional storage to the Node.
409  //The function obtains the Node's
410  //characteristics (e.g. its spatial dimension, the number of
411  //values to be stored, etc) from various virtual FiniteElement
412  //member functions, such as FiniteElement::required_nvalue(). 
413  //FiniteElement::construct_boundary_node(...) also
414  //stores a pointer to the newly created Node in the element's own
415  //Node_pt vector.
416  //Finally, the function returns the pointer to the
417  //newly created Node, so that it can be stored in the Mesh's Node_pt
418  //vector, as done here:
419  Node_pt.push_back(finite_element_pt(0)->construct_boundary_node(0));
420
421  //Find the number of nodes per element (N.B. all elements are identical
422  //so we can determine this value once and for all). 
423  unsigned n_node = finite_element_pt(0)->nnode();
424
425  //Loop over the remaning nodes of the first element
426  for(unsigned n=1;n<n_node;n++)
427   {
428    //Construct the next node and add it to the Mesh::Node_pt vector
429    //Note that these interior nodes need not (and should not)
430    //be boundary nodes, so they are created using the construct_node
431    //function, which has the same interface as
432    //construct_boundary_node()
433    Node_pt.push_back(finite_element_pt(0)->construct_node(n));
434   }
435
436  //Loop over the remaining elements apart from the last
437  for(unsigned e=1;e<(n_element-1);e++)
438   {
439    //Construct the e-th element
440    Element_pt[e] = new ELEMENT;
441
442    //The first local node of the e-th element is the last local node
443    //of the (e-1)-th element. We MUST NOT construct the node twice.
444    //Instead, we set the pointer in the e-th element to point to the
445    //previously created node in the (e-1)-th element.
446    finite_element_pt(e)->node_pt(0) = 
447      finite_element_pt(e-1)->node_pt(n_node-1);
448
449    //Loop over the remaining nodes of the e-th element
450    for(unsigned n=1;n<n_node;n++)
451     {
452      //Construct the next node and add it to the Mesh::Node_pt vector
453      //Note that these interior nodes need not (and should not)
454      //be boundary nodes, so they are created using the construct_node
455      //function, which has the same interface as
456      //construct_boundary_node()
457      Node_pt.push_back(finite_element_pt(e)->construct_node(n));
458     }
459   } //End of loop over elements
460  
461  
462  //Construct the final element
463  Element_pt[n_element-1] = new ELEMENT;
464  
465  //The first local node of the final element is the last local node
466  //of the penultimate element. We MUST NOT construct the node twice.
467  //Instead, we set the pointer in the final element to point to the
468  //previously created node in the penultimate element.
469  finite_element_pt(n_element-1)->node_pt(0) = 
470   finite_element_pt(n_element-2)->node_pt(n_node-1);
471
472  //Loop over the remaining central nodes of the final element
473  for(unsigned n=1;n<(n_node-1);n++)
474   {
475    //Construct the next node and add it to the Mesh::Node_pt vector
476    //Note that these interior nodes need not (and should not)
477    //be boundary nodes, so they are created using the construct_node
478    //function()
479    Node_pt.push_back(finite_element_pt(n_element-1)->construct_node(n));
480   }
481
482  //Construct the final node and add it to the Mesh::Node_pt vector.
483  //This node will be located on a boundary, and hence we use
484  //the construct_boundary_node function.
485  Node_pt.push_back(finite_element_pt(n_element-1)
486                    ->construct_boundary_node(n_node-1));
487
488  //We've now created all the nodes -- let's set their positions:
489
490  //Find the total number of nodes
491  unsigned n_global_node = nnode();
492
493  //Loop over all nodes
494  for(unsigned n=0;n<n_global_node;n++)
495   {
496    //Set the position of the node (equally spaced through the unit interval)
497    Node_pt[n]->x(0) = double(n)/double(n_global_node-1);
498   }
499 
500  //Set the boundary data:
501
502  //There are two boundaries in this mesh
503  set_nboundary(2);
504
505  //Boundary 0 contains the first node in the mesh:
506  add_boundary_node(0,Node_pt[0]);
507
508  //Boundary 1 contains the final node in the mesh:
509  add_boundary_node(1,Node_pt[n_global_node-1]); 
510
511 } // End of constructor
512
513}; // End of OneDimMesh class.
514
515//-----------------------PROBLEM CLASS-------------------------------------
516
517//==========================================================================
518/// Define the Problem which creates the mesh, applies the
519/// boundary conditions, and assigns equation numbers.
520//==========================================================================
521class DemoPoissonProblem : public Problem 
522 {
523  public:
524
525  /// Problem constructor: Pass the sign of the source function (default
526  /// is +1)
527  DemoPoissonProblem(const int& sign=1) : Sign(sign)
528   {
529    //Create a OneDimMesh Mesh object and set it to be the problem's mesh.
530    //The element type, TwoNodePoissonElement, is passed  as a template 
531    //parameter to the mesh. The argument to the constructor indicates
532    //the number of elements in the mesh.
533    Problem::mesh_pt() = new OneDimMesh<TwoNodePoissonElement>(10);
534
535    //Pin the unknowns at the ends of the 1D domain:
536    //The 1D mesh has 2 boundaries, each of which contains a single node;
537    //the nodes on the boundary are available from Mesh::boundary_node_pt(...)
538
539    //Pin the single nodal value at the single node on mesh 
540    //boundary 0 (= the left domain boundary at x=0)
541    mesh_pt()->boundary_node_pt(0,0)->pin(0);
542
543    //Pin the single nodal value at the single node on mesh 
544    //boundary 1 (= the right domain boundary at x=1)
545    mesh_pt()->boundary_node_pt(1,0)->pin(0);
546
547    // All values are initialised to zero. This is consistent with the
548    // boundary condition at x=0 and no further action is required
549    // at that node.
550
551    // Apply the boundary condition at x=1: u(x=1)=-/+1
552    mesh_pt()->boundary_node_pt(1,0)->set_value(0,-double(Sign));
553
554    
555    // Finish problem setup: Set the sign for the source function
556    // in all elements
557    
558    //Find number of elements in mesh
559    unsigned n_element = mesh_pt()->nelement();
560    
561    // Loop over the elements 
562    for(unsigned i=0;i<n_element;i++)
563     {
564      // The sign() member function is defined in the 
565      // TwoNodePoissonElement class not the base GeneralisedElement class.
566      // In order to use it, we must
567      // upcast from GeneralisedElement to the specific element type,
568      // which is achieved by a C++ dynamic_cast.
569      TwoNodePoissonElement *specific_element_pt 
570       = dynamic_cast<TwoNodePoissonElement*>(mesh_pt()->element_pt(i));
571      
572      // Set the sign of the source function
573      specific_element_pt->sign() = Sign;
574     }
575
576    //Assign the global and local equations numbers for the problem
577    cout << "Number of equations is " << assign_eqn_numbers() << std::endl;
578   }
579
580  
581  /// Check that everything has been set up properly
582  void actions_before_newton_solve()
583   {
584    if (0==self_test())
585     {
586      cout << "Problem has been set up correctly and can be solved." 
587           << std::endl;
588     }
589    else
590     {
591      throw 
592       OomphLibError("Trouble! Check error messages and fix the problems.\n",
593                     "DemoPoissonProblem::actions_before_newton_solve()",
594                     OOMPH_EXCEPTION_LOCATION);
595     }
596   }
597
598  
599  /// Print out the result after the solve
600  void actions_after_newton_solve() 
601   {
602    ofstream file ("result.dat");
603    mesh_pt()->output(file);
604   }
605
606 private:
607
608  /// The sign of the source function
609  int Sign;
610
611
612}; //End of problem definition
613
614//----------------------------MAIN FUNCTION-------------------------------
615
616int main()
617 {
618  //Build the problem 
619  DemoPoissonProblem problem;
620
621  //Solve the problem, using Newton's method
622  problem.newton_solve();
623
624 }