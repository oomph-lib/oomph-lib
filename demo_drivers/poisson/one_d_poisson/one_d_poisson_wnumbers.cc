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
26 //kruemelmonster
27 //Driver for a simple 1D poisson problem
28 
29 // Generic oomph-lib routines
30 #include "generic.h"
31 
32 // Include Poisson elements/equations
33 #include "poisson.h"
34 
35 // Include the mesh
36 #include "meshes/one_d_mesh.h"
37 
38 using namespace std;
39 
40 using namespace oomph;
41 
42 //==start_of_namespace================================================
43 /// Namespace for fish-shaped solution of 1D Poisson equation
44 //====================================================================
45 namespace FishSolnOneDPoisson
46 {
47 
48  /// Sign of the source function 
49  /// (- gives the upper half of the fish, + the lower half)
50  int Sign=-1;
51 
52 
53  /// Exact, fish-shaped solution as a 1D vector
54  void get_exact_u(const Vector<double>& x, Vector<double>& u)
55  {
56   u[0] = double(Sign)*((sin(sqrt(30.0))-1.0)*x[0]-sin(sqrt(30.0)*x[0]));
57  }
58 
59 
60  /// Source function required to make the fish shape an exact solution 
61  void source_function(const Vector<double>& x, double& source)
62  {
63   source = double(Sign)*30.0*sin(sqrt(30.0)*x[0]);
64  }
65 
66 } // end of namespace
67 
68 
69 
70 
71 
72 
73 
74 //==start_of_problem_class============================================
75 /// 1D Poisson problem in unit interval.
76 //====================================================================
77 template<class ELEMENT> 
78 class OneDPoissonProblem : public Problem
79 {
80 
81 public:
82 
83  /// Constructor: Pass number of elements and pointer to source function
84  OneDPoissonProblem(const unsigned& n_element, 
85                     PoissonEquations<1>::PoissonSourceFctPt source_fct_pt);
86 
87  /// Destructor (empty)
88  ~OneDPoissonProblem()
89   {
90    delete mesh_pt();
91   }
92 
93  /// Update the problem specs before solve: (Re)set boundary conditions
94  void actions_before_newton_solve();
95 
96  /// Update the problem specs after solve (empty)
97  void actions_after_newton_solve(){}
98 
99  /// Doc the solution, pass the number of the case considered,
100 /// so that output files can be distinguished.
101 void doc_solution(const unsigned& label);
102
103private:
104
105 /// Pointer to source function
106 PoissonEquations<1>::PoissonSourceFctPt Source_fct_pt;
107
108}; // end of problem class
109
110
111
112
113
114//=====start_of_constructor===============================================
115/// Constructor for 1D Poisson problem in unit interval.
116/// Discretise the 1D domain with n_element elements of type ELEMENT.
117/// Specify function pointer to source function. 
118//========================================================================
119template<class ELEMENT>
120OneDPoissonProblem<ELEMENT>::OneDPoissonProblem(const unsigned& n_element,
121 PoissonEquations<1>::PoissonSourceFctPt source_fct_pt) : 
122 Source_fct_pt(source_fct_pt)
123{ 
124 Problem::Sparse_assembly_method = Perform_assembly_using_two_arrays;
125
126// Problem::Problem_is_nonlinear = false;
127 // Set domain length 
128 double L=1.0;
129
130 // Build mesh and store pointer in Problem
131 Problem::mesh_pt() = new OneDMesh<ELEMENT>(n_element,L);
132
133 // Set the boundary conditions for this problem: By default, all nodal
134 // values are free -- we only need to pin the ones that have 
135 // Dirichlet conditions. 
136
137 // Pin the single nodal value at the single node on mesh 
138 // boundary 0 (= the left domain boundary at x=0)
139 mesh_pt()->boundary_node_pt(0,0)->pin(0);
140 
141 // Pin the single nodal value at the single node on mesh 
142 // boundary 1 (= the right domain boundary at x=1)
143 mesh_pt()->boundary_node_pt(1,0)->pin(0);
144
145 // Complete the setup of the 1D Poisson problem:
146
147 // Loop over elements and set pointers to source function
148 for(unsigned i=0;i<n_element;i++)
149  {
150   // Upcast from GeneralisedElement to the present element
151   ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
152   
153   //Set the source function pointer
154   elem_pt->source_fct_pt() = Source_fct_pt;
155  }
156
157 // Setup equation numbering scheme
158 assign_eqn_numbers();
159
160} // end of constructor
161
162
163
164
165//===start_of_actions_before_newton_solve========================================
166/// Update the problem specs before solve: (Re)set boundary values
167/// from the exact solution. 
168//========================================================================
169template<class ELEMENT>
170void OneDPoissonProblem<ELEMENT>::actions_before_newton_solve()
171{
172 
173 // Assign boundary values for this problem by reading them out
174 // from the exact solution.
175
176 // Left boundary is node 0 in the mesh:
177 Node* left_node_pt=mesh_pt()->node_pt(0);
178
179 // Determine the position of the boundary node (the exact solution
180 // requires the coordinate in a 1D vector!)
181 Vector<double> x(1);
182 x[0]=left_node_pt->x(0);
183 
184 // Boundary value (read in from exact solution which returns
185 // the solution in a 1D vector)
186 Vector<double> u(1);
187 FishSolnOneDPoisson::get_exact_u(x,u);
188 
189 // Assign the boundary condition to one (and only) nodal value
190 left_node_pt->set_value(0,u[0]);
191
192
193 // Right boundary is last node in the mesh:
194 unsigned last_node=mesh_pt()->nnode()-1;
195 Node* right_node_pt=mesh_pt()->node_pt(last_node);
196
197 // Determine the position of the boundary node
198 x[0]=right_node_pt->x(0);
199 
200 // Boundary value (read in from exact solution which returns
201 // the solution in a 1D vector)
202 FishSolnOneDPoisson::get_exact_u(x,u);
203 
204 // Assign the boundary condition to one (and only) nodal value
205 right_node_pt->set_value(0,u[0]);
206
207 
208} // end of actions before solve
209
210
211
212//===start_of_doc=========================================================
213/// Doc the solution in tecplot format. Label files with label.
214//========================================================================
215template<class ELEMENT>
216void OneDPoissonProblem<ELEMENT>::doc_solution(const unsigned& label)
217{ 
218 using namespace StringConversion;
219
220 // Number of plot points
221 unsigned npts;
222 npts=5; 
223
224 // Output solution with specified number of plot points per element
225 ofstream solution_file(("soln" + to_string(label) + ".dat").c_str());
226 mesh_pt()->output(solution_file,npts);
227 solution_file.close();
228
229 // Output exact solution at much higher resolution (so we can
230 // see how well the solutions agree between nodal points)
231 ofstream exact_file(("exact_soln" + to_string(label) + ".dat").c_str());
232 mesh_pt()->output_fct(exact_file,20*npts,FishSolnOneDPoisson::get_exact_u); 
233 exact_file.close();
234
235 // Doc pointwise error and compute norm of error and of the solution
236 double error,norm;
237 ofstream error_file(("error" + to_string(label) + ".dat").c_str());
238 mesh_pt()->compute_error(error_file,FishSolnOneDPoisson::get_exact_u,
239                          error,norm); 
240 error_file.close();
241
242 // Doc error norm:
243 cout << "\nNorm of error    : " << sqrt(error) << std::endl; 
244 cout << "Norm of solution : " << sqrt(norm) << std::endl << std::endl;
245 cout << std::endl;
246
247} // end of doc
248
249 
250
251/// /////////////////////////////////////////////////////////////////////
252/// /////////////////////////////////////////////////////////////////////
253/// /////////////////////////////////////////////////////////////////////
254
255
256//======start_of_main==================================================
257/// Driver for 1D Poisson problem
258//=====================================================================
259int main()
260{
261
262 // Set up the problem: 
263 // Solve a 1D Poisson problem using a source function that generates
264 // a fish shaped exact solution
265 unsigned n_element=40; //Number of elements
266 OneDPoissonProblem<QPoissonElement<1,4> > //Element type as template parameter
267  problem(n_element,FishSolnOneDPoisson::source_function);
268
269 // Check whether the problem can be solved
270 cout << "\n\n\nProblem self-test ";
271 if (problem.self_test()==0)  
272  {
273   cout << "passed: Problem can be solved." << std::endl;
274  }
275 else 
276  {
277   throw OomphLibError("failed!",
278                       OOMPH_CURRENT_FUNCTION,OOMPH_EXCEPTION_LOCATION);
279  }
280
281 // Set the sign of the source function:
282 cout << "\n\n\nSolving with negative sign:\n" << std::endl;
283 FishSolnOneDPoisson::Sign=-1;
284
285 // Solve the problem with this Sign
286 problem.newton_solve();
287
288 //Output solution for this case (label output files with "0")
289 problem.doc_solution(0);
290
291
292 // Change the sign of the source function:
293 cout << "\n\n\nSolving with positive sign:\n" << std::endl;
294 FishSolnOneDPoisson::Sign=1;
295
296 // Re-solve the problem with this Sign (boundary conditions get
297 // updated automatically when Problem::actions_before_newton_solve() is
298 // called.
299 problem.newton_solve();
300
301 //Output solution for this case (label output files with "1")
302 problem.doc_solution(1);
303
304} // end of main
305
306
307
308
309
310
311
312
313
314