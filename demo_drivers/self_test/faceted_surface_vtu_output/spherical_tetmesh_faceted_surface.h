#ifndef OOMPH_SPHERICAL_TETMESH_FACETED_SURFACES_HEADER
#define OOMPH_SPHERICAL_TETMESH_FACETED_SURFACES_HEADER

namespace oomph
{

  //=============================================================
  /// TetMeshFacetedSurface that defines inner boundary
  //=============================================================
  class SphericalTetMeshFacetedSurface
    : public virtual TetMeshFacetedClosedSurface
  {
  public:
    /// Constructor
    SphericalTetMeshFacetedSurface()
    {
      // Golden ratio
      const double phi = 0.5 * (1.0 + sqrt(5.0));

      // Make vertices
      unsigned n_vertex = 12;
      Vertex_pt.resize(n_vertex);

      // Set basic icosahedron points
      Vector<double> icosa_point(3);

      icosa_point[0] = 0.0;
      icosa_point[1] = 1.0;
      icosa_point[2] = phi;
      Vertex_pt[0] = new TetMeshVertex(icosa_point);

      icosa_point[0] = 0.0;
      icosa_point[1] = -1.0;
      icosa_point[2] = phi;
      Vertex_pt[1] = new TetMeshVertex(icosa_point);

      icosa_point[0] = 0.0;
      icosa_point[1] = 1.0;
      icosa_point[2] = -phi;
      Vertex_pt[2] = new TetMeshVertex(icosa_point);

      icosa_point[0] = 0.0;
      icosa_point[1] = -1.0;
      icosa_point[2] = -phi;
      Vertex_pt[3] = new TetMeshVertex(icosa_point);

      icosa_point[0] = 1.0;
      icosa_point[1] = phi;
      icosa_point[2] = 0.0;
      Vertex_pt[4] = new TetMeshVertex(icosa_point);

      icosa_point[0] = -1.0;
      icosa_point[1] = phi;
      icosa_point[2] = 0.0;
      Vertex_pt[5] = new TetMeshVertex(icosa_point);

      icosa_point[0] = 1.0;
      icosa_point[1] = -phi;
      icosa_point[2] = 0.0;
      Vertex_pt[6] = new TetMeshVertex(icosa_point);

      icosa_point[0] = -1.0;
      icosa_point[1] = -phi;
      icosa_point[2] = 0.0;
      Vertex_pt[7] = new TetMeshVertex(icosa_point);

      icosa_point[0] = phi;
      icosa_point[1] = 0.0;
      icosa_point[2] = 1.0;
      Vertex_pt[8] = new TetMeshVertex(icosa_point);

      icosa_point[0] = phi;
      icosa_point[1] = 0.0;
      icosa_point[2] = -1.0;
      Vertex_pt[9] = new TetMeshVertex(icosa_point);

      icosa_point[0] = -phi;
      icosa_point[1] = 0.0;
      icosa_point[2] = 1.0;
      Vertex_pt[10] = new TetMeshVertex(icosa_point);

      icosa_point[0] = -phi;
      icosa_point[1] = 0.0;
      icosa_point[2] = -1.0;
      Vertex_pt[11] = new TetMeshVertex(icosa_point);

      // Make facets
      unsigned n_facet = 20;
      Facet_pt.resize(n_facet);

      unsigned n_vertex_on_facet = 3;
      Facet_pt[0] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[0]->set_vertex_pt(0, Vertex_pt[0]);
      Facet_pt[0]->set_vertex_pt(1, Vertex_pt[1]);
      Facet_pt[0]->set_vertex_pt(2, Vertex_pt[8]);

      Facet_pt[1] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[1]->set_vertex_pt(0, Vertex_pt[0]);
      Facet_pt[1]->set_vertex_pt(1, Vertex_pt[10]);
      Facet_pt[1]->set_vertex_pt(2, Vertex_pt[1]);

      Facet_pt[2] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[2]->set_vertex_pt(0, Vertex_pt[0]);
      Facet_pt[2]->set_vertex_pt(1, Vertex_pt[5]);
      Facet_pt[2]->set_vertex_pt(2, Vertex_pt[10]);

      Facet_pt[3] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[3]->set_vertex_pt(0, Vertex_pt[0]);
      Facet_pt[3]->set_vertex_pt(1, Vertex_pt[4]);
      Facet_pt[3]->set_vertex_pt(2, Vertex_pt[5]);

      Facet_pt[4] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[4]->set_vertex_pt(0, Vertex_pt[0]);
      Facet_pt[4]->set_vertex_pt(1, Vertex_pt[8]);
      Facet_pt[4]->set_vertex_pt(2, Vertex_pt[4]);

      Facet_pt[5] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[5]->set_vertex_pt(0, Vertex_pt[5]);
      Facet_pt[5]->set_vertex_pt(1, Vertex_pt[11]);
      Facet_pt[5]->set_vertex_pt(2, Vertex_pt[10]);

      Facet_pt[6] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[6]->set_vertex_pt(0, Vertex_pt[5]);
      Facet_pt[6]->set_vertex_pt(1, Vertex_pt[2]);
      Facet_pt[6]->set_vertex_pt(2, Vertex_pt[11]);

      Facet_pt[7] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[7]->set_vertex_pt(0, Vertex_pt[4]);
      Facet_pt[7]->set_vertex_pt(1, Vertex_pt[2]);
      Facet_pt[7]->set_vertex_pt(2, Vertex_pt[5]);

      Facet_pt[8] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[8]->set_vertex_pt(0, Vertex_pt[4]);
      Facet_pt[8]->set_vertex_pt(1, Vertex_pt[9]);
      Facet_pt[8]->set_vertex_pt(2, Vertex_pt[2]);

      Facet_pt[9] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[9]->set_vertex_pt(0, Vertex_pt[8]);
      Facet_pt[9]->set_vertex_pt(1, Vertex_pt[9]);
      Facet_pt[9]->set_vertex_pt(2, Vertex_pt[4]);

      Facet_pt[10] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[10]->set_vertex_pt(0, Vertex_pt[6]);
      Facet_pt[10]->set_vertex_pt(1, Vertex_pt[9]);
      Facet_pt[10]->set_vertex_pt(2, Vertex_pt[8]);

      Facet_pt[11] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[11]->set_vertex_pt(0, Vertex_pt[1]);
      Facet_pt[11]->set_vertex_pt(1, Vertex_pt[6]);
      Facet_pt[11]->set_vertex_pt(2, Vertex_pt[8]);

      Facet_pt[12] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[12]->set_vertex_pt(0, Vertex_pt[1]);
      Facet_pt[12]->set_vertex_pt(1, Vertex_pt[7]);
      Facet_pt[12]->set_vertex_pt(2, Vertex_pt[6]);

      Facet_pt[13] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[13]->set_vertex_pt(0, Vertex_pt[10]);
      Facet_pt[13]->set_vertex_pt(1, Vertex_pt[7]);
      Facet_pt[13]->set_vertex_pt(2, Vertex_pt[1]);

      Facet_pt[14] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[14]->set_vertex_pt(0, Vertex_pt[10]);
      Facet_pt[14]->set_vertex_pt(1, Vertex_pt[11]);
      Facet_pt[14]->set_vertex_pt(2, Vertex_pt[7]);

      Facet_pt[15] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[15]->set_vertex_pt(0, Vertex_pt[11]);
      Facet_pt[15]->set_vertex_pt(1, Vertex_pt[3]);
      Facet_pt[15]->set_vertex_pt(2, Vertex_pt[7]);

      Facet_pt[16] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[16]->set_vertex_pt(0, Vertex_pt[7]);
      Facet_pt[16]->set_vertex_pt(1, Vertex_pt[3]);
      Facet_pt[16]->set_vertex_pt(2, Vertex_pt[6]);

      Facet_pt[17] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[17]->set_vertex_pt(0, Vertex_pt[6]);
      Facet_pt[17]->set_vertex_pt(1, Vertex_pt[3]);
      Facet_pt[17]->set_vertex_pt(2, Vertex_pt[9]);

      Facet_pt[18] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[18]->set_vertex_pt(0, Vertex_pt[9]);
      Facet_pt[18]->set_vertex_pt(1, Vertex_pt[3]);
      Facet_pt[18]->set_vertex_pt(2, Vertex_pt[2]);

      Facet_pt[19] = new TetMeshFacet(n_vertex_on_facet);
      Facet_pt[19]->set_vertex_pt(0, Vertex_pt[2]);
      Facet_pt[19]->set_vertex_pt(1, Vertex_pt[3]);
      Facet_pt[19]->set_vertex_pt(2, Vertex_pt[11]);


      // Set one-based boundary IDs
      unsigned one_based_boundary_id = 1;
      for (unsigned f = 0; f < n_facet; f++)
      {
        Facet_pt[f]->set_one_based_boundary_id(one_based_boundary_id);
      }

      // Identify point in hole
      Vector<double> inner_point(3, 0.0);
      set_hole_for_tetgen(inner_point);
    }
  };

} // namespace oomph

#endif