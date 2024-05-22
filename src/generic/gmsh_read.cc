//
// Created by iqraa on 20-1-24.
//
#include "gmsh_read.h"


#include <iostream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <set>



/*! Constructor*/
GMSH::Gmsh::Gmsh(const std::string &filename, bool verbose):verbose_(verbose) {
    mshFile.open(filename);
    std::string line;

    if (mshFile.is_open()) {
        while (std::getline(mshFile, line)) {
            if (line == "$MeshFormat") {
                double MeshVersion;
                int suffix1, suffix2;
                std::getline(mshFile, line);
                std::istringstream(line) >> MeshVersion >> suffix1 >> suffix2;

                // Check msh file format
                if (MeshVersion != 4.1 | suffix1 != 0 | suffix2 != 8) {
                    std::cerr << "Error: Msh file format must be 4.1 0 8" << std::endl;
                    exit(EXIT_FAILURE);
                }
            } else if (line == "$PhysicalNames") {
                int nPhysical;
                std::getline(mshFile, line);  // Read the number of physical groups
                std::istringstream(line) >> nPhysical;

                for (int i = 0; i < nPhysical; ++i) {
                    Boundary boundary;

                    std::getline(mshFile, line);
                    std::istringstream iss(line);
                    int dim, tag;
                    std::string name;

                    iss >> dim >> tag;  // Read dimension and tag
                    iss.ignore();  // Ignore the space before the name
                    std::getline(iss, name, '\"');  // Read the name within quotes
                    std::getline(iss, name, '\"');

                    boundary.dim = dim;
                    boundary.tag = tag;
                    boundary.physicalNames[name] = tag;

                    // save it
                    boundaries.push_back(boundary);
                }
            } else if (line == "$Entities") {
                int numPEntity, numCEntity, numSEntity, numVEntity;

                std::getline(mshFile, line);
                std::istringstream(line) >> numPEntity >> numCEntity >> numSEntity >> numVEntity;

                /// point Entities
                pEntities.reserve(numPEntity);
                for (int i = 0; i < numPEntity; ++i) {
                    Point p;
                    std::getline(mshFile, line);
                    std::istringstream iss(line);
                    iss >> p.pointTag >> p.x >> p.y >> p.z >> p.numPhysicalTags;

                    // fill the physical tags
                    if (p.numPhysicalTags > 0) {
                        for (int j = 0; j < p.numPhysicalTags; ++j) {
                            int dummy;
                            iss >> dummy;
                            p.physicalTags[j] = dummy;
                        }
                    }

                    pEntities.push_back(p);
                }

                /// curve Entities
                cEntities.reserve(numCEntity);
                for (int i = 0; i < numCEntity; ++i) {
                    double tmp = 0.0;
                    Curve curve;
                    std::getline(mshFile, line);

                    std::istringstream iss(line);
                    iss >> curve.curveTag;
                    // skip
                    for (int j = 0; j < 6; ++j) {
                        iss >> tmp;
                    }
                    iss >> curve.numPhysicalTags;

                    // if numPhysicalTags > 0 assign the Tags
                    if (curve.numPhysicalTags > 0) {
                        int dummy;
                        for (int j = 0; j < curve.numPhysicalTags; ++j) {
                            iss >> dummy;
                            curve.physicalTags[j] = dummy;
                        }
                    }
                    // if numBoundingPoints > 0 assign point tags
                    iss >> curve.numBoundingPoints;
                    if (curve.numBoundingPoints > 0) {
                        int dummy;
                        for (int j = 0; j < curve.numBoundingPoints; ++j) {
                            iss >> dummy;
                            curve.pointTags[j] = dummy;
                        }
                    }
                    cEntities.push_back(curve);
                }

                /// Surface Entities
                sEntities.reserve(numSEntity);
                for (int i = 0; i < numSEntity; ++i) {
                    double tmp;
                    Surface surface;
                    std::getline(mshFile, line);
                    std::istringstream iss(line);
                    iss >> surface.surfaceTag;
                    // skip
                    for (int j = 0; j < 6; ++j) {
                        iss >> tmp;
                    }
                    iss >> surface.numPhysicalTags;

                    surface.physicalTags.resize(surface.numPhysicalTags, 0);
                    if (surface.numPhysicalTags > 0) {
                        for (int j = 0; j < surface.numPhysicalTags; ++j) {
                            int dummy;
                            iss >> dummy;
                            surface.physicalTags[j] = dummy;
                        }
                    }

                    iss >> surface.numBoundingCurves;
                    surface.curveTags.resize(surface.numBoundingCurves, 0);
                    if (surface.numBoundingCurves > 0) {

                        for (int j = 0; j < surface.numBoundingCurves; ++j) {
                            int dummy;
                            iss >> dummy;
                            surface.curveTags[j] = dummy;
                        }
                    }

                    sEntities.push_back(surface);
                }

                /// Volume Entities
                vEntities.reserve(numVEntity);
                for (int i = 0; i < numVEntity; ++i) {
                    double tmp;
                    Volume volume;
                    std::getline(mshFile, line);
                    std::istringstream iss(line);
                    iss >> volume.volumeTag;
                    // skip
                    for (int j = 0; j < 6; ++j) {
                        iss >> tmp;
                    }

                    iss >> volume.numPhysicalTags;
                    volume.physicalTags.reserve(volume.numPhysicalTags);
                    if (volume.numPhysicalTags > 0) {
                        int dummy;
                        for (int j = 0; j < volume.numPhysicalTags; ++j) {
                            iss >> dummy;
                            volume.physicalTags.push_back(dummy);
                        }
                    }

                    iss >> volume.numBoundingSurfaces;
                    volume.surfaceTags.reserve(volume.numBoundingSurfaces);
                    if (volume.numBoundingSurfaces > 0) {
                        int dummy;
                        for (int j = 0; j < volume.numBoundingSurfaces; ++j) {
                            iss >> dummy;
                            volume.surfaceTags.push_back(dummy);
                        }
                    }
                    vEntities.push_back(volume);
                }
            } else if (line == "$Nodes") {
                int tmp;
                int numEntityBlocks, numNodes, maxNodeTag;
                std::getline(mshFile, line);

                std::istringstream(line) >> numEntityBlocks >> numNodes >> minNodeTag >> maxNodeTag;
                if (numNodes > numEntityBlocks) vertices.resize(numNodes);
                else vertices.resize(numEntityBlocks);

                // initialize the nodes vector [for later use in collectNodes() not here]
                nodes.resize(numNodes);

                int numNodesInBlock;
                int count = 0;
                for (int i = 0; i < numEntityBlocks; ++i) {
                    std::getline(mshFile, line);

                    std::istringstream(line)
                            >> vertices[count].entityDim >> vertices[count].entityTag
                            >> vertices[count].parametric_ >> vertices[count].numNodesInBlock;

                    numNodesInBlock = vertices[count].numNodesInBlock;

                    // Now read nodes ids
                    for (int j = 0; j < numNodesInBlock; ++j) {
                        std::getline(mshFile, line);

                        std::istringstream(line) >> tmp;
                        vertices[count + j].nodeTag =
                                tmp - minNodeTag; // to reset numbering to zero use [tmp - minNodeTag];
                    }
                    // Read the coord.
                    for (int j = 0; j < numNodesInBlock; ++j) {
                        std::getline(mshFile, line);

                        std::istringstream(line) >> vertices[count + j].x >> vertices[count + j].y
                                                 >> vertices[count + j].z;
                        if (j > 0) {
                            vertices[count + j].entityDim = vertices[count + j - 1].entityDim;
                            vertices[count + j].entityTag = vertices[count + j - 1].entityTag;
                            vertices[count + j].parametric_ = vertices[count + j - 1].parametric_;
                            vertices[count + j].numNodesInBlock = vertices[count + j - 1].numNodesInBlock;
                        }
                    }

                    // reset count
                    if (numNodesInBlock == 0)
                        count += 1;
                    else
                        count += numNodesInBlock;
                }
            } else if (line == "$Elements") {
                int tmp;
                int numEntityBlocks, numElem, minElementTag, maxElementTag;
                std::getline(mshFile, line);

                std::istringstream
                        (line) >> numEntityBlocks >> numElem >> minElementTag >> maxElementTag;
                elements.resize(numElem);

                // counter to maintain the number of elements
                int count = 0;
                for (int i = 0; i < numEntityBlocks; ++i) {
                    std::getline(mshFile, line);

                    std::istringstream
                            (line) >> elements[count].entityDim >> elements[count].entityTag
                                   >> elements[count].elementType >> elements[count].numElementsInBlock;

                    // Read elem tag and node tag
                    for (int j = 0; j < elements[count].numElementsInBlock; ++j) {
                        if (j > 0) {
                            elements[count + j].entityDim = elements[count + j - 1].entityDim;
                            elements[count + j].entityTag = elements[count + j - 1].entityTag;
                            elements[count + j].elementType = elements[count + j - 1].elementType;
                            elements[count + j].numElementsInBlock = elements[count + j - 1].numElementsInBlock;
                        }
                        std::getline(mshFile, line);
                        std::istringstream iss(line);

                        // Read the element tag and make start from zero
                        iss >> tmp;
                        elements[count + j].elementTag = tmp - minElementTag;

                        // if it's a node
                        if (elements[count + j].elementType == ElementType::Vertex_) {
                            iss >> tmp;
                            //shift it to start from zero
                            elements[count + j].nodeTag = tmp - minNodeTag;
                        }

                        // if it's an edge/line
                        if (elements[count + j].elementType == ElementType::Line_) {
                            for (int k = 0; k < 2; ++k) {
                                iss >> tmp;
                                // shift to start from zero
                                elements[count + j].edgeTags[k] = tmp - minNodeTag;
                            }
                        }

                        // if it's quad
                        if (elements[count + j].elementType == ElementType::Quadrilateral_) {
                            for (int k = 0; k < 4; ++k) {
                                iss >> tmp;
                                // shift it to start from zero
                                elements[count + j].quadTags[k] = tmp - minNodeTag;
                            }
                        }

                        // if it's hex
                        if (elements[count + j].elementType == ElementType::Hexahedron_) {
                            for (int k = 0; k < 8; ++k) {
                                iss >> tmp;
                                // shift to start from zero
                                elements[count + j].hexTags[k] = tmp - minNodeTag;
                            }
                        }
                    }

                    count += elements[count].numElementsInBlock;
                }
            }
        }
        mshFile.close();
    } else {
        std::cerr << "Unable to open " << filename << std::endl;
    }

    /// Section of collecting already complete elements
    entityTagToBoundaries();

    renumberBCondition();

    // collect Nodes
    collectNodes();

    // Find internals nodes and create new [Quads and/or Hexas]
    kernel();

    fillInTags();
    /// renumbering B.C. section
    renumber();

    // mesh info
    printInfo(verbose_);

}

/*! Get boundaries from entities. */
void GMSH::Gmsh::entityTagToBoundaries() {
    for (auto &boundary: boundaries) {
        Tuple tuple;
        // get the id and the dimension of the BC from PhysicalNames
        if (boundary.dim == EntityType::Point_) {
            // loop over the entity of similar dim
            for (auto &pEntity: pEntities) {
                tuple = {};
                if (pEntity.numPhysicalTags != 0) {
                    tuple.tag = pEntity.pointTag;
                    tuple.boundary = pEntity.physicalTags;
                    tuple.dim = EntityType::Point_;
                    // save it
                    tuples.push_back(tuple);
                }
            }
        } else if (boundary.dim == EntityType::Curve_) {
            // loop over the entity of similar dim
            for (auto &cEntity: cEntities) {
                tuple = {};
                if (cEntity.numPhysicalTags != 0) {
                    tuple.tag = cEntity.curveTag;
                    tuple.boundary = cEntity.physicalTags;
                    tuple.dim = EntityType::Curve_;
                    // save it
                    tuples.push_back(tuple);
                }
            }
        } else if (boundary.dim == EntityType::Surface_) {
            // loop over the entity of similar dim
            for (auto &sEntity: sEntities) {
                tuple = {};
                if (sEntity.numPhysicalTags != 0) {
                    tuple.tag = sEntity.surfaceTag;
                    tuple.boundary = sEntity.physicalTags;
                    tuple.dim = EntityType::Surface_;
                    // save it
                    tuples.push_back(tuple);
                }
            }
        } else if (boundary.dim == EntityType::Volume_) {
            // loop over the entity of similar dim
            for (auto &vEntity: vEntities) {
                tuple = {};
                if (vEntity.numPhysicalTags != 0) {
                    tuple.tag = vEntity.volumeTag;
                    tuple.boundary = vEntity.physicalTags;
                    tuple.dim = EntityType::Volume_;
                    // save it
                    tuples.push_back(tuple);
                }
            }
        }
    }
}

/*!@brief renumberBCondition the boundaries. Gmsh gives the user the flexibility to
 * choose the numbering, however oomph-lib require the numbering starting from 0.
 * */
void GMSH::Gmsh::renumberBCondition() {
    for (auto &b: boundaries) {
        renumberOld.push_back(b.tag);
    }
    sort(renumberOld.begin(), renumberOld.end());

    for (int i = 0; i < renumberOld.size(); ++i) {
        orderedBC[renumberOld[i]] = i + 1; // delete the one if you want 0 based b.c.
    }
}

/*! Fill in [nodes] vector */
void GMSH::Gmsh::collectNodes() {
    for (auto &vertex: vertices) {
        if (vertex.nodeTag == -1) continue; // do not go down but increment i
        nodes[vertex.nodeTag].nodeTag = vertex.nodeTag;
        nodes[vertex.nodeTag].coord = {vertex.x, vertex.y, vertex.z};

        if (vertex.numNodesInBlock != 0) {
            nodes[vertex.nodeTag].entityDim = vertex.entityDim;
            nodes[vertex.nodeTag].entityTag = vertex.entityTag;
        }

        // assign BC
        for (auto &tuple: tuples)
        {
            if (tuple.tag == nodes[vertex.nodeTag].entityTag &&
                tuple.dim == EntityType::Surface_) {
                nodes[vertex.nodeTag].boundaries = tuple.boundary;
            }
        }
    }
}

/*! Fill in internals based on the element type. */
void GMSH::Gmsh::kernel()
{
    /// Gmsh always put the highest element type last in order.
    if (elements.back().elementType == ElementType::Quadrilateral_) {
        // Set nodesPerElement
        nodesPerElement = 4;

        // Set edgesPerElement
        edgesPerElement = 4;

        // Set nFacePerElement keep the initial val. [No faces for Quads]
        facesPerElement = -1;

        // collect Edges
        collectEdges();

        // collect Quads
        collectQuads();

        // Attempt to reorder quads nodes and add missing edges
        for (auto &quad: quads)
        {
            correctOrientation(quad);

            addMissingEdge(quad);
        }

    }

    // if the last element is hexahedral
    if (elements.back().elementType == ElementType::Hexahedron_)
    {
        // Set nodesPerElement
        nodesPerElement = 8;

        // Set edgesPerElement
        edgesPerElement = 12;

        // Set nFacePerElement
        facesPerElement = 6;

        // Collect Edges
        collectEdges();

        // Collect Quads
        collectQuads();

        // Collect Hexahedral
        collectHexas();

        /// Attempt to reorder hexas nodes and add missing quads/faces and edges
        for (auto &hex: hexas)
        {
            if (!correctOrientation(hex)) {
                // copy the hexa element
                Hexa h = hex;

                if (verbose_)
                    std::cout << "Inverting hex " << hex.hexaTag << std::endl;

                // Reorder hex.
                h.nodesTag[0] = hex.nodesTag[4];
                h.nodesTag[1] = hex.nodesTag[5];
                h.nodesTag[2] = hex.nodesTag[6];
                h.nodesTag[3] = hex.nodesTag[7];
                h.nodesTag[4] = hex.nodesTag[0];
                h.nodesTag[5] = hex.nodesTag[1];
                h.nodesTag[6] = hex.nodesTag[2];
                h.nodesTag[7] = hex.nodesTag[3];
                hex = h;
            }

            // Add quad if it is not in the quads list
            addMissingQuad(hex);
        }
    }
}

/*! Fill in [edges] vector */
bool GMSH::Gmsh::collectEdges()
{
    bool isEdges = false;
    int counter = 0;

    for (auto const &element: elements) {
        Edge edge;
        if (element.elementType == ElementType::Line_) {
            edge.entityDim = element.entityDim;
            edge.entityTag = element.entityTag;

            //edge.edgeTag = element.elementTag;
            edge.edgeTag = counter; counter++;
            edge.nodesTag = element.edgeTags;

            // assign BC
            for (auto &tuple: tuples) {
                if (tuple.tag == edge.entityTag && tuple.dim == EntityType::Curve_) {
                    edge.boundaries = tuple.boundary;

                    // and then, assign BC to nodes
                    for (auto &nodeTag: edge.nodesTag) {
                        nodes[nodeTag].boundaries = tuple.boundary;
                    }
                }
            }

            // save the edges
            edges.push_back(edge);

            // I found edges element
            isEdges = true;
        }
    }

    return isEdges;
}

/*! Fill in [quads] vector */
bool GMSH::Gmsh::collectQuads()
{
    bool isQuad = false;
    int counter = 0;
    // if elementType = 3 add it to the list of quads
    for (auto const &element: elements) {
        if (element.elementType == ElementType::Quadrilateral_) {
            Quad quad;
            quad.entityTag = element.entityTag;
            quad.entityDim = element.entityDim;

            //quad.quadTag  = element.elementTag;
            quad.quadTag = counter; counter++;
            quad.nodesTag = element.quadTags;

            // assign  the BC
            for (auto &tuple: tuples) {
                if (tuple.tag == quad.entityTag && tuple.dim == quad.entityDim) {
                    // assign  the BC to the current quad
                    quad.boundaries = tuple.boundary;

                    // and then, loop over the current quad nodes and give them the
                    // same BC of the quad
                    for (auto &nodeTag: quad.nodesTag) {
                        nodes[nodeTag].boundaries = tuple.boundary;
                    }
                }
            }

            // save the quad
            quads.push_back(quad);

            // I found quad element
            isQuad = true;
        }
    }

    return isQuad;
}

/*! Fill in [hexas] vector */
bool GMSH::Gmsh::collectHexas() {
    bool isHexas = false;
    int counter = 0;
    for (auto &element: elements) {
        Hexa hexa;

        if (element.elementType == ElementType::Hexahedron_) {
            hexa.entityDim = element.entityDim;
            hexa.entityTag = element.entityTag;

            // hexa.hexaTag = element.elementTag; // original gmsh file numbering
            hexa.hexaTag = counter; counter++;
            hexa.nodesTag = element.hexTags;
            hexa.boundaries = element.boundaries;

            // assign BC
            for (auto &tuple: tuples) {
                if (tuple.tag == hexa.entityTag && tuple.dim == hexa.entityDim) {
                    hexa.boundaries = tuple.boundary;

                    // and then, loop over the current quad nodes and give them the
                    // same BC of the quad
                    for (auto &nodeTag: hexa.nodesTag) {
                        nodes[nodeTag].boundaries = tuple.boundary;
                    }
                }
            }

            // save the element
            hexas.push_back(hexa);

            // I found hexas element
            isHexas = true;
        }
    }

    return isHexas;
}

void GMSH::Gmsh::addMissingQuad(Hexa& hex)
{
    std::vector<Quad> faces = getFaces(hex);
    for (auto &face: faces)
    {
        if (!inList(face.nodesTag, quads))
            addInternalQuad(face.nodesTag);

        // get edges of the current face/quad
        addMissingEdge(face);
    }
}

void GMSH::Gmsh::addMissingEdge(Quad& quad)
{
    std::vector<Edge> lines = getEdges(quad);
    for (auto &line: lines) {
        if (!inList(line.nodesTag, edges))
            addInternalEdge(line.nodesTag);
    }
}

void GMSH::Gmsh::renumber()
{
    // renumber nodes bc
    for (auto &node: nodes) {
        for (int &bound: node.boundaries) {
            if (bound != -1)
                bound = orderedBC[bound];
        }
    }

    // renumber quads bc
    for (auto &quad: quads) {
        for (int &bound: quad.boundaries) {
            if (bound != -1)
                bound = orderedBC[bound];
        }
    }

    /// Assign b.c. for edges and hexas
    // Assign b.c. to edges
    for (auto &edge: edges)
    {
        int n0 = edge.nodesTag[0];
        int n1 = edge.nodesTag[1];
        if (nodes[n0].boundaries == nodes[n1].boundaries)
        {
            edge.boundaries = nodes[n0].boundaries;
        }
    }

    // assign b.c. for hexas
    for (auto& hexa: hexas) {
        for (int qTag: hexa.quadsTag) {
            if (quads[qTag].boundaries[0]>0) hexa.boundaries[0] = quads[qTag].boundaries[0];
            if (hexa.boundaries[0] == -1) hexa.boundaries[0] = 0;
        }
    }
}

void GMSH::Gmsh::addInternalEdge(std::vector<int> &line) {
    // from edges get the last edge
    Edge e = edges.back();
    int tag = e.edgeTag;

    Edge edge;
    edge.edgeTag = tag + 1;
    edge.nodesTag = line;
    edges.emplace_back(edge);
}

void GMSH::Gmsh::addInternalQuad(std::vector<int> &face) {
    Quad q = quads.back();
    int tag = q.quadTag;

    Quad quad;
    quad.quadTag = tag + 1;
    quad.nodesTag = face;
    quads.push_back(quad);
}

/*! Function to check vectors of the same length if they are similar.
 * @param[in] a
 * @param[in] b
 * @param[out] boolean
 * */
bool GMSH::Gmsh::areSimilar(std::vector<int> a, std::vector<int> b) {
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    return a == b;
}

/*! @brief Check the orientation. */
bool GMSH::Gmsh::correctOrientation(Quad &quad) {
    bool ordered = false;
    int num_ccw = 0;
    for (int i = 0; i < quad.nodesTag.size(); ++i) {
        const int tag1 = quad.nodesTag[i];
        const int tag2 = quad.nodesTag[(i + 1) % 4];
        const int tag3 = quad.nodesTag[(i + 2) % 4];

        std::vector<double> node1 = nodes[tag1].coord;
        std::vector<double> node2 = nodes[tag2].coord;
        std::vector<double> node3 = nodes[tag3].coord;
        const double dotprod =
                (node2[0] - node1[0]) * (node3[1] - node2[1]) - (node2[1] - node1[1]) * (node3[0] - node2[0]);
        if (dotprod > 0.0) ++num_ccw;
    }

    if (num_ccw == 0) {
        // flip the orientation, first the first 4 nodes
        std::reverse(quad.nodesTag.begin(), quad.nodesTag.begin() + 3);
        // then the rest of the nodes
        std::reverse(quad.nodesTag.begin() + 4, quad.nodesTag.end());

        if (verbose_) std::cout << "Inverting quad " << quad.quadTag << std::endl;
    } else if (num_ccw == 4) {
        ordered = true;
    } else {
        ordered = false;
        std::cerr << "quad " << quad.quadTag << " is not convex - errors may emerge!" << std::endl;
    }
    return ordered;
}

/*! @brief Check the orientation. */
bool GMSH::Gmsh::correctOrientation(const Hexa &hex)
{
    const std::vector<double> cc = getCenter(hex);

    // Get outwards pointing faces.
    std::vector<Quad> faces = getFaces(hex);

    for (auto &face: faces) {
        std::vector<double> aa = getArea(face.nodesTag);
        std::vector<double> pp = nodes[face.nodesTag[0]].coord;
        //std::cout << " Dot Product= " << ((pp - cc) & aa) << std::endl;

        // Check if vector from any point on face to cc points outwards
        if (((pp - cc) & aa) < 0) {
            // Incorrectly oriented
            return false;
        }
    }

    return true;
}

/*! @brief get faces form current hexa. */
std::vector<GMSH::Quad> GMSH::Gmsh::getFaces(const Hexa &hexa) {
    std::vector<Quad> faces(6);
    faces[0].nodesTag = {hexa.nodesTag[3], hexa.nodesTag[2], hexa.nodesTag[1], hexa.nodesTag[0]};
    faces[1].nodesTag = {hexa.nodesTag[4], hexa.nodesTag[5], hexa.nodesTag[6], hexa.nodesTag[7]};
    faces[2].nodesTag = {hexa.nodesTag[0], hexa.nodesTag[1], hexa.nodesTag[5], hexa.nodesTag[4]};
    faces[3].nodesTag = {hexa.nodesTag[2], hexa.nodesTag[3], hexa.nodesTag[7], hexa.nodesTag[6]};
    faces[4].nodesTag = {hexa.nodesTag[0], hexa.nodesTag[4], hexa.nodesTag[7], hexa.nodesTag[3]};
    faces[5].nodesTag = {hexa.nodesTag[1], hexa.nodesTag[2], hexa.nodesTag[6], hexa.nodesTag[5]};
    return faces;
}

std::vector<GMSH::Edge> GMSH::Gmsh::getEdges(const Quad &quad) {
    std::vector<Edge> lines(4);

    lines[0].nodesTag = {quad.nodesTag[0], quad.nodesTag[1]};
    lines[1].nodesTag = {quad.nodesTag[1], quad.nodesTag[2]};
    lines[2].nodesTag = {quad.nodesTag[2], quad.nodesTag[3]};
    lines[3].nodesTag = {quad.nodesTag[3], quad.nodesTag[0]};
    return lines;
}

std::vector<double> GMSH::Gmsh::getCenter(const std::vector<int> &nodesTag) {
    std::vector<double> center{0.0, 0.0, 0.0};
    for (auto &nodeTag: nodesTag) {
        center += nodes[nodeTag].coord;
    }
    return center / double(nodesTag.size());
}

std::vector<double> GMSH::Gmsh::getCenter(const Hexa &hex) {
    // estimation of the center of the hexa element
    std::vector<double> apex(3, 0);
    apex = getCenter(hex.nodesTag);

    // Get Hexa faces
    std::vector<Quad> faces = getFaces(hex);

    // Calculate the centre by breaking the hexa into pyramids and
    // volume-weighted averaging their centres
    double sumV = 0.0;
    std::vector<double> sumVc(3, 0.0);

    // loop over faces
    for (auto &face: faces) {
        std::vector<double> fCenter = getCenter(face.nodesTag);

        std::vector<double> pyrHeight = apex - fCenter;
        std::vector<double> fArea = getArea(face.nodesTag);

        double pyrVol = std::abs((1.0 / 3.0) * (fArea & pyrHeight));

        std::vector<double> pyrCenter(3, 0);
        pyrCenter = (3.0 / 4.0) * fCenter + (1.0 / 4.0) * apex;

        if (pyrVol < small) {
            std::cerr << "zero or negative pyramid volume: " << -pyrVol
                      << " for hexa " << hex.hexaTag
                      << std::endl;
        }

        sumVc -= pyrVol * pyrCenter;
        sumV -= pyrVol;
    }
    return sumVc / (sumV + vSmall);
}

std::vector<double> GMSH::Gmsh::getArea(const std::vector<int> &points) {
    // If the face is a triangle, do a direct calculation
    if (points.size() == 3) {
        return 0.5 *
               ((nodes[points[1]].coord - nodes[points[0]].coord) ^ (nodes[points[2]].coord - nodes[points[0]].coord));
    }

    // For more complex faces, decompose into triangles ...

    // Compute an estimate of the centre as the average of the points
    std::vector<double> pAvg = {0.0, 0.0, 0.0};
    for (auto point: points) {
        pAvg += nodes[point].coord;
    }
    pAvg /= double(points.size());

    // Compute the face area normal and unit normal by summing up the
    // normals of the triangles formed by connecting each edge to the
    // point average.
    std::vector<double> sumA = {0., 0., 0.};
    for (int pi = 0; pi < points.size(); ++pi) {
        const std::vector<double> &p = nodes[points[pi]].coord;
        // Get the next point, considering cyclic indexing
        const std::vector<double> &pNext = nodes[points[(pi + 1) % points.size()]].coord;

        const std::vector<double> a = (pNext - p) ^ (pAvg - p);

        sumA += a;
    }

    return 0.5 * sumA;
}

void GMSH::Gmsh::fillInTags()
{
    for (auto& hexa: hexas)
    {
        /// filling quads tags for hexas
        int qCount = 0;
        // getting quads/faces for current hexa
        for (auto &face: getFaces(hexa))
        {
            for (auto &quad: quads)
            {
                if (areSimilar(face.nodesTag, quad.nodesTag))
                {
                    hexa.quadsTag[qCount] = quad.quadTag;
                    qCount++;
                }

                /// filling edges tags for quads
                int eCount = 0;
                for (auto &line: getEdges(quad))
                {
                    for (auto & edge: edges)
                    {
                        if (areSimilar(line.nodesTag, edge.nodesTag))
                        {
                            quad.edgesTag[eCount] = edge.edgeTag;
                            eCount++;
                        }
                    }
                }

            }
        }

        /// filling edges tags for hexas
        std::set<int> edgeSet;
        for (auto& quad: hexa.quadsTag)
        {
            for (auto& edge: quads[quad].edgesTag)
            {
                edgeSet.insert(edge);
            }
        }

        // if the size of the edgeSet is not 12 exit
        if (edgeSet.size()!= 12)
        {
            std::cout <<"number of edges per hex is wrong!"<<std::endl;
            exit(-1);
        }

        std::vector<int> lines(edgeSet.begin(), edgeSet.end());

        for (int i = 0; i < lines.size(); ++i)
        {
            hexa.edgesTag[i] = lines[i];
        }

    }
}

/// Getters
unsigned GMSH::Gmsh::getBoundSize() const { return boundaries.size(); }

unsigned GMSH::Gmsh::getNodesSize() const { return nodes.size(); }

unsigned GMSH::Gmsh::getEdgesSize() const { return edges.size(); }

unsigned GMSH::Gmsh::getQuadsSize() const { return quads.size(); }

unsigned GMSH::Gmsh::getHexasSize() const { return hexas.size(); }

unsigned GMSH::Gmsh::getMeshDim() const {
    unsigned dim = 3;
    if (elements.back().elementType == ElementType::Hexahedron_) { dim = 3; }

    if (elements.back().elementType == ElementType::Quadrilateral_) { dim = 2; }
    return dim;
}

int GMSH::Gmsh::getNodesPerElement() const { return nodesPerElement; }

int GMSH::Gmsh::getEdgesPerElement() const { return edgesPerElement; }

int GMSH::Gmsh::getFacesPerElement() const { return facesPerElement; }

std::vector<GMSH::Node> &GMSH::Gmsh::getNodes() { return nodes; }

std::vector<GMSH::Edge> &GMSH::Gmsh::getEdges() { return edges; }

std::vector<GMSH::Quad> &GMSH::Gmsh::getQuads() { return quads; }

std::vector<GMSH::Hexa> &GMSH::Gmsh::getHexas() { return hexas; }

void GMSH::Gmsh::printInfo(bool info = false) {
    int width = 2, precision = 6;
    std::cout << "GMSH Reader Info ..." << std::endl;
    std::cout << "   Statistics:" << std::endl;
    std::cout << "      # of Nodes: " << nodes.size() << std::endl;
    std::cout << "      # of Edges: " << edges.size() << std::endl;
    std::cout << "      # of Quads: " << quads.size() << std::endl;
    std::cout << "      # of Hexas: " << hexas.size() << std::endl;
    std::cout << "      # of bc   : " << boundaries.size() << std::endl;
    std::cout << "   Renumbering: " << std::endl;

    for (int i: renumberOld) {
        std::cout << "      bc: " << i << " ---> " << std::setw(width)
                  << std::setfill('0') << std::setprecision(precision)
                  << orderedBC[i] << std::endl;
    }

    // print all the info
    if (info) {
        nodesInfo();
        edgesInfo();
        quadsInfo();
        hexasInfo();
    }
}

void GMSH::Gmsh::nodesInfo() {
    int width = 2, precision = 6;
    int i = 0;
    std::cout << "Nodes: \n";
    for (auto &node: nodes) {
        std::cout << " tag: "
                  << std::setw(width) << std::setprecision(precision) << node.nodeTag;


        std::cout << " BC: ";
        for (int boundary: node.boundaries) {
            std::cout << std::setw(2) << std::setfill('0')
                      << boundary << " ";
        }


        //std::cout << " E.Tag: " << std::setw(width) << std::setfill('0')<<node.entityTag;
        //std::cout << " E.Dim: " << std::setw(width) << std::setfill('0')<<node.entityDim;

        std::cout << " Coord: "
                  << "(" << std::setw(width) << std::setprecision(precision) << node.coord[0]
                  << " " << std::setw(width) << std::setprecision(precision) << node.coord[1]
                  << " " << std::setw(width) << std::setprecision(precision) << node.coord[2];
        std::cout << ")";

        std::cout << std::endl;
        i += 1;
    }
    std::cout << std::endl;
}

void GMSH::Gmsh::edgesInfo() {
    int i = 0;
    std::cout << "Edges: \n";
    for (auto &edge: edges) {
        std::cout << " tag: " << std::setw(2) << std::setfill('0') << edge.edgeTag;

        std::cout << " BC: ";
        for (int boundary: edge.boundaries) {
            std::cout << std::setw(2) << std::setfill('0')
                      << boundary << " ";
        }

        std::cout << " Nodes: ( ";
        for (int node: edge.nodesTag) {
            std::cout << std::setw(2) << std::setfill('0')
                      << node << " ";
        }
        std::cout << ")";

        std::cout << std::endl;
        i += 1;
    }
    std::cout << std::endl;
}

void GMSH::Gmsh::quadsInfo() {
    int i = 0;
    std::cout << "Quads: \n";
    for (auto &quad: quads) {
        std::cout << " tag: " << std::setw(2) << std::setfill('0') << quad.quadTag;

        std::cout << " BC: ";
        for (int boundary: quad.boundaries) {
            std::cout << std::setw(2) << std::setfill('0')
                      << boundary << " ";
        }

        std::cout <<" Edges: (" ;
        for (int edge : quad.edgesTag) {
            std::cout   << std::setw(2) << std::setfill('0')
                        << edge << " ";
        }
        std::cout <<")";

        std::cout << " Nodes: (";
        for (int node: quad.nodesTag) {
            std::cout << std::setw(2) << std::setfill('0')
                      << node << " ";
        }
        std::cout << ")" << std::endl;
        i += 1;
    }
    std::cout << std::endl;
}

void GMSH::Gmsh::hexasInfo() {
    if (hexas.empty()) {
        std::cout << "It is 2D mesh!";
        return;
    }

    int i = 0;
    std::cout << "Hexas: \n";
    for (auto &hexa: hexas) {
        std::cout <<" tag: "
                  << std::setw(2) << std::setfill('0') << hexa.hexaTag;

        std::cout << " BC: ";
        for (auto &boundary: hexa.boundaries) {
            std::cout << std::setw(2) << std::setfill('0')
                      << boundary << " ";
        }

        std::cout << " Faces: (";
        for (auto &quadTag: hexa.quadsTag) {
            std::cout << std::setw(2) << std::setfill('0')
                      << quadTag << " ";
        }
        std::cout << ")";

        std::cout << " Edges: (";
        for (int edgeTag: hexa.edgesTag) {
            std::cout << std::setw(2) << std::setfill('0')
                      << edgeTag << " ";
        }
        std::cout << ")";

        std::cout << " Nodes: (";
        for (int node: hexa.nodesTag) {
            std::cout << std::setw(2) << std::setfill('0')
                      << node << " ";
        }
        std::cout << ")" << std::endl;
        i += 1;
    }
    std::cout << std::endl;
}

/// Utilities
bool GMSH::Gmsh::inList(std::vector<int> &q, std::vector<Edge> &E) {
    // loop over the quads and see if face is there
    for (auto &e: E) {
        if (areSimilar(q, e.nodesTag)) return true;
    }
    return false;
}


bool GMSH::Gmsh::inList(std::vector<int> &q, std::vector<Quad> &Q) {
    // loop over the quads and see if face is there
    for (auto &quad: Q) {
        if (areSimilar(q, quad.nodesTag)) return true;
    }
    return false;
}