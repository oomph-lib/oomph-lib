//
// Created by iqraa on 20-1-24.
//

#ifndef READ_GMSH_H
#define READ_GMSH_H


#include "Vector.h"

#include <fstream>
#include <map>
#include <limits>
#include <unordered_map>

namespace GMSH
{

    static const double small = std::numeric_limits<double>::epsilon();
    static const double vSmall= std::numeric_limits<double>::min();

    struct Tuple
    {
        int tag = -1;
        std::vector<int> boundary;
        int dim = -1;
    };

    struct Boundary
    {
        int dim;
        int tag;
        std::map<std::string, int> physicalNames;
    };

    /*! \brief: This for point entities. */
    struct Point
    {
        int pointTag=-1, numPhysicalTags = -1;
        double x=0.0, y=0.0, z=0.0;
        std::vector<int> physicalTags = {-1};
    };

    /*! \brief: This for curve entities. */
    struct Curve
    {
        int curveTag=-1, numPhysicalTags=-1;
        std::vector<int> physicalTags={-1};
        int numBoundingPoints = -1;
        std::vector<int> pointTags ={-1,-1};
    };

    /*! \brief: This for surface entities. */
    struct Surface
    {
        int surfaceTag=-1, numPhysicalTags=-1;
        std::vector<int> physicalTags={-1};
        int numBoundingCurves=-1;
        std::vector<int> curveTags={-1,-1,-1,-1};
    };

    /*! \brief: This for volume entities. */
    struct Volume
    {
        int volumeTag=-1, numPhysicalTags=-1;
        std::vector<int> physicalTags={-1};
        int numBoundingSurfaces= 6;
        std::vector<int> surfaceTags = {-1,-1,-1,-1,-1,-1};
    };

    /*! \brief: This for the point entity of the entity mesh. */
    struct Vertex {
        int entityDim = -1, entityTag = -1, parametric_ = -1, numNodesInBlock = -1;
        int nodeTag = -1;
        double x = 0.0, y = 0.0, z = 0.0;
        std::vector<int> boundaries={-1};
    };

    /*! \brief: This for elements including points, edges, quads and hex elements. */
    struct Element
    {
        int entityDim = -1, entityTag = -1, elementType = -1, numElementsInBlock = -1;
        int elementTag = -1;
        int nodeTag;
        std::vector<int> edgeTags = {-1, -1};
        std::vector<int> quadTags = {-1, -1, -1, -1};
        std::vector<int> hexTags  = {-1, -1, -1, -1, -1, -1, -1, -1};
        std::vector<int> boundaries = {-1};
    };

    /*! \brief: This for actual node of the mesh. */
    struct Node
    {
        int entityDim = -1, nodeTag = -1, entityTag = -1;
        std::vector<double> coord {0.0,0.0,0.0};
        std::vector<int> boundaries = {0};
    };

    struct Edge
    {
        int entityDim = -1, entityTag = -1, edgeTag = -1;
        std::vector<int> nodesTag = {-1, -1};
        std::vector<int> boundaries = {0};
    };

    struct Quad
    {
        int entityDim = -1, entityTag = -1,  quadTag = -1;
        std::vector<int> nodesTag = {-1, -1, -1, -1};
        std::vector<int> edgesTag = {-1, -1, -1, -1};
        std::vector<int> boundaries={0};
    };

    struct Hexa
    {
        int entityDim = -1, entityTag = -1,  hexaTag = -1;
        std::vector<int> nodesTag = {-1, -1, -1, -1, -1, -1, -1, -1};
        std::vector<int> edgesTag = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
        std::vector<int> quadsTag = {-1, -1, -1, -1, -1, -1};
        std::vector<int> boundaries = {0};
    };

    template<typename T>
    std::vector<T> operator+(const std::vector<T>& lhs, const std::vector<T>& rhs) {
        if (lhs.size() != rhs.size()) {
            std::cerr << "Error: Vector sizes do not match for addition." << std::endl;
            return std::vector<T>(); // Return an empty vector
        }
        std::vector<T> result(lhs.size());
        for (size_t i = 0; i < lhs.size(); ++i) {
            result[i] = lhs[i] + rhs[i];
        }
        return result;
    }

    template<typename T>
    std::vector<T> operator-(const std::vector<T>& lhs, const std::vector<T>& rhs) {
        if (lhs.size() != rhs.size()) {
            std::cerr << "Error: Vector sizes do not match for subtraction." << std::endl;
            return std::vector<T>(); // Return an empty vector
        }
        std::vector<T> result(lhs.size());
        for (size_t i = 0; i < lhs.size(); ++i) {
            result[i] = lhs[i] - rhs[i];
        }
        return result;
    }

    template<typename T>
    std::vector<T> operator*(const std::vector<T>& vec, const T& scalar) {
        std::vector<T> result(vec.size());
        for (size_t i = 0; i < vec.size(); ++i) {
            result[i] = vec[i] * scalar;
        }
        return result;
    }

    template<typename T>
    std::vector<T> operator*(const T& scalar, const std::vector<T>& vec) {
        std::vector<T> result(vec.size());
        for (size_t i = 0; i < vec.size(); ++i) {
            result[i] = vec[i] * scalar;
        }
        return result;
    }

    template<typename T>
    std::vector<T> operator/(const std::vector<T>& vec, const T& scalar) {
        if (scalar == 0) {
            throw std::invalid_argument( "Error: Division by zero.");
        }
        std::vector<T> result(vec.size());
        for (size_t i = 0; i < vec.size(); ++i) {
            result[i] = vec[i] / scalar;
        }
        return result;
    }

    template<typename T>
    std::vector<T> operator/(const T& scalar, const std::vector<T>& vec) {
        if (scalar == 0) {
            throw std::invalid_argument( "Error: Division by zero.");
        }
        std::vector<T> result(vec.size());
        for (size_t i = 0; i < vec.size(); ++i) {
            result[i] = vec[i] / scalar;
        }
        return result;
    }

// Overloading compound assignment operators
    template<typename T>
    std::vector<T>& operator+=(std::vector<T>& lhs, const std::vector<T>& rhs) {
        if (lhs.size() != rhs.size()) {
            throw std::invalid_argument( "Error: Vector sizes do not match for addition.");
        }
        for (size_t i = 0; i < lhs.size(); ++i) {
            lhs[i] += rhs[i];
        }
        return lhs;
    }

    template<typename T>
    std::vector<T>& operator-=(std::vector<T>& lhs, const std::vector<T>& rhs) {
        if (lhs.size() != rhs.size()) {
            throw std::invalid_argument( "Error: Vector sizes do not match for subtraction." );
        }
        for (size_t i = 0; i < lhs.size(); ++i) {
            lhs[i] -= rhs[i];
        }
        return lhs;
    }

    template<typename T>
    std::vector<T>& operator/=(std::vector<T>& vec, const T& scalar) {
        if (scalar == 0) {
            throw std::invalid_argument("Error: Division by zero." );
        }
        for (size_t i = 0; i < vec.size(); ++i) {
            vec[i] /= scalar;
        }
        return vec;
    }


    template<typename T>
    std::vector<T> operator^(const std::vector<T>& lhs, const std::vector<T>& rhs)
    {
        if (lhs.size() != rhs.size()) {
            throw std::invalid_argument("Vectors must have the same dimensionality for cross product.");
        }

        std::vector<T> result(lhs.size());

        result[0] = lhs[1] * rhs[2] - lhs[2] * rhs[1];
        result[1] = lhs[2] * rhs[0] - lhs[0] * rhs[2];
        result[2] = lhs[0] * rhs[1] - lhs[1] * rhs[0];

        return result;
    }


    template<typename T>
    T operator&(const std::vector<T>& lhs, const std::vector<T>& rhs)
    {
        if (lhs.size() != rhs.size()) {
            throw std::invalid_argument("Vectors must have the same length for dot product.");
        }

        T result = 0;

        size_t n = lhs.size();

        for (size_t i = 0; i < n; ++i) {
            result += lhs[i] * rhs[i];
        }

        return result;
    }

    template<typename T>
    T mag(const std::vector<T>& vec)
    {
        return std::sqrt(vec & vec);
    }

    template<typename T>
    T mag(const std::vector<T>& veca, const std::vector<T>& vecb)
    {
        return std::sqrt(veca & vecb);
    }

    class Gmsh  {
    public:

        explicit Gmsh(const std::string &filename, bool verbose=false);

        /// No copy
        Gmsh(const Gmsh &) = delete;

        /// No move
        Gmsh(Gmsh &&) = delete;

        /// No move assignment
        Gmsh &operator=(Gmsh &&) = delete;

        /// Destructor
        ~Gmsh() = default;

        std::vector<Node>& getNodes();

        std::vector<Edge>& getEdges();

        std::vector<Quad>& getQuads();

        std::vector<Hexa>& getHexas();

        static std::vector<Quad> getFaces(const Hexa& hexa);

        /// Getters
        unsigned getBoundSize() const;

        unsigned getNodesSize() const;

        unsigned getEdgesSize() const;

        unsigned getQuadsSize() const;

        unsigned getHexasSize() const;

        unsigned getMeshDim () const;

        int getNodesPerElement() const;

        int getEdgesPerElement() const;

        int getFacesPerElement() const;

        /// \brief: INFO
        void printInfo(bool);

        void nodesInfo();

        void edgesInfo();

        void quadsInfo();

        void hexasInfo();

        int getBoundaryId(const std::string& bName);
    private:
        static bool inList(std::vector<int>& q, std::vector<Edge>& E);

        static bool inList(std::vector<int>& q, std::vector<Quad>& Q);

        /*! \brief: collect nodes */
        void collectNodes();

        /*! \brief: collect the edges into a vector of nodes tags*/
        bool collectEdges();

        /*! \brief: collect quadrilateral elements into vector of tags and ids*/
        bool collectQuads();

        /*! \brief: collect hexahedrals elements into vector of tags and ids*/
        bool collectHexas();

        void addMissingEdge(Quad& quad);

        void addMissingQuad(Hexa& hex);

        void kernel();

        void fillInTags();

        void addInternalEdge(std::vector<int>& line);

        void addInternalQuad(std::vector<int>& face);

        void renumber();

        static bool areSimilar(std::vector<int> a, std::vector<int> b);

        /*!@brief:
         * get the dimension and the tag of the boundary
         * @param [in] boundary
         * @param [out] DimsAndTags
         * */
        void entityTagToBoundaries();

        void renumberBCondition();

        bool correctOrientation( Quad &quad);

        bool correctOrientation(const Hexa& hex);

        static std::vector<Edge> getEdges(const Quad& quad);

        std::vector<double> getCenter(const Hexa& hex);

        std::vector<double> getCenter(const std::vector<int>& nodesTag);

        std::vector<double> getArea(const std::vector<int>& points);

    private:
        // flag to show all info
        bool verbose_ = false;

        std::ifstream mshFile;
        int minNodeTag = 0;
        std::vector<bool> found4;
        std::vector<bool> found2;
        std::vector<bool> found20;

        /// Entities
        std::vector<Point> pEntities;
        std::vector<Curve> cEntities;
        std::vector<Surface> sEntities;
        std::vector<Volume> vEntities;
        std::vector<Boundary> boundaries;


        /// Vertex, and Elements
        std::vector<Vertex> vertices;
        std::vector<Element> elements;

        ///
        std::vector<Node> nodes;
        std::vector<Edge> edges;
        std::vector<Quad> quads;
        std::vector<Hexa> hexas;
        std::vector<Quad> quadList;
        std::vector<Edge> edgeList;

        // internal added new quads and edges
        std::vector<std::vector<int>> fourNodesList;
        std::vector<std::vector<int>> commonInternalEdges;
        std::vector<std::vector<int>> twoNodesList;
        std::vector<std::vector<int>> internalQuads;

        /// Element types
        enum ElementType{
            Vertex_ = 15, // point
            Line_ = 1,    // edge
            Quadrilateral_ = 3, // face
            Hexahedron_ = 5 // brick
        };

        enum EntityType{
            Point_ = 0,
            Curve_ = 1,
            Surface_ = 2,
            Volume_ = 3
        };

        int nodesPerElement = -1;

        int edgesPerElement = -1;

        int facesPerElement = -1;

        //
        std::vector<Tuple> tuples;

        std::vector<int> renumberOld;
        std::unordered_map<int, int> orderedBC;

    };
}
#endif //GMSH_H
