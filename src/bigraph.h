#ifndef bigraph_h_
#define bigraph_h_

#include "coord.h"

#include <cassert>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

enum GraphColor {
    GC_NONE = -1, 
    GC_WHITE = 0, 
    GC_GRAY, 
    GC_BLACK, 
    GC_RED, 
    GC_GREEN, 
    GC_BLUE 
};

class Vertex;
class BigraphVisitor;

//
// Edge
//
class Edge {
public:
    // The directions an edge can take.
    // Do not change the values
    enum Dir {
        ED_SENSE = 0,
        ED_ANTISENSE = 1,
        ED_COUNT = 2
    };

    static Dir EDGE_DIRECTIONS[2];

    // Flag indicating whether the sequences linked by an edge
    // are from the same strand or not.
    // Do not change the values
    enum Comp {
        EC_SAME = 0,
        EC_REVERSE = 1, 
        EC_COUNT
    };

    Edge(Vertex* end, Dir dir, Comp comp, const SeqCoord& coord) : _end(end), _dir(dir), _comp(comp), _coord(coord), _color(GC_NONE) {
    }

    Vertex* start() const {
        assert(_twin != NULL);
        return _twin->end();
    }
    Vertex* end() const {
        return _end;
    }
    Dir dir() const {
        return _dir;
    }
    Comp comp() const {
        return _comp;
    }
    void twin(Edge* twin) {
        _twin = twin;
    }
    Edge* twin() const {
        return _twin;
    }
    const SeqCoord& coord() const {
        return _coord;
    }
    SeqCoord& coord() {
        return _coord;
    }
    Match match() const {
        return Match(_coord, _twin->_coord, _comp == EC_REVERSE, 0);
    }
    void color(GraphColor color) {
        _color = color;
    }
    GraphColor color() const {
        return _color;
    }
    bool isSelf() const {
        return start() == end();
    }
    bool operator==(const Edge& edge) const;

    std::string label() const;
    void join(Edge* edge);
    void extend(Edge* edge);
    void validate() const;

private:
    Edge* _twin;
    Vertex* _end;
    Dir _dir;
    Comp _comp;
    SeqCoord _coord;
    GraphColor _color;
};

typedef std::vector<Edge *> EdgePtrList;

//
// Vertex
//
class Vertex {
public:
    typedef std::string Id;
    typedef std::unordered_map<std::string, size_t> IndexTable;

    Vertex(const Id& id, const std::string& seq, bool contained=false, const std::string& index="");
    ~Vertex();

    const Vertex::Id& id() const {
        return _id;
    }
    const std::string& seq() const {
        return _seq;
    }
    std::string index() const;
    size_t coverage() const {
        return _coverage;
    }

    // Merge another vertex into this vertex, as specified by pEdge
    void merge(Edge* edge);

    // Edge list operations
    void addEdge(Edge* edge);
    void removeEdge(Edge* edge);
    bool hasEdge(Edge* edge) const;
    void deleteEdges();
    size_t sweepEdges(GraphColor c);

    const EdgePtrList& edges() const {
        return _edges;
    }
    EdgePtrList edges(Edge::Dir dir) const {
        EdgePtrList ev;
        for (auto i = _edges.begin(); i != _edges.end(); ++i) {
            if ((*i)->dir() == dir) {
                ev.push_back(*i);
            }
        }
        return ev;
    }

    size_t degrees() const {
        return _edges.size();
    }
    size_t degrees(Edge::Dir dir) const {
        EdgePtrList ev = edges(dir);
        return ev.size();
    }
    void color(GraphColor c) {
        _color = c;
    }
    GraphColor color() const {
        return _color;
    }
    void contained(bool c) {
        _contained = c;
    }
    bool contained() const {
        return _contained;
    }

    void validate() const;
private:
    Id _id;
    GraphColor _color;
    std::string _seq;
    IndexTable _indexTbl;
    size_t _coverage; // Number of vertices that have been merged into this one
    bool _contained;

    EdgePtrList _edges;
};

typedef std::unordered_map<Vertex::Id, Vertex *> VertexTable;

//
// Bigraph
//
class Bigraph {
public:
    Bigraph(size_t n = 0) : _containment(false) {
        if (n > 0) {
            _vertices.reserve(n);
        }
    }
    ~Bigraph();

    bool addVertex(Vertex* vertex);
    Vertex* getVertex(const Vertex::Id& id) const;
    void removeVertex(Vertex* vertex);
    size_t sweepVertices(GraphColor c);

    void addEdge(Vertex* vertex, Edge* edge);
    size_t sweepEdges(GraphColor c);

    // Merge vertices that are joined by the specified edge
    void merge(Vertex* vertex, Edge* edge);

    // Simplify the graph by removing transitive edges
    void simplify();

    // Validate the graph is sane
    void validate() const;

    // Visit each vertex in the graph and call the visit functor object 
    bool visit(BigraphVisitor* vistor);

    void color(GraphColor c);
    void containment(bool c) {
        _containment = c;
    }
    bool containment() const {
        return _containment;
    }

    static bool load(std::istream& stream, size_t minOverlap, bool allowContainments, size_t maxEdges, Bigraph* g);
    static bool load(const std::string& filename, size_t minOverlap, bool allowContainments, size_t maxEdges, Bigraph* g);
    static bool save(std::ostream& stream, const Bigraph* g);
    static bool save(const std::string& filename, const Bigraph* g);
private:
    // Simplify the graph by compacting edges in the given direction
    void simplify(Edge::Dir dir);

    VertexTable _vertices;
    bool _containment;
};

// helper class
class EdgeCreator {
public:
    EdgeCreator(Bigraph* g, bool allowContainments, size_t maxEdges) : _graph(g), _allowContainments(allowContainments), _maxEdges(maxEdges) {
    }
    bool create(const Overlap& overlap, GraphColor color=GC_NONE);
private:
    Bigraph* _graph;
    bool _allowContainments;
    size_t _maxEdges;
};

#endif // bigraph_h_
