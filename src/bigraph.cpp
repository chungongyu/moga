#include "bigraph.h"
#include "asqg.h"
#include "kseq.h"

#include <fstream>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <log4cxx/logger.h>

#define GZIP_EXT ".gz"

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.Bigraph"));

//
// Edge
//
std::string Edge::label() const {
    // get the unmatched coordinates in end vertex
    const SeqCoord& coord = _twin->coord();
    SeqCoord unmatched = coord.complement();
    const std::string& seq = _end->seq();
    std::string label = seq.substr(unmatched.interval.start, unmatched.length());
    if (comp() == EC_REVERSE) {
        make_reverse_complement_dna(label);
    }
    return label;
}

//
// Vertex
//
Vertex::~Vertex() {
    for (EdgePtrList::iterator i = _edges.begin(); i != _edges.end(); ++i) {
        delete *i;
    }
}

void Vertex::merge(Edge* edge) {
    // Merging two string vertices has two parts
    // First, the sequence of the vertex is extended
    // by the the content of the edge label
    // Then, all the edges that are pointing to this node
    // must be updated to contain the extension of the vertex
    Edge* twin = edge->twin();

    // Merge the sequence
    std::string label = edge->label();
}

void Vertex::addEdge(Edge* edge) {
    assert(edge->start() == this);
    _edges.push_back(edge);
}

//
// Bigraph
//
Bigraph::~Bigraph() {
    for (VertexTable::iterator i = _vertices.begin(); i != _vertices.end(); ++i) {
        delete i->second;
    }
}

bool Bigraph::addVertex(Vertex* vertex) {
    VertexTable::iterator i = _vertices.find(vertex->id());
    if (i != _vertices.end()) {
        return false;
    }
    _vertices[vertex->id()] = vertex;
    return true;
}

Vertex* Bigraph::getVertex(const Vertex::Id& id) const {
    VertexTable::const_iterator i = _vertices.find(id);
    if (i != _vertices.end()) {
        return i->second;
    }
    return NULL;
}

void Bigraph::addEdge(Vertex* vertex, Edge* edge) {
    vertex->addEdge(edge);
}

void Bigraph::simplify() {
    simplify(Edge::ED_SENSE);
    simplify(Edge::ED_ANTISENSE);
}

void Bigraph::merge(Vertex* start, Edge* edge) {
    Vertex* end = edge->end();

    // Merge the data
    start->merge(edge);
}

void Bigraph::simplify(Edge::Dir dir) {
    bool changed = true;
    while (changed) {
        changed = false;
        for (VertexTable::iterator i = _vertices.begin(); i != _vertices.end(); ++i) {
            // Get the edges for this direction
            EdgePtrList edges = i->second->edges(dir);

            // If there is a single edge in this direction, merge the vertices
            // Don't merge singular self edges though
            if (edges.size() == 1 && !edges[0]->isSelf()) {
                // Check that the edge back is singular as well
                Edge* single = edges[0];
                Edge* twin = single->twin();
                Vertex* end = single->end();
                if (end->degrees(twin->dir()) == 1) {
                    merge(i->second, single);
                    changed = true;
                }
            }
        }
    }
}

class EdgeCreator {
public:
    EdgeCreator(Bigraph* g, bool allowContainments, size_t maxEdges) : _graph(g), _allowContainments(allowContainments), _maxEdges(maxEdges) {
    }

    bool create(const Overlap& overlap) {
        // Initialize data and perform checks
        bool isContainment = overlap.match.isContainment();
        Edge::Comp comp = (overlap.match.isRC) ? Edge::EC_REVERSE : Edge::EC_SAME;

        Vertex* verts[2];
        for (size_t i = 0; i < 2; ++i) {
            verts[i] = _graph->getVertex(overlap.id[i]);

            // If one of the vertices is not in the graph, skip this edge
            // This can occur if one of the verts is a strict substring of some other vertex so it will
            // never be added to the graph
            if (verts[i] == NULL) {
                return false;
            }
        }

        // Check if this is a substring containment, if so mark the contained read
        // but do not create edges
        for (size_t i = 0; i < 2; ++i) {
            if (!overlap.match.coords[i].isExtreme()) {
                return false;
            }
        }

        // If either vertex has the maximum number of edges,
        // do not add any more. This is to protect against ultra-dense
        // regions of the graph inflating memory usage. The nodes that reach
        // this limit, and nodes connected to them are marked as super repeats.
        // After loading the graph, all edges to super repeats are cut to prevent
        // misassembly.
        {
            size_t degrees0 = verts[0]->degrees(), degrees1 = verts[1]->degrees();
            if (degrees0 > _maxEdges || degrees1 > _maxEdges) {
                return NULL;
            }
        }

        if (!isContainment) {
            Edge* edges[2];
            for (size_t i = 0; i < 2; ++i) {
                const SeqCoord& coord = overlap.match.coords[i];
                Edge::Dir dir = coord.isLeftExtreme() ? Edge::ED_ANTISENSE : Edge::ED_SENSE;
                edges[i] = new Edge(verts[1 - i], dir, comp, coord);
            }

            edges[0]->twin(edges[1]);
            edges[1]->twin(edges[0]);

            _graph->addEdge(verts[0], edges[0]);
            _graph->addEdge(verts[1], edges[1]);
        } else {
            // Contained edges don't have a direction, they can be travelled from
            // one vertex to the other in either direction. Hence, we 
            // add two edges per vertex. Later during the contain removal
            // algorithm this is important to determine transitivity
            Edge* edges[4];
            for (size_t i = 0; i < 2; ++i) {
                const SeqCoord& coord = overlap.match.coords[i];
                edges[i    ] = new Edge(verts[1 - i], Edge::ED_SENSE,     comp, coord);
                edges[i + 2] = new Edge(verts[1 - i], Edge::ED_ANTISENSE, comp, coord);
            }

            // Twin the edges and add them to the graph
            edges[0]->twin(edges[1]);
            edges[1]->twin(edges[0]);

            edges[2]->twin(edges[3]);
            edges[3]->twin(edges[2]);

            _graph->addEdge(verts[0], edges[0]);
            _graph->addEdge(verts[1], edges[1]);

            _graph->addEdge(verts[0], edges[2]);
            _graph->addEdge(verts[1], edges[3]);
        }

        return true;
    }
private:
    Bigraph* _graph;
    bool _allowContainments;
    size_t _maxEdges;
};

bool loadASQG(std::istream& stream, size_t minOverlap, bool allowContainments, size_t maxEdges, Bigraph* g) {
    enum {
        STAGE_HEAD, 
        STAGE_VERTEX, 
        STAGE_EDGE
    };
    
    int stage = STAGE_HEAD;
    std::string line;
    while (std::getline(stream, line)) {
        ASQG::RecordType rt = ASQG::recordtype(line);
        switch(rt) {
            case ASQG::RT_HEADER: {
                if (stage != STAGE_HEAD) {
                    LOG4CXX_ERROR(logger, boost::format("Error: Unexpected header record found at line %s") % line);
                    return false;
                }
                ASQG::HeaderRecord record;
                if (!ASQG::HeaderRecord::parse(line, record)) {
                    LOG4CXX_ERROR(logger, boost::format("Error: Unexpected header record found at line %s") % line);
                    return false;
                }
                break;
            }
            case ASQG::RT_VERTEX: {
                if (stage == STAGE_HEAD) {
                    stage = STAGE_VERTEX;
                }
                if (stage != STAGE_VERTEX) {
                    LOG4CXX_ERROR(logger, boost::format("Error: Unexpected vertex record found at line %s") % line);
                    return false;
                }
                ASQG::VertexRecord record;
                if (!ASQG::VertexRecord::parse(line, record)) {
                    LOG4CXX_ERROR(logger, boost::format("Error: Unexpected vertex record found at line %s") % line);
                    return false;
                }
                Vertex* vertex = new Vertex(record.id, record.seq);
                if (!g->addVertex(vertex)) {
                    LOG4CXX_ERROR(logger, boost::format("Error: Attempted to insert vertex into graph with a duplicate id: %s") % vertex->id());
                    LOG4CXX_ERROR(logger, "All reads must have a unique identifier");
                    delete vertex;
                    return false;
                }
                break;
            }
            case ASQG::RT_EDGE: {
                if (stage == STAGE_VERTEX) {
                    stage = STAGE_EDGE;
                }
                if (stage != STAGE_EDGE) {
                    LOG4CXX_ERROR(logger, boost::format("Error: Unexpected edge record found at line %s") % line);
                    return false;
                }
                ASQG::EdgeRecord record;
                if (!ASQG::EdgeRecord::parse(line, record)) {
                    LOG4CXX_ERROR(logger, boost::format("Error: Unexpected edge record found at line %s") % line);
                    return false;
                }
                const Overlap& ovr = record.overlap();
                // Add the edge to the graph
                if (ovr.match.length() >= minOverlap) {
                    EdgeCreator creator(g, allowContainments, maxEdges);
                    if (!creator.create(ovr)) {
                        return false;
                    }
                }

                break;
            }
            default: {
                return false;
            }
        }
    }

    return true;
}

bool loadASQG(const std::string& filename, size_t minOverlap, bool allowContainments, size_t maxEdges, Bigraph* g) {
    std::ifstream file(filename.c_str());

    boost::iostreams::filtering_istreambuf buf;
    if (boost::algorithm::ends_with(filename, GZIP_EXT)) {
        buf.push(boost::iostreams::gzip_decompressor());
    }
    buf.push(file);

    std::istream stream(&buf);
    return loadASQG(stream, minOverlap, allowContainments, maxEdges, g);
}

