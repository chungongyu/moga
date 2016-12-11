#include "bigraph_visitors.h"
#include "bigraph.h"
#include "kseq.h"
#include "reads.h"
#include "utils.h"

#include <boost/assign.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.BigraphVisitor"));

//
// ChimericVisitor
//
void ChimericVisitor::previsit(Bigraph* graph) {
    _chimeric = 0;
    graph->color(GC_WHITE);
}

bool ChimericVisitor::visit(Bigraph* graph, Vertex* vertex) {
    // Check if this node is chimeric
    const std::string& seq = vertex->seq();
    if (vertex->degrees(Edge::ED_SENSE) == 1 && vertex->degrees(Edge::ED_ANTISENSE) == 1 && seq.length() < _minLength) {
        Edge* prevEdge = vertex->edges(Edge::ED_ANTISENSE)[0];
        Edge* nextEdge = vertex->edges(Edge::ED_SENSE)[0];
        Vertex* prevVert = prevEdge->end();
        Vertex* nextVert = nextEdge->end();

        bool chimeric = true;
        if (chimeric) {
            chimeric &= (prevVert->degrees(Edge::ED_SENSE) >= 2);
        }
        if (chimeric) {
            chimeric &= (nextVert->degrees(Edge::ED_ANTISENSE) >= 2);
        }
        if (chimeric) {
            // smallest?
            bool smallest = false;
            {
                EdgePtrList edges = prevVert->edges(Edge::ED_SENSE);
                for (size_t k = 0; k < edges.size(); ++k) {
                    if (edges[k]->coord().length() > prevEdge->coord().length() && edges[k]->coord().length() - prevEdge->coord().length() >= _delta) {
                        smallest = true;
                    }
                }
            }
            {
                EdgePtrList edges = nextVert->edges(Edge::ED_ANTISENSE);
                for (size_t k = 0; k < edges.size(); ++k) {
                    if (edges[k]->coord().length() > nextEdge->coord().length() && edges[k]->coord().length() - nextEdge->coord().length() >= _delta) {
                        smallest = true;
                    }
                }
            }
            chimeric &= smallest;
        }
        if (chimeric) {
            vertex->color(GC_BLACK);
            ++_chimeric;
            return true;
        }
    }
    return false;
}

void ChimericVisitor::postvisit(Bigraph* graph) {
    graph->sweepVertices(GC_BLACK);
    LOG4CXX_INFO(logger, boost::format("[ChimericVisitor]: Removed %d chimeric") % _chimeric);
}
//
// ContainRemoveVisitor
//
void ContainRemoveVisitor::previsit(Bigraph* graph) {
    graph->color(GC_WHITE);

    // Clear the containment flag, if any containments are added
    // during this algorithm the flag will be reset and another
    // round must be re-run
    graph->containment(false);

    _contained = 0;
}

bool ContainRemoveVisitor::visit(Bigraph* graph, Vertex* vertex) {
    if (vertex->contained()) {
        // Add any new irreducible edges that exist when pToRemove is deleted
        // from the graph
        EdgePtrList edges = vertex->edges();

        // Delete the edges from the graph
        for (EdgePtrList::iterator i = edges.begin(); i != edges.end(); ++i) {
           Edge* edge = *i;
           Edge* twin = edge->twin();
           Vertex* end = edge->end();

           end->removeEdge(twin);
           vertex->removeEdge(edge);

           SAFE_DELETE(twin);
           SAFE_DELETE(edge);
        }

        vertex->color(GC_BLACK);
        ++_contained;

        return true;
    }
    return false;
}

void ContainRemoveVisitor::postvisit(Bigraph* graph) {
    graph->sweepVertices(GC_BLACK);
    LOG4CXX_INFO(logger, boost::format("[ContainRemoveVisitor] Removed %d containment vertices") % _contained);
}

//
// FastaVisitor
//
bool FastaVisitor::visit(Bigraph* graph, Vertex* vertex) {
    DNASeq seq(vertex->id(), vertex->seq());
    _stream << seq;
    return false;
}

//
// MaximalOverlapVisitor
//
void MaximalOverlapVisitor::previsit(Bigraph* graph) {
    // The graph must not have containments
    assert(!graph->containment());

    // Set all the vertices in the graph to "vacant"
    graph->color(GC_WHITE);

    _dummys = 0;
}

class OverlapCmp {
public:
    bool operator()(const Edge* x, const Edge* y) const {
        return x->coord().length() > y->coord().length();
    }
};

bool MaximalOverlapVisitor::visit(Bigraph* graph, Vertex* vertex) {
    bool modified = false;

    typedef bool(*PredicateEdge)(const Edge*);
    static PredicateEdge PredicateEdgeArray[Edge::ED_COUNT] = {
        MaximalOverlapVisitor::isSenseEdge, 
        MaximalOverlapVisitor::isAntiSenseEdge
    };

    for (size_t i = 0; i < Edge::ED_COUNT; ++i) {
        Edge::Dir dir = Edge::EDGE_DIRECTIONS[i];
        EdgePtrList edges = vertex->edges(dir);

        std::sort(edges.begin(), edges.end(), OverlapCmp());

        for (size_t j = 1; j < edges.size(); ++j) {
            if (edges[j]->color() == GC_BLACK) {
                continue;
            }

            if (edges[0]->coord().length() - edges[j]->coord().length() <= _delta) {
                continue;
            }

            EdgePtrList redges = edges[j]->end()->edges();
            EdgePtrList::iterator last = std::remove_if(redges.begin(), redges.end(), PredicateEdgeArray[i]);
            redges.resize(std::distance(redges.begin(), last));
            assert(!redges.empty());

            EdgePtrList::const_iterator largest = std::min_element(redges.begin(), redges.end(), OverlapCmp());

            if ((*largest)->coord().length() - edges[j]->coord().length() <= _delta) {
                continue;
            }

            if (dir == Edge::ED_SENSE) {
                LOG4CXX_INFO(logger, boost::format("[MaximalOverlapVisitor] remove edge %s->%s (%d)") % edges[j]->start()->id() % edges[j]->end()->id() % edges[j]->coord().length());
            } else {
                LOG4CXX_INFO(logger, boost::format("[MaximalOverlapVisitor] remove edge %s->%s (%d)") % edges[j]->end()->id() % edges[j]->start()->id() % edges[j]->coord().length());
            }
            edges[j]->color(GC_BLACK);
            edges[j]->twin()->color(GC_BLACK);
            ++_dummys;
            modified = true;
        }
    }

    return modified;
}

void MaximalOverlapVisitor::postvisit(Bigraph* graph) {
    graph->sweepEdges(GC_BLACK);
    LOG4CXX_INFO(logger, boost::format("[MaximalOverlapVisitor] Removed %d dummy edges") % _dummys);
}

bool MaximalOverlapVisitor::isSenseEdge(const Edge* edge) {
    return (!edge->match().isRC && edge->dir() == Edge::ED_SENSE) || (edge->match().isRC && edge->dir() == Edge::ED_ANTISENSE);
}

bool MaximalOverlapVisitor::isAntiSenseEdge(const Edge* edge) {
    return !isSenseEdge(edge);
}

// 
// BigraphSearchTree
//

template< typename DistanceT >
class BigraphSearchTree {
    class Node;
    typedef std::vector< Node* > NodePtrList;
    class Node {
        Node(Vertex* vertex, Edge* edge, size_t distance, Node* parent) : vertex(vertex), edge(edge), distance(distance), parent(parent) {
            if (parent != NULL) {
                ++parent->children;
            }
        }
        ~Node() {
            assert(children == 0);
            if (parent != NULL) {
                assert(parent->children > 0);
                --parent->children;
            }
        }

        Vertex* vertex;
        Edge* edge;
        size_t distance;
        Node* parent;

        size_t children;
    };
public:
    BigraphSearchTree(Vertex* start, Vertex* end, Edge::Dir searchDir, size_t minDistance, size_t maxDistance, size_t maxNodes) : _searchDir(searchDir), _minDistance(minDistance), _maxDistance(maxDistance), _maxNodes(maxNodes), _numNodes(1) {
        _leaves = boost::assign::list_of(new Node(start, NULL, 0, NULL));
    }
    ~BigraphSearchTree() {
        // Delete the tree: 
        // delete each leaf and recurse up the tree iteratively.
        for (typename NodePtrList::iterator i = _leaves.begin(); i != _leaves.end(); ++i) {
            Node* curr = *i;
            // loop invariant: curr is a deletable node
            do {
                assert(curr->children == 0);
                Node* next = curr->parent;
                delete curr;
                curr = next;
            } while (curr != NULL && curr->children == 0);
        }
    }

    void build() {
        typedef std::deque< Node* > NodePtrQueue;

        NodePtrQueue Q;
        std::copy(_leaves.begin(), _leaves.end(), std::back_inserter(Q));

        _leaves.clear();
        while (!Q.empty()) {
            if (_numNodes >= _maxNodes) {
                std::copy(Q.begin(), Q.end(), std::back_inserter(Q));
                break;
            }

            Node* curr = Q.front();
            Q.pop_front();

            assert(curr->children == 0);
            EdgePtrList edges = curr->vertex->getEdges(_searchDir);

            if (curr->distance >= _maxDistance || edges.empty()) {
                _leaves.push_back(curr);
            } else {
                for (EdgePtrList::const_iterator i = edges.begin(); i != edges.end(); ++i) {
                    Edge* edge = *i;
                    Node* child = new Node(edge->end(), edge, curr->distance + 1, curr);
                    _leaves.push_back(child);
                }
            }
        }
    }
private:
    NodePtrList _leaves;

    Vertex* _end;
    Edge::Dir _searchDir;
    size_t _numNodes;

    size_t _minDistance;
    size_t _maxDistance;
    size_t _maxNodes;
};

//
// PairedReadVisitor
//

void PairedReadVisitor::previsit(Bigraph* graph) {
    graph->color(GC_WHITE);
    _dummys = 0;
}

bool PairedReadVisitor::visit(Bigraph* graph, Vertex* vertex) {
    if (vertex->color() == GC_WHITE) {
        /*
        // forward extend
        while (true) {
        }
        // backward extend
        while (true) {
        }
        */
    }
    return false;
}

void PairedReadVisitor::postvisit(Bigraph* graph) {
    graph->sweepEdges(GC_BLACK);
    LOG4CXX_INFO(logger, boost::format("[PairedReadVisitor] Removed %d dummy edges") % _dummys);
}

//
// SmoothingVisitor
//
void SmoothingVisitor::previsit(Bigraph* graph) {
    graph->color(GC_WHITE);
    _simple = 0;
    _complex = 0;
}

bool SmoothingVisitor::visit(Bigraph* graph, Vertex* vertex) {
    return false;
}

void SmoothingVisitor::postvisit(Bigraph* graph) {
    graph->sweepVertices(GC_RED);
    LOG4CXX_INFO(logger, boost::format("[SmoothingVisitor] Removed %s simple and %d complex bubbles") % _simple % _complex);
}

// 
// StatisticsVisitor
//
void StatisticsVisitor::previsit(Bigraph*) {
    _terminal = 0;
    _island = 0;
    _monobranch = 0;
    _dibranch = 0;
    _simple = 0;
    _edges = 0;
    _vertics = 0;
}

bool StatisticsVisitor::visit(Bigraph* graph, Vertex* vertex) {
    size_t fdegrees = vertex->degrees(Edge::ED_SENSE);
    size_t rdegrees = vertex->degrees(Edge::ED_ANTISENSE);

    if (fdegrees == 0 && rdegrees == 0) {
        ++_island;
    } else if (fdegrees == 0 || rdegrees == 0) {
        ++_terminal;
    }

    if (fdegrees > 1 && rdegrees > 1) {
        ++_dibranch;
    } else if (fdegrees > 1 || rdegrees > 1) {
        ++_monobranch;
    }

    if (fdegrees == 1 || rdegrees == 1) {
        ++_simple;
    }

    _edges += fdegrees + rdegrees;
    ++_vertics;

    return false;
}

void StatisticsVisitor::postvisit(Bigraph*) {
    LOG4CXX_INFO(logger, boost::format("[StatisticsVisitor] Vertices: %d Edges: %d Islands: %d Tips: %d Monobranch: %d Dibranch: %d Simple: %d") % _vertics % _edges % _island % _terminal % _monobranch % _dibranch % _simple);
}

//
// TrimVisitor
//
void TrimVisitor::previsit(Bigraph* graph) {
    _island = 0;
    _terminal = 0;
    graph->color(GC_WHITE);
}

bool TrimVisitor::visit(Bigraph* graph, Vertex* vertex) {
    bool modified = false;

    const std::string& seq = vertex->seq();
    if (vertex->degrees() == 0) {
        // Is an island, remove if the sequence length is less than the threshold
        if (seq.length() < _minLength) {
            LOG4CXX_TRACE(logger, boost::format("[TrimVisitor] island %s %d") % vertex->id() % seq.length());
            vertex->color(GC_BLACK);
            ++_island;
            modified = true;
        }
    } else {
        // Check if this node is a dead-end
        for (size_t idx = 0; idx < Edge::ED_COUNT; idx++) {
            Edge::Dir dir = Edge::EDGE_DIRECTIONS[idx];
            if (vertex->degrees(dir) == 0 && seq.length() < _minLength) {
                LOG4CXX_TRACE(logger, boost::format("[TrimVisitor] terminal %s %d") % vertex->id() % seq.length());
                vertex->color(GC_BLACK);
                ++_terminal;
                modified = true;
                break;
            }
        }
    }

    return modified;
}

void TrimVisitor::postvisit(Bigraph* graph) {
    graph->sweepVertices(GC_BLACK);
    LOG4CXX_INFO(logger, boost::format("[TrimVisitor] Removed %d island and %d dead-end short vertices") % _island % _terminal);
}
