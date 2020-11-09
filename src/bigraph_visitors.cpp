#include "bigraph_visitors.h"
#include "asqg.h"
#include "bigraph_search.h"
#include "kseq.h"
#include "reads.h"
#include "sequence_process_framework.h"
#include "utils.h"

#include <functional>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.BigraphVisitor"));

//
// EdgeColorVisitor
//
class EdgeColorVisitor : public BigraphVisitor {
public:
    EdgeColorVisitor(GraphColor c, bool twin=false) : _color(c), _twin(twin) {
    }
    EdgeColorVisitor(GraphColor c, std::function<bool(const Vertex*, const Edge* edge)> filter, bool twin=false) : _color(c), _filter(filter), _twin(twin) {
    }

    void previsit(Bigraph* graph) {
    }
    bool visit(Bigraph* graph, Vertex* vertex) {
        bool modified = false;
        const EdgePtrList& edges = vertex->edges();
        for (auto& edge : edges) {
            if (!_filter || _filter(vertex, edge)) {
                edge->color(_color);
                if (_twin) {
                    edge->twin()->color(_color);
                }
                modified = true;
            }
        }
        return modified;
    }
    void postvisit(Bigraph* graph) {
    }
private:
    GraphColor _color;
    bool _twin;
    std::function<bool(const Vertex*, const Edge*)> _filter;
};

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
        const EdgePtrList& edges = vertex->edges();

        // Delete the edges from the graph
        for (auto i = edges.begin(); i != edges.end(); ++i) {
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
    const std::string& index = vertex->index();
    if (!index.empty()) {
        ASQG::StringTagValue bx(index);
        if (!seq.comment.empty()) {
            seq.comment += ' ';
        }
        seq.comment += bx.tostring(ASQG::BARCODE_TAG);
    }
    _stream << seq;
    return false;
}

//
// IdenticalRemoveVisitor
//
void IdenticalRemoveVisitor::previsit(Bigraph* graph) {
    graph->color(GC_WHITE);
    _count = 0;
}

bool IdenticalRemoveVisitor::visit(Bigraph* graph, Vertex* vertex) {
    bool modified = false;
    if (vertex->contained()) {
        // Check if this vertex is identical to any other vertex
        const EdgePtrList& edges = vertex->edges();
        for (const auto& edge : edges) {
            Vertex* other = edge->end();
            if (vertex->seq().length() != other->seq().length()) {
                continue;
            }
            Overlap ovr(edge->start()->id(), edge->end()->id(), edge->match());
            if (!ovr.isContainment() || ovr.containedIdx() != 0) {
                continue;
            }
            if (vertex->seq() == other->seq()) {
                vertex->color(GC_BLACK);
                LOG4CXX_DEBUG(logger, boost::format("[IdenticalRemoveVisitor] Remove identical vertex: %s/%s") % vertex->id() % other->id());
                ++_count;
                break;
            }
        }
    }
    return modified;
}

void IdenticalRemoveVisitor::postvisit(Bigraph* graph) {
    graph->sweepVertices(GC_BLACK);
    LOG4CXX_INFO(logger, boost::format("[IdenticalRemoveVisitorRemoveVisitor] Removed %lu identical vertices") % _count);
}

//
// LoopRemoveVisitor
//
void LoopRemoveVisitor::previsit(Bigraph* graph) {
    _loops.clear();
}

bool LoopRemoveVisitor::visit(Bigraph* graph, Vertex* vertex) {
    bool modified = false;
    if (vertex->degrees(Edge::ED_SENSE) == 1 && vertex->degrees(Edge::ED_ANTISENSE) == 1) {
        Edge* prevEdge = vertex->edges(Edge::ED_ANTISENSE)[0];
        Edge* nextEdge = vertex->edges(Edge::ED_SENSE)[0];
        Vertex* prevVert = prevEdge->end();
        Vertex* nextVert = nextEdge->end();
        if (!prevEdge->isSelf() && !nextEdge->isSelf() && prevVert == nextVert) {
            //vertex->color(GC_BLACK);
            _loops.push_back(vertex);
            modified = true;
        }
    }
    return modified;
}

void LoopRemoveVisitor::postvisit(Bigraph* graph) {
    /////////////////////////////////////////
    // R4 is a Loop node
    //
    //   R1-->R2-->R3
    //      /  \
    //      \  /
    //       R4
    //
    // We change it to:
    //   R1-->R2-->R4-->R2-->R3
    //
    /////////////////////////////////////////
    for (auto& vertex : _loops) {
        assert(vertex->degrees(Edge::ED_SENSE) == 1 && vertex->degrees(Edge::ED_ANTISENSE) == 1);
        Edge* prevEdge = vertex->edges(Edge::ED_ANTISENSE)[0];
        Edge* nextEdge = vertex->edges(Edge::ED_SENSE)[0];
        Vertex* prevVert = prevEdge->end();
        Vertex* nextVert = nextEdge->end();
        assert(prevVert == nextVert);

        LOG4CXX_INFO(logger, boost::format("[LoopRemoveVisitor] Remove loop vertex: id=%s, coverage=%d, seq=%d, prev vertex: id=%s, coverage=%d, seq=%d") % vertex->id() % vertex->coverage() % vertex->seq().length() % prevVert->id() % prevVert->coverage() % prevVert->seq().length());

        /////////////////////////////////////////
        //
        // R1-->R2-->R3
        //       \
        //        R4-->R2
        //
        ////////////////////////////////////////
        Edge* nextTwin = nextEdge->twin();
        vertex->merge(nextEdge);
        vertex->removeEdge(nextEdge);
        SAFE_DELETE(nextEdge);
        nextVert->removeEdge(nextTwin);
        SAFE_DELETE(nextTwin);

        /////////////////////////////////////////
        //
        //   R1-->R2-->R4-->R2-->R3
        //
        /////////////////////////////////////////
        Edge* prevTwin = prevEdge->twin();
        std::string label = prevTwin->label();
        bool prepend = false;
        if (prevTwin->dir() == Edge::ED_ANTISENSE) {
            prepend = true;
        }
        prevVert->merge(prevTwin);
        {
            EdgePtrList transEdges = prevVert->edges(Edge::EDGE_DIRECTIONS[Edge::ED_COUNT - prevEdge->dir() - 1]);
            for (auto& transEdge : transEdges) {
                if (transEdge != prevTwin && !prepend) {
                    SeqCoord& coord = transEdge->coord();
                    coord.interval.offset(label.length());
                }
            }
        }
        prevVert->removeEdge(prevTwin);
        SAFE_DELETE(prevTwin);
        vertex->removeEdge(prevEdge);
        SAFE_DELETE(prevEdge);

        graph->removeVertex(vertex);
        SAFE_DELETE(vertex);
    }
    LOG4CXX_INFO(logger, boost::format("[LoopRemoveVisitor] Removed %lu loop vertices") % _loops.size());
}

//
// MaximalOverlapVisitor
//
void MaximalOverlapVisitor::previsit(Bigraph* graph) {
    // The graph must not have containments
    assert(!graph->containment());

    // Set all the edges in the graph to "vacant"
    EdgeColorVisitor ecVisit(GC_WHITE);
    graph->visit(&ecVisit);

    _dummys = 0;
}

class OverlapCmp {
public:
    bool operator()(const Edge* x, const Edge* y) const {
        return x->coord().length() > y->coord().length();
    }
};

class EdgeDirCmp {
public:
    EdgeDirCmp(const Edge* edge) : _edge(edge) {
    }
    bool operator()(const Edge* edge) const {
        return _edge->twin()->dir() == (edge->dir() + 1) % Edge::ED_COUNT;
    }
private:
    const Edge* _edge;
};

bool MaximalOverlapVisitor::visit(Bigraph* graph, Vertex* vertex) {
    bool modified = false;

    for (size_t i = 0; i < Edge::ED_COUNT; ++i) {
        Edge::Dir dir = Edge::EDGE_DIRECTIONS[i];
        EdgePtrList fwdlist = vertex->edges(dir);

        std::sort(fwdlist.begin(), fwdlist.end(), OverlapCmp());

        for (size_t j = 1; j < fwdlist.size(); ++j) {
            if (fwdlist[j]->color() == GC_BLACK) {
                continue;
            }

            if (fwdlist[0]->coord().length() - fwdlist[j]->coord().length() < _delta) {
                continue;
            }

            EdgePtrList revlist = fwdlist[j]->end()->edges();
            auto last = std::remove_if(revlist.begin(), revlist.end(), EdgeDirCmp(fwdlist[j]));
            if (last != revlist.end()) {
                revlist.resize(std::distance(revlist.begin(), last));
            }
            assert(!revlist.empty());

            std::sort(revlist.begin(), revlist.end(), OverlapCmp());

            bool largest = revlist[0]->end() == vertex;
            for (size_t k = 1; k < revlist.size() && revlist[k]->coord().length() - revlist[0]->coord().length() < _delta && !largest; ++k) {
                largest = revlist[k]->end() == vertex; 
            }
            if (largest) {
                continue;
            }

            if (dir == Edge::ED_SENSE) {
                LOG4CXX_INFO(logger, boost::format("[MaximalOverlapVisitor] remove edge %s->%s (%d)") % fwdlist[j]->start()->id() % fwdlist[j]->end()->id() % fwdlist[j]->coord().length());
            } else {
                LOG4CXX_INFO(logger, boost::format("[MaximalOverlapVisitor] remove edge %s->%s (%d)") % fwdlist[j]->end()->id() % fwdlist[j]->start()->id() % fwdlist[j]->coord().length());
            }
            fwdlist[j]->color(GC_BLACK);
            fwdlist[j]->twin()->color(GC_BLACK);
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

// 
// InsertSizeEstimateVisitor
//
void InsertSizeEstimateVisitor::previsit(Bigraph* graph) {
    graph->color(GC_GREEN);

    _samples.clear();
}

bool InsertSizeEstimateVisitor::visit(Bigraph* graph, Vertex* vertex) {
    if (vertex->color() == GC_GREEN) {

        std::unordered_map<Vertex::Id, int> distancelist = {{vertex->id(), 0}};
        vertex->color(GC_RED);

        ////////////////////////////////////////
        // search forward & backward
        //
        // cases:
        //              c
        //      f    fc | rc   f
        //    R1-->R2-->R3-->R4-->R5
        //    ^    |^   |^   |^   |
        //     \__/ \__/ \__/ \__/
        //      r    fc   rc   r

        //              c    c
        //      f    fc | r  | rc   f
        //    R1-->R2-->R3-->R4-->R5-->R6
        //    ^    |^   |^   |^   |^   |
        //     \__/ \__/ \__/ \__/ \__/
        //      r    fc   f    rc   r
        //
        // the state machine:
        //                     ---------    
        //                    |         |
        //                    V         | f
        //    -----  f/rc  ---------    |
        //   |start| ---> | forward | --
        //    -----        ---------   
        //      |          ^     |fc   
        //      |          |rc   V       
        //       \  fc/r  ----------   
        //         ----> | backword | --
        //                ----------    |
        //                    ^         | r
        //                    |         |
        //                     ---------
        //
        ////////////////////////////////////////
        for (size_t idx = 0; idx < Edge::ED_COUNT; idx++) {
            Edge::Dir searchDir = Edge::EDGE_DIRECTIONS[idx];
            int distance = 0, flag = (searchDir == Edge::ED_SENSE ? 1 : -1);
            Vertex* p = vertex;

            while (true) {
                const EdgePtrList& straight = this->edges(p, searchDir); // forward & backword edges

                // Don't walk thru singular self edges
                if (straight.size() != 1 || straight[0]->isSelf() || straight[0]->end()->color() == GC_RED) {
                    break;
                }
                Edge* single = straight[0];
                Edge* twin = single->twin();
                Vertex* end = single->end();
                const EdgePtrList& opposite = this->edges(end, twin->dir());
                if (opposite.size() != 1) {
                    break;
                }

                p = end;
                if (searchDir == Edge::ED_SENSE) {      // forward
                    const SeqCoord& coord = single->coord();
                    distance += flag * (coord.seqlen - coord.length());
                    //distance += flag * single->coord().complement().length();
                } else {                                // backward
                    const SeqCoord& coord = twin->coord();
                    distance += flag * (coord.seqlen - coord.length());
                    //distance += flag * twin->coord().complement().length();
                }
                if (single->comp() == Edge::EC_REVERSE) {
                    searchDir = Edge::EDGE_DIRECTIONS[Edge::EC_COUNT - searchDir - 1];
                }
                distancelist[p->id()] = distance;
                LOG4CXX_DEBUG(logger, boost::format("InsertSizeEstimateVisitor::visit addVertex(%s,%s,%d)") % vertex->id() % p->id() % distance);

                p->color(GC_RED);
            }
        }

        for (auto i = distancelist.begin(); i != distancelist.end(); ++i) {
            Vertex::Id  pairId = PairEnd::id(i->first);
            if (i->first < pairId) {
                auto j = distancelist.find(pairId);
                if (j != distancelist.end()) {
                    size_t distance = std::abs(j->second - i->second);
                    _samples.push_back(distance);
                    LOG4CXX_DEBUG(logger, boost::format("InsertSizeEstimateVisitor::visit addLink(%s,%s,%d)") % i->first % j->first % distance);
                }
            }
        }
    }
    return false;
}

void InsertSizeEstimateVisitor::postvisit(Bigraph* graph) {
    LOG4CXX_INFO(logger, boost::format("InsertSizeEstimateVisitor::samples=%d") % _samples.size());

    typedef boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::count, boost::accumulators::tag::mean, boost::accumulators::tag::moment<2> > > Accumulator;
    Accumulator acc;
    std::for_each(_samples.begin(), _samples.end(), std::ref(acc));
    if (boost::accumulators::count(acc) > 0) {
        _average = (size_t)boost::accumulators::mean(acc);
        _delta = std::sqrt(
                boost::accumulators::moment<2>(acc) - std::pow(boost::accumulators::mean(acc), 2)
                );

        LOG4CXX_INFO(logger, boost::format("InsertSizeEstimateVisitor::average=%d, delta=%d") % _average % _delta);
    }

    graph->color(GC_GREEN);
}

EdgePtrList InsertSizeEstimateVisitor::edges(const Vertex* vertex, Edge::Dir dir) { // reduced edges
    EdgePtrList edges = vertex->edges(dir);
    std::sort(edges.begin(), edges.end(), OverlapCmp());
    // remove if the edge links to the same seq with `vertex` or it's subseq
    {
        size_t k = 0;
        for (size_t i = 0; i < edges.size(); ++i) {
            if (!edges[i]->coord().isContained() || !edges[i]->coord().isExtreme()) { // containment
                edges[k++] = edges[i];
            }
        }
        edges.resize(k);
    }
    // remove the other edges if they links to the same seq
    {
        size_t k = 0;
        for (size_t i = 1; i < edges.size(); ++i) {
            if (edges[i]->coord().length() != edges[k]->coord().length() || edges[i]->label() != edges[k]->label()) {
                edges[++k] = edges[i];
            }
        }
        if (!edges.empty()) {
            edges.resize(k + 1);
        }
    }
    return edges;
}

//
// PairedReadVisitor
//

void PairedReadVisitor::previsit(Bigraph* graph) {
    _vertices.clear();
}

bool PairedReadVisitor::visit(Bigraph* graph, Vertex* vertex1) {
    const Vertex* paired_v1 = graph->getVertex(PairEnd::id(vertex1->id()));
    assert(paired_v1 != NULL);

    //if (PairEnd::id(paired_v1->id()) < PairEnd::id(vertex1->id())) {
        _vertices.push_back(vertex1);
    //}
    return false;
}

typedef std::unordered_map<Vertex::Id, BigraphWalk::DistanceAttr> PairedDistanceMap;
typedef std::unordered_map<Vertex::Id, PairedDistanceMap> PairedLinkList;

class PairedVertexProcess {
public:
    PairedVertexProcess(Bigraph* graph, PairedReadVisitor* vistor) : _graph(graph), _visitor(vistor) {
    }
    BigraphWalk::NodePtrList process(const Vertex* vertex1) {
        BigraphWalk::NodePtrList linklist;

        const Vertex* paired_v1 = _graph->getVertex(PairEnd::id(vertex1->id()));
        assert(paired_v1 != NULL);

        // BFS the adjacents of vertex1
        BigraphWalk::NodePtrList adjacents;
        if (vertex1->seq().length() > _visitor->_maxDistance) {
            std::function<bool(const Edge* edge)> filter = [](const Edge* edge) -> bool {
                    if (edge->dir() == Edge::ED_SENSE || edge->comp() == Edge::EC_REVERSE) {
                        const Edge* e = NULL;
                        if (edge->dir() == Edge::ED_SENSE) {
                            e = edge;
                        } else {
                            e = edge->twin();
                        }
                        const SeqCoord& coord = e->coord();
                        return coord.seqlen > coord.length();
                    }
                    return false;
                };
            BigraphWalk::build(vertex1, filter, NULL, 0, _visitor->_maxDistance, _visitor->_maxNodes, &adjacents);
        }
        std::sort(adjacents.begin(), adjacents.end(), [](const BigraphWalk::NodePtr& x, const BigraphWalk::NodePtr& y) -> bool {
                    return std::abs(x->attr.distance) < std::abs(y->attr.distance);
                }); // sort by distance

        // Match each virtual read with paired vertex1
        for (const auto& node1 : adjacents) {
            const Vertex* paired_v2 = _graph->getVertex(PairEnd::id(node1->vertex->id()));
            assert(paired_v2 != NULL);
            LOG4CXX_DEBUG(logger, boost::format("vertex1: %s<->%s, vertex2: %s<->%s") % vertex1->id() % paired_v1->id() % node1->vertex->id() % paired_v2->id());

            BigraphWalk::NodePtrList faraways;
            for (size_t i = 0; i < Edge::ED_COUNT && faraways.empty(); ++i) {
                Edge::Dir dir = Edge::EDGE_DIRECTIONS[i];
                BigraphWalk::build(paired_v1, [dir](const Edge* edge)->bool {
                            return edge->dir() == dir;
                        }, paired_v2, 0, std::abs(node1->attr.distance) + _visitor->_insertDelta*4, 1, &faraways);
            }
            for (const auto& node2 : faraways) {
                linklist.push_back(node1);
                LOG4CXX_DEBUG(logger, boost::format("paired_read_all\t%s\t%s\t%d\t%s\t%s\t%d") % vertex1->id() % node1->vertex->id() % node1->attr.distance % paired_v1->id() % node2->vertex->id() % node2->attr.distance);
                break;
            }
        }

        return linklist;
    }
private:
    Bigraph* _graph;
    PairedReadVisitor* _visitor;
};

class PairedVertexPostProcess {
public:
    PairedVertexPostProcess(PairedLinkList* links) : _links(links) {
    }
    void process(const Vertex* vertex1, BigraphWalk::NodePtrList& linklist) {
        std::sort(linklist.begin(), linklist.end(), [](const BigraphWalk::NodePtr& x, const BigraphWalk::NodePtr& y) -> bool {
                    return std::abs(x->attr.distance) < std::abs(y->attr.distance);
                });
        for (size_t i = 0; i < linklist.size(); ++i) {
            const BigraphWalk::NodePtr& xi = linklist[i];

            BigraphWalk::DistanceAttr e = BigraphWalk::attrLink(xi->attr);
            addLink(vertex1->id(), xi->vertex->id(), e, _links);

            for (size_t j = i + 1; j < linklist.size(); ++j) {
                const BigraphWalk::NodePtr& xj = linklist[j];
                if (BigraphWalk::diffDir(xi->attr, xj->attr) || xi->attr.distance == xj->attr.distance) { // different search dir or equal
                    continue;
                }

                assert(std::abs(xi->attr.distance) < std::abs(xj->attr.distance));
                BigraphWalk::DistanceAttr e = BigraphWalk::attrLink(xi->attr, xj->attr);
                if (BigraphWalk::hasLink(xi->vertex, xj->vertex, e)) {
                    addLink(xi->vertex->id(), xj->vertex->id(), e, _links);
                } else {
                    LOG4CXX_DEBUG(logger, boost::format("paired_read_branch\t%s\t%s\t%d\t%s\t%d") % vertex1->id() % xi->vertex->id() % xi->attr.distance % xj->vertex->id() % xj->attr.distance);
                }
            }
        }
    }
private:
    void addLink(const Vertex::Id& v1, const Vertex::Id& v2, const BigraphWalk::DistanceAttr& e, PairedLinkList* links) {
        //LOG4CXX_TRACE(logger, boost::format("PairedLinkList::addLink\t%s\t%s\t%d\t%d\t%d") % v1 % v2 % e.distance % e.dir % e.comp);
        if (e.distance < 0) {
            BigraphWalk::DistanceAttr t = e.twin();
            t.distance = -t.distance;
            addLink(v2, v1, t, links);
        } else {
            auto i = links->find(v1);
            if (i == links->end()) {
                PairedDistanceMap tbl = {{v2, e}}; 
                (*links)[v1] = tbl;
            } else {
                auto k = i->second.find(v2);
                if (k != i->second.end() && k->second.distance != e.distance) {
                    LOG4CXX_WARN(logger, boost::format("PairedReadVisitor::addLink(%s, %s, %d) != %d") % v1 % v2 % e.distance % k->second.distance);
                }
                if (k == i->second.end() || k->second.distance > e.distance) {
                    i->second[v2] = e;
                }
            }
        }
    }

    PairedLinkList* _links;
};

class PairedEdgeCreator {
public:
    PairedEdgeCreator(Bigraph* graph) : _graph(graph) {
    }
    void create(const Vertex::Id& v1, const Vertex::Id& v2, const BigraphWalk::DistanceAttr& attr, GraphColor color) {
        assert(attr.distance > 0);

        Vertex* vertex[2] = {
            _graph->getVertex(v1), 
            _graph->getVertex(v2), 
        };

        bool found = false;
        const EdgePtrList& edges = vertex[0]->edges();
        for (auto& edge : edges) {
            //if (edge->dir() == attr.dir && edge->comp() == attr.comp && edge->end() == vertex[1]) {
            if (edge->dir() == attr.dir && edge->end() == vertex[1]) {
                if (edge->comp() == attr.comp && edge->coord().complement().length() == attr.distance) {
                    edge->color(color);
                    edge->twin()->color(color);
                    found = true;
                    break;
                }
                LOG4CXX_WARN(logger, boost::format("PairedReadVisitor::addEdge(%s, %s, %d) != %d") % v1 % v2 % attr.distance % edge->coord().interval.start);
            }
        }
        if (!found) {
            // sequnces
            const std::string& seq1 = vertex[0]->seq();
            const std::string& seq2 = vertex[1]->seq();

            // coordinates
            SeqCoord coord[2] = {
                SeqCoord(attr.distance, seq1.length() - 1, seq1.length()),
                SeqCoord(0, seq1.length() - attr.distance - 1, seq2.length())
            };
            if (attr.dir == Edge::ED_ANTISENSE) {
                coord[0] = SeqCoord(0, seq2.length() - attr.distance - 1, seq1.length());
                coord[1] = SeqCoord(attr.distance, seq2.length() - 1, seq2.length());
            }
            if (attr.comp == Edge::EC_REVERSE) {
                coord[1].flip();
            }

            EdgeCreator creator(_graph, true, -1);
            Overlap ovr(v1, coord[0], v2, coord[1], attr.comp == Edge::EC_REVERSE, 0);
            creator.create(ovr, color);
        }
    }
private:
    Bigraph* _graph;
};

struct PairedEdgeFilter {
    PairedEdgeFilter(GraphColor c) : _color(c), _vertex(NULL) {
        for (size_t i = 0; i < Edge::ED_COUNT; ++i) {
            _hasColor[i] = false;
        }
    }
    bool operator()(const Vertex* vertex, const Edge* edge) {
        //return edge->color() != _color;
        assert(vertex != NULL);
        if (vertex != _vertex) {
            _vertex = vertex;
            for (size_t i = 0; i < Edge::ED_COUNT; ++i) {
                _hasColor[i] = false;
            }
            const EdgePtrList& edges = vertex->edges();
            for (const auto& e : edges) {
                if (e->color() == _color) {
                    _hasColor[e->dir()] = true;
                }
            }
        }

        return (_hasColor[edge->dir()] && edge->color() != _color) || edge->coord().isFull();
    }
private:
    GraphColor _color;
    const Vertex* _vertex;
    bool _hasColor[Edge::ED_COUNT];
};

void PairedReadVisitor::postvisit(Bigraph* graph) {
    PairedLinkList links;

    SequenceProcessFramework::VectorWorkItemGenerator<const Vertex *> generator(_vertices);
    std::vector<PairedVertexProcess *> proclist(_threads);
    PairedVertexPostProcess postproc(&links);
    for (size_t i = 0; i < proclist.size(); ++i) {
        proclist[i] = new PairedVertexProcess(graph, this);
    }
#ifdef _OPENMP
    if (_threads > 1) {
        SequenceProcessFramework::ParallelWorker<
            const Vertex*, 
            BigraphWalk::NodePtrList, 
            SequenceProcessFramework::VectorWorkItemGenerator<const Vertex *>, 
            PairedVertexProcess, 
            PairedVertexPostProcess
            > worker;
        size_t num = worker.run(generator, &proclist, &postproc, _batch);
    } else { // single thread
#endif // _OPENMP
        SequenceProcessFramework::SerialWorker<
            const Vertex*, 
            BigraphWalk::NodePtrList, 
            SequenceProcessFramework::VectorWorkItemGenerator<const Vertex *>, 
            PairedVertexProcess, 
            PairedVertexPostProcess
            > worker;
        size_t num = worker.run(generator, proclist[0], &postproc);
#ifdef _OPENMP
    }
#endif // _OPENMP
    for (size_t i = 0; i < proclist.size(); ++i) {
        delete proclist[i];
    }

    ///////////////////////////////////////////
    // simplify
    //
    // NOTE that all the distance>=0
    //
    // GC_GRAY  => unkown
    // GC_WHITE => true positives
    // GC_BLACK => false positives
    //
    //////////////////////////////////////////
    {
        EdgeColorVisitor ecVist(GC_GRAY);
        graph->visit(&ecVist);
    }

    // create true positive edges.
    PairedEdgeCreator creator(graph);
    for (auto i = links.begin(); i != links.end(); ++i) {
        std::vector<std::pair<Vertex::Id, BigraphWalk::DistanceAttr> > nodelist;

        // sorted by distance
        std::copy(i->second.begin(), i->second.end(), std::back_inserter(nodelist));
        std::sort(nodelist.begin(), nodelist.end(), [](const std::pair<Vertex::Id, BigraphWalk::DistanceAttr>& x, const std::pair<Vertex::Id, BigraphWalk::DistanceAttr>& y) -> bool {
                    return x.second.distance < y.second.distance;
                });

        for (size_t j = 0; j < nodelist.size(); ++j) {
            const std::pair<Vertex::Id, BigraphWalk::DistanceAttr>& xj = nodelist[j];

            bool hasLink = false;
            for (size_t k = 0; k < j; ++k) {
                const std::pair<Vertex::Id, BigraphWalk::DistanceAttr>& xk = nodelist[k];
                //fixme: check correctness
                if (xk.second.dir == xj.second.dir && BigraphWalk::hasLink(graph->getVertex(xk.first), xk.second, graph->getVertex(xj.first), xj.second)) {
                    hasLink = true;
                    break;
                }
            }
            if (!hasLink) {
                creator.create(i->first, xj.first, xj.second, GC_WHITE);
                LOG4CXX_DEBUG(logger, boost::format("paired_read_simplify\t%s\t%s\t%d") % i->first % xj.first % xj.second.distance);
            }
        }
    }

    // cleanup false positive edges.
    {
        PairedEdgeFilter peFilter(GC_WHITE);
        EdgeColorVisitor ecVisit(GC_BLACK, peFilter, true);
        graph->visit(&ecVisit);
    }

    graph->sweepEdges(GC_BLACK);
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
// TransitiveReductionVisitor
//
void TransitiveReductionVisitor::previsit(Bigraph* graph) {
}

bool TransitiveReductionVisitor::visit(Bigraph* graph, Vertex* vertex) {
    bool modified = false;
    return modified;
}

void TransitiveReductionVisitor::postvisit(Bigraph* graph) {
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
        if (seq.length() <= _minLength) {
            LOG4CXX_TRACE(logger, boost::format("[TrimVisitor] island %s %d") % vertex->id() % seq.length());
            vertex->color(GC_BLACK);
            ++_island;
            modified = true;
        }
    } else {
        // Check if this node is a dead-end
        for (size_t idx = 0; idx < Edge::ED_COUNT; idx++) {
            Edge::Dir dir = Edge::EDGE_DIRECTIONS[idx];
            if (vertex->degrees(dir) == 0 && seq.length() <= _minLength) {
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
