#ifndef bigraph_visitors_h_
#define bigraph_visitors_h_

#include "bigraph.h"

#include <iostream>
#include <unordered_map>

class BigraphVisitor {
public:
    virtual void previsit(Bigraph* graph) {
    }
    virtual bool visit(Bigraph* graph, Vertex* vertex) = 0;
    virtual void postvisit(Bigraph* graph) {
    }
};

// Detects and removes small "chimeric" vertices from the graph
// when they are less than minLength in size
class ChimericVisitor : public BigraphVisitor {
public:
    ChimericVisitor(size_t minLength, size_t delta) : _minLength(minLength), _delta(delta) {
    }
    void previsit(Bigraph* graph);
    bool visit(Bigraph* graph, Vertex* vertex);
    void postvisit(Bigraph* graph);
private:
    size_t _minLength;

    size_t _delta;
    size_t _chimeric;
};

// Remove contained vertices from the graph
class ContainRemoveVisitor : public BigraphVisitor {
public:
    ContainRemoveVisitor() : _contained(0) {
    }
    void previsit(Bigraph* graph);
    bool visit(Bigraph* graph, Vertex* vertex);
    void postvisit(Bigraph* graph);
private:
    size_t _contained;
};

// Visit each node, writing it to a file as a fasta record
class FastaVisitor : public BigraphVisitor {
public:
    FastaVisitor(std::ostream& stream) : _stream(stream) {
    }
    bool visit(Bigraph* graph, Vertex* vertex);
private:
    std::ostream& _stream;
};

// Remove identical vertices from the graph
class IdenticalRemoveVisitor : public BigraphVisitor {
public:
    IdenticalRemoveVisitor() : _count(0) {
    }
    void previsit(Bigraph* graph);
    bool visit(Bigraph* graph, Vertex* vertex);
    void postvisit(Bigraph* graph);
private:
    size_t _count;
};

// Remove loops from the graph
class LoopRemoveVisitor : public BigraphVisitor {
public:
    LoopRemoveVisitor() {
    }
    void previsit(Bigraph* graph);
    bool visit(Bigraph* graph, Vertex* vertex);
    void postvisit(Bigraph* graph);
private:
    std::vector<Vertex *> _loops;
};

// Run the YU LIN's maximal overlap algorithm on each node
class MaximalOverlapVisitor : public BigraphVisitor {
public:
    MaximalOverlapVisitor(size_t delta) : _delta(delta), _dummys(0) {
    }
    void previsit(Bigraph* graph);
    bool visit(Bigraph* graph, Vertex* vertex);
    void postvisit(Bigraph* graph);
private:
    size_t _delta;
    size_t _dummys;
};

// Visit each paired node to estimate insert size variant under normal distribution assumption.
class InsertSizeEstimateVisitor : public BigraphVisitor {
public:
    InsertSizeEstimateVisitor(size_t& average, size_t& delta) : _average(average), _delta(delta) {
    }
    void previsit(Bigraph* graph);
    bool visit(Bigraph* graph, Vertex* vertex);
    void postvisit(Bigraph* graph);
public:
    EdgePtrList edges(const Vertex* vertex, Edge::Dir dir);
    std::vector<size_t> _samples;

    size_t& _average;
    size_t& _delta;
};

// Visit each paired node via zigzag, sweep false positive edges.
class PairedReadVisitor : public BigraphVisitor {
public:
    PairedReadVisitor(size_t maxDistance, size_t insertSize, size_t insertDelta, size_t maxNodes, size_t threads=1, size_t batch=1000) : _maxDistance(maxDistance), _insertSize(insertSize), _insertDelta(insertDelta), _maxNodes(maxNodes), _threads(threads), _batch(batch) {
    }
    void previsit(Bigraph* graph);
    bool visit(Bigraph* graph, Vertex* vertex);
    void postvisit(Bigraph* graph);
private:
    size_t _maxDistance;
    size_t _insertSize;
    size_t _insertDelta;
    size_t _maxNodes;
    size_t _threads;
    size_t _batch;

    std::vector<const Vertex *> _vertices;
    friend class PairedVertexProcess;
    friend class PairedVertexPostProcess;
};

// Visit each node via 10x linked read
class LinkedReadVisitor : public BigraphVisitor {
public:
    LinkedReadVisitor() {
    }
    void previsit(Bigraph* graph);
    bool visit(Bigraph* graph, Vertex* vertex);
    void postvisit(Bigraph* graph);
};

// Smooth out variation in the graph
class SmoothingVisitor : public BigraphVisitor {
public:
    SmoothingVisitor() {
    }
    void previsit(Bigraph* graph);
    bool visit(Bigraph* graph, Vertex* vertex);
    void postvisit(Bigraph* graph);
private:
    size_t _simple;
    size_t _complex;
};

// Compile summary statistics for the graph 
class StatisticsVisitor : public BigraphVisitor {
public:
    StatisticsVisitor() {
    }
    void previsit(Bigraph* graph);
    bool visit(Bigraph* graph, Vertex* vertex);
    void postvisit(Bigraph* graph);
private:
    size_t _terminal;
    size_t _island;
    size_t _monobranch;
    size_t _dibranch;
    size_t _simple;
    size_t _edges;
    size_t _vertics;
};

// Run the Myers transitive reduction algorithm on each node
class TransitiveReductionVisitor : public BigraphVisitor {
public:
    TransitiveReductionVisitor() {
    }
    void previsit(Bigraph* graph);
    bool visit(Bigraph* graph, Vertex* vertex);
    void postvisit(Bigraph* graph);
private:
};

// Detects and removes small "tip" vertices from the graph
// when they are less than minLength in size
class TrimVisitor : public BigraphVisitor {
public:
    TrimVisitor(size_t minLength) : _minLength(minLength) {
    }
    void previsit(Bigraph* graph);
    bool visit(Bigraph* graph, Vertex* vertex);
    void postvisit(Bigraph* graph);
private:
    size_t _minLength;

    size_t _island;
    size_t _terminal;
};

#endif // bigraph_visitors_h_
