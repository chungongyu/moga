#ifndef suffix_array_builder_h_
#define suffix_array_builder_h_

#include "kseq.h"

#include <string>

class SuffixArray;

class SuffixArrayBuilder {
public:
    static SuffixArrayBuilder* create(const std::string& algorithm);

    virtual SuffixArray* build(const DNASeqList& sequences, size_t threads = 1) = 0;
};

#endif // suffix_array_builder_h_
