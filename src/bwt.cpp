#include "bwt.h"
#include "suffix_array.h"

BWT::BWT(const SuffixArray& sa, const DNASeqList& sequences) {
}

std::ostream& operator<<(std::ostream& stream, const BWT& bwt) {
    return stream;
}

std::istream& operator>>(std::istream& stream, BWT& bwt) {
    return stream;
}

const uint16_t BWT_FILE_MAGIC = 0xCACA;

bool BWTReader::read(BWT& bwt) {
    BWFlag flag;
    if (!readHeader(bwt._strings, bwt._suffixes, flag)) {
        return false;
    }
    if (!readRuns(bwt._runs, _numRuns)) {
        return false;
    }
    return true;
}

bool BWTReader::readHeader(size_t& num_strings, size_t& num_suffixes, BWFlag& flag) {
    uint16_t magic;
    if (!_stream.read((char *)&magic, sizeof(magic)) || magic != BWT_FILE_MAGIC) {
        return false;
    }
    if (!_stream.read((char *)&num_strings, sizeof(num_strings))) {
        return false;
    }
    if (!_stream.read((char *)&num_suffixes, sizeof(num_suffixes))) {
        return false;
    }
    if (!_stream.read((char *)&_numRuns, sizeof(_numRuns))) {
        return false;
    }
    if (!_stream.read((char *)&flag, sizeof(flag))) {
        return false;
    }
    return true;
}

bool BWTReader::readRuns(RLList& runs, size_t numRuns) {
    runs.resize(numRuns);
    if (!runs.empty()) {
        if (!_stream.read((char *)&runs[0], numRuns * sizeof(runs[0]))) {
            return false;
        }
    }
    return true;
}

bool BWTWriter::write(const SuffixArray& sa, const DNASeqList& sequences) {
    size_t num_strings = sa.strings(), num_suffixes = sa.size();

    if (!writeHeader(num_strings, num_suffixes, BWF_NOFMI)) {
        return false;
    }
    for (size_t i = 0; i < num_suffixes; ++i) {
        const SuffixArray::Elem& elem = sa[i];
        const DNASeq& read = sequences[elem.i];
        char c = (elem.j == 0 ? '$' : read.seq[elem.j - 1]);
        if (!writeChar(c)) {
            return false;
        }
    }
    if (!finalize()) {
        return false;
    }
    return true;
}

bool BWTWriter::writeHeader(size_t num_strings, size_t num_suffixes, BWFlag flag) {
    if (!_stream.write((const char *)&BWT_FILE_MAGIC, sizeof(BWT_FILE_MAGIC))) {
        return false;
    }
    if (!_stream.write((const char *)&num_strings, sizeof(num_strings))) {
        return false;
    }
    if (!_stream.write((const char *)&num_suffixes, sizeof(num_suffixes))) {
        return false;
    }

    // Here we do not know the number of runs that are going to be written to the file
    // so we save the offset in the file and write a dummy value. After the bwt string
    // has been written, we return here and fill in the correct value
    _posRun = _stream.tellp();
    _numRuns = 0;
    if (!_stream.write((const char *)&_numRuns, sizeof(_numRuns))) {
        return false;
    }

    assert(flag == BWF_NOFMI);
    if (!_stream.write((const char *)&flag, sizeof(flag))) {
        return false;
    }

    return true;
}

bool BWTWriter::writeChar(char c) {
    if (_currRun.initialized()) {
        if (_currRun == c && !_currRun.full()) {
            ++_currRun;
        } else {
            // Write out the old run and start a new one
            if (!writeRun(_currRun)) {
                return false;
            }
            _currRun = RLUnit(c);
        }
    } else {
        // Start a new run 
        _currRun = RLUnit(c);
    }
    return true;
}

bool BWTWriter::finalize() {
    if (_currRun.initialized()) {
        if (!writeRun(_currRun)) {
            return false;
        }
    }
    _stream.seekp(_posRun);
    _stream.write((const char *)&_numRuns, sizeof(_numRuns));
    _stream.seekp(std::ios_base::end);
    return _stream;
}

bool BWTWriter::writeRun(const RLUnit& run) {
    if (!_stream.write((const char *)&run.data, sizeof(run.data))) {
        return false;
    }
    ++_numRuns;
    return true;
}
