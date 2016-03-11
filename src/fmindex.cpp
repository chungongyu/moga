#include "fmindex.h"
#include "bwt.h"
#include "utils.h"

#include <iterator>

#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.FMIndex"));

//
// MarkerFill
//
template< class MarkerList >
class MarkerFill {
public:
    MarkerFill(MarkerList& markers, size_t symbols, size_t sampleRate) : _markers(markers), _sampleRate(sampleRate) {
        initialize(symbols);
    }
    ~MarkerFill() {
        assert(_currIdx == _markers.size());
    }

    virtual void fill(size_t* counts, size_t total, size_t unitIndex, bool lastOne) = 0;
protected:
    void initialize(size_t n) {
        // we place a marker at the beginning (with no accumulated counts), every sampleRate
        // bases and one at the very end (with the total counts)
        size_t required_markers = (n % _sampleRate == 0) ? (n / _sampleRate) + 1 : (n / _sampleRate) + 2;
        _markers.resize(required_markers);

        // Place a blank markers at the start of the data
        if (!_markers.empty()) {
            memset(_markers[0].counts, 0, sizeof(_markers[0].counts));
            _markers[0].unitIndex = 0;
        }
        _currIdx = 1;
        _nextPos = _sampleRate;
    }

    MarkerList& _markers;
    size_t _sampleRate;

    size_t _currIdx;
    size_t _nextPos;
};

class LargeMarkerFill : public MarkerFill< LargeMarkerList > {
public:
    LargeMarkerFill(LargeMarkerList& markers, size_t symbols, size_t sampleRate) : MarkerFill< LargeMarkerList >(markers, symbols, sampleRate) {
        _currIdx = 1;
        _nextPos = _sampleRate;
    }

    void fill(size_t* counts, size_t total, size_t unitIndex, bool lastOne) {
        // Check whether to place a new marker
        while (total >= _nextPos || lastOne) {
            // Sanity checks
            // The marker position should always be less than the running total unless 
            // the number of symbols is smaller than the sample rate
            size_t expectedPos = _currIdx * _sampleRate;
            assert(expectedPos <= total || lastOne);
            assert(_currIdx < _markers.size());

            LargeMarker& marker = _markers[_currIdx++];
            marker.unitIndex = unitIndex;
            memcpy(marker.counts, counts, sizeof(marker.counts));

            _nextPos += _sampleRate;
            lastOne = lastOne && total >= _nextPos;
        }
    }
};

class SmallMarkerFill : public MarkerFill< SmallMarkerList > {
public:
    SmallMarkerFill(const LargeMarkerList& lmarkers, SmallMarkerList& smarkers, size_t symbols, size_t sampleRate) : MarkerFill< SmallMarkerList >(smarkers, symbols, sampleRate), _lmarkers(lmarkers) {
    }

    void fill(size_t* counts, size_t total, size_t unitIndex, bool lastOne) {
        while (total >= _nextPos || lastOne) {
            // Sanity checks
            // The marker position should always be less than the running total unless 
            // the number of symbols is smaller than the sample rate
            size_t expectedPos = _currIdx * _sampleRate;
            assert(expectedPos <= total || lastOne);
            assert(_currIdx < _markers.size());

            // Calculate the large marker to set the relative count from
            // This is generally the most previously placed large block except it might 
            // be the second-previous in the case that we placed the last large marker.
            size_t largeIndex = expectedPos / DEFAULT_SAMPLE_RATE_LARGE;
            const LargeMarker& lmarker = _lmarkers[largeIndex];
            SmallMarker& smarker = _markers[_currIdx++];

            for (size_t i = 0; i < SIZEOF_ARRAY(lmarker.counts); ++i) {
                //assert(counts[i] - lmarker.counts);
                smarker.counts[i] = counts[i] - lmarker.counts[i];
            }
            smarker.unitIndex = unitIndex - lmarker.unitIndex;

            _nextPos += _sampleRate;
            lastOne = lastOne && total >= _nextPos;
        }
    };

    const LargeMarkerList& _lmarkers;
};

void FMIndex::initialize() {
    assert(IS_POWER_OF_2(_sampleRate));

    size_t counts[DNAAlphabet::ALL_SIZE];
    memset(counts, 0, sizeof(counts));
    size_t total = 0;

    // Fill in the marker values
    // We wish to place markers every sampleRate symbols however since a run may
    // not end exactly on sampleRate boundaries, we place the markers AFTER
    // the run crossing the boundary ends

    LargeMarkerFill f1(_lmarkers, _bwt.length(), DEFAULT_SAMPLE_RATE_LARGE);
    SmallMarkerFill f2(_lmarkers, _smarkers, _bwt.length(), _sampleRate);

    const RLString& runs = _bwt.str();
    for (size_t i = 0; i < runs.size(); ++i) {
        const RLUnit& run = runs[i];
        char c = (char)run;
        size_t len = run.count();

        // Update the count and advance the running total
        counts[DNAAlphabet::torank(c)] += len;
        total += len;

        // Check whether to place a new large marker
        f1.fill(counts, total, i + 1, i == runs.size() - 1);

        // Check whether to place a new small marker
        f2.fill(counts, total, i + 1, i == runs.size() - 1);
    }

    // Initialize C(a)
    /*
    memset(_pred, 0, sizeof(_pred));
    _pred[DNAAlphabet::torank(DNAAlphabet::DNA_ALL[0])] = 0;
    for (size_t i = 1; i < DNAAlphabet::ALL_SIZE; ++i) {
        size_t curr = DNAAlphabet::torank(DNAAlphabet::DNA_ALL[i]), prev = DNAAlphabet::torank(DNAAlphabet::DNA_ALL[i - 1]);
        _pred[curr] = _pred[prev] + counts[prev];
    }
    */
    _pred[DNAAlphabet::DNA_ALL[0]] = 0;
    for (size_t i = 1; i < DNAAlphabet::ALL_SIZE; ++i) {
        char curr = DNAAlphabet::DNA_ALL[i], prev = DNAAlphabet::DNA_ALL[i - 1];
        _pred[curr] = _pred[prev] + counts[DNAAlphabet::torank(prev)];
    }
}

void FMIndex::info() const {
    const RLString& runs = _bwt.str();
    LOG4CXX_INFO(logger, "FMIndex info:");
    LOG4CXX_INFO(logger, boost::format("Large Sample rate: %d") % DEFAULT_SAMPLE_RATE_LARGE);
    LOG4CXX_INFO(logger, boost::format("Small Sample rate: %d") % _sampleRate);
    LOG4CXX_INFO(logger, boost::format("Contains %d symbols in %d runs (%1.4lf symbols per run)") % _bwt.length() % runs.size() % (runs.empty() ? 0 : ((double)_bwt.length() / runs.size())));
    LOG4CXX_INFO(logger, boost::format("Marker Memory -- Small Markers: %d (%.1lf MB) Large Markers: %d (%.1lf MB)") % _smarkers.size() % 0 % _lmarkers.size() % 0);
}

//
// MarkerFind
//
class MarkerFind {
public:
    MarkerFind(const RLString& runs, const LargeMarkerList& lmarkers, const SmallMarkerList& smarkers, size_t sampleRate) : _runs(runs), _lmarkers(lmarkers), _smarkers(smarkers), _sampleRate(sampleRate) {
    }
    size_t find(char c, size_t i) const {
        // The counts in the marker are not inclusive (unlike the Occurrence class)
        // so we increment the index by 1.
        ++i;

        LargeMarker lmarker = nearest(i);
        size_t position = lmarker.total();

        size_t r = lmarker.counts[DNAAlphabet::torank(c)];
        size_t currIdx = lmarker.unitIndex;

        // Search forwards (towards 0) until idx is found
        while (position < i) {
            size_t delta = i - position;

            assert(currIdx < _runs.size());

            const RLUnit& run = _runs[currIdx++];
            size_t n = run.count();
            if (n > delta) {
                n = delta;
            }
            if (run == c) {
                r += n;
            }
            position += n;
        }
        // Search backwards (towards 0) until idx is found
        while (position > i) {
            size_t delta = position - i;

            assert(currIdx < _runs.size());

            const RLUnit& run = _runs[currIdx--];
            size_t n = run.count();
            if (n > delta) {
                n = delta;
            }
            if (run == c) {
                r -= n;
            }
            position -= n;
        }
        
        assert(position == i);

        return r;
    }
private:
    LargeMarker nearest(size_t i) const {
        size_t baseIdx = i / _sampleRate; 
        size_t offset = MOD_POWER_2(i, _sampleRate); // equivalent to position % sampleRate
        if (offset >= _sampleRate>>1) {
            ++baseIdx;
        }
        return interpolated(baseIdx);
    }

    // Return a LargeMarker with values that are interpolated by adding
    // the relative count nearest to the requested position to the last
    // LargeMarker
    LargeMarker interpolated(size_t smallIdx) const {
        // Calculate the position of the LargeMarker closest to the target SmallMarker
        size_t largeIdx = smallIdx * _sampleRate / DEFAULT_SAMPLE_RATE_LARGE;

        LargeMarker absolute = _lmarkers[largeIdx];
        const SmallMarker& relative = _smarkers[smallIdx];
        for (size_t i = 0; i < SIZEOF_ARRAY(absolute.counts); ++i) {
            absolute.counts[i] += relative.counts[i];
        }
        absolute.unitIndex += relative.unitIndex;

        return absolute;
    }

    const RLString _runs;
    const LargeMarkerList _lmarkers;
    const SmallMarkerList _smarkers;
    size_t _sampleRate;
};

size_t FMIndex::getOcc(char c, size_t i) const {
    MarkerFind finder(_bwt.str(), _lmarkers, _smarkers, _sampleRate);
    return finder.find(c, i);
}

std::ostream& operator<<(std::ostream& stream, const FMIndex& index) {
    //stream << index._bwt;
    stream << "lmarkers" << std::endl;
    for (size_t i = 0; i < index._lmarkers.size(); ++i) {
        const LargeMarker& marker = index._lmarkers[i];
        stream << "--------------" << std::endl;
        stream << marker.unitIndex << std::endl;
        std::copy(marker.counts, marker.counts + SIZEOF_ARRAY(marker.counts), std::ostream_iterator< size_t >(stream, " "));
        stream << std::endl;
    }
    stream << "smarkers" << std::endl;
    for (size_t i = 0; i < index._smarkers.size(); ++i) {
        const SmallMarker& marker = index._smarkers[i];
        stream << "--------------" << std::endl;
        stream << marker.unitIndex << std::endl;
        std::copy(marker.counts, marker.counts + SIZEOF_ARRAY(marker.counts), std::ostream_iterator< uint16_t >(stream, " "));
        stream << std::endl;
    }
    //std::copy(index._pred, index._pred + DNAAlphabet::ALL_SIZE, std::ostream_iterator< size_t >(stream, " "));
    stream << std::endl;

    return stream;
}

std::istream& operator>>(std::istream& stream, FMIndex& index) {
    stream >> index._bwt;
    index.initialize();
    return stream;
}
