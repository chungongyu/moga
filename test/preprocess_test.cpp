#define BOOST_TEST_MODULE "siga.preprocess"
#include <boost/test/included/unit_test.hpp>

#include <sstream>

#include "kseq.h"
#include "reads.h"
#include "primer_screen.h"

BOOST_AUTO_TEST_SUITE(preprocess);

BOOST_AUTO_TEST_CASE(PrimerScreen_contains) {
    {
        std::string pcr_free_a = "AATGATACGGCGACCACCGAGATCTACA";
        BOOST_CHECK(PrimerScreen::containsPrimer(pcr_free_a));
    }
    {
        std::string pcr_free_b = "GATCGGAAGAGCGGTTCAGCAGGAATGC";
        BOOST_CHECK(PrimerScreen::containsPrimer(pcr_free_b));
    }
    {
        std::string pcr_free_b = "AGATCGGAAGAGCGGTTCAGCAGGAATGC";
        BOOST_CHECK(!PrimerScreen::containsPrimer(pcr_free_b));
    }
}

BOOST_AUTO_TEST_CASE(KSeq_transform) {
    DNASeq seq("test", "ACGTGAC");
    BOOST_CHECK_EQUAL(seq.name, "test");
    BOOST_CHECK_EQUAL(seq.seq, "ACGTGAC");
    BOOST_CHECK(seq.quality.empty());
    seq.make_reverse();
    BOOST_CHECK_EQUAL(seq.seq, "CAGTGCA");
    seq.make_complement();
    BOOST_CHECK_EQUAL(seq.seq, "GTCACGT");
}

BOOST_AUTO_TEST_CASE(KSeq_read) {
    // FASTA
    {
        std::stringstream stream(">test\tcomment\nACGTGAC\n");
        DNASeqList sequences;
        BOOST_CHECK(ReadDNASequences(stream, sequences));
        BOOST_CHECK_EQUAL(sequences.size(), 1);
    }
    // FASTQ
}

BOOST_AUTO_TEST_CASE(PairEnd_test) {
    BOOST_CHECK_EQUAL("R", PairEnd::basename("R/1"));
    BOOST_CHECK_EQUAL("R", PairEnd::basename("R/A"));
    BOOST_CHECK_EQUAL("R", PairEnd::basename("R/f"));

    BOOST_CHECK_EQUAL("R/2", PairEnd::id("R/1"));
    BOOST_CHECK_EQUAL("R/1", PairEnd::id("R/2"));
    BOOST_CHECK_EQUAL("R/B", PairEnd::id("R/A"));
    BOOST_CHECK_EQUAL("R/A", PairEnd::id("R/B"));
    BOOST_CHECK_EQUAL("R/r", PairEnd::id("R/f"));
    BOOST_CHECK_EQUAL("R/f", PairEnd::id("R/r"));
}

BOOST_AUTO_TEST_SUITE_END();
