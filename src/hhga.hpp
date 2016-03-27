#ifndef HHGA_H
#define HHGA_H

#include <map>
#include <vector>
#include <string>
#include <getopt.h>
#include <stdlib.h>
#include <iostream>
#include "bamtools/api/BamMultiReader.h"
#include "bamtools/api/BamWriter.h"
#include "bamtools/api/SamReadGroup.h"
#include "Fasta.h"
#include "Variant.h"
#include <cmath>

namespace hhga {

//using namespace BamTools;
//using namespace vcflib;
using namespace std;

#define PHRED_MAX 1000

typedef double prob_t;
typedef int32_t pos_t;

class allele_t {
public:
    friend ostream& operator<<(ostream& out, allele_t& var);
    friend bool operator<(const allele_t& a, const allele_t& b);
    friend allele_t operator+(const allele_t& a, const allele_t& b);
    prob_t prob;
    string ref;
    string alt;
    string repr;
    int32_t position;
    allele_t(string r, string a, long p, prob_t t)
        : ref(r), alt(a), position(p), prob(t)
    {
        // a string representation is made and used for sorting
        stringstream s;
        s << position << ":" << ref << "/" << alt;
        repr = s.str();
    }
};

short qualityChar2ShortInt(char c);
long double qualityChar2LongDouble(char c);
long double lnqualityChar2ShortInt(char c);
char qualityInt2Char(short i);
long double ln2log10(long double prob);
long double log102ln(long double prob);
long double phred2ln(int qual);
long double ln2phred(long double prob);
long double phred2float(int qual);
long double float2phred(long double prob);
void parse_region(const string& region,
                  string& startSeq,
                  int32_t& startPos,
                  int32_t& stopPos);
void set_region(BamTools::BamMultiReader& reader,
                const string& startSeq,
                int startPos,
                int stopPos);
void set_region(vcflib::VariantCallFile& vcffile, const string& region_str);
vector<prob_t> deletion_probs(const vector<prob_t>& quals, size_t sp, size_t l);
vector<prob_t> insertion_probs(const vector<prob_t>& quals, size_t sp, size_t l);

class HHGA {
public:
    string chrom_name;
    int32_t begin_pos;
    int32_t end_pos;

    typedef BamTools::BamAlignment alignment_t;
    //set<allele_t> alleles;
    vector<alignment_t> alignments;
    map<alignment_t*, vector<allele_t> > alignment_alleles;

    // the class label for the example
    string label;
    
    // the feature model
    vector<allele_t> reference;
    vector<vector<allele_t> > haplotypes; // haps/genotypes
    vector<vector<allele_t> > alleles;
    vector<prob_t> mapping_qualities;

    // helpers for construction 
    vector<allele_t> pad_alleles(vector<allele_t> aln_alleles,
                                 pos_t bal_min, pos_t bal_max);
    void project_positions(vector<allele_t>& aln_alleles,
                           map<pair<int32_t, size_t>, size_t>& pos_proj);

    // construct the hhga of a particular region
    HHGA(size_t window_size,
         BamTools::BamMultiReader& bam_reader,
         FastaReference& fasta_ref,
         vcflib::Variant& var,
         const string& class_label);

    const string str(void);
    const string vw(void);
};

}

#endif
