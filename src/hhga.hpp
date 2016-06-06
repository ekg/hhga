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
#include <iomanip>
#include "vg.hpp"
#include "multichoose.h"
#include "join.h"

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
    const string str(void) const { stringstream s; s << position << ":" << ref << "/" << alt; return s.str(); }
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

std::vector<std::string> &split_delims(const std::string &s,
                                       const std::string& delims,
                                       std::vector<std::string> &elems);
std::vector<std::string> split_delims(const std::string &s,
                                      const std::string& delims);

map<int, double> test_labels(int alt_count);
// fraction of times where both are non-missing where they agree
double pairwise_identity(const vector<allele_t>& h1, const vector<allele_t>& h2);
double pairwise_qualsum(const vector<allele_t>& h1, const vector<allele_t>& h2);
int missing_count(const vector<allele_t>& hap);
double entropy(const string& st);
bool is_repeat_unit(const string& seq, const string& unit);
string repeat(const string& s, int n);
map<string, int> repeat_counts(int position, const string& sequence, int maxsize);
pair<int, int> callable_window(int pos,
                               string sequence,
                               string alleleseq,
                               int min_repeat_size,
                               double min_repeat_entropy);
vector<vector<int> > possible_genotypes(int allele_count, int ploidy);
string string_for_genotype(const vector<int>& gt);
int label_for_genotype(const string& gt, const vector<vector<int> >& genotypes);
string genotype_for_label(int label, const vector<vector<int> >& genotypes);
pair<int, int> pair_for_gt_class(int gt);

class HHGA {
public:
    string chrom_name;
    int32_t begin_pos;
    int32_t end_pos;
    string repr; // for representing the site and variants

    typedef BamTools::BamAlignment alignment_t;

    // graph
    vg::VG graph;
    vector<vg::Alignment> graph_alns;
    map<int, double> graph_weights;
    map<int, double> graph_coverage;

    //set<allele_t> alleles;
    int alignment_count;
    vector<alignment_t> alignments;
    set<alignment_t*> unitigs;
    map<alignment_t*, vector<allele_t> > alignment_alleles;
    map<alignment_t*, map<int, double> > matches;
    map<alignment_t*, map<int, double> > qualsum;
    // handling genotype likelihoods
    double prob_aln_gt(alignment_t* aln, int gt);
    map<alignment_t*, map<int, double> > prob_aln_given_genotype;;
    map<int, double> likelihoods;
    map<int, map<double, vector<alignment_t*> > > allele_support;
    map<int, vector<alignment_t*> > allele_examples; // same as support but limited number
    vector<alignment_t*> softclipped;
    map<alignment_t*, vector<string> > alignment_groups;
    vector<pair<string, alignment_t*> > grouped_normal_alignments;
    vector<pair<string, alignment_t*> > grouped_unitig_alignments;
    map<string, map<int, int> > allele_counts;
    map<alignment_t*, int> missing_counts;

    // the class label for the example
    string label;

    // do we express things in exponentiated or phred form
    bool exponentiate;

    // the feature model
    vector<allele_t> reference;
    vector<vector<allele_t> > haplotypes; // genotypes
    vector<vector<allele_t> > genotypes; // genotypes
    map<int, int> sample_id;
    vector<vector<allele_t> > alleles;
    vector<prob_t> mapping_qualities;
    map<string, double> call_info_num; // from input VCFs, numbers
    map<string, string> call_info_str; // from input VCFs, strings

    // helpers for construction 
    vector<allele_t> pad_alleles(vector<allele_t> aln_alleles,
                                 pos_t bal_min, pos_t bal_max);
    void project_positions(vector<allele_t>& aln_alleles,
                           map<pair<int32_t, size_t>, size_t>& pos_proj);
    void flatten_to_ref(vector<allele_t>& alleles);
    void strandify(vector<allele_t>& alleles, bool is_rev);
    void missing_to_ref(vector<vector<allele_t> >& obs);

    // construct the hhga of a particular region
    HHGA(size_t window_size,
         BamTools::BamMultiReader& bam_reader,
         BamTools::BamMultiReader& unitig_reader,
         FastaReference& fasta_ref,
         vcflib::VariantCallFile& graph_vcf,
         size_t graph_window,
         vcflib::Variant& var,
         const string& input_name,
         const string& class_label,
         const string& gt_class,
         const vector<vector<int> >& all_genotypes,
         int max_depth = 0,
         int min_allele_count = 0,
         double min_repeat_entropy = 0,
         bool full_overlap = false,
         int max_node_size = 0,
         bool expon = false,
         bool show_bases = false,
         bool assume_ref = true);

    const string str(void);
    const string vw(void);
};

}

#endif
