#include "hhga.hpp"

namespace hhga {

short qualityChar2ShortInt(char c) {
    return static_cast<short>(c) - 33;
}

long double qualityChar2LongDouble(char c) {
    return static_cast<long double>(c) - 33;
}

long double lnqualityChar2ShortInt(char c) {
    return log(static_cast<short>(c) - 33);
}

char qualityInt2Char(short i) {
    return static_cast<char>(i + 33);
}

long double ln2log10(long double prob) {
    return M_LOG10E * prob;
}

long double log102ln(long double prob) {
    return M_LN10 * prob;
}

long double phred2ln(int qual) {
    return M_LN10 * qual * -.1;
}

long double ln2phred(long double prob) {
    return -10 * M_LOG10E * prob;
}

long double phred2float(int qual) {
    return pow(10, qual * -.1);
}

long double float2phred(long double prob) {
    if (prob == 1)
        return PHRED_MAX;  // guards against "-0"
    long double p = -10 * (long double) log10(prob);
    if (p < 0 || p > PHRED_MAX) // int overflow guard
        return PHRED_MAX;
    else
        return p;
}

std::vector<std::string> &split_delims(const std::string &s,
                                       const std::string& delims,
                                       std::vector<std::string> &elems) {
    char* tok;
    char cchars [s.size()+1];
    char* cstr = &cchars[0];
    strcpy(cstr, s.c_str());
    tok = strtok(cstr, delims.c_str());
    while (tok != NULL) {
        elems.push_back(tok);
        tok = strtok(NULL, delims.c_str());
    }
    return elems;
}

std::vector<std::string> split_delims(const std::string &s,
                                      const std::string& delims) {
    std::vector<std::string> elems;
    return split_delims(s, delims, elems);
}


void parse_region(
    const string& region,
    string& startSeq,
    int32_t& startPos,
    int32_t& stopPos) {

    size_t foundFirstColon = region.find(":");

    // we only have a single string, use the whole sequence as the target
    if (foundFirstColon == string::npos) {
        startSeq = region;
        startPos = 0;
        stopPos = -1;
    } else {
        startSeq = region.substr(0, foundFirstColon);
        string sep = "..";
        size_t foundRangeSep = region.find(sep, foundFirstColon);
        if (foundRangeSep == string::npos) {
            sep = "-";
            foundRangeSep = region.find("-", foundFirstColon);
        }
        if (foundRangeSep == string::npos) {
            startPos = atoi(region.substr(foundFirstColon + 1).c_str());
            // differ from bamtools in this regard, in that we process only
            // the specified position if a range isn't given
            stopPos = startPos + 1;
        } else {
            startPos = atoi(region.substr(foundFirstColon + 1, foundRangeSep - foundFirstColon).c_str());
            // if we have range sep specified, but no second number, read to the end of sequence
            if (foundRangeSep + sep.size() != region.size()) {
                stopPos = atoi(region.substr(foundRangeSep + sep.size()).c_str()); // end-exclusive, bed-format
            } else {
                //stopPos = reference.sequenceLength(startSeq);
                stopPos = -1;
            }
        }
    }
}

void set_region(BamTools::BamMultiReader& reader,
                const string& startSeq,
                int startPos,
                int stopPos) {

    map<string, int> refLength;
    map<string, int> refID;

    int id = 0;
    BamTools::RefVector references = reader.GetReferenceData();
    for (BamTools::RefVector::iterator r = references.begin(); r != references.end(); ++r) {
        refLength[r->RefName] = r->RefLength;
        refID[r->RefName] = id++;
    }

    // parse the region string

    if (stopPos == -1) {
        stopPos = refLength[startSeq];
    }

    int startSeqRefID = refID[startSeq];

    if (!reader.LocateIndexes()) {
        cerr << "[hhga] could not load BAM index" << endl;
        exit(1);
    } else {
        reader.SetRegion(startSeqRefID, startPos, startSeqRefID, stopPos);
    }
}

void set_region(vcflib::VariantCallFile& vcffile, const string& region_str) {
    if (vcffile.is_open()) {
        vcffile.setRegion(region_str);
    }
}

vector<vector<int> > possible_genotypes(int allele_count, int ploidy) {
    vector<int> alleles;
    for (int i = 0; i < allele_count; ++i) alleles.push_back(i);
    return multichoose(ploidy, alleles);
}

vector<int> genotype_for_string(const string& gt) {
    auto gtvs = split(gt, "/|");
    vector<int> gtvi;
    for (auto& s : gtvs) {
        gtvi.push_back(atoi(s.c_str()));
    }
    std::sort(gtvi.begin(), gtvi.end());
    return gtvi;
}

string string_for_genotype(const vector<int>& gt) {
    vector<string> gs;
    for (auto& i : gt) {
        gs.push_back(convert(i));
    }
    return join(gs, "/");
}

int label_for_genotype(const string& gt, const vector<vector<int> >& genotypes) {
    //auto genotypes = possible_genotypes(allele_count, ploidy);
    // split the string and turn it into a vector
    auto genotype = genotype_for_string(gt);
    int i = 0;
    for (auto& s : genotypes) {
        ++i;
        if (s == genotype) {
            return i;
        }
    }
    cerr << "warning: unknown genotype '" << gt << "'" << endl;
    return genotypes.size()+1;
}

string genotype_for_label(int label, const vector<vector<int> >& genotypes) {
    if (label > genotypes.size()) {
        cerr << "warning: unknown label '" << label << "'" << endl;
        return "./.";        
    } else {
        return string_for_genotype(genotypes.at(label));
    }
}

map<int, double> test_labels(int alt_count) {
    map<int, double> labels;
    for (int i = 1; i < alt_count+1; ++i) {
        labels[i] = 0;
    }
    return labels;
}

pair<int, int> pair_for_gt_class(int gt) {
    switch (gt) {
    case 1: // 0/0
        return make_pair(0,0);
        break;
    case 2: // 0/1
        return make_pair(0,1);
        break;
    case 3: // 0/2
        return make_pair(0,2);
        break;
    case 4: // 1/1
        return make_pair(1,1);
        break;
    case 5: // 1/2
        return make_pair(1,2);
        break;
    case 6: // 2/2
        return make_pair(2,2);
        break;
    case 7: // ./.
        return make_pair(8,8);
        break;
    default:
        break;
    }
}

double HHGA::prob_aln_gt(alignment_t* aln, int gt) {
    auto& match = matches[aln];
    auto& qsum = qualsum[aln];
    double prob = 0;
    auto gtp = pair_for_gt_class(gt);
    int a = gtp.first;;
    int b = gtp.second;
    double prob_sample = match[a]/2 + match[b]/2;
    double prob_in_gt = 1-(phred2float(qsum[a])/2 + phred2float(qsum[b])/2);
    double qsumsum = 0;
    for (auto& q : qsum) {
        qsumsum += q.second;
    }
    double prob_out_gt = 1-phred2float(qsumsum - (qsum[a]/2 + qsum[b]/2));
    //double prob_error = 1-(phred2float(qsum[a])/2 + phred2float(qsum[b])/2);
    //double prob_error = phred2float(qsum[a])/2 + phred2float(qsum[b])/2;
    return prob_sample * prob_in_gt;
}

double pairwise_qualsum(const vector<allele_t>& h1, const vector<allele_t>& h2) {
    // assert they are normalized
    // that given, we can compare directly
    map<int, const allele_t*> p1;
    map<int, const allele_t*> p2;
    for (size_t i = 0; i < h1.size(); ++i) {
        p1[i] = &h1[i];
    }
    for (size_t i = 0; i < h2.size(); ++i) {
        p2[i] = &h2[i];
    }
    int possible = 0;
    for (auto& p : p2) {
        auto a2 = p.second->alt;
        if (a2 != "M") ++possible;
    }
    double qualsum = 0;
    for (auto& p : p1) {
        auto a1 = p1[p.first]->alt;
        if (p2.find(p.first) == p2.end()) continue;
        auto a2 = p2[p.first]->alt;
        //if (possible > 1 && a1 == "R" && a1 == "R") continue; // avoid counting ref bases in indels
        if (a1 == "M" || a2 == "M") continue;
        if (a1 == a2) {
            qualsum += p.second->prob;
        }
    }
    return (possible ? (double) qualsum / (double) possible : 0);
}

double pairwise_identity(const vector<allele_t>& h1, const vector<allele_t>& h2) {
    int count = 0;
    // assert they are normalized
    // that given, we can compare directly
    // if we have a pure reference pair of alleles
    // the match is 1
    // otherwise, we measure identity in terms of the number of non-ref positions that match
    map<int, const allele_t*> p1;
    map<int, const allele_t*> p2;
    for (size_t i = 0; i < h1.size(); ++i) {
        if (h1[i].alt == "M") continue;
        p1[i] = &h1[i];
    }
    for (size_t i = 0; i < h2.size(); ++i) {
        if (h2[i].alt == "M") continue;
        p2[i] = &h2[i];
    }
    int possible = p2.size();
    int seen = 0;
    for (auto& p : p2) {
        auto a2 = p2[p.first]->alt;
        if (p1.find(p.first) == p1.end()) continue;
        auto a1 = p1[p.first]->alt;
        if (a1 == "M" || a2 == "M") continue;
        ++seen;
        if (a1 == a2) {
            ++count;
        }
    }
    return (possible ? (double) count / (double) possible : 0);
}

bool has_softclip(const vector<allele_t>& aln_alleles) {
    for (auto& allele : aln_alleles) {
        if (allele.alt == "S") return true;
    }
    return false;
}

int missing_count(const vector<allele_t>& hap) {
    int m = 0;
    for (auto& a : hap) if (a.alt == "M") ++m;
    return m;
}

void HHGA::missing_to_ref(vector<vector<allele_t> >& obs) {
    for (auto& hap : obs) {
        size_t i = 0;
        for (vector<allele_t>::iterator a = hap.begin(); a != hap.end(); ++a, ++i) {
            if (a->alt == "M") *a = reference[i];
        }
    }
}

map<string, int> repeat_counts(int position, const string& sequence, int maxsize) {
    map<string, int> counts;
    for (int i = 1; i <= maxsize; ++i) {
        // subseq here i bases
        string seq = sequence.substr(position, i);
        // go left.

        int j = position - i;
        int leftsteps = 0;
        while (j >= 0 && seq == sequence.substr(j, i)) {
            j -= i;
            ++leftsteps;
        }

        // go right.
        j = position;

        int rightsteps = 0;
        while (j + i <= sequence.size() && seq == sequence.substr(j, i)) {
            j += i;
            ++rightsteps;
        }
        // if we went left and right a non-zero number of times, 
        if (leftsteps + rightsteps > 1) {
            counts[seq] = leftsteps + rightsteps;
        }
    }

    // filter out redundant repeat information
    if (counts.size() > 1) {
        map<string, int> filteredcounts;
        map<string, int>::iterator c = counts.begin();
        string prev = c->first;
        filteredcounts[prev] = c->second;  // shortest sequence
        ++c;
        for (; c != counts.end(); ++c) {
            int i = 0;
            string seq = c->first;
            while (i + prev.length() <= seq.length() && seq.substr(i, prev.length()) == prev) {
                i += prev.length();
            }
            if (i < seq.length()) {
                filteredcounts[seq] = c->second;
                prev = seq;
            }
        }
        return filteredcounts;
    } else {
        return counts;
    }
}

double entropy(const string& st) {
    vector<char> stvec(st.begin(), st.end());
    set<char> alphabet(stvec.begin(), stvec.end());
    vector<double> freqs;
    for (set<char>::iterator c = alphabet.begin(); c != alphabet.end(); ++c) {
        int ctr = 0;
        for (vector<char>::iterator s = stvec.begin(); s != stvec.end(); ++s) {
            if (*s == *c) {
                ++ctr;
            }
        }
        freqs.push_back((double)ctr / (double)stvec.size());
    }
    double ent = 0;
    double ln2 = log(2);
    for (vector<double>::iterator f = freqs.begin(); f != freqs.end(); ++f) {
        ent += *f * log(*f)/ln2;
    }
    ent = -ent;
    return ent;
}

bool is_repeat_unit(const string& seq, const string& unit) {
    if (seq.size() % unit.size() != 0) {
        return false;
    } else {
        int maxrepeats = seq.size() / unit.size();
        for (int i = 0; i < maxrepeats; ++i) {
            if (seq.substr(i * unit.size(), unit.size()) != unit) {
                return false;
            }
        }
        return true;
    }
}

string repeat(const string& s, int n) {
    ostringstream os;
    for(int i = 0; i < n; i++)
        os << s;
    return os.str();
}

pair<int, int> callable_window(int pos,
                               string sequence,
                               string alleleseq,
                               int min_repeat_size,
                               double min_entropy) {

    // force the callable window to extend across
    // tandem repeats and homopolymers when indels are present
    int right_boundary = pos+1;
    int left_boundary = pos;
    // check for repeats of up to length 16
    for (auto& r : repeat_counts(pos, sequence, 16)) {
        auto repeatunit = r.first;
        int rptcount = r.second;
        //cerr << repeatunit << ":" << rptcount << endl;
        string repeatstr = repeat(repeatunit, rptcount);
        // assumption of left-alignment may be problematic... so this should be updated
        if (repeatstr.size() >= min_repeat_size) { // && is_repeat_unit(alleleseq, repeatunit)) {
            // determine the boundaries of the repeat
            // adjust to ensure we hit the first of the repeatstr
            size_t startpos = sequence.find(repeatstr, pos-1-repeatstr.size());
            left_boundary = min((int)startpos, left_boundary);
            if (startpos == string::npos) {
                /*
                cerr << "could not find repeat sequence?" << endl;
                cerr << "repeat sequence: " << repeatstr << endl;
                cerr << sequence << endl;
                cerr << "matched repeats:" << endl;
                */
                break; // ignore right-repeat boundary in this case
            }
            right_boundary = max(right_boundary, (int)(left_boundary + repeatstr.size() + 1)); // 1 past edge of repeat
        }
    }

    while (min_entropy > 0 && // ignore if turned off
           left_boundary > 0 &&
           right_boundary < sequence.size()-1 && //guard
           entropy(sequence.substr(left_boundary, right_boundary - left_boundary)) < min_entropy) {
        --left_boundary;
        ++right_boundary;
    }

    return make_pair(left_boundary, right_boundary);
    // edge case, the indel is an insertion and matches the reference to the right
    // this means there is a repeat structure in the read, but not the ref
    /*
    if (currentSequence.substr(pos - currentSequenceStart, length) == readSequence) {
        repeatRightBoundary = max(repeatRightBoundary, pos + length + 1);
    }
    */
}

HHGA::HHGA(size_t window_length,
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
           int max_depth,
           int min_allele_count,
           double min_repeat_entropy,
           bool full_overlap,
           int max_node_size,
           bool expon,
           bool show_bases,
           bool assume_ref) {

    exponentiate = expon;

    if (!gt_class.empty()) {
        // convert the genotype into
        // require that it be in
        // 0/0, 0/1, 1/1, 0/2, 1/2, 2/2
        auto gt = var.samples[var.sampleNames.front()][gt_class].front();
        label = convert(label_for_genotype(gt, all_genotypes));
    } else {
        label = class_label;
    }

    // store the names of all the reference sequences in the BAM file
    map<int, string> referenceIDToName;
    vector<BamTools::RefData> referenceSequences = bam_reader.GetReferenceData();
    int i = 0;
    for (BamTools::RefVector::iterator r = referenceSequences.begin();
         r != referenceSequences.end(); ++r) {
        referenceIDToName[i] = r->RefName;
        ++i;
    }


    int32_t begin_pos = var.position-1 - window_length/2;
    int32_t end_pos = begin_pos + window_length;
    string seq_name = var.sequenceName;
    int32_t center_pos = var.position-1;//begin_pos + (end_pos - begin_pos) / 2;
    //int32_t center_pos = var.position-1 + var.ref.size()/2;

    // we'll use this later to cut and pad the matrix
    string window_ref_seq = fasta_ref.getSubSequence(seq_name, begin_pos, window_length);

    int repeat_window_length = window_length * 8;
    int repeat_window_start = var.position-1 - repeat_window_length/2;
    string repeat_window = fasta_ref.getSubSequence(seq_name, repeat_window_start, repeat_window_length);;
    int callable_begin_pos = repeat_window_start + repeat_window_length/2 +1;
    int callable_end_pos = repeat_window_start + repeat_window_length/2 +1;
    for (auto& allele_seq : var.alleles) {
        auto f = callable_window(repeat_window_length/2 +2,
                                 repeat_window,
                                 allele_seq,
                                 3, min_repeat_entropy);
        //cerr << "callable window for " << allele_seq << " " << f.first << "-" << f.second << endl;
        callable_begin_pos = min((int)(repeat_window_start + f.first), callable_begin_pos);
        callable_end_pos = max((int)(repeat_window_start + f.second), callable_end_pos);
    }

    //cerr << "callable window for site " << callable_begin_pos << "-" << callable_end_pos << endl;
    
    //vcflib::VariantCallFile& graph_vcf;
    stringstream targetss;
    auto graph_begin_pos = var.position-1 - graph_window/2;
    auto graph_end_pos = var.position-1 + var.ref.size() + graph_window/2;
    targetss << seq_name << ":" << graph_begin_pos << "-" << graph_end_pos;
    auto target = targetss.str();
    set<string> allowed_variants;
    //cerr << target << " " << var.vrepr() << endl;
    allowed_variants.insert(var.vrepr());
    graph = vg::VG(graph_vcf,
                   fasta_ref,
                   target,
                   false,
                   1000,
                   0,
                   true, // don't manipulate the VCF records
                   false, false, false,
                   &allowed_variants);
    if (max_node_size > 0) {
        graph.dice_nodes(max_node_size); // force nodes to be 1bp
    }
    graph.compact_ids();
    int head_tail_id = 100;
    for (auto& n : graph.head_nodes()) {
        graph.swap_node_id(n, head_tail_id++);
    }
    for (auto& n : graph.tail_nodes()) {
        graph.swap_node_id(n, head_tail_id++);
    }
    // now set the ids of the ref and alts to keep sort invariance
    map<string, int> allele_seq_to_id;
    {
        int k = 0;
        for (auto& a : var.alleles) {
            allele_seq_to_id[a] = ++k;
        }
    }

    set<vg::id_t> allele_nodes;
    /// todo ... switch k to use fraction mapping to each allele
    graph.for_each_node([&](vg::Node* n) {
            if (!graph.is_head_node(n)
                && !graph.is_tail_node(n)) {
                // find out which allele it is
                // we're using a literal graphification of the VCF here due to options to vg construction
                // so we should be able to map from non-head, non-tail node to VCF allele
                auto f = allele_seq_to_id.find(n->sequence());
                if (f == allele_seq_to_id.end()) {
                    cerr << "could not find allele sequence for " << pb2json(*n) << endl;
                    cerr << var << endl;
                } else {
                    graph.swap_node_id(n, f->second + 200);
                    allele_nodes.insert(f->second + 200);
                }
            }
        });
    
    graph.rebuild_indexes();
    //graph.serialize_to_file("graphs/"+target+ ".vg");
    graph.for_each_node([&](vg::Node* n) {
            graph_coverage[n->id()] = 0;
            graph_weights[n->id()] = 0;
        });

    long int lowestReferenceBase = 0;
    long unsigned int referenceBases = 0;
    unsigned int currentRefSeqID = 0;
    bool biallelic_snp = var.alleles.size() == 2 && var.ref.size() == 1 && var.alleles.back().size() == 1;
    bool use_repeat_window = min_repeat_entropy && !biallelic_snp;

    // set up our readers
    set_region(bam_reader, seq_name, begin_pos, end_pos);
    // get the alignments at the locus
    // aligning them to the graph
    BamTools::BamAlignment aln;
    while (bam_reader.GetNextAlignment(aln)) {
        if (aln.IsMapped()) {
            if (!use_repeat_window) {
                alignments.push_back(aln);
            } else if (aln.Position <= callable_begin_pos
                       && aln.GetEndPosition() > callable_end_pos) {
                alignments.push_back(aln);
            }
        }
    }

    // now handle the graph region, which can be bigger
    set_region(bam_reader, seq_name, graph_begin_pos, graph_end_pos);
    //stringstream t; t << "graphs/" << target << ".reads";
    //ofstream oalns(t.str());
    while (bam_reader.GetNextAlignment(aln)) {
        auto vgaln = graph.align(aln.QueryBases);
        vgaln.set_quality(aln.Qualities);
        graph_alns.push_back(vgaln);
        //oalns << aln.QueryBases << endl;
    }
    //oalns.close();

    // handle the unitigs
    set_region(unitig_reader, seq_name, begin_pos, end_pos);
    //BamTools::BamAlignment aln;
    int unitig_count = 0;
    while (unitig_reader.GetNextAlignment(aln)) {
        if (aln.IsMapped()) {
            if (!use_repeat_window) {
                alignments.push_back(aln);
                ++unitig_count;
                unitigs.insert(&alignments.back());
            } else if (aln.Position <= callable_begin_pos
                       && aln.GetEndPosition() > callable_end_pos) {
                alignments.push_back(aln);
                ++unitig_count;
                unitigs.insert(&alignments.back());
            }
        }
    }

    // compress the alignment information into the graph
    for (auto& vgaln : graph_alns) {
        auto& path = vgaln.path();
        for (int i = 0; i < path.mapping_size(); ++i) {
            auto mapping = path.mapping(i);
            graph_coverage[mapping.position().node_id()]++;
        }
        auto qual_per_node = alignment_quality_per_node(vgaln);
        for (auto& n : qual_per_node) {
            graph_weights[n.first] += (double)n.second;
        }
    }

    double mapping_to_alleles = 0;
    for (auto& c : graph_coverage) {
        if (allele_nodes.count(c.first)) {
            mapping_to_alleles += c.second;
        }
    }

    if (mapping_to_alleles > 0) { // avoid -nan
        for (auto& c : graph_coverage) {
            if (allele_nodes.count(c.first)) {
                c.second /= (double) mapping_to_alleles;
            } else {
                c.second = 0;
            }
        }
    }

    // highest position
    int32_t min_pos = alignments.front().Position;
    int32_t max_pos = min_pos;

    for (auto& aln : alignments) {
        //cerr << "on alignment " << aln.Name << endl;
        int32_t endpos = aln.GetEndPosition();

        min_pos = min(aln.Position, min_pos);
        max_pos = max(endpos, max_pos);
        // iterate through the alignment
        // converting it into a series of alleles

        // record the qualities
        vector<prob_t> quals;
        assert(aln.Qualities.size() == aln.QueryBases.size());
        for (string::iterator c = aln.Qualities.begin(); c != aln.Qualities.end(); ++c) {
            if (exponentiate) {
                quals.push_back(
                    1-phred2float(
                        qualityChar2ShortInt(*c)));
            } else {
                quals.push_back(
                    qualityChar2ShortInt(*c));
            }
        }

        string refseq = fasta_ref.getSubSequence(referenceIDToName[aln.RefID],
                                                 aln.Position,
                                                 aln.GetEndPosition() - (aln.Position - 1));
        const string& readseq = aln.QueryBases;
        int rel_pos = aln.Position - this->begin_pos;

        int rp = 0; int sp = 0;

        vector<allele_t>& aln_alleles = alignment_alleles[&aln];

        vector<BamTools::CigarOp>::const_iterator cigarIter = aln.CigarData.begin();
        vector<BamTools::CigarOp>::const_iterator cigarEnd  = aln.CigarData.end();
        for ( ; cigarIter != cigarEnd; ++cigarIter ) {
            unsigned int len = cigarIter->Length;
            char t = cigarIter->Type;        
            switch (t) {
            case 'I':
            {
                auto iprobs = insertion_probs(quals, sp, len);
                for (int i = 0; i < len; ++i) {
                    aln_alleles.push_back(
                        allele_t("U",
                                 readseq.substr(sp + i, 1),
                                 rp + aln.Position-1,
                                 iprobs[i]));

                }
                sp += len;
            }
            break;
            case 'D':
            {
                auto dprobs = deletion_probs(quals, sp, len);
                for (int i = 0; i < len; ++i) {
                        aln_alleles.push_back(
                            allele_t(refseq.substr(rp + i, 1),
                                     "U",
                                     rp + i + aln.Position,
                                     dprobs[i]));
                }
                rp += len;
            }
            break;
            case 'X':
            case 'M':
            {
                for (int i = 0; i < len; ++i) {
                    aln_alleles.push_back(
                        allele_t(refseq.substr(rp + i, 1),
                                 readseq.substr(sp + i, 1),
                                 rp + i + aln.Position,
                                 quals[sp+i]));
                }
            }
            rp += len;
            sp += len;
            break;
            case 'S':
                // position is -1 if at the beginning
                // or +1 if at the end
                if (cigarIter == aln.CigarData.begin()) {
                    aln_alleles.push_back(allele_t("", "S", rp + i + aln.Position-2, len));
                } else {
                    aln_alleles.push_back(allele_t("", "S", rp + i + aln.Position+2, len));
                }
                sp += len;
                break;
            case 'H':
                // clipped sequence not present in the read
                break;
            default:
                cerr << "do not recognize cigar element " << t <<":"<< len << endl;
                break;
            }
        }
    }

    // make the reference haplotype
    for (size_t i = 0; i < window_length; ++i) {
        string base = window_ref_seq.substr(i, 1);
        reference.push_back(allele_t(base, base, begin_pos + i, 1));
    }

    // make each alt into a haplotype
    map<string, vector<allele_t> > vhaps;
    /*
    auto& vref = vhaps[var.ref];
    for (size_t i = 0; i < var.ref.size(); ++i) {
        string base = var.ref.substr(i, 1);
        vref.push_back(allele_t(base, base, var.position-1 + i, 1));
    }
    */
    int max_haplotype_length = 0;
    for (auto& allele : var.alleles) {
        max_haplotype_length = max((int)allele.size(), max_haplotype_length);
    }

    // handle out input haplotypes
    // note that parsedalternates is giving us 1-based positions
    bool has_insertion = false;
    bool has_deletion = false;
    for (auto& p : var.parsedAlternates()) {
        auto& valleles = vhaps[p.first];
        for (auto& a : p.second) {
            if (a.ref == a.alt && a.alt.size() > 1) {
                // break it apart
                for (size_t i = 0; i < a.ref.size(); ++i) {
                    valleles.push_back(allele_t(a.ref.substr(i,1),
                                                a.alt.substr(i,1),
                                                a.position+i-1, 1));
                }
            } else {
                // cluster insertions behind the previous base
                if (a.ref.empty()) {
                    has_insertion = true;
                    for (size_t i = 0; i < a.alt.size(); ++i) {
                        valleles.push_back(allele_t("U",
                                                    a.alt.substr(i,1),
                                                    a.position-2, 1));
                    }
                } else if (a.alt.empty()) {
                    has_deletion = true;
                    // deletions get broken into individual bases
                    for (size_t i = 0; i < a.ref.size(); ++i) {
                        valleles.push_back(allele_t(a.ref.substr(i,1),
                                                    "U",
                                                    a.position+i-1, 1));
                    }
                } else {
                    valleles.push_back(allele_t(a.ref,
                                                a.alt,
                                                a.position-1, 1));
                }
            }
        }
        while (valleles.size() < max_haplotype_length) {
            valleles.push_back(allele_t("",
                                        "U",
                                        valleles.back().position,
                                        1));
        }
    }

    // normalize away the VCF funk
    // by removing any reference-matching bases
    if (has_insertion || has_deletion) {
        set<string> first_bases;
        for (auto& v : vhaps) {
            auto& valleles = v.second;
            first_bases.insert(valleles[0].alt);
        }
        if (first_bases.size() == 1) {
            for (auto& v : vhaps) {
                auto& valleles = v.second;
                valleles.front().alt = "M";
            }
            // TODO it might be nice to re-center
            // but i'm not sure the right way to do it for insertions and deletions
        }
    }

    // for each sample
    // get the genotype
    vector<string> haplotype_seqs;
    for (auto& allele : var.alleles) {
        haplotypes.push_back(vhaps[allele]);
        haplotype_seqs.push_back(allele);
    }
    vector<string> genotype_seqs;
    int sid = 0; int gid = 0;
    for (auto& s : var.samples) {
        ++sid;
        auto& gtstr = s.second["GT"].front();
        auto gt = vcflib::decomposeGenotype(gtstr);
        for (auto& g : gt) {
            if (g.first != vcflib::NULL_ALLELE) {
                for (size_t i = 0; i < g.second; ++i){
                    genotypes.push_back(vhaps[var.alleles[g.first]]);
                    sample_id[gid++] = sid;
                    genotype_seqs.push_back(var.alleles[g.first]);
                }
            }
        }
    }

    // for all the info fields
    for (auto& f : var.info) {
        // what kind of field is this?
        auto field_name = f.first;
        auto field_type = var.infoType(field_name);
        auto& fields = f.second;
        int i = 0;
        for (auto& field : fields) {
            ++i;
            stringstream k;
            k << input_name << field_name << "_" << i;
            string key = k.str();
            try {
                if (field_type == vcflib::FIELD_FLOAT
                    || field_type == vcflib::FIELD_INTEGER) {
                    call_info_num[key] = stod(field);
                } else if (field_type == vcflib::FIELD_BOOL
                           || field_type == vcflib::FIELD_STRING) {
                    call_info_str[key] = field;
                }
            } catch (...) {
                // do nothing if the field is invalid
                // wtf -- only VCF would be impossible to get right
            }
        }
    }

    // do the same for QUAL
    call_info_num[input_name + "QUAL"] = var.quality;

    // find alleles above a threshold rate of incidence
    for (auto a = alignment_alleles.begin(); a != alignment_alleles.end(); ++a) {
        vector<allele_t>& aln_alleles = a->second;
        map<int, int> pos_count;
        for (auto& allele : aln_alleles) {
            allele_counts[allele.str()][pos_count[allele.position]++]++;
        }
    }

    // keep only those alleles > our threshold
    for (auto a = alignment_alleles.begin(); a != alignment_alleles.end(); ++a) {
        vector<allele_t>& aln_alleles = a->second;
        vector<allele_t> filtered_alleles;
        int last_pos = 0;
        map<int, int> pos_count;
        for (auto& allele : aln_alleles) {
            if (allele_counts[allele.str()][pos_count[allele.position]++] < min_allele_count) {
                if (last_pos && last_pos != allele.position) {
                    filtered_alleles.push_back(allele_t("", "M", allele.position, 1));
                }
            } else {
                filtered_alleles.push_back(allele);
            }
            last_pos = allele.position;
        }
        aln_alleles = filtered_alleles;
    }

    map<int32_t, size_t> pos_max_length;

    
    // trim the reads to the right size and determine the maximum indel length at each reference position
    for (auto a = alignment_alleles.begin(); a != alignment_alleles.end(); ++a) {
        vector<allele_t>& aln_alleles = a->second;
        map<int32_t, size_t> pos_counts;
        for (auto& allele : aln_alleles) {
            ++pos_counts[allele.position];
        }
        // now record the indels
        for (auto p : pos_counts) {
            pos_max_length[p.first] = max(pos_max_length[p.first], p.second);
        }
    }

    // do for the alternate haps too
    for (auto& v : vhaps) {
        map<int32_t, size_t> pos_counts;
        for (auto& allele : v.second) {
            ++pos_counts[allele.position];
        }
        // now record the indels
        for (auto p : pos_counts) {
            pos_max_length[p.first] = max(pos_max_length[p.first], p.second);
        }
    }

    // maps position/indels into offsets
    // pair<i, 0> -> reference
    // pari<i, j> -> jth insertion after base
    map<pair<int32_t, size_t>, size_t> pos_proj;
    size_t j = 0;
    for (auto p : pos_max_length) {
        int32_t pos = p.first;
        for (int32_t i = 0; i < p.second; ++i) {
            pos_proj[make_pair(pos, i)] = j++;
        }
    }

    /*
    map<size_t, pair<int32_t, size_t> > inv_pos_proj;
    for (auto p : pos_proj) {
        inv_pos_proj[p.second] = p.first;
        cerr << p.first.first << ":" << p.first.second << " " << p.second << endl;
    }
    */

    // convert positions into the new frame
    for (auto a = alignment_alleles.begin(); a != alignment_alleles.end(); ++a) {
        project_positions(a->second, pos_proj);
    }
    // same for ref
    project_positions(reference, pos_proj);
    // and genotype/haps
    for (auto& hap : haplotypes) {
        project_positions(hap, pos_proj);
    }
    for (auto& hap : genotypes) {
        project_positions(hap, pos_proj);
    }

    // get the min/max of the vector
    // the min should be 0
    // the max tells us how wide the MSA should be
    pos_t msav_min = pos_proj.begin()->second;
    pos_t msav_max = pos_proj.rbegin()->second;
    
    // where is the new center
    pos_t center = pos_proj[make_pair(center_pos, 0)];// + shift_center;
    //cerr << "center is " << center << endl;
    pos_t bal_min = max(center - window_length/2, (size_t)0);
    pos_t bal_max = bal_min + window_length;

    // re-center
    // re-strip out our limits
    // add the missing bases
    // add the gap bases
    vector<alignment_t*> to_erase;
    for (auto a = alignment_alleles.begin(); a != alignment_alleles.end(); ++a) {
        //cerr << a->first->Name << endl;
        a->second = pad_alleles(a->second, bal_min, bal_max);
        if (a->second.empty()) to_erase.push_back(a->first);
    }
    for (auto e : to_erase) alignment_alleles.erase(e);

    reference = pad_alleles(reference, bal_min, bal_max);
    
    for (auto& hap : haplotypes) {
        hap = pad_alleles(hap, bal_min, bal_max);
    }
    for (auto& hap : genotypes) {
        hap = pad_alleles(hap, bal_min, bal_max);
    }

    if (assume_ref) {
        // determine limits
        // and pad within them
        missing_to_ref(haplotypes);
        missing_to_ref(genotypes);
    }

    // optionally force the reference matching alleles to be R
    if (!show_bases) {
        for (auto a = alignment_alleles.begin(); a != alignment_alleles.end(); ++a) {
            flatten_to_ref(a->second);
        }
        for (auto& hap : haplotypes) {
            flatten_to_ref(hap);
        }
        for (auto& hap : genotypes) {
            flatten_to_ref(hap);
        }
    }

    for (auto& aln : alignments) {
        if (alignment_alleles.find(&aln) == alignment_alleles.end()) continue;
        auto h = alignment_alleles[&aln];
        missing_counts[&aln] = missing_count(h);
    }
    
    alignment_count = 0;
    for (auto& aln : alignments) {
        if (alignment_alleles.find(&aln) == alignment_alleles.end()) continue;
        if (full_overlap && missing_counts[&aln] > 0) continue;
        ++alignment_count;
    }

    // establish the allele/hap/ref matches
    for (auto& aln : alignments) {
        if (alignment_alleles.find(&aln) == alignment_alleles.end()) continue;
        if (full_overlap && missing_counts[&aln] > 0) continue;
        int i = 0;
        auto& weight = matches[&aln];
        for (auto& hap : haplotypes) {
            weight[i] = pairwise_identity(
                alignment_alleles[&aln],
                hap);
            ++i;
        }
    }

    // sum up the quality support for the allele/hap/ref matches
    for (auto& aln : alignments) {
        if (alignment_alleles.find(&aln) == alignment_alleles.end()) continue;
        if (full_overlap && missing_counts[&aln] > 0) continue;
        int i = 0;
        auto& weight = qualsum[&aln];
        for (auto& hap : haplotypes) {
            weight[i] = pairwise_qualsum(
                alignment_alleles[&aln],
                hap);
            ++i;
        }
    }

    for (int i = 1; i <= 7; ++i) {
        likelihoods[i] = 1;
    }

    for (auto& aln : alignments) {
        if (alignment_alleles.find(&aln) == alignment_alleles.end()) continue;
        if (full_overlap && missing_counts[&aln] > 0) continue;
        int i = 0;
        auto& prob = prob_aln_given_genotype[&aln];
        // for all possible genotype
        // estimate prob(aln | gentoype)
        bool supports_anything = false;
        for (int i = 1; i <= 7; ++i) {
            prob[i] = prob_aln_gt(&aln, i);
            if (prob[i] > 0) supports_anything = true;
        }
        if (supports_anything) {
            for (int i = 1; i <= 7; ++i) {
                likelihoods[i] *= prob[i];
            }
        }
    }

    double maxlikelihood = 0;
    for (int i = 1; i <= 7; ++i) {
        //likelihoods[i] = float2phred(1-likelihoods[i]);
        maxlikelihood = max(maxlikelihood, likelihoods[i]);
    }
    if (maxlikelihood > 0) {
        for (int i = 1; i <= 7; ++i) {
            //likelihoods[i] = float2phred(1-likelihoods[i]);
            auto l = likelihoods[i]/maxlikelihood;
            if (l < 1e-3) l = 0;
            likelihoods[i] = l;
        }
    }

    auto aln_sort = [&](alignment_t* a1, alignment_t* a2) {
        auto m1 = missing_counts[a1];
        auto m2 = missing_counts[a2];
        if (m1 < m2) {
            return true;
        } else if (m1 == m2) {
            return a1->Position < a2->Position;
        } else {
            return false;
        }
    };
    
    // collect allele supports
    for (auto& aln : alignments) {
        if (alignment_alleles.find(&aln) == alignment_alleles.end()) continue;
        if (full_overlap && missing_counts[&aln] > 0) continue;
        // use all the best supports
        map<double, vector<int> > bests;
        for (int i = 0; i < haplotypes.size(); ++i) {
            bests[matches[&aln][i]].push_back(i);
        }
        if (bests.rbegin()->first) {
            for (auto i : bests.rbegin()->second) {
                allele_support[i][1-matches[&aln][i]].push_back(&aln);
            }
        } else {
            for (auto i : bests.rbegin()->second) {
                // 8 is the magic "extra" namespace
                allele_support[8][1-matches[&aln][i]].push_back(&aln);
            }
        }
    }
    int x = 0, y = 0;
    for (auto& supp : allele_support) {
        for (auto& hsup : supp.second) {
            auto& sup = hsup.second;
            // sort it baby
            std::sort(sup.begin(), sup.end(), aln_sort);
            sup.erase(std::unique(sup.begin(), sup.end()), sup.end());
        }
    }

    // collect soft clips
    //vector<alignment_t*> softclipped;
    for (auto& aln : alignments) {
        if (alignment_alleles.find(&aln) == alignment_alleles.end()) continue;
        if (full_overlap && missing_counts[&aln] > 0) continue;
        // if we have a softclip
        if (has_softclip(alignment_alleles[&aln])) {
            softclipped.push_back(&aln);
        }
    }
    std::sort(softclipped.begin(), softclipped.end(), aln_sort);
    softclipped.erase(std::unique(softclipped.begin(), softclipped.end()),
                      softclipped.end());
    if (max_depth && softclipped.size() > max_depth) {
        softclipped.erase(softclipped.begin() + max_depth, softclipped.end());
    }

    // organize the ordered alignments
    for (auto& supp : allele_support) {
        int i = 0;
        int j = 0;
        int u = 0; // unitigs
        for (auto& hsup : supp.second) {
            auto& sup = hsup.second;
            for (auto& aln : sup) {
                stringstream ss;
                // we limit ourselves to only 7 alleles, 1 softclip (=9) and 1 degenerate (OB=8)
                if (unitigs.count(aln)) {
                    ss << min(supp.first, 8) << "u" << u++;
                    if (max_depth && u+i > max_depth) break;
                } else {
                    if (aln->IsReverseStrand()) {
                        if (max_depth && i >= max_depth) continue;
                        ss << min(supp.first, 8) << "-" << i++;
                    } else {
                        if (max_depth && j >= max_depth) continue;
                        ss << min(supp.first, 8) << "+" << j++;
                    }
                }
                alignment_groups[aln].push_back(ss.str());
                if (!unitigs.count(aln)) {
                    grouped_normal_alignments.push_back(make_pair(ss.str(), aln));
                } else {
                    grouped_unitig_alignments.push_back(make_pair(ss.str(), aln));
                }
            }
        }
    }

    // and the soft clips
    {
        int i = 0; // forward strand
        int j = 0; // reverse strand
        int u = 0; // unitigs
        for (auto& aln : softclipped) {
            stringstream ss;
            // we keep soft clips in the special namespace 9
            if (unitigs.count(aln)) {
                ss << 9 << "u" << u++;
            } else {
                if (aln->IsReverseStrand()) {
                    ss << 9 << "-" << i++;
                } else {
                    ss << 9 << "+" << j++;
                }
            }
            alignment_groups[aln].push_back(ss.str());
            if (!unitigs.count(aln)) {
                grouped_normal_alignments.push_back(make_pair(ss.str(), aln));
            } else {
                grouped_unitig_alignments.push_back(make_pair(ss.str(), aln));
            }
        }
    }

    // make the label that represents our hhga site
    // and which we will later use to project back into VCF
    stringstream vrep;
    vrep << var.sequenceName << "_" << var.position;
    vrep << "_" << var.ref << "_";
    vrep << join(var.alt, ",");
    repr = vrep.str();

}

void HHGA::flatten_to_ref(vector<allele_t>& alleles) {
    for (auto& allele : alleles) {
        if (allele.alt != "U"
            && allele.alt != "M"
            && allele.alt == allele.ref) {
            allele.alt = "R";
        }
    }
}

void HHGA::strandify(vector<allele_t>& alleles, bool is_rev) {
    if (is_rev) {
        for (auto& allele : alleles) {
            // set to lower case
            std::transform(allele.alt.begin(), allele.alt.end(), allele.alt.begin(), ::tolower);
        }
    }
}

void HHGA::project_positions(vector<allele_t>& aln_alleles,
                             map<pair<int32_t, size_t>, size_t>& pos_proj) {
    // adjust the allele positions
    // if the new position is not the same as the last
    // set j = 0
    size_t j = 0;
    pos_t last = aln_alleles.begin()->position;
    for (auto& allele : aln_alleles) {
        if (last != allele.position) j = 0;
        last = allele.position;
        allele.position = pos_proj[make_pair(allele.position, j++)];
    }
}

vector<allele_t> HHGA::pad_alleles(vector<allele_t> aln_alleles,
                                   pos_t bal_min, pos_t bal_max) {
    vector<allele_t> padded;
    if (aln_alleles.empty()) return padded;

    // remove the bits outside the window
    pos_t aln_start = aln_alleles.front().position;
    pos_t aln_end = aln_alleles.back().position;

    // pad the sides
    if (aln_start > bal_max) return padded;
    // pad the beginning with "missing" features
    for (int32_t q = bal_min; q < aln_start; ++q) {
        padded.push_back(allele_t("", "M", q, 1));
    }
    // pad the gaps
    bool first = true;
    pos_t last = aln_alleles.front().position;
    for (auto& allele : aln_alleles) {
        if (!first &&
            last+1 != allele.position) {
            for (int32_t j = 0; j < allele.position - (last + 1); ++j) {
                padded.push_back(allele_t("", "U", j + last + 1, 1));
            }
        }
        last = allele.position;
        padded.push_back(allele);
        first = false;
    }
    // pad the end with "missing" features
    for (int32_t q = aln_end+1; q < bal_max; ++q) {
        padded.push_back(allele_t("", "M", q, 1));
    }

    int i = padded.front().position;
    for (auto& p : padded) {
        p.position = i++;
    }

    padded.erase(std::remove_if(padded.begin(), padded.end(),
                                [&](const allele_t& allele) {
                                    return allele.position < bal_min || allele.position >= bal_max;
                                }),
                      padded.end());

    return padded;
}

const string HHGA::str(void) {
    //return std::to_string(alleles.size());
    stringstream out;
    //out << std::fixed << std::setprecision(1);
    out << repr << endl;
    out << "reference          ";
    for (auto& allele : reference) {
        if (allele.alt == "M") out << " ";
        else if (allele.alt == "U") out << "-";
        else if (allele.alt == "R") out << ".";
        else out << allele.alt;
    }
    out << endl;
    for (auto& hap : haplotypes) {
        out << "hap                ";
        for (auto& allele : hap) {
            if (allele.alt == "M") out << " ";
            else if (allele.alt == "U") out << "-";
            else if (allele.alt == "R") out << ".";
            else out << allele.alt;
        }
        out << endl;
    }
    for (auto& hap : genotypes) {
        out << "geno               ";
        for (auto& allele : hap) {
            if (allele.alt == "M") out << " ";
            else if (allele.alt == "U") out << "-";
            else if (allele.alt == "R") out << ".";
            else out << allele.alt;
        }
        out << endl;
    }

    auto do_alignment = [&](alignment_t* aln,
                            const string& name,
                            const string& prepend) {
        out << prepend << " ";
        // print out the stuff
        if (aln->IsReverseStrand())     out << "S"; else out << "s";
        if (aln->IsMateReverseStrand()) out << "O"; else out << "o";
        if (aln->IsDuplicate())         out << "D"; else out << "d";
        if (aln->IsFailedQC())          out << "Q"; else out << "q";
        if (aln->IsFirstMate())         out << "F"; else out << "f";
        if (aln->IsSecondMate())        out << "X"; else out << "x";
        if (aln->IsMateMapped())        out << "Y"; else out << "y";
        if (aln->IsPaired())            out << "P"; else out << "p";
        if (aln->IsPrimaryAlignment())  out << "Z"; else out << "z";
        if (aln->IsProperPair())        out << "I"; else out << "i";
        out << "  ";
        for (auto& allele : alignment_alleles[aln]) {
            if (allele.alt == "M") out << " ";
            else if (allele.alt == "U") out << "-";
            else if (allele.alt == "R") out << ".";
            else out << allele.alt;
        }
        // the grouping
        out << " " << name;

        // now the matches
        out << " : ";
        for (auto w : matches[aln]) {
            out << w.second << " ";
        }
        out << ": ";
        // now the qualsums
        for (auto w : qualsum[aln]) {
            out << w.second << " ";
        }
        out << ": ";
        for (auto w : prob_aln_given_genotype[aln]) {
            out << w.second << " ";
        }
        out << ": ";
        out << aln->MapQuality;
        out << " " << aln->Name;
        out << endl;
    };
    
    for (auto g : grouped_normal_alignments) {
        auto& name = g.first;
        auto& aln = g.second;
        do_alignment(aln, name, "aln   ");
    }

    for (auto g : grouped_unitig_alignments) {
        auto& name = g.first;
        auto& aln = g.second;
        do_alignment(aln, name, "unitig");
    }

    out << "p(obs|genotype)s ";
    for (auto& l : likelihoods) {
        out << l.first << ":" << l.second << " ";
    }
    out << endl;

    out << "graph_coverage ";
    for (auto& w : graph_coverage) {
        out << w.first << "N:" << w.second << " ";
    }
    out << endl;

    /*
    out << "graph_weight ";
    for (auto& w : graph_weights) {
        out << w.first << "N:" << w.second << " ";
    }
    out << endl;
    */

    // now handle caller input features
    for (auto& f : call_info_num) {
        out << f.first << ":" << f.second << " ";
    }
    for (auto& f : call_info_str) {
        out << f.first << ":" << f.second << " ";
    }
    out << endl;
    
    return out.str();
}

const string HHGA::vw(void) {
    stringstream out;
    // write the class of the example
    out << label << " ";
    out << "'" << repr << " ";
    // do the ref
    out << "|ref ";
    size_t idx = 0;
    for (auto& allele : reference) {
        out << ++idx << allele.alt << ":" << allele.prob << " ";;
    }
    // do the haps
    size_t i = 1;
    for (auto& hap : haplotypes) {
        out << "|hap" << i << " ";
        ++i;
        idx = 0;
        for (auto& allele : hap) {
            if (allele.alt != "M") {
                out << ++idx << allele.alt << ":" << allele.prob << " ";;
            }
        }
    }
    i = 1;
    int gid = 0;
    for (auto& geno : genotypes) {
        out << "|geno" << i << " ";
        ++i;
        idx = 0;
        for (auto& allele : geno) {
            if (allele.alt != "M") {
                out << ++idx << allele.alt << ":" << allele.prob << " ";
            }
        }
    }

    // do the row wise alignment features
    for (auto g : grouped_normal_alignments) {
        auto& name = g.first;
        auto& aln = g.second;
        
        out << "|aln" << name << " ";
        idx = 0;
        for (auto& allele : alignment_alleles[aln]) {
            out << ++idx << allele.alt << ":" << allele.prob << " ";
        }
    }

    // tranposed into colum wise
    int coln = 0;
    for (auto& allele : reference) { //this coud just be the lenght not sure where to get it from
        out << "|col" << coln << " ";
        for (auto g : grouped_normal_alignments) {
            auto&  alle = alignment_alleles[g.second][coln];
            out << alle.alt << ":" << alle.prob << " ";
        }
        ++coln;
    }
    
    for (auto g : grouped_normal_alignments) {
        auto& name = g.first;
        auto& aln = g.second;
        out << "|match" << name << " ";
        // match properties
        for (auto w : matches[aln]) {
            out << w.first+1 << "H:" << w.second << " ";
        }
    }

    for (auto g : grouped_normal_alignments) {
        auto& name = g.first;
        auto& aln = g.second;
        out << "|qual" << name << " ";
        // match properties
        for (auto w : qualsum[aln]) {
            out << w.first+1 << "H:" << w.second << " ";
        }
    }

    // do the row wise unitig features
    for (auto g : grouped_unitig_alignments) {
        auto& name = g.first;
        auto& aln = g.second;
        out << "|unitig" << name << " ";
        idx = 0;
        for (auto& allele : alignment_alleles[aln]) {
            out << ++idx << allele.alt << ":" << 1 << " "; //allele.prob << " ";
        }
    }

    for (auto g : grouped_unitig_alignments) {
        auto& name = g.first;
        auto& aln = g.second;
        out << "|xmatch" << name << " ";
        // match properties
        for (auto w : matches[aln]) {
            out << w.first+1 << "H:" << w.second << " ";
        }
    }

    out << "|depth ";
    out << "bam:" << alignment_count << " ";
    out << "graph:" << graph_alns.size() << " ";

    out << "|likelihood ";
    for (auto& l : likelihoods) {
        out << l.first << "G:" << l.second << " ";
    }

    for (auto g : grouped_normal_alignments) {
        auto& name = g.first;
        auto& aln = g.second;
        out << "|properties" << name << " ";
        if (exponentiate) {
            out << "mapqual:" << 1-phred2float(min(aln->MapQuality, (uint16_t)60)) << " ";
        } else {
            out << "mapqual:" << aln->MapQuality << " ";
        }
        // handle flags
        if (aln->IsReverseStrand())     out << "strand:1"; else out << "strand:0"; out << " ";
        if (aln->IsMateReverseStrand()) out << "ostrand:1"; else out << "ostrand:0"; out << " ";
        if (aln->IsDuplicate())         out << "dup:1"; else out << "dup:0"; out << " ";
        if (aln->IsFailedQC())          out << "qcfail:1"; else out << "qcfail:0"; out << " ";
        if (aln->IsFirstMate())         out << "fmate:1"; else out << "fmate:0"; out << " ";
        if (aln->IsSecondMate())        out << "xmate:1"; else out << "xmate:0"; out << " ";
        if (aln->IsMateMapped())        out << "ymap:1"; else out << "ymap:0"; out << " ";
        if (aln->IsPaired())            out << "paired:1"; else out << "paired:0"; out << " ";
        if (aln->IsPrimaryAlignment())  out << "zprimary:1"; else out << "zprimary:0"; out << " ";
        if (aln->IsProperPair())        out << "iproper:1"; else out << "iproper:0"; out << " ";
    }

    out << "|vgraph ";
    graph.for_each_node([&](vg::Node* n) {
            auto& seq = n->sequence();
            for (int j = 0; j < seq.size(); ++j) {
                out << n->id() << "_" << j << "_" << seq[j] << " ";
            }
        });

    out << "|kgraph ";
    for (auto& w : graph_coverage) {
        out << w.first << "C:" << w.second << " ";
    }

    out << "|software ";
    // now handle caller input features
    for (auto& f : call_info_num) {
        out << f.first << ":" << f.second << " ";
    }

    // and the alignment supports
    
    
    /*
    for (auto& f : call_info_str) {
        out << f.first << "_" << f.second << " ";
    }
    */
    return out.str();
}


vector<prob_t> deletion_probs(const vector<prob_t>& quals, size_t sp, size_t l) {
    
    // because deletions have no quality information,
    // use the surrounding sequence quality as a proxy
    // to provide quality scores of equivalent magnitude to insertions,
    // take N bp, right-centered on the position of the deletion
    // this function ensures that the window is fully contained within the read

    int spanstart = 0;

    // this is used to calculate the quality string adding 2bp grounds
    // the indel in the surrounding sequence, which it is dependent
    // upon
    int L = l + 2;

    // if the event is somehow longer than the read (???)
    // then we need to bound it at the read length
    if (L > quals.size()) {
        L = quals.size();
        spanstart = 0;
    } else {
        if (sp < (L / 2)) {
            // if the read pointer is less than half of the deletion length
            // we need to avoid running past the end of the read, so bound the spanstart
            spanstart = 0;
        } else {
            spanstart = sp - (L / 2);
        }
        // set upper bound to the string length
        if (spanstart + L > quals.size()) {
            spanstart = quals.size() - L;
        }
    }

    auto qual_begin = (quals.begin() + spanstart);
    auto qual_end = qual_begin + L;

    vector<prob_t> del_quals;
    while (qual_begin != qual_end) {
        del_quals.push_back(*qual_begin++);
    }

    return del_quals;
}

vector<prob_t> insertion_probs(const vector<prob_t>& quals, size_t sp, size_t l) {

    // insertion quality is taken as the minimum of
    // the inserted bases and the two nearest flanking ones
    // this function ensures that the window is fully contained within the read

    int spanstart = 0;
        
    // this is used to calculate the quality string adding 2bp grounds
    // the indel in the surrounding sequence, which it is dependent
    // upon
    int L = l + 2;

    // if the event is somehow longer than the read (???)
    // then we need to bound it at the read length        
    if (L > quals.size()) {
        L = quals.size();
        spanstart = 0;
    } else {
        // set lower bound to 0
        if (sp < 1) {
            spanstart = 0;
        } else {
            // otherwise set it to one back
            spanstart = sp - 1;
        }
        // set upper bound to the string length
        if (spanstart + L > quals.size()) {
            spanstart = quals.size() - L;
        }
    }

    auto qual_begin = (quals.begin() + spanstart);
    auto qual_end = qual_begin + L;

    vector<prob_t> ins_quals;
    while (qual_begin != qual_end) {
        ins_quals.push_back(*qual_begin++);
    }

    return ins_quals;
}

ostream& operator<<(ostream& out, allele_t& var) {
    out << var.position << ":" << var.ref << "/" << var.alt << ":" << var.prob;
    return out;
}

allele_t operator+(const allele_t& a, const allele_t& b) {
    return allele_t(a.ref + b.ref, a.alt + b.alt, a.position, a.prob * b.prob);
}

bool operator<(const allele_t& a, const allele_t& b) {
    return a.repr < b.repr;
}

}
