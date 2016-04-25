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

string label_for_genotype(const string& gt) {
    if        (gt == "0/0") {
        return "1";
    } else if (gt == "0/1") {
        return "2";
    } else if (gt == "0/2") {
        return "3";
    } else if (gt == "1/1") {
        return "4";
    } else if (gt == "1/2") {
        return "5";
    } else if (gt == "2/2") {
        return "6";
    } else {
        cerr << "warning: unknown genotype '" << gt << "'" << endl;
        return "7";
    }
}

string genotype_for_label(const string& gt, int alt_count) {
    if        (gt == "1") {
        return "0/0";
    } else if (gt == "2") {
        return "0/1";
    } else if (gt == "3") {
        if (alt_count > 1) {
            return "0/2";
        } else {
            // emit a valid het if we don't have 2 alts
            return "0/1";
        }
    } else if (gt == "4") {
        return "1/1";
    } else if (gt == "5") {
        if (alt_count > 1) {
            return "1/2";
        } else {
            // emit a valid het if we don't have 2 alts
            return "0/1";
        }
    } else if (gt == "6") {
        if (alt_count > 1) {
            return "2/2";
        } else {
            // emit a valid hom if we don't have 2 alts
            return "1/1";
        }
    } else {
        //cerr << "warning: unknown genotype '" << gt << "'" << " expected one of 0/0, 0/1, 0/2, 1/1, 1/2, 2/2" << endl;
        return "./.";
    }
}

map<int, double> test_labels(int alt_count) {
    map<int, double> labels;
    for (int i = 1; i < alt_count+1; ++i) {
        labels[i] = 0;
    }
    return labels;
}

string multiclass_label_for_genotype(const string& gt) {
    vector<string> label;
    auto labg = labels_for_genotype(gt);
    for (auto m : labg) {
        stringstream lss;
        lss << m.first << ":" << m.second;
        label.push_back(lss.str());
    }
    return join(label, " ");
}

map<int, double> labels_for_genotype(const string& gt) {
    map<int, double> labels;
    for (auto allele : vcflib::decomposeGenotype(gt)) {
        labels[allele.first+1] = allele.second;
    }
    return labels;
}

string genotype_for_labels(const map<int, double>& gt,
                           int alt_count) {
    // count up weight until we get to 2
    map<double, vector<int> > allele_by_weight;
    for (auto& g : gt) {
        allele_by_weight[g.second].push_back(g.first);
    }
    vector<pair<double, int> > flattened_by_weight;
    for (auto& a : allele_by_weight) {
        for (auto& v : a.second) {
            flattened_by_weight.push_back(make_pair(a.first, v));
        }
    }
    double weight = 0;
    vector<string> gts;
    for (auto& w : flattened_by_weight) {
        int count = round(w.first);
        for (int i = 0; i < count; ++i) {
            gts.push_back(std::to_string(w.second-1));
        }
        weight += count;
        if (weight >= 2) break; // diploid
    }
    return join(gts, "/");
}

double pairwise_identity(const vector<allele_t>& h1, const vector<allele_t>& h2) {
    int count = 0;
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
    int covered = 0;
    for (auto& p : p1) {
        auto a1 = p1[p.first]->alt;
        if (p2.find(p.first) == p2.end()) continue;
        auto a2 = p2[p.first]->alt;
        if (a1 == "M" || a2 == "M") continue;
        ++covered;
        if (a1 == a2) {
            ++count;
        }
    }
    return (possible ? (double) count / (double) possible : 0);
}

int missing_count(const vector<allele_t>& hap) {
    int m = 0;
    for (auto& a : hap) if (a.alt == "M") ++m;
    return m;
}

HHGA::HHGA(size_t window_length,
           BamTools::BamMultiReader& bam_reader,
           FastaReference& fasta_ref,
           vcflib::Variant& var,
           const string& input_name,
           const string& class_label,
           const string& gt_class,
           bool multiclass,
           bool expon,
           bool show_bases,
           bool assume_ref) { // assumes the haplotypes are ref everywhere

    exponentiate = expon;

    if (!gt_class.empty()) {
        // convert the genotype into
        // require that it be in
        // 0/0, 0/1, 1/1, 0/2, 1/2, 2/2
        auto gt = var.samples[var.sampleNames.front()][gt_class].front();
        if (!multiclass) {
            label = label_for_genotype(gt);
        } else {
            // handles generic case with many alleles
            label = multiclass_label_for_genotype(gt);
        }
    } else {
        label = class_label;
    }

    // store the names of all the reference sequences in the BAM file
    map<int, string> referenceIDToName;
    vector<BamTools::RefData> referenceSequences = bam_reader.GetReferenceData();
    int i = 0;
    for (BamTools::RefVector::iterator r = referenceSequences.begin(); r != referenceSequences.end(); ++r) {
        referenceIDToName[i] = r->RefName;
        ++i;
    }

    int32_t begin_pos = var.position-1 - window_length/2;
    int32_t end_pos = begin_pos + window_length;
    string seq_name = var.sequenceName;
    int32_t center_pos = var.position-1;//begin_pos + (end_pos - begin_pos) / 2;

    // we'll use this later to cut and pad the matrix
    string window_ref_seq = fasta_ref.getSubSequence(seq_name, begin_pos, window_length);

    // set up our readers
    set_region(bam_reader, seq_name, begin_pos, end_pos);

    long int lowestReferenceBase = 0;
    long unsigned int referenceBases = 0;
    unsigned int currentRefSeqID = 0;

    // get the alignments at the locus
    BamTools::BamAlignment aln;
    while (bam_reader.GetNextAlignment(aln)) {
        if (aln.IsMapped()) {
            alignments.push_back(aln);
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

    // handle out input haplotypes
    // note that parsedalternates is giving us 1-based positions
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
                    for (size_t i = 0; i < a.alt.size(); ++i) {
                        valleles.push_back(allele_t("U",
                                                    a.alt.substr(i,1),
                                                    a.position-2, 1));
                    }
                } else if (a.alt.empty()) {
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

    map<int32_t, size_t> pos_max_length;

    // trim the reads to the right size and determine the maximum indel length at each reference position
    for (auto a = alignment_alleles.begin(); a != alignment_alleles.end(); ++a) {
        vector<allele_t>& aln_alleles = a->second;
        aln_alleles.erase(std::remove_if(aln_alleles.begin(), aln_alleles.end(),
                                         [&](const allele_t& allele) {
                                             return allele.position < begin_pos || allele.position >= end_pos;
                                         }),
                          aln_alleles.end());
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
    pos_t center = pos_proj[make_pair(center_pos, 0)];
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
        for (auto& hap : haplotypes) {
            size_t i = 0;
            for (vector<allele_t>::iterator a = hap.begin(); a != hap.end(); ++a, ++i) {
                if (a->alt == "M") *a = reference[i];
            }
        }
        for (auto& hap : genotypes) {
            size_t i = 0;
            for (vector<allele_t>::iterator a = hap.begin(); a != hap.end(); ++a, ++i) {
                if (a->alt == "M") *a = reference[i];
            }
        }
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

    // establish the allele/hap/ref matches
    for (auto& aln : alignments) {
        if (alignment_alleles.find(&aln) == alignment_alleles.end()) continue;
        int i = 0;
        auto& weight = matches[&aln];
        for (auto& hap : haplotypes) {
            weight[i] = pairwise_identity(
                alignment_alleles[&aln],
                hap);
            ++i;
        }
    }

    // sort the alignments by the fraction of match and then by position
    // TODO record an alignment sort order
    // and use it on output
    for (auto& aln : alignments) {
        if (alignment_alleles.find(&aln) == alignment_alleles.end()) continue;
        ordered_alignments.push_back(&aln);
    }
    // sort by mismatch count
    // then by position
    std::sort(ordered_alignments.begin(), ordered_alignments.end(),
              [&](alignment_t* a1, alignment_t* a2) {
                  auto& h1 = alignment_alleles[a1];
                  auto& h2 = alignment_alleles[a2];
                  auto m1 = missing_count(h1);
                  auto m2 = missing_count(h2);
                  if (m1 < m2) {
                      return true;
                  } else if (m1 == m2) {
                      return a1->Position < a2->Position;
                  } else {
                      return false;
                  }
              });

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
    // remove the bits outside the window
    /*
    aln_alleles.erase(std::remove_if(aln_alleles.begin(), aln_alleles.end(),
                                     [&](const allele_t& allele) {
                                         return allele.position < bal_min || allele.position >= bal_max;
                                     }),
                      aln_alleles.end());
    if (aln_alleles.empty()) return padded;
    */
    /*
    cerr << "pre:    ";
    for (auto& a : aln_alleles) cerr << a << " ";
    cerr << endl;
    */
    // pad the sides
    pos_t aln_start = aln_alleles.front().position;
    pos_t aln_end = aln_alleles.back().position;
    // pad the beginning with "missing" features
    for (int32_t q = bal_min; q != aln_start; ++q) {
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

    padded.erase(std::remove_if(padded.begin(), padded.end(),
                                [&](const allele_t& allele) {
                                    return allele.position < bal_min || allele.position >= bal_max;
                                }),
                      padded.end());
    //if (aln_alleles.empty()) return padded;
    
    /*
    cerr << "padded: ";
    for (auto& a : padded) cerr << a << " ";
    cerr << endl;
    */
    return padded;
}

const string HHGA::str(void) {
    //return std::to_string(alleles.size());
    stringstream out;
    out << std::fixed << std::setprecision(1);
    out << repr << endl;
    out << "reference   ";
    for (auto& allele : reference) {
        if (allele.alt == "M") out << " ";
        else if (allele.alt == "U") out << "-";
        else if (allele.alt == "R") out << ".";
        else out << allele.alt;
    }
    out << endl;
    for (auto& hap : haplotypes) {
        out << "hap         ";
        for (auto& allele : hap) {
            if (allele.alt == "M") out << " ";
            else if (allele.alt == "U") out << "-";
            else if (allele.alt == "R") out << ".";
            else out << allele.alt;
        }
        out << endl;
    }
    for (auto& hap : genotypes) {
        out << "geno        ";
        for (auto& allele : hap) {
            if (allele.alt == "M") out << " ";
            else if (allele.alt == "U") out << "-";
            else if (allele.alt == "R") out << ".";
            else out << allele.alt;
        }
        out << endl;
    }

    size_t i = 0;
    for (auto a : ordered_alignments) {
        auto& aln = *a;
        // print out the stuff
        if (aln.IsReverseStrand())     out << "S"; else out << "s";
        if (aln.IsMateReverseStrand()) out << "O"; else out << "o";
        if (aln.IsDuplicate())         out << "D"; else out << "d";
        if (aln.IsFailedQC())          out << "Q"; else out << "q";
        if (aln.IsFirstMate())         out << "F"; else out << "f";
        if (aln.IsSecondMate())        out << "X"; else out << "x";
        if (aln.IsMateMapped())        out << "Y"; else out << "y";
        if (aln.IsPaired())            out << "P"; else out << "p";
        if (aln.IsPrimaryAlignment())  out << "Z"; else out << "z";
        if (aln.IsProperPair())        out << "I"; else out << "i";
        out << "  ";
        for (auto& allele : alignment_alleles[&aln]) {
            if (allele.alt == "M") out << " ";
            else if (allele.alt == "U") out << "-";
            else if (allele.alt == "R") out << ".";
            else out << allele.alt;
        }
        // now the matches
        out << " ";
        for (auto w : matches[&aln]) {
            out << w.second << " ";
        }
        out << aln.MapQuality;
        out << " " << aln.Name;
        out << endl;
    }
    // now handle caller input features
    for (auto& f : call_info_num) {
        out << f.first << ":" << f.second << " ";
    }
    for (auto& f : call_info_str) {
        out << f.first << ":" << f.second << " ";
    }
    
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
            out << ++idx << allele.alt << ":" << allele.prob << " ";;
        }
    }
    i = 1;
    int gid = 0;
    for (auto& geno : genotypes) {
        out << "|geno" << i << " ";
        ++i;
        idx = 0;
        for (auto& allele : geno) {
            out << ++idx << allele.alt << ":" << allele.prob << " ";;
        }
    }

    i = 0;
    // do the alignments
    for (auto a : ordered_alignments) {
        auto& aln = *a;
        auto name = "aln" + std::to_string(i++);
        out << "|" << name << " ";
        // now alleles
        idx = 0;
        for (auto& allele : alignment_alleles[&aln]) {
            out << ++idx << allele.alt << ":" << allele.prob << " ";
        }
    }

    i = 0;
    // do the strands
    // do the mapping probs
    for (auto a : ordered_alignments) {
        auto& aln = *a;
        out << "|properties" << i++ << " ";
        if (exponentiate) {
            out << "mapqual:" << 1-phred2float(min(aln.MapQuality, (uint16_t)60)) << " ";
        } else {
            out << "mapqual:" << aln.MapQuality << " ";
        }
        // handle flags
        if (aln.IsReverseStrand())     out << "strand:1"; else out << "strand:0"; out << " ";
        if (aln.IsMateReverseStrand()) out << "ostrand:1"; else out << "ostrand:0"; out << " ";
        if (aln.IsDuplicate())         out << "dup:1"; else out << "dup:0"; out << " ";
        if (aln.IsFailedQC())          out << "qcfail:1"; else out << "qcfail:0"; out << " ";
        if (aln.IsFirstMate())         out << "fmate:1"; else out << "fmate:0"; out << " ";
        if (aln.IsSecondMate())        out << "xmate:1"; else out << "xmate:0"; out << " ";
        if (aln.IsMateMapped())        out << "ymap:1"; else out << "ymap:0"; out << " ";
        if (aln.IsPaired())            out << "paired:1"; else out << "paired:0"; out << " ";
        if (aln.IsPrimaryAlignment())  out << "zprimary:1"; else out << "zprimary:0"; out << " ";
        if (aln.IsProperPair())        out << "iproper:1"; else out << "iproper:0"; out << " ";
    }
    
    i = 0;
    for (auto a : ordered_alignments) {
        auto& aln = *a;
        out << "|match" << i++ << " ";
        for (auto w : matches[&aln]) {
            out << w.first+1 << "H:" << w.second << " ";
        }
    }

    out << "|software ";
    // now handle caller input features
    for (auto& f : call_info_num) {
        out << f.first << ":" << f.second << " ";
    }
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
