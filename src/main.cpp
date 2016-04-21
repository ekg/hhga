#include "hhga.hpp"

using namespace std;
using namespace hhga;

void printUsage(int argc, char** argv) {

    cerr << "usage: " << argv[0] << " [-b FILE]" << endl
         << endl
         << "options:" << endl
         << "    -h, --help            this dialog" << endl
         << "    -f, --fasta-reference FILE  the reference sequence" << endl
         << "    -b, --bam FILE        use this BAM as input (multiple allowed)" << endl
         << "    -v, --vcf FILE        derive an example from every record in this file" << endl
         << "    -n, --name NAME       apply NAME as the prefix for the annotations in --vcf" << endl
         << "    -w, --window-size N   use a fixed window of this size in the MSA matrix" << endl
         << "    -r, --region REGION   limit variants to those in this region (chr:start-end)" << endl
         << "    -t, --text-viz        make a human-readible, compact output" << endl
         << "    -c, --class-label X   add this label (e.g. -1 for false, 1 for true)" << endl
         << "    -g, --gt-class        use the GT field for the sample to make diploid genotype class labels" << endl
         << "    -e, --exponentiate    convert features that come PHRED-scaled to [0,1]" << endl
         << "    -s, --show-bases      show all the bases in the alignments instead of R ref match symbol" << endl
         << "    -p, --predictions-in  stream in predictions and output an annotated VCF" << endl
         << "    -N, --sample NAME     name of the sample in the output VCF (use with --gt-class)" << endl
         << "    -d, --debug           print useful debugging information to stderr" << endl
         << endl
         << "Generates examples for vw using a VCF file and BAM file." << endl
         << "May optionally convert vw predictions into an annotated VCF file for downstream integration." << endl
         << endl
         << "authors: Erik Garrison <erik.garrison@gmail.com> and Nicolás Della Penna <nikete@gmail.com>" << endl;

}

int main(int argc, char** argv) {

    vector<string> inputFilenames;
    string vcf_file_name;
    string vcf_feature_prefix;
    string region_string;
    string fastaFile;
    string output_format = "vw";
    string class_label;
    size_t window_size = 50;
    bool debug = false;
    bool exponentiate = false;
    bool show_bases = false;
    bool assume_ref = true; // assume haps are ref when not given
    bool predictions_in = false;
    bool gt_class = false;
    string sample_name;

    // parse command-line options
    int c;

    while (true) {

        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"bam",  required_argument, 0, 'b'},
            {"vcf", required_argument, 0, 'v'},
            {"name", required_argument, 0, 'n'},
            {"region", required_argument, 0, 'r'},
            {"fasta-reference", required_argument, 0, 'f'},
            {"text-viz", no_argument, 0, 't'},
            {"class-label", no_argument, 0, 'c'},
            {"gt-class", no_argument, 0, 'g'},
            {"window-size", required_argument, 0, 'w'},
            {"exponentiate", no_argument, 0, 'e'},
            {"show-bases", no_argument, 0, 's'},
            {"predictions-in", no_argument, 0, 'p'},
            {"sample-name", required_argument, 0, 'S'},
            {"debug", no_argument, 0, 'd'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hb:r:f:v:tc:w:dn:espgS:",
                         long_options, &option_index);

        if (c == -1)
            break;
 
        switch (c) {

        case '?':
            printUsage(argc, argv);
            return 0;
            break;

        case 'h':
            printUsage(argc, argv);
            return 0;
            break;

        case 'b':
            inputFilenames.push_back(optarg);
            break;

        case 'v':
            vcf_file_name = optarg;
            break;

        case 'n':
            vcf_feature_prefix = optarg;
            break;

        case 'r':
            region_string = optarg;
            break;

        case 'f':
            fastaFile = optarg;
            break;

        case 't':
            output_format = "text-viz";
            break;

        case 'c':
            class_label = optarg;
            break;

        case 'g':
            gt_class = true;
            break;

        case 'w':
            window_size = atoi(optarg);
            break;

        case 'e':
            exponentiate = true;
            break;

        case 's':
            show_bases = true;
            break;

        case 'p':
            predictions_in = true;
            break;

        case 'S':
            sample_name = optarg;
            break;

        case 'd':
            debug = true;
            break;

        default:
            return 1;
            break;
        }
    }

    if (predictions_in) {

        stringstream headerss;
        if (sample_name.empty()) sample_name = "unknown";
        headerss 
            << "##fileformat=VCFv4.3" << endl
            << "##source=hhga" << endl
            << "##INFO=<ID=prediction,Number=1,Type=Integer,Description=\"hhga+vw prediction for site\">" << endl
            << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl
            << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
        headerss << sample_name;
        vcflib::VariantCallFile vcf_file;
        string header = headerss.str();
        vcf_file.openForOutput(header);
        // add sample
        //vcf_file.sampleNames.push_back(sample_name);
        cout << vcf_file.header << endl;
        
        
        // stream in predictions, use the annotation we apply to the vw input
        // to reconstruct a VCF file with our model's predictions as an INFO field
        for (std::string line; std::getline(std::cin, line); ) {
            auto fields = split_delims(line, " \t");
            // we want the second field
            auto& prediction = fields[0];
            auto comment = fields[1];
            // strip off the leading ' if present
            if (comment.find("'") != string::npos) {
                comment = comment.substr(1);
            }
            auto vcf_fields = split_delims(comment, "_");
            auto& seqname = vcf_fields[0];
            auto pos = stol(vcf_fields[1].c_str());
            auto& ref = vcf_fields[2];
            auto& alt = vcf_fields[3];

            vcflib::Variant var(vcf_file);
            var.sequenceName = seqname;
            var.position = pos;
            var.quality = 0;
            var.ref = ref;
            var.alt.push_back(alt);
            var.id = ".";
            var.filter = ".";
            var.info["prediction"].push_back(convert(prediction));
            if (gt_class) {
                var.samples[sample_name]["GT"].clear();
                var.samples[sample_name]["GT"].push_back(genotype_for_label(prediction));
                var.format.push_back("GT");
            }
            cout << var << endl;
        }
        return 0;
    }
    
    if (fastaFile.empty()) {
        cerr << "no FASTA reference specified" << endl;
        printUsage(argc, argv);
        return 1;
    }

    if (inputFilenames.empty()) {
        cerr << "no input files specified" << endl;
        printUsage(argc, argv);
        return 1;
    }

    BamTools::BamMultiReader bam_reader;
    if (!bam_reader.Open(inputFilenames)) {
        cerr << "could not open input BAM files" << endl;
        return 1;
    }

    vcflib::VariantCallFile vcf_file;
    if (!vcf_file_name.empty()) {
        vcf_file.open(vcf_file_name);
        if (!vcf_file.is_open()) {
            cerr << "could not open " << vcf_file_name << endl;
            return 1;
        }
    }

    FastaReference fasta_ref;
    fasta_ref.open(fastaFile);

    // if we've got a limiting region, use it
    if (!region_string.empty()) {
        set_region(vcf_file, region_string);
    }
    // iterate through all the vcf records, building one hhga matrix for each
    vcflib::Variant var(vcf_file);
    while (vcf_file.getNextVariant(var)) {
        if (debug) { cerr << "Got variant " << var << endl; }
        HHGA hhga(window_size,
                  bam_reader,
                  fasta_ref,
                  var,
                  vcf_feature_prefix,
                  class_label,
                  gt_class,
                  exponentiate,
                  show_bases,
                  assume_ref);
        if (output_format == "vw") {
            cout << hhga.vw() << endl;
        } else if (output_format == "text-viz") {
            cout << hhga.str() << endl;
        }
    }

    return 0;

}
