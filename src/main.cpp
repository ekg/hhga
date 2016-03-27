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
         << "    -w, --window-size N   use a fixed window of this size in the MSA matrix" << endl
         << "    -r, --region REGION   limit variants to those in this region (chr:start-end)" << endl
         << "    -t, --text-viz        make a human-readible, compact output" << endl
         << "    -c, --class-label X   add this label (e.g. -1 for false, 1 for true)" << endl
         << endl
         << "Generates reports on the rate of putative mutations or errors in the input alignment data." << endl
         << "Alignments are read from the specified files, or stdin if none are specified" << endl
         << endl
         << "authors: Erik Garrison <erik.garrison@gmail.com> and Nicolás Della Penna <nikete@gmail.com>" << endl;

}

int main(int argc, char** argv) {

    vector<string> inputFilenames;
    string vcf_file_name;
    string region_string;
    string fastaFile;
    string output_format = "vw";
    string class_label;
    size_t window_size = 50;

    // parse command-line options
    int c;

    while (true) {

        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"bam",  required_argument, 0, 'b'},
            {"vcf", required_argument, 0, 'v'},
            {"region", required_argument, 0, 'r'},
            {"fasta-reference", required_argument, 0, 'f'},
            {"text-viz", no_argument, 0, 't'},
            {"class-label", no_argument, 0, 'c'},
            {"window-size", required_argument, 0, 'w'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hb:r:f:v:tc:w:",
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

        case 'w':
            window_size = atoi(optarg);
            break;

        default:
            return 1;
            break;
        }
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
        //cerr << "Got variant " << var << endl;
        HHGA hhga(window_size,
                  bam_reader,
                  fasta_ref,
                  var,
                  class_label);
        if (output_format == "vw") {
            cout << hhga.vw() << endl;
        } else if (output_format == "text-viz") {
            cout << hhga.str() << endl;
        }
    }

    return 0;

}
