#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include "command_line_parser.hpp" 

using namespace std;

enum Mode {w,r};
enum Zygosity {homo_ref, homo_var, hetero, undefined};

Zygosity zygosity(char * sample){
    if(sample[0]=='.' || sample[1]=='.') return undefined;
    if(sample[0]=='0' && sample[1]=='0') return homo_ref;
    if(sample[0]=='1' && sample[1]=='1') return homo_var;
    if(sample[0]=='1' && sample[1]=='0') return hetero;
    if(sample[0]=='0' && sample[1]=='1') return hetero;
    else throw std::runtime_error("Genotype "+string(1,sample[0])+":"+string(1,sample[1]) + " unknown.");
}

void checkFile(string path,Mode mode){  

            if(mode==r){
                std::ifstream stream(path);
                if(!stream.good()){
                    throw std::runtime_error(path + " could not be opened!");
                }
                stream.close(); 
            } else{
                std::ofstream stream(path);
                if(!stream.good()){
                    throw std::runtime_error(path + " could not be opened!");
                }
                stream.close();
            }       
}

void parse_line(vector<string>& result, string line, char sep) {
	string token;
	istringstream iss (line);
	while (getline(iss, token, sep)) {
		result.push_back(token);
	}
}

bool is_header(vector<string> tokens){
    if (tokens[0].substr(0,2) == "##") return true;
    if (tokens[0].at(0) == '#') {
        // check number of samples/paths given
        if (tokens.size() < 9) {
            throw runtime_error("VariantReader::VariantReader: not a proper VCF-file.");
        }
        if (tokens.size() < 10) {
            throw runtime_error("VariantReader::VariantReader: no haplotype paths given.");
        }
        // validate header line
        map<unsigned int, string> fields = { {0, "#CHROM"}, {1, "POS"}, {2, "ID"}, {3, "REF"}, {4, "ALT"}, {5, "QUAL"}, {6, "FILTER"}, {7, "INFO"}, {8, "FORMAT"}, {9, "sample"}};
        for (unsigned int i = 0; i < 9; ++i) {
            if (tokens[i] != fields[i]) {
                throw runtime_error("VariantReader::VariantReader: VCF header line is malformed.");
            }
        }
        return true;
    }
    if (tokens.size() < 10) {
        throw runtime_error("VariantReader::VariantReader: malformed VCF-file, or no haplotype paths given in VCF.");
    }
    return false;
}

bool next_result(vector<string>& result, ifstream &stream){
    
    bool header = true;
    string line;
    
    while(header){
        if(getline(stream,line)){
            if (line.size() == 0) continue;

            vector<string> tokens;
		    parse_line(tokens, line, '\t');
            if(!is_header(tokens)){
                header = false;
                result = tokens;
            }

        } else { return false;   }   
    }
    return true;
}

//Case insensitive string comparison
bool iequals(const string& a, const string& b)
{
    unsigned int sz = a.size();
    if (b.size() != sz)
        return false;
    for (unsigned int i = 0; i < sz; ++i)
        if (tolower(a[i]) != tolower(b[i]))
            return false;
    return true;
}

int main (int argc, char* argv[])
{
    //Initialise command line parser
    CommandLineParser parser;
    parser.addAuthorInformation("GenotypeConcordance,  author: Hartmut Haentze \nCompare two vcf files of the same individual with same structure and analyze differences.");
    
    parser.addRequiredArgument('c',"Call set vcf file");
    parser.addRequiredArgument('t',"True vcf file");
    //parser.addOptionalArgument('o',"","Output file");
    parser.addFlagArgument('O',"Print output in  machine readable format to console");

    //Parse command line
    parser.parse(argc,argv);
    string callset_vcf = parser.getArgument('c');
    string true_vcf = parser.getArgument('t');
    //string output_file = parser.getArgument('o');
    bool printOutput = parser.getFlag('O');

    //test files
    checkFile(callset_vcf,r);
    checkFile(true_vcf,r);
    /*if(output_file!=""){
        checkFile(output_file,w);
    }*/
    
    int true_homo_var = 0;
    int true_homo_ref = 0;
    int true_hetero = 0;
    int tp_homo_var= 0;
    int tp_homo_ref = 0;
    int tp_hetero = 0;
    int fp_homo_var= 0;
    int fp_homo_ref = 0;
    int fp_hetero = 0;
    int fn_homo_var= 0;
    int fn_hetero = 0;
    int true_undefined = 0;
    int sample_undefined = 0;

    //predicted hetero, but true homo
    int lenient_predicted_hetero = 0;


    /** Code partially taken from: https://github.com/eblerjana/pangenie/blob/master/src/variantreader.cpp **/    
    ifstream callset_stream(callset_vcf);
    ifstream true_stream(true_vcf);

    while(true){
        vector<string> callset_tokens;
        vector<string> true_tokens;

        //Get next line and break if EOF
        if(!next_result(callset_tokens,callset_stream)) break;
        if(!next_result(true_tokens,true_stream)) break;

        // compare chromosomes
        if(callset_tokens[0] != true_tokens [0]){
            string s("");
            string t("");
            for(size_t i = 0; i < callset_tokens.size(); i++){
                s += callset_tokens[i];
                s += " ";
                t += true_tokens[i];
                t += " ";
            }
            cerr<< "Sample: " << s << endl;
            cerr<< "True: " << t << endl;
            throw runtime_error("true and sample vcf files have different structures!");
        }
        //compare positions
        if(callset_tokens[1] != true_tokens [1]){
            string s("");
            string t("");
            for(size_t i = 0; i < callset_tokens.size(); i++){
                s += callset_tokens[i];
                s += " ";
                t += true_tokens[i];
                t += " ";
            }
            cerr << "Sample: " << s << endl;
            cerr << "True: " << t << endl;
            throw runtime_error("true and sample vcf files have different structures!");
        }
        //compare reference
        if(!iequals(callset_tokens[3],true_tokens[3])){
            throw runtime_error("reference of true and sample vcf files different!");
        }
        //compare alternative allele
        if(!iequals(callset_tokens[4],true_tokens[4])){
            throw runtime_error("alternative allele of true and sample vcf files different!");
        }

        //get alleles
        vector<string> sample;
        parse_line(sample,callset_tokens[9],':');
        char callset_alleles[2] = {sample[0][0],sample[0][2]};
        char true_alleles[2] = {true_tokens[9][0],true_tokens[9][2]};

        
        if(zygosity(true_alleles) == undefined){
            true_undefined++;
            continue;
        }

        //homo-ref does not count into fn
        if(zygosity(true_alleles) == homo_ref){
            true_homo_ref++;
            if(zygosity(callset_alleles) == undefined) sample_undefined++;
            else if(zygosity(callset_alleles) == homo_ref) tp_homo_ref++;
            else if (zygosity(callset_alleles) == homo_var) fp_homo_var++;
            else fp_hetero++;
        }

        else if(zygosity(true_alleles) == homo_var){
            true_homo_var++;
            if(zygosity(callset_alleles) == undefined){
                sample_undefined++;
                fn_homo_var++;
            } 
            else if(zygosity(callset_alleles) == homo_ref) {
                fp_homo_ref++;
                fn_homo_var++;
            }
            else if (zygosity(callset_alleles) == homo_var) tp_homo_var++;
            else {
                fp_hetero++;
                fn_homo_var++;
            }
        }
        else if(zygosity(true_alleles) == hetero){
            true_hetero++;
            if(zygosity(callset_alleles) == undefined){
                sample_undefined++;
                fn_hetero++;
            } 
            else if(zygosity(callset_alleles) == homo_ref) {
                fp_homo_ref++;
                fn_hetero++;
            }
            else if (zygosity(callset_alleles) == homo_var) {
                fp_homo_var++;
                fn_hetero++;
                lenient_predicted_hetero++;
            }
            else tp_hetero++;
        }
    } 

    cerr<< "\nHomo-Var: " << true_homo_var << "\tCorrectly Predicted: " << tp_homo_var<< "\tWrongly Predicted: " << fp_homo_var<< endl;
    cerr<< "Homo-Ref: " << true_homo_ref << "\tCorrectly Predicted: " << tp_homo_ref << "\tWrongly Predicted: " << fp_homo_ref << endl;
    cerr<< "Hetero: " << true_hetero << "\tCorrectly Predicted: " << tp_hetero << "\tWrongly Predicted: " << fp_hetero << endl;
    cerr<< "Lenient Hetero: " << lenient_predicted_hetero << endl;
    cerr<< "Undefined in True vcf: " << true_undefined << " (should be 0)\tUndefined in sample vcf: " << sample_undefined << endl;

    //Quality control
    int true_variants = true_homo_ref + true_homo_var + true_hetero;
    int correct_sample_variants = tp_homo_ref + tp_homo_var+ tp_hetero;
    int wrong_sample_variants = fp_homo_ref + fp_homo_var+ fp_hetero;
    int sample_variants = correct_sample_variants + wrong_sample_variants + sample_undefined;
    cerr<< "True variants: " << true_variants <<  "\tSample variants: " << sample_variants << "\t(Should be equal)" << endl;

    //Sensitivity
    float sens_homo_var = (float)tp_homo_var/ (float) true_homo_var;
    float sens_homo_ref = (float)tp_homo_ref/ (float) true_homo_ref;
    float sens_hetero = (float)tp_hetero / (float)true_hetero;
    float sens_lenient_hetero = (float) (tp_hetero + lenient_predicted_hetero) / (float)true_hetero;
    float sens_lenient = (float) (tp_hetero + lenient_predicted_hetero + tp_homo_var) / (float)(true_hetero + true_homo_var);
    cerr<< "\nSensitivity homo_var: \t" << sens_homo_var << endl;
    cerr<< "Sensitivity hetero (strict): \t" << sens_hetero << endl;
    cerr<< "Sensitivity hetero (lenient): \t" << sens_lenient_hetero << endl;
    cerr<< "Sensitivity variant (lenient): \t" << sens_lenient << endl;     

    //Genotype Concordance
    float concordance = (float)correct_sample_variants / (float)sample_variants;
    cerr<< "\nGenotype concordance: \t" << concordance << endl;

    //Weighted Genotype Concordance (as specified in PanGenie supplements)
    float weighted_concordance = (sens_homo_var + sens_hetero + sens_homo_ref) / 3;
    cerr<< "Weighted genotype concordance: \t" << weighted_concordance << endl;

    //Meassures (as specified in PanGenie supplements)
    float tp = tp_homo_var + tp_hetero;
    float fn = fn_homo_var + fn_hetero;
    float fp = fp_homo_var + fp_hetero;

    float prec = tp / (tp+fp);
    float recall = tp / (tp+fn);
    float f_score = 2 * (prec*recall) / (prec+recall);

    cerr<< "\nPrecision: \t" << prec << endl;
    cerr<< "Recall: \t" << recall << endl;
    cerr<< "F-score: \t" << f_score << endl;

    //Print machine readable format
    if(printOutput){
        vector<pair<string,float>> scores{make_pair("GC",concordance),
                                make_pair("wGC", weighted_concordance),
                                make_pair("Precision",prec),
                                make_pair("Recall",recall),
                                make_pair("F-score", f_score)};
        string line1("");
        string line2("");
        for(auto s : scores){
            line1+= s.first + "\t";
            line2+= to_string(s.second) + "\t";
        }
        cout << line1 << endl;
        cout << line2 << endl;
    }

    return 0;
}