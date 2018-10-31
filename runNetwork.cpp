#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <sstream>
#include <cstring>

#include <math.h>

#include <algorithm>

#include <vector>
#include <map>

#include <string>

#include <getopt.h>

#include <dirent.h>

#include <omp.h>

using namespace std;

static string VERSION = "1.8";

static int DEBUG = 0;

//Integer to Identifier mapping
vector<string> int_id;
map<string, unsigned int> id_int;

//pValues index synced to identifier mapping
vector<double> int_pValue;

//keeps record of which p-Values have been accessed
vector<int> visited;

//Holds bin with highest index
unsigned int maxSize = 0;

//Holds max index of proteins in analysis
unsigned int proteinsInAnalysis;
//Holds number of proteins in network ==> NMB will be run for
unsigned int proteinsUnderConsideration;

//Bins
vector<vector<unsigned> > bin_prot;
vector<vector<unsigned> > prot_bin;
vector<unsigned> prot_con;

//Protein protein network
vector<vector<unsigned> > prot_prot;

//Random networks
vector<vector<vector<unsigned> > > randomNetworks;

//Keep genes to run NMB for
map<string, int> gene_compute;

void usage(){
  cout << "runNetwork, Version: " << VERSION << endl << endl;
  cout << "Usage:" << endl <<
  "runNetwork " << "[arguments]" << endl << endl <<
  "Required arguments:" << endl << 
  "   -i | --iterations           number of iterations" << endl <<
  "   -j | --job                  jobnumber for partitioning input set (>0)" << endl <<
  "   -J | --jobs                 total number of jobs" << endl <<
  "   -n | --network              tab-delimited file for interaction network, has to be symetric" << endl <<
  "   -N | --name                 jobname, defines output filename" << endl <<
  "   -p | --pvalues              file to read p-values, 1st column is identifier" << endl <<
  "Optional arguments:" << endl <<
  "   -b | --binaryNetworkFile    flag for indicating binary network ('short' fields). Requires mapping file" << endl <<
  "   -c | --correctInflation     sets correction factor for genomic inflation" << endl <<
  "   -C | --calculateOnlyPhi     set to only calculate Phi and exit" << endl <<
  "   -d | --printDistributions   print distributions of phi-values" << endl <<
  "   -D | --networkDirectory     directory for getting the random networks" << endl <<
  "   -e | --estimatePerformance  just prints out runtime for NMB calculation, exluding loading/initialization" << endl <<
  "   -f | --flipped              inverted analysis, look for genes with non-mutated neighborhood" << endl <<
  "   -g | --ignoreInputNodes     set to include interactors when resampling" << endl <<
  "   -G | --ignoreCentralNode    set to include query node when resampling" << endl <<
  "   -h | --help                 print this message" << endl <<
  "   -l | --list                 file with list of genes for which to compute NMB" << endl <<
  "   -m | --minimumPvalues       sets the minimum neighbours with p-value < 1 for a node to be considered" << endl <<
  "   -M | --mappingFile          mapping file for p-Value IDs -> network IDs" << endl <<
  "   -o | --outputDirectory      sets output directory" << endl <<
  "   -P | --pValColumn           column index to take the p-value (counting starts at 0 ;) ) [default=3]" << endl <<
  "   -r | --ignoreConnectivity   set to ignore connectivity when permuting" << endl <<
  "   -R | --randomize            read in the p-values, return randomized and quit" << endl <<
  "   -s | --separator            set separator for network file ('\t' is default, ',' if a .csv file is provided)" << endl <<
  "   -t | --threads              set number of parallel processes run (default = 1), needs to be compiled using -fopenmp" << endl <<
  "   -v | --version              print version number" << endl <<
  "   -V | --verbose              extensive output" << endl <<
  "   -I | --included             list genes with significant p-Values not found in the network (and their p-Values)" << endl <<
  cout << endl << "Note: Genes with no assigned p-Value will be removed from the network" << endl << endl;
}

int myrandom (int i) {return rand()%i; }

long double calculateChiSquare(vector<unsigned>* nodes){
  double score = 0;

  for(int i=0;i<nodes->size();i++){
    score += log(int_pValue[(*nodes)[i]]);
    if(DEBUG) cerr << (*nodes)[i] << "\t" << int_id[(*nodes)[i]] << "\t" << int_pValue[(*nodes)[i]] << endl;
  }

  if(score != 0){score *= -2;}

  return score;
}
long double calculateChiSquareFromValues(vector<double>* values){
  double score = 0;

  for(int i=0;i<values->size();i++){
    score += log((*values)[i]);
  }

  if(score != 0){score *= -2;}

  return score;
}

unsigned int getNodesFromBin(unsigned int binIndex, unsigned int nodeCounter, vector<unsigned int>* newNodes, vector<int>* selectedNodes, unsigned int* binsProtsSelected){
  if(bin_prot[binIndex].size() == 0){return 0;}

  unsigned int selectedNodesCounter = 0;
  unsigned int selectedBinMax = bin_prot[binIndex].size();

  unsigned int minSelectForPermute = 10;

  if(binsProtsSelected[binIndex] < 1){return 0;}

  //TODO: CHECK!!!!!
  //if(nodeCounter<minSelectForPermute){
    unsigned int randomIndex = rand()%selectedBinMax;
    while((*selectedNodes)[bin_prot[binIndex][randomIndex]]){
      randomIndex = rand()%selectedBinMax;
    }
    unsigned int prot = bin_prot[binIndex][randomIndex];
    newNodes->push_back(prot);
    (*selectedNodes)[prot] = true;
    selectedNodesCounter++;
    //save bin status
    for(unsigned int k=0;k<prot_bin[prot].size();k++){binsProtsSelected[prot_bin[prot][k]]--;}
  //}else{
  //  random_shuffle(bin_prot[binIndex].begin(), bin_prot[binIndex].end(), myrandom);
  //
  //  for(unsigned int j=0;j<nodeCounter && j<selectedBinMax;j++){
  //    unsigned int prot = bin_prot[binIndex][j];
  //    if(!(*selectedNodes)[prot]){
  //      newNodes->push_back(prot);
  //      (*selectedNodes)[prot] = true;
  //      selectedNodesCounter++;
  //      //save bin status
  //      for(unsigned int k=0;k<prot_bin[prot].size();k++){binsProtsSelected[prot_bin[prot][k]]--;}
  //    }
  //  }
  //}

  return selectedNodesCounter;
}

//void getRandomNodes(vector<unsigned int>* nodes, vector<unsigned int>* newNodes, bool ignoreInputNodes){
void getRandomNodes(int index, vector<unsigned int>* newNodes, bool ignoreInputNodes, bool ignoreCentralNode){
  vector<int> selectedNodes(proteinsInAnalysis);

  if(DEBUG) cerr << "Index: " << index << endl;

  unsigned int binsProtsSelected[maxSize];
  unsigned int binsToSelectFrom[maxSize];

  for(int i=0;i<maxSize;i++){
    binsToSelectFrom[i] = 0;
    binsProtsSelected[i] = bin_prot[i].size();
  }
  for(int i=0;i<proteinsInAnalysis;i++){selectedNodes[i] = false;}

  if(ignoreInputNodes && !(prot_prot[index].size() > proteinsUnderConsideration/2.0)){
    for(unsigned int j=0;j<prot_prot[index].size();j++){
      unsigned int prot = prot_prot[index][j];
      selectedNodes[prot] = true;
      for(unsigned int k=0;k<prot_bin[prot].size();k++){
        binsProtsSelected[prot_bin[prot][k]]--;
      }
    }
  }

  if(ignoreCentralNode){
    selectedNodes[index] = true;
    for(unsigned int k=0;k<prot_bin[index].size();k++){
      binsProtsSelected[prot_bin[index][k]]--;
    }
  }

  //create vector of bins to select from
  for(unsigned int i=0;i<prot_prot[index].size();i++){
    unsigned int size = prot_con[prot_prot[index][i]];
    binsToSelectFrom[size]++;

    if(DEBUG) cerr << "BinToSelectFrom: " << size << "\t->\t" << binsToSelectFrom[size] << endl;
  }

  //select proteins from corresponding bins
  for(unsigned int i=0;i<maxSize;i++){
    if(binsToSelectFrom[i] > 0){
      unsigned int selectedNodesCounter = 0;
      unsigned int binCounter = 1;
      while(selectedNodesCounter < binsToSelectFrom[i]){
        int dir = binCounter % 2;
        if(dir != 1){dir = binCounter/2;}else{dir = -1*((binCounter-1)/2);}
        binCounter++;
        int newBin = i+dir;
        if(newBin >= maxSize){continue;}
        if(DEBUG) cerr << "GetNodeFrom: " << newBin << endl;
        selectedNodesCounter += getNodesFromBin(newBin, binsToSelectFrom[i]-selectedNodesCounter, newNodes, &selectedNodes, binsProtsSelected);
      }
    }
    
  }
}

vector<string> getFilesFromDirectory(string path = ".") {
  DIR*    dir;
  dirent* pdir;
  vector<string> files;

  dir = opendir(path.c_str());
  while ( (pdir = readdir(dir)) ) {
    //skip hidden files and special directories
    if(pdir->d_name[0] == '.'){continue;}
    //take only the _InWeb3_HUGO_sym files
    files.push_back(pdir->d_name);
  }

  return files;
}

int loadNetwork(string networkFile, vector<vector<unsigned> >* prot_prot_, char separator_='\t', int verbose=0, int mappingFile=0, int binaryNetworkFile=0, int symetric=1){
  visited.resize(25000);
  
  string line_;
  string type = networkFile.substr(networkFile.size()-4, networkFile.size()-1);

  ifstream IN_;
  // From now on by definition: mapping file = binary
  if(binaryNetworkFile){
    if (verbose) cout << "Loading PPI data using mapping file [binary]..." << endl;
    IN_.open(networkFile.c_str(), ios::binary);
  }else{
    if (verbose) cout << "Loading PPI data..." << endl;
    IN_.open(networkFile.c_str());
  }

  if(!IN_){
    cerr << "Network file not found: " << endl << networkFile.c_str() << endl;
    return 0;
  }

  // From now on by definition: mapping file = binary
  bool csv = false;
  if(strncmp(type.c_str(), ".csv", 4) == 0) csv = true;
  
  if(binaryNetworkFile){

    IN_.seekg (0, IN_.end);
    long length = (long)IN_.tellg();
    IN_.seekg (0, IN_.beg);

    if (verbose) cerr << "Filesize: " << length << endl;

    short int1, int2;

    while((long)IN_.tellg() != length){
      IN_.read((char*)&int1, 2);
      IN_.read((char*)&int2, 2);

      if(int_pValue[int1] == 0){continue;}
      if(int_pValue[int2] == 0){continue;}

      //Remove self interactions
      if (int1 == int2) {continue;}

      if(int1 >= prot_prot_->size()){
        prot_prot_->resize(int1+1);
      }
      (*prot_prot_)[int1].push_back(int2);

    }

    if (verbose) cerr << "Done reading" << endl;
  }else{
    while(getline(IN_, line_)){
      vector<string> tokens;
      if(line_[0] == '#'){continue;}
      istringstream iss(line_);
      string token;
      if(csv){ 
         while(getline(iss, token, ',')){
         if(token[0] == '"') token = token.substr(1, token.size()-2);
         if(strncmp(token.c_str(), "gene1", 5) == 0){continue;} 
         tokens.push_back(token);
         }
      }else{
        while(getline(iss, token, separator_)){
          tokens.push_back(token);
        }
      }

      short int1, int2;

      if(!mappingFile){
        if(id_int.find(tokens[0]) == id_int.end()){
          int_pValue.push_back(1);
          int_id.push_back(tokens[0]);
          id_int.insert(pair<string, unsigned int>(tokens[0], int_id.size()-1));
        }
        if(id_int.find(tokens[1]) == id_int.end()){
          int_pValue.push_back(1);
          int_id.push_back(tokens[1]);
          id_int.insert(pair<string, unsigned int>(tokens[1], int_id.size()-1));
        }
        int1 = id_int[tokens[0]];
        int2 = id_int[tokens[1]];
      }else{
        int1 = atoi(tokens[0].c_str());
        int2 = atoi(tokens[1].c_str());
      }

      //Remove self interactions
      if (int1 == int2) {continue;}
      
      if(int1 >= visited.size()){
         visited.resize(int1+1);
      }
      if(int2 >= visited.size()){
         visited.resize(int2+1);
      }
      //mark genes as visited
      visited[int1] = 1;
      visited[int2] = 1;
      
      if(int1 >= prot_prot_->size()){ 
        prot_prot_->resize(int1+1);
      }
      (*prot_prot_)[int1].push_back(int2);
    }
  }
  IN_.close();

  //if not symetric
  //just fill up networks
  unsigned i = 0;
  unsigned sym = 1;
  while((*prot_prot_)[i].size()<1){ i++; }
  for(unsigned j=0;j<(*prot_prot_)[i].size();j++){
    if(find((*prot_prot_)[(*prot_prot_)[i][j]].begin(), (*prot_prot_)[(*prot_prot_)[i][j]].end(), i)==(*prot_prot_)[(*prot_prot_)[i][j]].end()){ sym=0; }
  }
  if(!sym){
    // foreach prot_prot_, search and fill up if not present
    for(unsigned i=0;i<(*prot_prot_).size();i++){
      for(unsigned j=0;j<(*prot_prot_)[i].size();j++){
        unsigned newId = (*prot_prot_)[i][j];
        if(newId >= prot_prot_->size()){
          prot_prot_->resize(newId+1);
          (*prot_prot_)[newId].push_back(i);
        }else if(find((*prot_prot_)[newId].begin(), (*prot_prot_)[newId].end(), i)==(*prot_prot_)[newId].end()){
          (*prot_prot_)[newId].push_back(i);
        }
      }
    }
  }

  //Check if we have a pValue for the key
  //TODO: Set pValue to 1
  for(unsigned i=0;i<(*prot_prot_).size();i++){
  
  }

  return 1;
}

int loadMappingFile(string mapping, vector<string>* int_id_, map<string, unsigned int>* id_int_){
  ifstream in_mapping(mapping.c_str());
  if(!in_mapping){
    cerr << "Mapping file not found: " << endl << mapping.c_str() << endl;
    return 0;
  }

  unsigned i = 0;
  string line;
  while(getline(in_mapping, line)){
    vector<string> tokens;
    if(line[0] == '#'){continue;}
    istringstream iss(line);
    string token;
    while(getline(iss, token, '\t')){
      tokens.push_back(token);
    }
    //token[0]: pVal ID
    //token[1]: InWeb IDs
    unsigned int newInt = (unsigned int)atoi(tokens[1].c_str());
    if(newInt >= int_id_->size()){
      int_id_->resize(newInt+1);
    }
    (*int_id_)[newInt] = tokens[0];
    id_int_->insert(pair<string, unsigned int>(tokens[0], newInt));
  }
  return 1;
}

long double calculatePfromPhis(double chi_, long unsigned int iterations_, map<long double, unsigned int>* phis_, int highChiIsGood_){
  long posCounter_ = 0;
  if(highChiIsGood_ == 1){
    for (map<long double, unsigned int>::iterator itr=(*phis_).begin(); itr!=(*phis_).end(); itr++){
      if(chi_>itr->first){posCounter_+=itr->second;}
    }
  }else{
    for (map<long double, unsigned int>::iterator itr=(*phis_).begin(); itr!=(*phis_).end(); itr++){
      if(chi_<itr->first){posCounter_+=itr->second;}
    }
  }
  long double p_ = 1-posCounter_/(iterations_+1.0);
  return p_;
}

int main (int argc, char **argv) { 
  srand ( time(NULL) ); //initialize the random seed

  clock_t t1,t2;
  t1=clock();

  int ignoreInputNodes = 1;
  int ignoreCentralNode = 1;

  int highChiIsGood = 1;

  unsigned int job = 1;
  unsigned int jobs = 1;

  unsigned int pValColumn = 3;

  long unsigned int iterations = 10000;

  bool verbose = 0;
  bool included = 0;
  bool version = 0;
  bool help = 0;
  bool flipped = 0;
  bool randomize = 0;
  bool ignoreConnectivity = 0;
  int minimumPvalues = 1;
  bool loadRandomNetworks = 0;
  bool onlyChi = 0;
  bool mappingFile = 0;
  bool estimatePerformance = 0;
  bool printPhis = 0;
  bool limit = 0;
  bool binaryNetworkFile = false;
  int nthreads = 1;
  #if defined(_OPENMP)
  nthreads = omp_get_thread_num();
  #endif

  float lambda = 1;

  string pValueData, ppiData, jobName, limitList;
  string OUTDIR = "./";
  string NETWORKDIR = ".";
  string mapping = "";
  char separator = '\t';

  const struct option longopts[] =
  {
    {"binaryNetworkFile",  required_argument, 0,  'b'},
    {"correctInflation",   required_argument, 0,  'c'},
    {"onlyChi",            no_argument,       0,  'C'},
    {"printDistributions", no_argument,       0,  'd'},
    {"networkDirectory",   required_argument, 0,  'D'},
    {"estimatePerformance",no_argument,       0,  'e'},
    {"flipped",            no_argument,       0,  'f'},
    {"ignoreInputNodes",   no_argument,       0,  'g'},
    {"ignoreCentralNode",  no_argument,       0,  'G'},
    {"help",               no_argument,       0,  'h'},
    {"iterations",         required_argument, 0,  'i'},
    {"included",           no_argument,       0,  'I'},
    {"job",                required_argument, 0,  'j'},
    {"jobs",               required_argument, 0,  'J'},
    {"list",               required_argument, 0,  'l'},
    {"minimumPvalues",     required_argument, 0,  'm'},
    {"mappingFile",        required_argument, 0,  'M'},
    {"network",            required_argument, 0,  'n'},
    {"name",               required_argument, 0,  'N'},
    {"outputDirectory",    required_argument, 0,  'o'},
    {"pvalues",            required_argument, 0,  'p'},
    {"pValColumn",         required_argument, 0,  'P'},
    {"ignoreConnectivity", no_argument,       0,  'r'},
    {"randomize",          no_argument,       0,  'R'},
    {"separator",          required_argument, 0,  's'},
    {"nthreads",           required_argument, 0,  't'},
    {"version",            no_argument,       0,  'v'},
    {"verbose",            no_argument,       0,  'V'},
    {0,0,0,0}
  };

  opterr=1;

  int c;
  //while ((c = getopt_long(argc, argv, "vhVfRgGrCdeIi:j:J:p:P:n:N:o:m:c:D:M:", longopts, NULL)) != -1) {
  while ((c = getopt_long(argc, argv, "bc:CdD:efgGhi:I:j:J:l:m:M:n:N:o:p:P:rRs:t:vV", longopts, NULL)) != -1) {
    switch (c) {
      case 'b':
        binaryNetworkFile = true;
        break;
      case 'c':
        lambda = atof(optarg);
        break;
      case 'C':
        onlyChi = 1;
        break;
      case 'd':
        printPhis = 1;
        break;
      case 'D':
        NETWORKDIR = optarg;
        loadRandomNetworks = 1;
        break;
      case 'e':
        estimatePerformance = 1;
        break;
      case 'f':
        flipped = 1;
        highChiIsGood = (flipped==1) ? 0 : 1;
        break;
      case 'g':
        ignoreInputNodes = 0;
        break;
      case 'G':
        ignoreCentralNode = 0;
        break;
      case 'h':
        help = 1;
        usage();
        return 0;
      case 'i':
        iterations = atoi(optarg);
        break;
      case 'I':
         included = 1;
         break;
      case 'j':
        job = atoi(optarg);
        break;
      case 'J':
        jobs = atoi(optarg);
        break;
      case 'l':
        limit = 1;
        limitList = optarg;
        break;
      case 'm':
        minimumPvalues = atoi(optarg);
        break;
      case 'M':
        mapping = optarg;
        mappingFile = 1;
        break;
      case 'n':
        ppiData = optarg;
        break;
      case 'N':
        jobName = optarg;
        break;
      case 'o':
        OUTDIR = optarg;
        break;
      case 'p':
        pValueData = optarg;
        break;
      case 'P':
        pValColumn = atoi(optarg);
        break;
      case 'r':
        ignoreConnectivity = 1;
        break;
      case 'R':
        randomize = 1;
        break;
      case 's':
        if(optarg[0]=='\\' && optarg[1]=='t') {
          separator = '\t';
          break;
        }else{
          separator = *optarg;
          break;
        }
      case 't':
        nthreads = atoi(optarg);
        #if !defined(_OPENMP)
        if(nthreads > 1) cerr << "Warning: More than one thread requested, but compiled without openMP" << endl;
        nthreads = 1;
        #endif
        break;
      case 'v':
        cerr << "runNetwork, Version: " << VERSION << endl;
        return 0;;
      case 'V':
        verbose = 1;
        cerr << "Verbose mode, printing a lot..." << endl;
        break;
      case ':':
        break;
      case '?':
        break;
      default:
        abort ();
    }
  }


  if(!job){
    cerr << "ERROR: 'job' parameter missing" << endl << endl;
    usage();
    return 0;
  }

  if(!jobs){
    cerr << "ERROR: 'jobs' parameter missing" << endl << endl;
    usage();
    return 0;
  }
  job--;

  if(pValueData.size() < 1){
    cerr << "ERROR: no p-Value data available" << endl;
    usage();
    return 0;
  }

  if(ppiData.size() < 1){
    cerr << "ERROR: no network data available" << endl;
    usage();
    return 0;
  }

  if(jobName.size() < 1){
    cerr << "ERROR: 'jobName' parameter missing" << endl;
    usage();
    return 0;
  }

  //load mapping file
  // int_id
  // id_int
  if(mappingFile){
    if (verbose) cerr << "Loading mapping data..." << endl;
    if(!loadMappingFile(mapping, &int_id, &id_int)){
      return 0;
    }
    int_pValue.resize(int_id.size());
    for(int i=0;i<int_pValue.size();i++){int_pValue[i]=1;};
  }

  //Load list of genes for compute (skip computation for genes NOT on this list)
  if(limit){
    if (verbose) cerr << "Loading genes to compute for..." << endl;
    ifstream in_limit(limitList.c_str());
    if(!in_limit){
      cerr << "Limit file not found: " << endl << limitList.c_str() << endl;
      return 0;
    }
    string line;
    while(getline(in_limit, line)){
      gene_compute[line] = 1;
    }
  }

  //Load p-Values
  ifstream in(pValueData.c_str());
  if(!in){
    cerr << "ERROR: p-Value file not found:" << endl << pValueData.c_str() << endl;
    return 0;
  }

  //creates id_int mapping
  if (verbose) cerr << "Loading p-Values data..." << endl;
  unsigned i = 0;
  string line;
  while(getline(in, line)){
    vector<string> tokens;
    if(line[0] == '#'){continue;}
    istringstream iss(line);
    string token;
    while(getline(iss, token, '\t')){
      tokens.push_back(token);
    }
    if (pValColumn >= tokens.size()) { 
      cerr << "ERROR: invalid column reference entered" << endl; 
      return 0;
      }
    //if token casts inaccurately to 0 (meaning an invalid string input was found), skip
    if (((double)atof(tokens[pValColumn].c_str()) == 0) && (tokens[pValColumn].compare("0") != 0)){continue;}
    if(mappingFile){
      if(id_int.find(tokens[0]) == id_int.end()) continue;
      int_pValue[id_int[tokens[0]]] = (double)atof(tokens[pValColumn].c_str());
      //if(tokens[0] == "ABCF1") cerr <<  id_int[tokens[0]] << "\t" << int_pValue[id_int[tokens[0]]] << endl;
    }else{
      int_pValue.push_back((double)atof(tokens[pValColumn].c_str()));
      int_id.push_back(tokens[0]);
      id_int.insert(pair<string, unsigned int>(tokens[0], i++)); // i gets incremented here!!!
    }
  }
  in.close();

  if (verbose) {
    unsigned int count = 0;
    for(unsigned i=0;i<=int_pValue.size();i++){
      if(int_pValue[i] > 0) count++;
    }
    cout << count << " genes with p-values" << endl;
  }
   
  if(randomize == 1){
    random_shuffle(int_pValue.begin(), int_pValue.end(), myrandom);
    for(int i=0;i<int_pValue.size();i++){
      cout << int_id[i] << "\tx\tx\t" << int_pValue[i] << endl;
    }
    return 1;
  }

  if(int_pValue.size() < 1){
    cerr << "ERROR: no p-values loaded" << endl;
    return 0;
  }

  //Load network as int -> array of interactors
  if(!loadNetwork(ppiData, &prot_prot, separator, verbose, mappingFile, binaryNetworkFile)){
    return 0;
  }

  proteinsInAnalysis = prot_prot.size();

  if (verbose) cerr << "Proteins loaded: " << proteinsInAnalysis << endl;

  prot_con.reserve(prot_prot.size());
  for (int i=0; i<prot_prot.size(); i++) {
    prot_con.push_back(prot_prot[i].size());
    if(prot_prot[i].size() > 0) proteinsUnderConsideration++;
  }

  if (verbose) cerr << "Proteins in analysis: " << proteinsUnderConsideration << endl;
  if(proteinsUnderConsideration < 1){
    cerr << "ERROR: No proteins from network would be considered" << endl;
    cerr << "No protein from the input has at least one neighbour with a p-Value?" << endl;
    return 0;
  }

  if(!onlyChi){
    // Put proteins into bins
    for (int i=0; i<prot_prot.size(); i++) {
      unsigned size = prot_con[i];
      if(size > maxSize){
        maxSize = size;
        if(DEBUG) cerr << "MaxSize: " << maxSize << endl;
      }
    }
    maxSize++;
    if(DEBUG) cerr << "Done estimating protein degrees" << endl;

    bin_prot.resize(maxSize+1);
    prot_bin.resize(prot_prot.size()+1);

    //Fill initial bins
    for (int i=0; i<prot_prot.size(); i++) {
      unsigned size = prot_prot[i].size();
      if(size == 0){continue;}
      bin_prot[size].push_back(i);
      prot_bin[i].push_back(size);
    }
    if(DEBUG) cerr << "Done filling bins" << endl;

    vector<vector<unsigned> > ORG_bin_prot(bin_prot);

    //Now fill up to at least 20
    for (int i=0; i<bin_prot.size(); i++) {
      unsigned size = bin_prot[i].size();
      if(size == 0 || size >= 20){continue;}

      unsigned count = 2;
      while(bin_prot[i].size() < 20){
        int dir = count % 2;
        if(dir != 1){dir = count/2;}else{dir = -1*((count-1)/2);}
        int newBin = i+dir;      
        if (newBin <= ORG_bin_prot.size()) {
          if(newBin < maxSize && newBin > 0 && ORG_bin_prot[newBin].size() > 0){ //data in bin
            for(unsigned j=0;j<ORG_bin_prot[newBin].size();j++){ //foreach prot j in new bin
              bool proteinInBin = false;
              for(unsigned k=0;k<prot_bin[bin_prot[newBin][j]].size();k++){ //foreach bin k for protein
                if(prot_bin[bin_prot[newBin][j]][k] == i){proteinInBin = true;}
              }
              if(!proteinInBin){
                bin_prot[i].push_back(bin_prot[newBin][j]);
                prot_bin[bin_prot[newBin][j]].push_back(i);
              }
            }
          }
        }
        count++;
      }
    }
    if(DEBUG) cerr << "Done completing bins" << endl;
  }

  //if a directory is specified, load random networks
  if(loadRandomNetworks){
    if(verbose) cerr << "Loading random networks from '" << NETWORKDIR << "'" << endl;
    vector<string> files;
    files = getFilesFromDirectory(NETWORKDIR);
    randomNetworks.reserve(files.size());
    for(unsigned i=0;i<files.size();i++){
      ostringstream file;
      file << NETWORKDIR << "/" << files[i];
      loadNetwork(file.str(), &randomNetworks[i], verbose, mappingFile);
      if(verbose) cerr << i << "/" << files.size() << NETWORKDIR << files[i] << endl;
    }
  }

  t2=clock();
  float diff (((float)t2-(float)t1)/CLOCKS_PER_SEC);
  if(verbose) cerr << "Everything initialized, starting calculations... " << endl;
  if(verbose) cerr << "Initialization took " << diff << "s" << endl;

  /////////////////////////////////////////////////////////////////////////////
  //
  // Calculations start from here
  //
  /////////////////////////////////////////////////////////////////////////////


  // Initialize timer for computation
  t1=clock();

  //get list of indexes to compute
  int protCounter = 1;
  vector<unsigned int> indexes;
  indexes.reserve((int)prot_prot.size()/jobs+1);
  for (int index=0; index<prot_prot.size(); index++) {
    if(protCounter++%jobs != job){continue;}
    //skip if no interactions known
    if(prot_prot[index].size() <= 0){continue;}
    //skip if not in gene_compute
    if(limit){if(gene_compute.find(int_id[index]) == gene_compute.end()){continue;}}
    indexes.push_back(index);
  }

  vector<string> resultLines;
  //resultLines.reserve((int)prot_prot.size()/jobs+1);
  resultLines.resize((int)indexes.size());
  if(verbose){cerr << (int)indexes.size() << " genes to calculate" << endl;}

  int COUNTER = 0;

  #if defined(_OPENMP)
  if(verbose){cerr << nthreads << " thread(s) started" << endl;}
  clock_t clock_timer;
  double wall_timer = omp_get_wtime();
  #pragma omp parallel num_threads(nthreads)
  {
  #pragma omp for
  #endif
  for (int proteinLoopCounter=0; proteinLoopCounter<indexes.size(); proteinLoopCounter++) {
    //unsigned int tid = omp_get_thread_num();
    unsigned int index = indexes[proteinLoopCounter];

    if (DEBUG){cerr << index << "\t" << int_id[index] << endl;}

    map<long double, unsigned int> chis_count;
    vector<long double> chis;
    chis.reserve(iterations);
    int actualIterations = 0;
    unsigned int pCount = 0;
    for(unsigned int i=0;i<prot_prot[index].size();i++){
      if(int_pValue[prot_prot[index][i]] < 1){pCount++;}
    }
    double chi = calculateChiSquare(&prot_prot[index]);

    if(onlyChi){                              //calculate only Chi2
      //nothing
    }else if(pCount < minimumPvalues){        //minimum number of affected neighbors not matched
      chi = 0;                                //also takes care if 0 neighbors are affected
      chis_count[0] = iterations;
    }else if(ignoreConnectivity == 0){        //replace taking network connectivity into consideration
      for(int i=0;i<iterations;i++){
        vector<unsigned int> newNodes;
        newNodes.reserve(prot_prot[index].size());
        getRandomNodes(index, &newNodes, ignoreInputNodes, ignoreCentralNode);
        if (DEBUG) cerr << prot_prot[index].size() << " " << newNodes[0] << endl;
        chis_count[calculateChiSquare(&newNodes)]++;
        //estimate p value convergence
        if(i%1000==1){
          if(DEBUG) {cout << int_id[index] << "\t" << 1.0/i << "\t" << 1.0/(i*0.1) << "\t" << calculatePfromPhis(chi, i, &chis_count, highChiIsGood) << endl;};
          if(calculatePfromPhis(chi, i, &chis_count, highChiIsGood) >= 1.0/(i*0.9)){
            actualIterations=i;
            break;
          }
        }
      }
    } else if(ignoreConnectivity == 1){      //just select random replacment nodes
      //get list of values and remove central node or all incoming ones
      vector<double> T_int_pValue;
      vector<bool> selectedNodes(proteinsInAnalysis);

      for(int i=0;i<proteinsInAnalysis;i++){selectedNodes[i] = false;}
      if(ignoreCentralNode){selectedNodes[index] = true;}
      if(ignoreInputNodes){
        for(unsigned int j=0;j<prot_prot[index].size();j++){
          selectedNodes[prot_prot[index][j]] = true;
        }
      }
      for(unsigned int i=0;i<prot_prot.size();i++){
        if(!((prot_prot[i].size()<=0) | selectedNodes[i])){
          T_int_pValue.push_back(int_pValue[i]);
        }
      }
      int T_int_pValue_max = T_int_pValue.size();
      if(verbose){cerr << index << "\t" << T_int_pValue.size() << endl;}
      //for i iterations: randomize list and pick X new values
      vector<bool> selectedPs(T_int_pValue_max);
      for(int i=0;i<iterations;i++){
        vector<double> ps;
        for(int j=0;j<T_int_pValue_max;j++){selectedPs[j] = false;}
        for(unsigned int k=0;k<prot_prot[index].size();k++){
          unsigned int randomIndex = rand()%T_int_pValue_max;
          while(selectedPs[randomIndex]){
            randomIndex = rand()%T_int_pValue_max;
          }
          ps.push_back(T_int_pValue[randomIndex]);
        }
        chis_count[calculateChiSquareFromValues(&ps)]++;
        //estimate p value convergence
        if(i%1000==1){
          if(calculatePfromPhis(chi, i, &chis_count, highChiIsGood)>=1/(i*0.9)/i){
            actualIterations = i;
            break;
          }
        }
      }
    }

    unsigned int posCounter = 0;

    if(actualIterations==0){actualIterations=iterations;}
    long double nmbScore = calculatePfromPhis(chi, actualIterations, &chis_count, highChiIsGood);
    ostringstream oss;
    oss << int_id[index] << "\t" << chi << "\t" << nmbScore;

    //Add phi values for distribution plotting
    if(printPhis==1){
      vector<double> phis;
      phis.reserve(iterations);
      for (map<long double, unsigned int>::iterator itr=chis_count.begin(); itr!=chis_count.end(); itr++){
        for(int i=0;i<itr->second;++i){
          phis.push_back(itr->first);
        }
      }
      stringstream ss;
      for(int i=0; i<phis.size();++i){
        if(i != 0) ss << ",";
        ss << phis[i];
      }
      oss << "\t" << ss.str();
    }
    resultLines[proteinLoopCounter] = oss.str();
  } //for index end
  #if defined(_OPENMP)
  } //omp end
  #endif

  if(verbose){cerr << "Done calculating all p-values" << endl;}

  if(!estimatePerformance){
    ostringstream filename;
    job++;
    filename << OUTDIR << "/" << jobName << "_j" << job << "_J" << jobs <<
                "_i" << iterations << "_f" << flipped << "_g" << ignoreInputNodes << ignoreCentralNode <<
                "_P" << pValColumn << "_r" << ignoreConnectivity << "_m" << minimumPvalues << ".txt";
    ofstream myfile (filename.str().c_str());
    if (myfile.is_open()){
      for(int i=0;i<resultLines.size();i++){
        myfile << resultLines[i] << endl;
      }
      t2=clock();
      diff = (((float)t2-(float)t1)/CLOCKS_PER_SEC);
      myfile << "#" << diff << endl;
    }
    myfile.close();
  }

  if(included){
      cout << "Printing genes with significant p-Values not found in the network to file \"notIncluded.txt\"..." << endl;
      ofstream notIncluded("notIncluded.txt");
      notIncluded << "#gene   q" << endl;
      for(int i=0; i<int_pValue.size(); i++){
         if(visited[i] != 1){
          notIncluded << int_id[i] << "   " << int_pValue[id_int[int_id[i]]] << endl;  
          }       
      }
      notIncluded.close();
  }

  //Calculate runtime and print
  t2=clock();
  diff = (((float)t2-(float)t1)/CLOCKS_PER_SEC);
  #if defined(_OPENMP)
  if(verbose) cerr << "Total runtime for " << nthreads << " thread(s): " << diff << "s" << endl;
  if(verbose) cerr << "Real runtime: " << omp_get_wtime() - wall_timer << "s" << endl;
  if(estimatePerformance) cout << omp_get_wtime()-wall_timer << endl;
  #else
  if(verbose) cerr << "Runtime for computation: " << diff << "s" << endl;
  if(estimatePerformance) cout << diff << endl;
  #endif
  return 0;
}
