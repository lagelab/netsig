#include <iostream>
#include <fstream>
#include <sstream>

#include <string>

#include <vector>
#include <map>

#include <iterator>

#include <dirent.h>
#include <stdint.h>
#include <stdlib.h>

using namespace std;

vector<string> getFilesFromDirectory(string path = ".") {
  unsigned char isFile = 0x8;
  unsigned char isFolder = 0x4;

  DIR*    dir;
  dirent* pdir;
  vector<string> files;

  dir = opendir(path.c_str());
  while (pdir = readdir(dir)) {
    //skip hidden files and special directories
    if(pdir->d_name[0] == '.'){continue;}
    //skip directories
    if(pdir->d_type != isFile){continue;}
    //take only the _InWeb3_HUGO_sym files
    //string filename(pdir->d_name);
    //string search("_dataLocAnnot");
    //size_t found = filename.find(search);
    //if(found != string::npos){continue;}
    //search = "sif";
    //found = filename.find(search);
    //
    //if(found == string::npos){continue;}
    files.push_back(pdir->d_name);
  }

  return files;
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
    unsigned int newInt = (unsigned int)atoi(tokens[0].c_str());
    if(newInt >= int_id_->size()){
      int_id_->resize(newInt+1);
    }
    (*int_id_)[newInt] = tokens[1];
    id_int_->insert(pair<string, unsigned int>(tokens[1], newInt));
  }
  return 1;
}

int main (int argc, char* argv[]){
  vector<string> int_id;
  map<string, unsigned int> id_int;

  if(argc<5){
    cout << "Usage: " << endl << "./inwebBinaryNetwork SOURCEDIR TARGETDIR MAPPING_FILE ACTION" << endl <<
    "ACTION: 1 = create binary network" << endl <<
    "        2 = create plain text network" << endl;
    return 0;
  }

  string NETWORKDIR = argv[1];
  string OUTDIR = argv[2];
  string MAPPINGFILE = argv[3];
  int ACTION = atoi(argv[4]);

  vector<string> files;
  files = getFilesFromDirectory(NETWORKDIR);

  cerr << files.size() << " networks to transform" << endl;
  int size = files.size();

  static uint16_t prot_prot [23000][23000];

  unsigned long OFFSET = 0;

  if(!loadMappingFile(MAPPINGFILE, &int_id, &id_int)){
    return 0;
  }
  //set iterator
  for(unsigned long i=OFFSET;i<size+OFFSET;i++){

    //====================\\
    //                    \\
    // create binary file \\
    //                    \\
    //====================\\

    if(ACTION == 1){
      for(unsigned j=0;j<23000;j++){
        for(unsigned k=0;k<23000;k++){
          prot_prot[j][k]=0;
        }
      }
      ostringstream file;
      file << NETWORKDIR << files[i-OFFSET];

      ifstream in(file.str().c_str());
      
      if(!in){
        cerr << "Networkfile not found: " << endl << file.str() << endl;
        return 0;
      }

      ostringstream filename;
      filename << OUTDIR << "/" << files[i-OFFSET] << ".bin";
      ofstream o(filename.str().c_str(),ios::binary);

      string line;
      while(getline(in, line)){
        vector<string> tokens;
        if(line[0] == '#'){continue;}
        istringstream iss(line);
        string token;
        while(getline(iss, token, '-')){
          tokens.push_back(token);
        }
        uint16_t id1 = (uint16_t)atof(tokens[0].c_str());
        uint16_t id2 = (uint16_t)atof(tokens[1].c_str());
        if(prot_prot[id1][id2] == 1){ continue; }
        prot_prot[id1][id2] = 1;
        prot_prot[id2][id1] = 1;
        o.write((char*)&id1,sizeof(id1));
        o.write((char*)&id2,sizeof(id2));
      }
      in.close();
      o.close();

    //========================\\
    //                        \\
    // create plain text file \\
    //                        \\
    //========================\\

    }else{
      ostringstream file;
      file << NETWORKDIR << files[i-OFFSET];
      ifstream in(file.str().c_str());

      cerr << file.str().c_str() << endl;

      ostringstream filename;
      filename << OUTDIR << "/" << files[i-OFFSET] << ".txt";
      ofstream o(filename.str().c_str(),ios::binary);

      in.seekg (0, in.end);
      long length = (long)in.tellg();
      in.seekg (0, in.beg);

      short int1, int2;
      while((long)in.tellg() != length){
        in.read((char*)&int1, 2);
        in.read((char*)&int2, 2);
        //Check for non-existing mappings
        //cout << int1 << "\t" << int_id[int1] << "\t" << int2 << "\t" << int_id[int2] << endl;
        o << int_id[int1] << "\t" << int_id[int2] << endl;
      }
      in.close();
      o.close();
    }
  }
  return 0;
}
