/* Author: Lucijan Peš */
#include <iostream>
#include <fstream>
#include <string>

#include <jellyfish/mer_dna.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/hash_counter.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/mer_overlap_sequence_parser.hpp>
#include <jellyfish/mer_iterator.hpp>
#include <vector>
#include <set>
#include "Util.hpp"
#include "readerWriter.hpp"
#include "bloom_filter.hpp"

using namespace std;


typedef jellyfish::cooperative::hash_counter<jellyfish::mer_dna>                  mer_hash_type;
typedef jellyfish::mer_overlap_sequence_parser<jellyfish::stream_manager<char**>> sequence_parser_type;
typedef jellyfish::mer_iterator<sequence_parser_type, jellyfish::mer_dna>         mer_iterator_type;


class mer_counter : public jellyfish::thread_exec {
  mer_hash_type&                    mer_hash_;
  jellyfish::stream_manager<char**> streams_;
  sequence_parser_type              parser_;
  const bool                        canonical_;



public:
  mer_counter(int nb_threads, mer_hash_type& mer_hash,
              char** file_begin, char** file_end,
              bool canonical)
    : mer_hash_(mer_hash)
    , streams_(file_begin, file_end)
    , parser_(jellyfish::mer_dna::k(), streams_.nb_streams(), 3 * nb_threads, 4096, streams_)
    , canonical_(canonical)
  { }

  virtual void start(int thid) {
    mer_iterator_type mers(parser_, canonical_);

    for( ; mers; ++mers)
      mer_hash_.add(*mers, 1);
    mer_hash_.done();
  }
};


//Counts all k-mers and saves only ones that are repeated k or more times.
int countKmers(char** readFilename, const std::string& dumpFilename, int k, int discardThreshold) {


  jellyfish::mer_dna::k(k); // Set length of mers
  const uint64_t hash_size    = 10000000; // Initial size of hash.
  const uint32_t num_reprobes = 126;
  const uint32_t num_threads  = 16; // Number of concurrent threads
  const uint32_t counter_len  = 4;  //was 7, Minimum length of counting field
  const bool     canonical    = false;//true; // Use canonical representation

  // create the hash
  mer_hash_type mer_hash(hash_size, jellyfish::mer_dna::k()*2, counter_len, num_threads, num_reprobes);

  // count the kmers
  //first argument of all arguments
  mer_counter counter(num_threads, mer_hash, readFilename, readFilename + 1, canonical);

  counter.exec_join(num_threads);

  const auto jf_ary = mer_hash.ary();

  //get number of saved kmers later to be used for bloom aproximation
  int numSavedKmers = 0;

  // Open file and dump all the k-mers to file which satisfy the threshold
  Kmer kmer;
  kmer.setSize(k);
  ReaderWriter rw = ReaderWriter(fopen(dumpFilename.c_str(), "wb"));
  const auto end = jf_ary->end();
  for(auto it = jf_ary->begin(); it != end; ++it) {
    auto& key_val = *it;
    if(key_val.second >= discardThreshold){
      //GHETTOHACK
      std::stringstream buffer;
      buffer << key_val.first;
      //cout << buffer.str() << "\n";
      kmer.setValue(kmer.toInt128(buffer.str()));
      rw.writeKmer(kmer.value);
      numSavedKmers++;
    }
  }
  rw.close();
  return numSavedKmers;
}

bloom_filter constructFilter(int numSavedKmers, double falsePositiveProbability){
  bloom_parameters parameters;
  // How many elements roughly do we expect to insert?
  parameters.projected_element_count = numSavedKmers;
  // Maximum tolerable false positive probability? (0,1)
  parameters.false_positive_probability = falsePositiveProbability;
  // Simple randomizer (optional)
  parameters.random_seed = 0xA5A5A5A5;
  if (!parameters)
  {
     std::cout << "Error - Invalid set of bloom filter parameters!" << std::endl;
     exit(1);
  }
  parameters.compute_optimal_parameters();
  bloom_filter filter(parameters);
  return filter;
}

cFP constructCFP(const std::string& dumpFilename, bloom_filter& filter, int maximumMBUsage, int k, int numSavedKmers){
  cFP falsePositives;
  Kmer kmer;
  unsigned __int128 rawKmer;
  std::vector<unsigned __int128> kmerExtensions;
  kmer.setSize(k);


  ReaderWriter rwS = ReaderWriter(fopen(dumpFilename.c_str(), "rb"));
  ReaderWriter rwD = ReaderWriter(fopen("D.temp", "wb"));
  ReaderWriter rwDnew;


  int numWrittenKmers = 0;
  for (int i=0; i < numSavedKmers; i++){
    rawKmer = rwS.readOneKmer();
    kmer.setValue(rawKmer);
    kmerExtensions = kmer.getAllExtensions();
    for (auto &extension : kmerExtensions) // access by reference to avoid copying
    {
        if(filter.contains(extension)){
          rwD.writeKmer(extension);
          numWrittenKmers++;
        }
    }
  }
  cout << "Pass 0: wrote " << numWrittenKmers << " extension k-mers to file.\n";
  rwD.close();
  rwS.revert();

  int i = 0;
  int j = 0;
  int unitSize = sizeof(unsigned __int128)*4;
  int numKmersLeft;
  cFP P;


  while(j < numSavedKmers){

    rwD = ReaderWriter(fopen("D.temp", "rb"));
    rwDnew = ReaderWriter(fopen("Dnew.temp", "wb"));

    numKmersLeft = 0;
    P = cFP();

    //fill the P set with the maximum size allowed
    while((unitSize*P.numFP)/(1024.0*1024.0) < (double)maximumMBUsage && j < numSavedKmers){
      rawKmer = rwS.readOneKmer();
      P.add(rawKmer);
      j++;
    }

    //read i times filtered D
    for(int k = 0; k < numWrittenKmers; k++){
      rawKmer = rwD.readOneKmer();
      if(!P.contains(rawKmer)){
        rwDnew.writeKmer(rawKmer);
        numKmersLeft++;
      }
    }

    numWrittenKmers = numKmersLeft;
    i++;
    cout << "Pass " << i << ": written " << numWrittenKmers << " k-mers.\n";
    rwD.close();
    rwDnew.close();
    //remove old D.temp, replace with new
    std::remove("D.temp");
    std::rename("Dnew.temp", "D.temp");
    //free P from memory
    P.set.clear();
  }

  //put the final D in the structure
  rwD = ReaderWriter(fopen("D.temp", "rb"));
  for(int k = 0; k < numWrittenKmers; k++){
    rawKmer = rwD.readOneKmer();
    falsePositives.add(rawKmer);
  }

  rwD.close();
  rwS.close();
  return falsePositives;
}


//Traverses de Bruijn graph starting from a set of branching nodes and finds longest contigs.
void GraphTraversal(const std::string& outputFilename, const std::string& countingFilename, bloom_filter& filter, cFP& falsePositives, int k, int numSavedKmers, int maxDepth){

  //count branching kmers
  std::set<unsigned __int128> branchingKmers;
  Kmer kmer;
  kmer.setSize(k);
  unsigned __int128 simpleKmer;
  ReaderWriter rwS = ReaderWriter(fopen(countingFilename.c_str(), "rb"));
  for(int i=0; i<numSavedKmers; i++){
    simpleKmer = rwS.readOneKmer();
    kmer.setValue(simpleKmer);
    if(kmer.isBranching(filter, falsePositives)){
      branchingKmers.insert(simpleKmer);
    }
  }
  cout << "Inserted " << branchingKmers.size() << " branching k-mers.\n";

  set<unsigned __int128>::iterator bkmer;
  Kmer startingKmer;
  Kmer extensionKmer;
  startingKmer.setSize(k);
  extensionKmer.setSize(k);
  std::vector<unsigned __int128> extensions;
  std::set<unsigned __int128> markedKmers; //used to mark so we don't double traverse
  ofstream contigFile;
  contigFile.open(outputFilename);
  int contigNo = 0;
  int longestContig = 0;
  int len;
  int maxBP = maxDepth;

  for (bkmer = branchingKmers.begin(); bkmer != branchingKmers.end(); ++bkmer) {
      //find all extensions of a branching k-mer
      startingKmer.setValue(*bkmer);
      extensions = startingKmer.getValidExtensions(filter, falsePositives);
      for (auto &extension : extensions) // access by reference to avoid copying
      {
        extensionKmer.setValue(extension);
        //skip if branching or already marked
        if(extensionKmer.isBranching(filter, falsePositives) || markedKmers.find(extension) != markedKmers.end()){
          continue;
        }

        //mark the k-mer
        markedKmers.insert(extension);

        std::vector<unsigned __int128> pExtensions;
        //traverse starting from that extension
        string contig = extensionKmer.toString();
        contig += extensionKmer.traverse(filter, falsePositives, pExtensions);
        //TODO depth algo here
        if(pExtensions.size() > 1){
          //recursive breadth first search of longest contig
          int length;
          contig += findLongest(pExtensions, maxBP, 0, k, filter, falsePositives, markedKmers, length);
        }

        //only save longer ones
        len = contig.length();
        if(len > 2*(k+1)){
          string info = ">contig" + std::to_string(contigNo) + "__len__" + std::to_string(contig.length());
          contigNo+=1;
          if(len > longestContig)
            longestContig = len;
          contigFile << info << "\n" << contig << "\n";
        }
      }
  }

  cout << "Written " << contigNo << " contig(s) to a file.\n";
  cout << "Length of longest contig written: " <<  longestContig << "\n";
}


//recursively finds longest path starting from extensions of a branching node
std::string findLongest(std::vector<unsigned __int128> extensions, int maxBP, int depth,  int k, bloom_filter& filter, cFP& falsePositives, std::set<unsigned __int128>& markedKmers,  int& length){
  if(depth >= maxBP){
    //cout << "Warning! Max depth exceeded. Returning 0 string.\n";
    length = 0;
    return "";
  }

  Kmer kmer;
  kmer.setSize(k);
  vector <int> extensionLengths;
  vector <std::string> extensionStrings;
  bool passed = false;
  for (auto &extension : extensions) // access by reference to avoid copying
  {
    if(markedKmers.find(extension) != markedKmers.end()){
      continue;
    }
    int slength;
    kmer.setValue(extension);
    markedKmers.insert(extension);
    std::vector<unsigned __int128> pExtensions;
    std::string contig = "";
    contig += kmer.getExtensionChar();
    contig += kmer.traverse(filter, falsePositives, pExtensions);
    slength = contig.length();
    if(pExtensions.size() > 1){
      //recursive breadth first search of longest contig
      int extendedLength;
      //recursion
      contig += findLongest(pExtensions, maxBP, depth+1, k, filter, falsePositives, markedKmers, extendedLength);
      slength += extendedLength;
    }
    extensionLengths.push_back(slength);
    extensionStrings.push_back(contig);
    passed = true;
  }

  //return one with max length
  if(!passed){
    length = 0;
    return "";
  }

  int maxLength = extensionLengths[0];
  int maxInd = 0;
  for(int i=1; i < extensionLengths.size(); i++){
    if(extensionLengths[i] > maxLength){
      maxInd = i;
      maxLength = extensionLengths[i];
    }
  }

  length = maxLength;
  return extensionStrings[maxInd];
}
