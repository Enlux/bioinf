#include <iostream>
#include <string>
#include "Util.h"
#include "readerWriter.h"
#include <cstdlib>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

using namespace std;

void timeIt(std::chrono::_V2::system_clock::time_point& start, std::string action){
  std::chrono::duration<double> elapsed = Clock::now() - start;
  cout << "\033[92m" << action << " done in: " << elapsed.count() << " seconds.\033[0m\n\n";
  start = Clock::now();
}

int main(int argc, char *argv[]) {
  if(argc < 4 || argc > 7){
    cout << "Usage: ./main inputFilename kMerSize outputFilename [discardThreshold] [maximumMBUsage] [maxDepth]\n";
    exit(1);
  }

  int k = atoi(argv[2]); //27
  unsigned int discardThreshold;
  if(argc >= 5)
    discardThreshold = atoi(argv[4]);
  else
    discardThreshold = 3;
  string outputFilename = argv[3];

  unsigned int maximumMBUsage;
  if(argc >= 6)
    maximumMBUsage = atoi(argv[5]);
  else
    maximumMBUsage = 4096;//4 gb, for now

  unsigned int maxDepth;
  if(argc >= 7)
    maxDepth = atoi(argv[6]);
  else
    maxDepth = 5;


  cout << "Maximum MB usage for cFP construction is set to " << maximumMBUsage << ".\n";
  cout << "Discard threshold set to " << discardThreshold << ".\n";
  cout << "k-mer size set to " << k << ".\n";
  cout << "Max recursion depth set to "<< maxDepth << ".\n";
  //Clock start;
  double falsePositiveProbability = 2.98/16000; //optimal size from the paper
  string countingFilename = "counting.temp";
  Kmer kmer;
  kmer.setSize(k);

  cout << "--------------------------\n";
  auto globalStart = Clock::now();

  //count k-mers and get number of saved kmers for bloom
  cout << ">Counting k-mers, saving all that are repeated " << discardThreshold << " times or more.\n";
  auto start = Clock::now();
  int numSavedKmers = countKmers(argv + 1, countingFilename, k, discardThreshold);
  cout << "Wrote " << numSavedKmers << " k-mers to " << countingFilename << "\n";
  timeIt(start, "Counting k-mers");

  //construct optimal bloom filter
  bloom_filter filter = constructFilter(numSavedKmers, falsePositiveProbability);

  cout << ">Inserting k-mers into bloom...\n";
  //read and insert k-mers into bloom filter
  unsigned __int128 simpleKmer;
  ReaderWriter rw = ReaderWriter(fopen(countingFilename.c_str(), "rb"));
  for (int i=0; i < numSavedKmers; i++){
    simpleKmer = rw.readOneKmer();
    filter.insert(simpleKmer);
  }
  rw.close();
  timeIt(start, "Bloom insertion");

  //construct cFP from bloom + file
  cout << ">Building cFP...\n";
  cFP falsePositives = constructCFP(countingFilename, filter, maximumMBUsage, k, numSavedKmers);
  timeIt(start, "Building cFP set");

  cout << ">Finding contigs...\n";
  GraphTraversal(outputFilename, countingFilename, filter, falsePositives, k, numSavedKmers, maxDepth);
  timeIt(start, "Graph traversal");


  //remove temp files
  std::remove("D.temp");
  std::remove(countingFilename.c_str());

  timeIt(globalStart, "Algorithm");
  return 0;
}
