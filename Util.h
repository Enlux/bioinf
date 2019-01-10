#include "bloom_filter.hpp"
#include "kmer.h"
int countKmers(char** readFilename, const std::string& dumpFilename, int k, int discardThreshold);
bloom_filter constructFilter(int numSavedKmers, double falsePositiveProbability);
cFP constructCFP(const std::string& dumpFilename, bloom_filter& filter, int maximumMBUsage, int k, int numSavedKmers);
void GraphTraversal(const std::string& outputFilename, const std::string& countingFilename, bloom_filter& filter, cFP& falsePositives, int k, int numSavedKmers, int maxDepth);
std::string findLongest(std::vector<unsigned __int128> extensions, int maxBP, int depth,  int k, bloom_filter& filter, cFP& falsePositives, std::set<unsigned __int128>& markedKmers, int& length);
