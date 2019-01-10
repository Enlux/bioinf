#include <vector>
#include "bloom_filter.hpp"
#include "cFP.h"
#include "string"
//#include "cFP.h"
class Kmer {
    unsigned __int128 mask;
    int k;
  public:
    unsigned __int128 value;
    void setValue (unsigned __int128);
    void setSize (int);
    char getExtensionChar();
    __int128 revComp();
    __int128 revComp(__int128);
    std::vector<unsigned __int128> getAllExtensions();
    std::vector<unsigned __int128> getValidExtensions(bloom_filter&, cFP&);
    bool isBranching(bloom_filter&, cFP&);
    unsigned __int128 toInt128(std::string);
    std::string toString();
    std::string traverse(bloom_filter&, cFP&, std::vector<unsigned __int128>&);
};
