#include "cFP.h"

//checks if cFP set contains a k-mer
bool cFP::contains(unsigned __int128 value2){
  return set.find(value2) != set.end();
}

//adds false positive k-mer to the set
void cFP::add(unsigned __int128 value2){
  set.insert(value2);
  numFP += 1;
}
