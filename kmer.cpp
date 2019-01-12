#include <iostream>
#include "kmer.hpp"
#include "bloom_filter.hpp"
using namespace std;

char int2char[4] = {'T', 'C', 'A', 'G'};
char complement[4] = {2, 3, 0, 1};

//converts a TCAG character to his 2-bit representation
static int char2int(char c)
 {
     for(int i=0; i<4; i++){
       if(int2char[i] == c)
         return i;
     }
     return -1;
 }

//sets value of the k-mer
void Kmer::setValue (unsigned __int128 bin) {
  value = bin;
}

//sets k-mer length
void Kmer::setSize (int s) {
  k = s;
  mask = 3;
  for(int i=1; i<k; i++){
    mask *= 4;
    mask += 3;
  }
}

//converts k-mer string to __int128
unsigned __int128 Kmer::toInt128(std::string str){
  unsigned __int128 r = 0;
  r = char2int(str[0]);
  for(int i = 1; i < k; i++){
    r = r*4;
    r += char2int(str[i]);
  }
  return r;
}

//converts k-mer from __int128 to string
std::string Kmer::toString(){
  std::string str = "";
  unsigned __int128 temp = value;
  for(int i=0; i < k; i++){
    str += int2char[temp%4];
    temp /= 4;
  }

  //reverse string
  int n = k;
  for (int i = 0; i < n / 2; i++)
      swap(str[i], str[n - i - 1]);

  return str;
}

//gets all possible extensions of the k-mer
std::vector<unsigned __int128> Kmer::getAllExtensions(){
  std::vector<unsigned __int128> extensions;
  for(int i=0; i<4; i++){
    extensions.push_back((value * 4 + i) & mask);
  }
  return extensions;
}

//gets all valid extensions of the k-mer
std::vector<unsigned __int128> Kmer::getValidExtensions(bloom_filter& filter, cFP& falsePositives){
  std::vector<unsigned __int128> extensions;
  unsigned __int128 simpleKmer;
  for(int i=0; i<4; i++){
    simpleKmer = (value * 4 + i) & mask;
    if((filter.contains(simpleKmer) && !falsePositives.contains(simpleKmer)) /*||
  (filter.contains(revComp(simpleKmer)) && !falsePositives.contains(revComp(simpleKmer)))*/){
      extensions.push_back(simpleKmer);
    }
  }
  return extensions;
}

//returns true if it's a branching k-mer, otherwise false
bool Kmer::isBranching(bloom_filter& filter, cFP& falsePositives){
  std::vector<unsigned __int128> extensions;
  extensions = getValidExtensions(filter, falsePositives);
  if(extensions.size() > 1)
    return true;
  return false;
}

//creates reverse composition of a k-mer
__int128 Kmer::revComp(__int128 oldValue){
  __int128 newValue = 0;
  //__int128 oldValue = value;
  for(int i=0; i<k; i++){
    newValue += complement[oldValue%4];
    newValue *= 4;
    oldValue /= 4;
  }
  newValue /= 4;
  return newValue;
}

__int128 Kmer::revComp(){
  __int128 newValue = 0;
  __int128 oldValue = value;
  for(int i=0; i<k; i++){
    newValue += complement[oldValue%4];
    newValue *= 4;
    oldValue /= 4;
  }
  newValue /= 4;
  return newValue;
}

//traverse starting from given k-mer
std::string Kmer::traverse(bloom_filter& filter, cFP& falsePositives, std::vector<unsigned __int128>& pExtensions){
  unsigned __int128 startingKmer = value;
  std::vector<unsigned __int128> extensions;
  //TODO change
  std::string contig = "";//toString(); //start from current k-mer
  while(true){
    extensions = getValidExtensions(filter, falsePositives);
    if(extensions.size() != 1){
      //its a branch or a dead-end, return
      pExtensions = extensions;
      return contig;
    }
    value = extensions[0];
    contig += int2char[value%4];
  }
  return contig;
}

//converts binary representation of the extension to string
char Kmer::getExtensionChar(){
  return int2char[value%4];
}
