/* Author: Ante Šušak */
#include <set>
//#include <vector>
//#include <algorithm>
#include <unordered_set>



class cFP {
  public:
    std::set<unsigned __int128> set;
    //std::vector<unsigned __int128> vector;
    int numFP;
    bool contains(unsigned __int128 value);
    void add(unsigned __int128 value);
    cFP(){
      std::set<unsigned __int128> set;
      //std::vector<unsigned __int128> vector;
      numFP = 0;
    }
};
