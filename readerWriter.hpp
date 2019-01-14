/* Author: Ante Šušak */
class ReaderWriter {
  FILE * f;
  public:
    void writeKmer(unsigned __int128);
    void close();
    void revert();
    unsigned __int128 readOneKmer();
    ReaderWriter(FILE * f1){
      f = f1;
    }
    ReaderWriter(){
    }
};
