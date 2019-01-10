#include <iostream>
#include "readerWriter.h"
#include "kmer.h"

//writes one k-mer to file
void ReaderWriter::writeKmer(unsigned __int128 i) {
  fwrite(&i, 1, sizeof(unsigned __int128), f);
}

//closes the file
void ReaderWriter::close() {
  fclose(f);
}

//puts pointer of the file to the beginning
void ReaderWriter::revert() {
  rewind(f);
}

//reads one k-mer from file and moves pointer by the size of one k-mer
unsigned __int128 ReaderWriter::readOneKmer(){
  unsigned __int128 simpleKmer;
  fread(&simpleKmer, sizeof(unsigned __int128), 1, f);
  return simpleKmer;
}
