#Author: Lucijan Peš
#!/usr/bin/python3

#this is all very simplified, kmers are in memory all the time
#not dynamically loaded, trivial structures are used
class GraphTraversal():

    def __init__(self, bloom, cFP, kmers, kmerSize, outputFilename):
        self.bloom = bloom
        self.cFP = cFP
        self.kmers = kmers
        self.kmerSize = kmerSize
        self.branchingKmers = []
        self.markedKmers = {}
        self.bases = ['A', 'C', 'T', 'G']
        self.reverseDict = {
            'A' : 'T', 'C' : 'G', 'T' : 'A' , 'G' : 'C'
        }
        self.output = open(outputFilename, "w+")
        self.run()


    def revComp(self, kmer):
        newKmer = kmer[::-1]
        newKmerRev = ""
        for i in range(self.kmerSize):
            newKmerRev += self.reverseDict[newKmer[i]]
        return newKmerRev

    def isBranching(self, kmer):
        #print(self.revComp(kmer))
        #print(self.getNumExtensions(self.revComp(kmer)))
        return self.getNumExtensions(kmer) > 1 #\
        #and self.getNumExtensions(self.revComp(kmer)) == 1)

    def isExtension(self, kmer):
        if kmer.encode() in self.bloom and kmer not in self.cFP:
            return True
        return False

    def getNumExtensions(self, kmer):
        numExt = 0
        for extension in self.bases:
            newKmer = kmer[1:] + extension
            if self.isExtension(newKmer):
                numExt += 1
                #save it for retrieveing later
                self.extCache = extension
        return numExt

    def isMarked(self, kmer):
        return kmer in self.markedKmers

    def getStartingKmer(self, branchingKmer):
        for extension in self.bases:
            currentKmer = branchingKmer[1:] + extension
            if self.isExtension(currentKmer):
                if self.isMarked(currentKmer) or self.isBranching(currentKmer):
                    continue
                self.markedKmers[currentKmer] = True
                return currentKmer
        return None

    def avance(self, kmer, reverse = False):
        numExt = self.getNumExtensions(kmer)
        if numExt == 1:

            #check for inbranching
            newKmer = kmer[1:] + self.extCache
            revKmer = self.revComp(newKmer)
            if self.getNumExtensions(revKmer) > 1:
                print("Num Extensions:", self.getNumExtensions(revKmer))
                return False
            #true if no out or inbraching or deadend
            self.extCache = newKmer
            return True
        return False

    #flip k-mer before everything if right
    def traverse(self, startingKmer, right = False):
        resultingSeq = ""

        if right:
            currentKmer = self.revComp(startingKmer)
        else:
            currentKmer = startingKmer

        while self.avance(currentKmer):
            #retrieve current k-mer
            extension = self.extCache[-1]
            resultingSeq += extension
            currentKmer = currentKmer[1:] + extension
            self.markedKmers[currentKmer] = True
            if currentKmer == startingKmer:
                break

        return resultingSeq

    def run(self):

        #find branching k-mers
        for kmer in self.kmers:
            if self.isBranching(kmer):
                self.branchingKmers.append(kmer)

        contigNo = 0

        print("Broj branši:", len(self.branchingKmers))
        for kmer in self.branchingKmers:
            while True:
                #find all extensions of branching k-mer
                startingKmer = self.getStartingKmer(kmer)
                if startingKmer is None:
                    break

                #traverse the graph starting from the startingKmer
                leftTraversal = self.traverse(startingKmer)
                #rightTraversal = self.traverse(startingKmer, right = True)

                #create the contig
                print(leftTraversal)
                print(startingKmer)
                #print(rightTraversal)
                contig = startingKmer + leftTraversal #+ startingKmer + rightTraversal

                if(len(contig) < 100):
                    continue
                #write contig to file
                self.output.write(">" + str(contigNo) + "__len__" + str(len(contig)) + "\n")
                self.output.write(contig + "\n")
                contigNo += 1

        self.output.close()
