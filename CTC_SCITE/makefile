CXX	= g++
CXXFLAGS = -std=c++11 -O2


CTC-SCITE: findBestTrees.o CTC_treeScoring.o binTree_output.o doublets.o enum.o matrices.o mcmc.o mcmcBinTreeMove.o output.o rand.o recMut.o scoreBinTree.o scoreTree.o treelist.o trees.o
	${CXX} ${CXXFLAGS} -o CTC-SCITE findBestTrees.o CTC_treeScoring.o binTree_output.o doublets.o enum.o matrices.o mcmc.o mcmcBinTreeMove.o output.o rand.o recMut.o scoreBinTree.o scoreTree.o treelist.o trees.o



findBestTrees.o: doublets.h matrices.h treelist.h trees.h output.h mcmc.h rand.h scoreTree.h binTree_output.h CTC_treeScoring.h
	${CXX} ${CXXFLAGS} -c findBestTrees.cpp

CTC_treeScoring.o: binTree_output.h scoreTree.h matrices.h rand.h trees.h matrices.h mcmcBinTreeMove.h mcmc.h CTC_treeScoring.h
	${CXX} ${CXXFLAGS} -c CTC_treeScoring.cpp

binTree_output.o: output.h scoreTree.h matrices.h trees.h binTree_output.h
	${CXX} ${CXXFLAGS} -c binTree_output.cpp

doublets.o: matrices.h treeList.h trees.h mcmc.h scoreTree.h scoreBinTree.h rand.h output.h mcmcBinTreeMove.h mcmcTreeMove.h recMut.h doublets.h
	${CXX} ${CXXFLAGS} -c doublets.cpp

enum.o: treelist.h trees.h matrices.h rand.h scoreTree.h enum.h
	${CXX} ${CXXFLAGS} -c enum.cpp

matrices.o: matrices.h
	${CXX} ${CXXFLAGS} -c matrices.cpp

mcmc.o: matrices.h treelist.h trees.h mcmc.h scoreTree.h scoreBinTree.h rand.h output.h mcmcBinTreeMove.h mcmcTreeMove.h
	${CXX} ${CXXFLAGS} -c mcmc.cpp

mcmcBinTreeMove.o: matrices.h treelist.h trees.h scoreTree.h rand.h output.h mcmcBinTreeMove.h
	${CXX} ${CXXFLAGS} -c mcmcBinTreeMove.cpp

output.o: output.h scoreTree.h trees.h recMut.h matrices.h
	${CXX} ${CXXFLAGS} -c output.cpp

rand.o: rand.h matrices.h
	${CXX} ${CXXFLAGS} -c rand.cpp

recMut.o: matrices.h trees.h scoreTree.o recMut.h
	${CXX} ${CXXFLAGS} -c recMut.cpp

scoreBinTree.o: scoreBinTree.h matrices.h trees.h
	${CXX} ${CXXFLAGS} -c scoreBinTree.cpp

scoreTree.o: scoreTree.h matrices.h treelist.h trees.h
	${CXX} ${CXXFLAGS} -c scoreTree.cpp

treelist.o: matrices.h rand.h treelist.h
	${CXX} ${CXXFLAGS} -c treelist.cpp	

trees.o: matrices.h treeList.o trees.h rand.h output.h
	${CXX} ${CXXFLAGS} -c trees.cpp

clean:
	rm *.o