#ifndef AHO_CORASICK_H
#define AHO_CORASICK_H
#include"Vertex.h"
#include <vector>
#include <string>

class AhoCorasick {
private:
	vector<Vertex*> trie;
	vector<int> wordsLength;
	int size;
	int root;

	void calcSuffLink(int vertex);

public:
	AhoCorasick():size(0), root(0) {
		trie.push_back(new Vertex());
		size++;
	}
	
	void addString(const string &str, int wordID);
	void prepare();
	void proccesString(const string &text, vector<unsigned> &answer);

};



#endif AHO_CORASICK_H
