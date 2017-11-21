#include "AhoCorasick.h"
#include <queue>

void AhoCorasick::addString(const string &str, int wordID) {

	int curVertex = root;
	for (int i = 0; i < str.length(); i++) {
		char ch = str[i];

		auto children = trie[curVertex]->children.find(ch);
		if (children == trie[curVertex]->children.end()) {
			trie.push_back(new Vertex());
			trie[size]->suffixLink = -1;
			trie[size]->parent = curVertex;
			trie[size]->parentChar = ch;
			trie[curVertex]->children[ch] = size;
			size++;
		}
		curVertex = trie[curVertex]->children[ch];
	}

	trie[curVertex]->leaf = true;
	trie[curVertex]->wordID = wordID;
	wordsLength.push_back(str.length());
}


void AhoCorasick::prepare() {
	queue<int> vertexQueue;
	vertexQueue.push(root);
	while (vertexQueue.size() > 0) {
		int curVertex = vertexQueue.front();
		vertexQueue.pop();
		calcSuffLink(curVertex);

		for (auto && children : trie[curVertex]->children) {
			vertexQueue.push(trie[curVertex]->children[children.first]);
		}
	}
}

void AhoCorasick::calcSuffLink(int vertex) {
	
	// empty string
	if (vertex == root) {
		trie[vertex]->suffixLink = root;
		trie[vertex]->endWordLink = root;
		return;
	}

	// one character substring
	if (trie[vertex]->parent == root) {
		trie[vertex]->suffixLink = root;
		if (trie[vertex]->leaf)
			trie[vertex]->endWordLink = vertex;
		else
			trie[vertex]->endWordLink = trie[trie[vertex]->suffixLink]->endWordLink;
		return;
	}

	int curBetterVertex = trie[trie[vertex]->parent]->suffixLink;
	char chVertex = trie[vertex]->parentChar;

	while (true){
		auto children = trie[curBetterVertex]->children.find(chVertex);
		if (children != trie[curBetterVertex]->children.end()) {
			trie[vertex]->suffixLink = trie[curBetterVertex]->children[chVertex];
			break;
		}
		if (curBetterVertex == root) {
			trie[vertex]->suffixLink = root;
			break;
		}
		curBetterVertex = trie[curBetterVertex]->suffixLink;
	}
	
	if (trie[vertex]->leaf)
		trie[vertex]->endWordLink = vertex;
	else
		trie[vertex]->endWordLink = trie[trie[vertex]->suffixLink]->endWordLink;
}

void AhoCorasick::proccesString(const string &text, vector<unsigned> &answer) {

	int currentState = root;

	for (int i = 0; i < text.length(); i++) {
		while (true) {
			auto children = trie[currentState]->children.find(text[i]);
			if (children != trie[currentState]->children.end()) {
				currentState = trie[currentState]->children[text[i]];
				break;
			}
		
			if (currentState == root)
				break;
			currentState = trie[currentState]->suffixLink;
		}

		int checkState = currentState;

		while (true){
			checkState = trie[checkState]->endWordLink;
			if (checkState == root)
				break;
			answer.push_back(i + 1 - wordsLength[trie[checkState]->wordID]);
			checkState = trie[checkState]->suffixLink;
		}
	}
}