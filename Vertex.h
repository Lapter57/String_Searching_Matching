#ifndef VERTEX_H
#define VERTEX_H
#include <map>
using namespace std;

class Vertex {
public:
	map<char, int> children;
	bool leaf;
	int parent;
	char parentChar;
	char myChar;
	int suffixLink;
	int endWordLink;
	int wordID;

	Vertex():leaf(false), parent(-1), suffixLink(-1), wordID(-1), endWordLink(-1) {}
};

#endif // VERTEX_H