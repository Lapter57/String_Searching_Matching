#ifndef VERTEX_H
#define VERTEX_H
#include <unordered_map>
using namespace std;

class Vertex {
public:
	unordered_map<char, int> children;
	bool leaf;
	int parent;
	char parentChar;
	int suffixLink;
	int endWordLink;
	int wordID;

	Vertex():leaf(false), parent(-1), suffixLink(-1), wordID(-1), endWordLink(-1) {}
	~Vertex() {
		children.clear();
	}
};

#endif // VERTEX_H