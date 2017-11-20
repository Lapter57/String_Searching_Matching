#include <iostream>
#include <time.h>
#include "StringSearching.h"
using namespace std;

int main() {
	srand((unsigned)time(0));
	ofstream out;
	vector<char> alphabet;
	alphabet.push_back('0');
	alphabet.push_back('1');
	/*alphabet.push_back('2');
	alphabet.push_back('3');
	alphabet.push_back('4');
	alphabet.push_back('5');
	alphabet.push_back('6');
	alphabet.push_back('7');
	alphabet.push_back('8');
	alphabet.push_back('9');*/
	StringSearching str_search(alphabet, 1000, 5);

	out.open("dataNoFixN.txt");
	try {
		str_search.investigate("text", 500, 1000, 100, 10000, out);

		out.open("dataNoFixM.txt");
		str_search.investigate("pattern", 500, 1000, 100, 1000, out);

		out.open("dataNoFixNP.txt");
		str_search.investigate("num_pattern", 500, 1000, 100, 1000, out);
	}
	catch (const char* e) {
		cout << e;
	}




	out.close();
	return 0;
}