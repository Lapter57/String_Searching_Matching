#include <iostream>
#include <time.h>
#include "StringSearching.h"
using namespace std;

int main() {
	srand((unsigned)time(0));
	ofstream out;
	ofstream out2;
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

	unsigned max_len = 1000000;
	unsigned max_num_patterns = 50;
	StringSearching str_search(alphabet, max_len, max_num_patterns);

	try {
		unsigned change_value_start = 10000;
		unsigned change_value_end = 100000;
		unsigned step = 10000;
		unsigned num_iter = 10000;
		out.open("dataNoFixN.txt");
		out2.open("OccurN.txt");
		str_search.investigate("text", change_value_start, change_value_end, step, num_iter, out, out2);
		out.close();

		change_value_start = 90000;
		change_value_end = 100000;
		step = 10000;
		num_iter = 1000;
		out.open("dataNoFixM.txt");
		out2.open("OccurM.txt");
		str_search.investigate("pattern", change_value_start, change_value_end, step, num_iter, out, out2);
		out.close();
		out2.close();

		change_value_start = 10;
		change_value_end = 50;
		step = 10;
		num_iter = 1000;
		out2.open("OccurNP.txt");
		out.open("dataNoFixNP.txt");
		str_search.investigate("num_pattern", change_value_start, change_value_end, step, num_iter, out, out2);
		out.close();
		out2.close();
	}
	catch (const char* e) {
		cout << e;
	}

	return 0;
}