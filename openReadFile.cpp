/* ofstream: Stream class to write on files
ifstream: Stream class to read from files
fstream: Stream class to both read and write from/to files.*/

// writing on a text file
#include <iostream>
#include <fstream>
using namespace std;

int main() {
	ofstream myfile("example.txt");
	if (myfile.is_open())
	{
		myfile << "This is a line.\n";
		myfile << "This is another line.\n";
		myfile.close();
	}
	else cout << "Unable to open file";
	return 0;
}

// reading a text file
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main () {
  string line;
  ifstream myfile ("example.txt");
  if (myfile.is_open())
  {
    while ( getline (myfile,line) )
    {
      cout << line << '\n';
    }
    myfile.close();
  }

  else cout << "Unable to open file";

  return 0;
}
