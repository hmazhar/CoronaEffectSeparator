#include <iostream>
#include <sstream>
#include <algorithm>
#include <cstring>
#include <zlib.h>
#include "core/ChVector.h"
#include "core/ChQuaternion.h"
#include "collision/ChCModelBullet.h"

using namespace std;

int main(int argc, char* argv[]) {

  stringstream input_file_ss;
  input_file_ss << argv[1] << ".txt";

  gzFile gz_file = gzopen(input_file_ss.str().c_str(), "rb");
  unsigned long int size;
  gzread(gz_file, (void*)&size, sizeof(size));
  std::string data;
  data.resize(size / sizeof(char));
  gzread(gz_file, (void*)data.data(), size);
  gzclose(gz_file);

  stringstream output_file_ss;
  output_file_ss<<argv[1]<<"_chrono.txt";

  ofstream ofile(output_file_ss.str());
  ofile<<data;
  ofile.close();

  return 0;
}
