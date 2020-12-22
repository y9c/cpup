/*  mpileup2readcounts.cc -- Get base counts from mpileup output

    Copyright (c) 2016, Avinash Ramu

    Author: Avinash Ramu <aramu@genome.wustl.edu>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

void usage() {
  std::cerr << "samtools mpileup -f ref.fa -l regions.bed"
            << " alignments.bam | mpileup2readcounts";
}

// Convert a number to a string
inline int str_to_num(string num) {
  stringstream ss;
  int num_uint;
  ss << num;
  ss >> num_uint;
  return num_uint;
}

// DS to hold the pertinent information
class mpileup_line {
public:
  string chr;
  int pos;
  string ref_base;
  int depth;
  string bases;
  string qual;
  // Counts for different bases
  int refcount;
  int fwdcount, revcount;
  int Acount, Ccount, Gcount, Tcount, Ncount;
  int acount, ccount, gcount, tcount, ncount;
  int Insertcount, Deletecount;
  int insertcount, deletecount;

  mpileup_line() {
    chr = ref_base = bases = qual = "NA";
    depth = pos = 0;
    refcount = fwdcount = revcount = 0;
    Acount = Ccount = Gcount = Tcount = Ncount = 0;
    acount = ccount = gcount = tcount = ncount = 0;
    Insertcount = Deletecount = 0;
    insertcount = deletecount = 0;
  }
  // Set the appropriate count for ref nucleotide
  void set_ref_nuc_count() {
    switch (ref_base[0]) {
    case 'A':
    case 'a':
      Acount = fwdcount;
      acount = revcount;
      break;
    case 'C':
    case 'c':
      Ccount = fwdcount;
      ccount = revcount;
      break;
    case 'G':
    case 'g':
      Gcount = fwdcount;
      gcount = revcount;
      break;
    case 'T':
    case 't':
      Tcount = fwdcount;
      tcount = revcount;
      break;
    case 'N':
    case 'n':
      Ncount = fwdcount;
      ncount = revcount;
      break;
    // Deal with -,R,Y,K,M,S,W etc
    default:
      break;
    }
  }
  static void print_header(ostream &out = cout) {
    out << "chr"
        << "\t"
        << "pos"
        << "\t"
        << "depth"
        << "\t"
        << "ref_base"
        << "\t"
        << "refcount"
        << "\t"
        << "Acount"
        << "\t"
        << "Ccount"
        << "\t"
        << "Gcount"
        << "\t"
        << "Tcount"
        << "\t"
        << "Ncount"
        << "\t"
        << "acount"
        << "\t"
        << "ccount"
        << "\t"
        << "gcount"
        << "\t"
        << "tcount"
        << "\t"
        << "ncount"
        << "\t"
        << "Deletecount"
        << "\t"
        << "deletecount" << endl;
  }
  void print(ostream &out = cout) {
    out << chr << "\t" << pos << "\t" << depth << "\t" << ref_base << "\t"
        << refcount << "\t" << Acount << "\t" << Ccount << "\t" << Gcount
        << "\t" << Tcount << "\t" << Ncount << "\t" << acount << "\t" << ccount
        << "\t" << gcount << "\t" << tcount << "\t" << ncount << "\t"
        << Deletecount << "\t" << deletecount << endl;
  }
};

// Parse the pileup string
void parse_bases_to_readcounts(mpileup_line &ml1) {
  for (int i = 0; i < ml1.bases.length(); i++) {
    char base = ml1.bases[i];
    string indelsize_string;
    int indelsize_int = 0;
    switch (base) {
    // Match to reference
    case '.':
      ml1.fwdcount += 1;
      break;
    case ',':
      ml1.revcount += 1;
      break;
    case 'a':
      ml1.acount += 1;
      break;
    case 'A':
      ml1.Acount += 1;
      break;
    case 'c':
      ml1.ccount += 1;
      break;
    case 'C':
      ml1.Ccount += 1;
      break;
    case 'g':
      ml1.gcount += 1;
      break;
    case 'G':
      ml1.Gcount += 1;
      break;
    case 't':
      ml1.tcount += 1;
      break;
    case 'T':
      ml1.Tcount += 1;
      break;
    case 'n':
      ml1.ncount += 1;
      break;
    case 'N':
      ml1.Ncount += 1;
      break;
    // This base is deleted (--reverse-del suport)
    case '*':
      ml1.Deletecount += 1;
      break;
    case '#':
      ml1.deletecount += 1;
      break;
    // Insertion or deletion
    case '-':
    case '+':
      if (base == '+') {
        ml1.insertcount += 1;
      }
      i++;
      while (ml1.bases[i] >= 48 && ml1.bases[i] <= 57) {
        indelsize_string = indelsize_string + ml1.bases[i];
        i = i + 1;
      }
      indelsize_int = str_to_num(indelsize_string);
      i += indelsize_int - 1;
      break;
    // Reference skips
    case '<':
    case '>':
      break;
    // End of read segment
    case '$':
      break;
    // Beginning of read segment
    case '^':
      i = i + 1; // Skip quality
      break;
    default:
      string err = "Unknown ref base: ";
      err += base;
      throw runtime_error(err);
    }
  }

  ml1.refcount = ml1.fwdcount + ml1.revcount;
  ml1.set_ref_nuc_count();
}

// Split the line into the required fields and parse
void process_mpileup_line(std::string line) {
  stringstream ss(line);
  mpileup_line ml1;
  string pos, depth;

  // get chrosome ID
  getline(ss, ml1.chr, '\t');
  // get pos
  getline(ss, pos, '\t');
  ml1.pos = str_to_num(pos);
  // get ref_base
  getline(ss, ml1.ref_base, '\t');
  // get depth
  getline(ss, depth, '\t');
  ml1.depth = str_to_num(depth);
  // get bases
  getline(ss, ml1.bases, '\t');
  // get quals
  getline(ss, ml1.qual, '\t');

  parse_bases_to_readcounts(ml1);

  if (ml1.Deletecount > 0 || ml1.deletecount > 0) {
    ml1.print();
  }
}

int main(int argc, char *argv[]) {
  string line;
  mpileup_line::print_header();
  getline(cin, line);
  while (cin) {
    try {
      process_mpileup_line(line);
    } catch (const std::runtime_error &e) {
      cerr << e.what() << endl;
      cerr << "\nError parsing line " << line;
      break;
    }
    getline(cin, line);
  }
}
