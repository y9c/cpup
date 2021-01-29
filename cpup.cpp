/*
 * cpup.cpp
 * Copyright (C) 2021 Ye Chang <yech1990@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include <cstring>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

using namespace std;

void usage() {
  std::cerr << "Usage: "
            << "samtools mpileup -d 0 -Q 0 --reverse-del "
               "-f <.fa> <.bam> "
               "| cpup"
            << endl;
}

// global variables
vector<string> names = {
    "depth",
    "ref",
    // forward strand
    //"fwd",
    "A",
    "C",
    "G",
    "T",
    "N",
    "Gap",
    "Insert",
    "Delete",
    // reverse strand
    //"rev",
    "a",
    "c",
    "g",
    "t",
    "n",
    "gap",
    "insert",
    "delete"};

string count_sep = ",";
string sample_sep = "\t";

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
  int pos, nsample;
  string chr, ref_base;
  // Counts for different bases
  vector<map<string, int>> counts;

  mpileup_line() {
    chr = ref_base = "NA";
    pos = 0;
  }

  static void print_header(int nsample, ostream& out = cout) {
    out << "chr"
        << "\t"
        << "pos"
        << "\t"
        << "ref_base"
        << "\t";
    for (int i = 0; i < nsample; i++) {
      for (int j = 0; j < names.size(); j++) {
        out << names[j];
        if (j < names.size() - 1) {
          out << count_sep;
        }
      }
      if (i < nsample - 1) {
        out << sample_sep;
      }
    }
    out << endl;
  }

  void print_mutation(ostream& out = cout) {
    out << chr << "\t" << pos << "\t" << ref_base << "\t";
    for (int i = 0; i < nsample; i++) {
      map<string, int> m = counts[i];
      for (int j = 0; j < names.size(); j++) {
        out << m[names[j]];
        if (j < names.size() - 1) {
          out << count_sep;
        }
      }
      if (i < nsample - 1) {
        out << sample_sep;
      }
    }
    out << endl;
  }
};

// Parse the pileup string
map<string, int> parse_counts(string& bases, string& qual, int depth) {
  map<string, int> m{
      {"depth", depth}, {"ref", 0},    {"fwd", 0},    {"rev", 0},
      {"A", 0},         {"a", 0},      {"C", 0},      {"c", 0},
      {"G", 0},         {"g", 0},      {"T", 0},      {"t", 0},
      {"N", 0},         {"n", 0},      {"Gap", 0},    {"gap", 0},
      {"Insert", 0},    {"insert", 0}, {"Delete", 0}, {"delete", 0},
  };

  for (int i = 0; i < bases.length(); i++) {
    char base = bases[i];
    string indelsize_string;
    int indelsize_int = 0;
    switch (base) {
      // Match to reference
      case '.':
        m["fwd"] += 1;
        break;
      case ',':
        m["rev"] += 1;
        break;
      case 'a':
        m["a"] += 1;
        break;
      case 'A':
        m["A"] += 1;
        break;
      case 'c':
        m["c"] += 1;
        break;
      case 'C':
        m["C"] += 1;
        break;
      case 'g':
        m["g"] += 1;
        break;
      case 'G':
        m["G"] += 1;
        break;
      case 't':
        m["t"] += 1;
        break;
      case 'T':
        m["T"] += 1;
        break;
      case 'n':
        m["n"] += 1;
        break;
      case 'N':
        m["N"] += 1;
        break;
      // This base is a gap (--reverse-del suport)
      // similar with Deletecount and deletecount, but with some difference
      case '*':
        m["Gap"] += 1;
        break;
      case '#':
        m["gap"] += 1;
        break;
      // Insertion
      case '+':
        i++;
        // 48 is number '0', 57 is number '9'
        while (bases[i] >= 48 && bases[i] <= 57) {
          indelsize_string = indelsize_string + bases[i];
          i = i + 1;
        }
        if (isupper(bases[i])) {
          m["Insert"] += 1;
        } else {
          m["insert"] += 1;
        }
        indelsize_int = str_to_num(indelsize_string);
        i += indelsize_int - 1;
        break;
      // Deletion
      case '-':
        i++;
        // 48 is number '0', 57 is number '9'
        while (bases[i] >= 48 && bases[i] <= 57) {
          indelsize_string = indelsize_string + bases[i];
          i = i + 1;
        }
        if (isupper(bases[i])) {
          m["Delete"] += 1;
        } else {
          m["delete"] += 1;
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
      // Beginning of read segment, Skip
      case '^':
        i = i + 1;
        break;
      default:
        string err = "Unknown ref base: ";
        err += base;
        throw runtime_error(err);
    }
  }

  m["ref"] = m["fwd"] + m["rev"];

  return m;
}

// Set the appropriate count for ref nucleotide
map<string, int> adjust_counts(map<string, int> m, string& ref_base) {
  switch (ref_base[0]) {
    case 'A':
    case 'a':
      m["A"] = m["fwd"];
      m["a"] = m["rev"];
      break;
    case 'C':
    case 'c':
      m["C"] = m["fwd"];
      m["c"] = m["rev"];
      break;
    case 'G':
    case 'g':
      m["G"] = m["fwd"];
      m["g"] = m["rev"];
      break;
    case 'T':
    case 't':
      m["T"] = m["fwd"];
      m["t"] = m["rev"];
      break;
    case 'N':
    case 'n':
      m["N"] = m["fwd"];
      m["n"] = m["rev"];
      break;
    // Deal with -,R,Y,K,M,S,W etc
    default:
      break;
  }
  return m;
}

// Split the line into the required fields and parse
mpileup_line process_mpileup_line(string line) {
  stringstream ss(line);
  mpileup_line ml;

  int ncol = 0;
  int nsample = 0;
  string col;
  while (getline(ss, col, '\t')) {
    if (ncol == 0) {
      // get chrosome ID
      ml.chr = col;
    } else if (ncol == 1) {
      // get pos
      ml.pos = str_to_num(col);
    } else if (ncol == 2) {
      // get ref_base
      ml.ref_base = col;
    } else {
      int depth;
      string bases, quals;
      // get depth
      depth = str_to_num(col);
      // get bases
      getline(ss, bases, '\t');
      // get quals
      getline(ss, quals, '\t');
      ml.counts.push_back(
          adjust_counts(parse_counts(bases, quals, depth), ml.ref_base));
      nsample++;
    }
    ncol++;
  }
  ml.nsample = nsample;
  // throw runtime_error("ERROR!");

  return ml;
}

int main(int argc, char* argv[]) {
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
      usage();
      return 0;
    }
  }
  string line;
  getline(cin, line);

  // print header
  try {
    mpileup_line ml = process_mpileup_line(line);
    mpileup_line::print_header(ml.nsample);
  } catch (const std::runtime_error& e) {
    cerr << e.what() << endl;
    cerr << "\nError parsing line " << line;
  }

  // parse and print each line
  while (cin) {
    try {
      mpileup_line ml = process_mpileup_line(line);
      ml.print_mutation();
    } catch (const std::runtime_error& e) {
      cerr << e.what() << endl;
      cerr << "\nError parsing line " << line;
      break;
    }
    getline(cin, line);
  }
}
