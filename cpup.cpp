/*
 * cpup.cpp
 * Copyright (C) 2021 Ye Chang <yech1990@gmail.com>
 *
 * Distributed under terms of the MIT license.
 */

#include <algorithm>
#include <cctype>
#include <cstring>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using namespace std;

void usage() {
  std::cerr
      << "Usage: " << endl
      << "  samtools mpileup -d 0 -Q 0 --reverse-del -l <.bed> -f <.fa> <.bam>"
         " | cpup"
      << endl;
}

// global variables
vector<string> names = {
    "depth",
    //"ref",
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

vector<string> names_no_strand = {
    "depth",
    //"ref",
    "A",
    "C",
    "G",
    "T",
    "N",
    "Gap",
    "Insert",
    "Delete"};

string count_sep = ",";
string indel_sep = "|";
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
  vector<map<string, int>> counts, istats, dstats;

  mpileup_line() {
    chr = ref_base = "NA";
    pos = 0;
  }

  static void print_header(
      int nsample,
      ostream& out = cout,
      bool stat_indel = false,
      bool hide_strand = false) {
    out << "chr"
        << "\t"
        << "pos"
        << "\t"
        << "ref_base"
        << "\t";
    for (int i = 0; i < nsample; i++) {
      if (!hide_strand) {
        for (int j = 0; j < names.size(); j++) {
          out << names[j];
          if (j < names.size() - 1) {
            out << count_sep;
          }
        }
      } else {
        for (int j = 0; j < names_no_strand.size(); j++) {
          out << names_no_strand[j];
          if (j < names_no_strand.size() - 1) {
            out << count_sep;
          }
        }
      }
      if (stat_indel) {
        out << count_sep << "istat";
        out << count_sep << "dstat";
      }
      if (i < nsample - 1) {
        out << sample_sep;
      }
    }
    out << endl;
  }

  void print_counter(
      ostream& out = cout,
      bool stat_indel = false,
      bool hide_strand = false) {
    out << chr << "\t" << pos << "\t" << ref_base << "\t";
    for (int i = 0; i < nsample; i++) {
      // count
      map<string, int> m = counts[i];
      if (!hide_strand) {
        for (int j = 0; j < names.size(); j++) {
          out << m[names[j]];
          if (j < names.size() - 1) {
            out << count_sep;
          }
        }
      } else {
        for (int j = 0; j < names_no_strand.size(); j++) {
          string name1 = names_no_strand[j];
          string name2 = name1;
          std::transform(
              name2.begin(),
              name2.end(),
              name2.begin(),
              [](unsigned char c) { return ::tolower(c); });
          if (name2 == "depth") {
            out << m[name1];
          } else {
            out << m[name1] + m[name2];
          }
          if (j < names_no_strand.size() - 1) {
            out << count_sep;
          }
        }
      }
      if (stat_indel) {
        // istat
        map<string, int> istat = istats[i];
        out << count_sep;
        if (hide_strand) {
          map<string, int> istat_no_strand;
          for (auto iter = istat.begin(); iter != istat.end(); ++iter) {
            string motif = iter->first;
            std::transform(
                motif.begin(),
                motif.end(),
                motif.begin(),
                [](unsigned char c) { return ::toupper(c); });
            istat_no_strand[motif] += iter->second;
          }
          for (auto iter = istat_no_strand.begin();
               iter != istat_no_strand.end();
               ++iter) {
            if (std::next(iter) != istat_no_strand.end()) {
              out << iter->first << ':' << iter->second << indel_sep;
            } else {
              out << iter->first << ':' << iter->second;
            }
          }
        } else {
          for (auto iter = istat.begin(); iter != istat.end(); ++iter) {
            if (std::next(iter) != istat.end()) {
              out << iter->first << ':' << iter->second << indel_sep;
            } else {
              out << iter->first << ':' << iter->second;
            }
          }
        }
        // dstat
        map<string, int> dstat = dstats[i];
        out << count_sep;
        if (hide_strand) {
          map<string, int> dstat_no_strand;
          for (auto iter = dstat.begin(); iter != dstat.end(); ++iter) {
            string motif = iter->first;
            std::transform(
                motif.begin(),
                motif.end(),
                motif.begin(),
                [](unsigned char c) { return ::toupper(c); });
            dstat_no_strand[motif] += iter->second;
          }
          for (auto iter = dstat_no_strand.begin();
               iter != dstat_no_strand.end();
               ++iter) {
            if (std::next(iter) != dstat_no_strand.end()) {
              out << iter->first << ':' << iter->second << indel_sep;
            } else {
              out << iter->first << ':' << iter->second;
            }
          }
        } else {
          for (auto iter = dstat.begin(); iter != dstat.end(); ++iter) {
            if (std::next(iter) != dstat.end()) {
              out << iter->first << ':' << iter->second << indel_sep;
            } else {
              out << iter->first << ':' << iter->second;
            }
          }
        }
      }
      // end
      if (i < nsample - 1) {
        out << sample_sep;
      }
    }
    out << endl;
  }
};

// Parse the pileup string
tuple<map<string, int>, map<string, int>, map<string, int>>
parse_counts(string& bases, string& qual, int depth) {
  map<string, int> m{
      {"depth", depth}, {"coverage", 0}, {"ref", 0},    {"mut", 0},
      {"Gap+gap", 0},   {"fwd", 0},      {"rev", 0},    {"A", 0},
      {"a", 0},         {"C", 0},        {"c", 0},      {"G", 0},
      {"g", 0},         {"T", 0},        {"t", 0},      {"N", 0},
      {"n", 0},         {"Gap", 0},      {"gap", 0},    {"Insert", 0},
      {"insert", 0},    {"Delete", 0},   {"delete", 0},
  };
  map<string, int> istat;
  map<string, int> dstat;

  // check if site is a empty (depth == 0)
  if (bases == "*") {
    return make_tuple(m, istat, dstat);
  }
  for (int i = 0; i < bases.length(); i++) {
    char base = bases[i];
    string indelsize_string;
    string indelseq;
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
        // stat_indel
        indelseq = bases.substr(i, indelsize_int);
        istat[indelseq]++;
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
        // stat_indel
        indelseq = bases.substr(i, indelsize_int);
        dstat[indelseq]++;
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
  m["mut"] =
      m["A"] + m["a"] + m["C"] + m["c"] + m["G"] + m["g"] + m["T"] + m["t"];
  m["Gap+gap"] = m["Gap"] + m["gap"];
  m["coverage"] = m["ref"] + m["mut"];

  return make_tuple(m, istat, dstat);
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

      map<string, int> m, istat, dstat;
      tie(m, istat, dstat) = parse_counts(bases, quals, depth);
      ml.counts.push_back(adjust_counts(m, ml.ref_base));
      ml.istats.push_back(istat);
      ml.dstats.push_back(dstat);
      nsample++;
    }
    ncol++;
  }
  ml.nsample = nsample;
  return ml;
}

vector<string> split_string(const string& i_str, const string& i_delim) {
  vector<string> result;

  size_t found = i_str.find(i_delim);
  size_t startIndex = 0;

  while (found != string::npos) {
    result.push_back(string(i_str.begin() + startIndex, i_str.begin() + found));
    startIndex = found + i_delim.size();
    found = i_str.find(i_delim, startIndex);
  }
  if (startIndex != i_str.size())
    result.push_back(string(i_str.begin() + startIndex, i_str.end()));
  return result;
}

int main(int argc, char* argv[]) {
  bool hide_header = false;
  bool hide_strand = false;
  bool stat_indel = false;
  map<string, int> min_cutoffs;
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
      usage();
      return 0;
    } else if (!strcmp(argv[i], "-H") || !strcmp(argv[i], "--header")) {
      hide_header = true;
    } else if (!strcmp(argv[i], "-S") || !strcmp(argv[i], "--strandless")) {
      hide_strand = true;
    } else if (!strcmp(argv[i], "-i") || !strcmp(argv[i], "--indel")) {
      stat_indel = true;
    } else if (!strcmp(argv[i], "-f") || !strcmp(argv[i], "--filter")) {
      if (i + 1 != argc) {
        vector<string> filters = split_string(argv[i + 1], ",");
        for (auto& s : filters) {
          vector<string> filter = split_string(s, ":");
          // mut, ref, coverage, ...
          // Gap + gap is a special key
          string filter_name = filter[0];
          int min_cutoff = std::stoi(filter[1]);
          min_cutoffs[filter_name] = min_cutoff;
        }
      }
      i++;
    }
  }

  string line;
  getline(cin, line);

  // print header
  if (!hide_header) {
    try {
      mpileup_line ml = process_mpileup_line(line);
      mpileup_line::print_header(ml.nsample, cout, stat_indel, hide_strand);
    } catch (const std::runtime_error& e) {
      cerr << e.what() << endl;
      cerr << "\nError parsing line " << line;
    }
  }

  // parse and print each line
  while (cin) {
    try {
      mpileup_line ml = process_mpileup_line(line);
      map<string, bool> filters_results;
      bool is_passed = true;
      for (auto iter = min_cutoffs.begin(); iter != min_cutoffs.end(); ++iter) {
        string filter_name = iter->first;
        int min_cutoff = iter->second;
        filters_results[filter_name] = false;
        for (int i = 0; i < ml.nsample; i++) {
          if (ml.counts[i][filter_name] >= min_cutoff) {
            filters_results[filter_name] = true;
            break;
          }
        }
        is_passed = is_passed && filters_results[filter_name];
      }
      if (is_passed) {
        ml.print_counter(cout, stat_indel, hide_strand);
      }
    } catch (const std::runtime_error& e) {
      cerr << e.what() << endl;
      cerr << "\nError parsing line " << line;
      break;
    }
    getline(cin, line);
  }
}
