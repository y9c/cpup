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

static const unsigned char basemap[256] = {
    0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  14,
    15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,
    30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,
    45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,
    60,  61,  62,  63,  64,  'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J',
    'M', 'L', 'K', 'N', 'O', 'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R',
    'Z', 91,  92,  93,  94,  95,  96,  't', 'v', 'g', 'h', 'e', 'f', 'c', 'd',
    'i', 'j', 'm', 'l', 'k', 'n', 'o', 'p', 'q', 'y', 's', 'a', 'a', 'b', 'w',
    'x', 'r', 'z', 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134,
    135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149,
    150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164,
    165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179,
    180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194,
    195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209,
    210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224,
    225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
    240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254,
    255};

void usage() {
  std::cerr
      << "Usage: " << endl
      << "  samtools mpileup -d 0 -Q 0 --reverse-del -l <.bed> -f <.fa> <.bam>"
         " | cpup"
      << endl;
}

// global variables
vector<string> names = {"a", "c", "g", "t", "n", "gap", "insert", "delete"};
vector<string> names_upper =
    {"A", "C", "G", "T", "N", "Gap", "Insert", "Delete"};

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
  vector<map<string, int>> Counts, counts, Istats, istats, Dstats, dstats;
  vector<int> depths;

  mpileup_line() {
    chr = ref_base = "NA";
    pos = 0;
  }

  static void print_header(
      int nsample,
      ostream& out = cout,
      bool stat_indel = false,
      bool hide_strand = false,
      bool by_strand = false) {
    out << "chr" << sample_sep << "pos" << sample_sep << "ref_base";
    if (by_strand) {
      out << sample_sep << "strand";
    }
    for (int i = 0; i < nsample; i++) {
      out << sample_sep << "depth";
      if (!hide_strand && !by_strand) {
        for (int j = 0; j < names_upper.size(); j++) {
          out << count_sep << names_upper[j];
        }
        if (stat_indel) {
          out << count_sep << "Istat";
          out << count_sep << "Dstat";
        }
      }
      for (int j = 0; j < names.size(); j++) {
        out << count_sep << names[j];
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
      bool hide_strand = false,
      bool by_strand = false,
      char strands = '*') {
    // forward strand
    if (by_strand) {
      if (strands == '*' || strands == '+') {
        out << chr << sample_sep << pos << sample_sep << ref_base << sample_sep
            << "+";
        for (int i = 0; i < nsample; i++) {
          map<string, int> M = Counts[i];
          out << sample_sep << depths[i];
          for (int j = 0; j < names.size(); j++) {
            out << count_sep << M[names[j]];
          }
          // indel stat
          if (stat_indel) {
            vector<map<string, int>> indel_stats = {Istats[i], Dstats[i]};
            for (auto ids : indel_stats) {
              out << count_sep;
              for (auto iter = ids.begin(); iter != ids.end(); ++iter) {
                if (std::next(iter) != ids.end()) {
                  out << iter->first << ':' << iter->second << indel_sep;
                } else {
                  out << iter->first << ':' << iter->second;
                }
              }
            }
          }
        }
        out << endl;
      }
      // reverse strand
      if (strands == '*' || strands == '-') {
        out << chr << sample_sep << pos << sample_sep << basemap[ref_base[0]]
            << sample_sep << "-";
        for (int i = 0; i < nsample; i++) {
          map<string, int> m = counts[i];
          out << sample_sep << depths[i];
          for (int j = 0; j < names.size(); j++) {
            out << count_sep << m[names[j]];
          }
          // indel stat
          if (stat_indel) {
            vector<map<string, int>> indel_stats = {istats[i], dstats[i]};
            for (auto ids : indel_stats) {
              out << count_sep;
              for (auto iter = ids.begin(); iter != ids.end(); ++iter) {
                if (std::next(iter) != ids.end()) {
                  out << iter->first << ':' << iter->second << indel_sep;
                } else {
                  out << iter->first << ':' << iter->second;
                }
              }
            }
          }
        }
        out << endl;
      }
    }  // end by_strand

    else if (hide_strand) {
      out << chr << sample_sep << pos << sample_sep << ref_base;
      for (int i = 0; i < nsample; i++) {
        map<string, int> M = Counts[i];
        map<string, int> m = counts[i];
        out << sample_sep << depths[i];
        for (int j = 0; j < names.size(); j++) {
          out << count_sep << M[names[j]] + m[names[j]];
        }
        // indel stat
        if (stat_indel) {
          vector<map<string, int>> indel_stats = {
              Istats[i],
              istats[i],
              Dstats[i],
              dstats[i]};
          for (auto ids : indel_stats) {
            map<string, int> ids_no_strand;
            for (auto iter = ids.begin(); iter != ids.end(); ++iter) {
              string motif = iter->first;
              std::transform(
                  motif.begin(),
                  motif.end(),
                  motif.begin(),
                  [](unsigned char c) { return ::toupper(c); });
              ids_no_strand[motif] += iter->second;
            }
            out << count_sep;
            for (auto iter = ids_no_strand.begin(); iter != ids_no_strand.end();
                 ++iter) {
              if (std::next(iter) != ids_no_strand.end()) {
                out << iter->first << ':' << iter->second << indel_sep;
              } else {
                out << iter->first << ':' << iter->second;
              }
            }
          }
        }
      }
      out << endl;
    }  // end by_strand

    else {
      out << chr << sample_sep << pos << sample_sep << ref_base;
      for (int i = 0; i < nsample; i++) {
        out << sample_sep << depths[i];
        map<string, int> M = Counts[i];
        for (int j = 0; j < names.size(); j++) {
          out << count_sep << M[names[j]];
        }
        // indel stat
        if (stat_indel) {
          vector<map<string, int>> indel_stats = {Istats[i], Dstats[i]};
          for (auto ids : indel_stats) {
            out << count_sep;
            for (auto iter = ids.begin(); iter != ids.end(); ++iter) {
              if (std::next(iter) != ids.end()) {
                out << iter->first << ':' << iter->second << indel_sep;
              } else {
                out << iter->first << ':' << iter->second;
              }
            }
          }
        }
        map<string, int> m = counts[i];
        for (int j = 0; j < names.size(); j++) {
          out << count_sep << m[names[j]];
        }
        // indel stat
        if (stat_indel) {
          vector<map<string, int>> indel_stats = {istats[i], dstats[i]};
          for (auto ids : indel_stats) {
            out << count_sep;
            for (auto iter = ids.begin(); iter != ids.end(); ++iter) {
              if (std::next(iter) != ids.end()) {
                out << iter->first << ':' << iter->second << indel_sep;
              } else {
                out << iter->first << ':' << iter->second;
              }
            }
          }
        }
      }
      out << endl;
    }  // end
  }
};

// Parse the pileup string
tuple<
    map<string, int>,
    map<string, int>,
    map<string, int>,
    map<string, int>,
    map<string, int>,
    map<string, int>>
parse_counts(string& bases, string& qual) {
  // forward strand
  map<string, int> M{
      {"coverage", 0},
      {"ref", 0},
      {"mut", 0},
      {"gap", 0},
      {"a", 0},
      {"c", 0},
      {"g", 0},
      {"t", 0},
      {"n", 0},
      {"gap", 0},
      {"insert", 0},
      {"delete", 0},
  };
  // reverse strand
  map<string, int> m{
      {"coverage", 0},
      {"ref", 0},
      {"mut", 0},
      {"gap", 0},
      {"a", 0},
      {"c", 0},
      {"g", 0},
      {"t", 0},
      {"n", 0},
      {"gap", 0},
      {"insert", 0},
      {"delete", 0},
  };
  map<string, int> istat;
  map<string, int> dstat;
  map<string, int> Istat;
  map<string, int> Dstat;

  // check if site is a empty (depth == 0)
  if (bases == "*") {
    return make_tuple(M, m, Istat, istat, Dstat, dstat);
  }
  for (int i = 0; i < bases.length(); i++) {
    char base = bases[i];
    string indelsize_string;
    string indelseq;
    int indelsize_int = 0;
    switch (base) {
      // Match to reference
      case '.':
        M["ref"] += 1;
        break;
      case ',':
        m["ref"] += 1;
        break;
      case 'A':
        M["a"] += 1;
        break;
      case 'a':
        m["a"] += 1;
        break;
      case 'C':
        M["c"] += 1;
        break;
      case 'c':
        m["c"] += 1;
        break;
      case 'G':
        M["g"] += 1;
        break;
      case 'g':
        m["g"] += 1;
        break;
      case 'T':
        M["t"] += 1;
        break;
      case 't':
        m["T"] += 1;
        break;
      case 'N':
        M["n"] += 1;
        break;
      case 'n':
        m["N"] += 1;
        break;
      // This base is a gap (--reverse-del suport)
      // similar with Deletecount and deletecount, but with some difference
      case '*':
        M["gap"] += 1;
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
        indelsize_int = str_to_num(indelsize_string);
        indelseq = bases.substr(i, indelsize_int);
        if (isupper(bases[i])) {
          M["insert"] += 1;
          Istat[indelseq]++;
        } else {
          m["insert"] += 1;
          istat[indelseq]++;
        }
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
        indelsize_int = str_to_num(indelsize_string);
        indelseq = bases.substr(i, indelsize_int);
        if (isupper(bases[i])) {
          M["delete"] += 1;
          Dstat[indelseq]++;
        } else {
          m["delete"] += 1;
          dstat[indelseq]++;
        }
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

  m["mut"] = m["a"] + m["c"] + m["g"] + m["t"];
  m["coverage"] = m["ref"] + m["mut"];

  return make_tuple(M, m, Istat, istat, Dstat, dstat);
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

      map<string, int> M, m, Istat, istat, Dstat, dstat;
      tie(M, m, Istat, istat, Dstat, dstat) = parse_counts(bases, quals);
      ml.depths.push_back(depth);
      ml.Counts.push_back(adjust_counts(M, ml.ref_base));
      ml.counts.push_back(adjust_counts(m, ml.ref_base));
      ml.Istats.push_back(Istat);
      ml.istats.push_back(istat);
      ml.Dstats.push_back(Dstat);
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
  bool by_strand = false;
  bool stat_indel = false;
  map<string, int> min_cutoffs;
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
      usage();
      return 0;
    } else if (!strcmp(argv[i], "-H") || !strcmp(argv[i], "--headerless")) {
      hide_header = true;
    } else if (!strcmp(argv[i], "-S") || !strcmp(argv[i], "--strandless")) {
      hide_strand = true;
    } else if (!strcmp(argv[i], "-s") || !strcmp(argv[i], "--by-strand")) {
      if (hide_strand) {
        cerr << "\n"
                "Can not use the `--by_strand (-s)` parameter together with "
                "the `--strandless (-S)` parameter"
             << endl;
        return 1;
      }
      by_strand = true;
    } else if (!strcmp(argv[i], "-i") || !strcmp(argv[i], "--indel")) {
      stat_indel = true;
    } else if (!strcmp(argv[i], "-f") || !strcmp(argv[i], "--filter")) {
      if (i + 1 != argc) {
        vector<string> filters = split_string(argv[i + 1], ",");
        for (auto& s : filters) {
          vector<string> filter = split_string(s, ":");
          // mut, ref, coverage, gap...
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
      mpileup_line::print_header(
          ml.nsample,
          cout,
          stat_indel,
          hide_strand,
          by_strand);
    } catch (const std::runtime_error& e) {
      cerr << e.what() << endl;
      cerr << "\nError parsing line " << line;
    }
  }

  // parse and print each line
  while (cin) {
    try {
      mpileup_line ml = process_mpileup_line(line);
      if (by_strand) {
        vector<bool> is_passed = {true, true};
        for (auto iter = min_cutoffs.begin(); iter != min_cutoffs.end();
             ++iter) {
          bool filters_result_fwd = false;
          bool filters_result_rev = false;
          string filter_name = iter->first;
          int min_cutoff = iter->second;
          for (int i = 0; i < ml.nsample; i++) {
            if (ml.Counts[i][filter_name] >= min_cutoff) {
              filters_result_fwd = true;
              break;
            }
          }
          for (int i = 0; i < ml.nsample; i++) {
            if (ml.counts[i][filter_name] >= min_cutoff) {
              filters_result_rev = true;
              break;
            }
          }
          is_passed[0] = is_passed[0] && filters_result_fwd;
          is_passed[1] = is_passed[1] && filters_result_rev;
        }

        if (is_passed[0] and is_passed[1]) {
          ml.print_counter(cout, stat_indel, hide_strand, by_strand, '*');
        } else if (is_passed[0] and !is_passed[1]) {
          ml.print_counter(cout, stat_indel, hide_strand, by_strand, '+');
        } else if (!is_passed[0] and is_passed[1]) {
          ml.print_counter(cout, stat_indel, hide_strand, by_strand, '-');
        }
      } else {
        bool are_passed = true;
        for (auto iter = min_cutoffs.begin(); iter != min_cutoffs.end();
             ++iter) {
          bool filters_result = false;
          string filter_name = iter->first;
          int min_cutoff = iter->second;
          for (int i = 0; i < ml.nsample; i++) {
            if (ml.Counts[i][filter_name] + ml.counts[i][filter_name] >=
                min_cutoff) {
              filters_result = true;
              break;
            }
          }
          are_passed = are_passed && filters_result;
        }
        if (are_passed) {
          ml.print_counter(cout, stat_indel, hide_strand, by_strand);
        }
      }
    } catch (const std::runtime_error& e) {
      cerr << e.what() << endl;
      cerr << "\nError parsing line " << line;
      break;
    }
    getline(cin, line);
  }
}
