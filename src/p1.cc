#include <iostream>
#include <libgen.h>
#include <string>
#include <unordered_map>
#include <map>
#include <cmath>
#include <unistd.h>
#include <vector>
#include <chrono>
#include <tuple>
#include <atomic>
#include <mutex>


using namespace std;

struct degree {
    int in;
    int out;

    degree(int in, int out) : in(in), out(out) {}
  };

struct edge_node {
  int e;
  std::atomic<bool> used;

  edge_node(int e) : e(e) {
    used.store(false);
  }
};

  unordered_map<int, unordered_map<int, vector<edge_node *>>> node_map = {};
  unordered_map<int, degree *> degree_map = {};
  int map_size = 0;
  int count_idx = 0;

/// arg_t represents the command-line arguments to the server
struct arg_t {
  string filename;      // fasta file name
  string output;        //output file name
  int overlap;          //number of chars of overlap of each read

  /// Construct an arg_t from the command-line arguments to the program
  ///
  /// @param argc The number of command-line arguments passed to the program
  /// @param argv The list of command-line arguments
  ///
  /// @throw An intmd5eger exception (1) if an invalid argument is given, or if
  ///        `-h` is passed in
  arg_t(int argc, char **argv) {
    long opt;
    while ((opt = getopt(argc, argv, "f:o:m:")) != -1) {
      switch (opt) {
      case 'f':
        filename = string(optarg);
        break;
      case 'o':
        output = string(optarg);
        break;
      case 'm':
        overlap = atoi(optarg);
        break;
      default: // on any error, print a help message.  This case subsumes `-h`
        throw 1;
        return;
      }
    }
  }

  /// Display a help message to explain how the command-line parameters for this
  /// program work
  ///
  /// @progname The name of the program
  static void usage(char *progname) {
    cout << basename(progname) << ": company user directory server\n"
         << "  -f [string] File of the fasta file\n"
         << "  -o [string] File of the output file\n"
         << "  -m [int]    # of overlap for each read\n"
         << "  -h          Print help (this message)\n";
  }
};

int getPrefixIndex(const char* read) {
  int sum = 0;
  for (int idx=0; idx<15; ++idx) {
    switch(read[idx]) {
      case 'A':
        sum += 0*pow(4, idx);
        break;
      case 'C':
        sum += 1*pow(4, idx);
        break;
      case 'G':
        sum += 2*pow(4, idx);
        break;
      case 'T':
        sum += 3*pow(4, idx);
        break;
      default:
        return -1;
    }
  }
  return sum;
}

int getSufixIndex(const char* read) {
  int sum = 0;
  for (int idx=15; idx<30; ++idx) {
    switch(read[idx]) {
      case 'A':
        sum += 0*pow(4, idx-15);
        break;
      case 'C':
        sum += 1*pow(4, idx-15);
        break;
      case 'G':
        sum += 2*pow(4, idx-15);
        break;
      case 'T':
        sum += 3*pow(4, idx-15);
        break;
      default:
        return -1;
    }
  }
  return sum;
}


string getRead(int val) {
  string read;
  int div, i;
  for (int idx = 14; idx >= 0; idx--) {
    div = pow(4, idx);
    i = val / div;
    val -= div * i;
    switch (i) {
      case 0:
        read.append("A");
        break;
      case 1:
        read.append("C");
        break;
      case 2:
        read.append("G");
        break;
      case 3:
        read.append("T");
        break;
      default:
        return "ERR";
    }
  }
  reverse(read.begin(), read.end());
  return read;
}




void insert(int pre, int suf, int edge) {
  if (!node_map.contains(pre)){
    node_map.insert({pre, unordered_map<int, vector<edge_node *>>()});
    node_map.at(pre).insert({suf, vector<edge_node *>()});
  } else {
    if(!node_map.at(pre).contains(suf)) {
      node_map.at(pre).insert({suf, vector<edge_node *>()});
    }
  }
  node_map.at(pre).at(suf).push_back(new edge_node(edge));

  if (!degree_map.contains(pre)) {
    degree_map.insert({pre, new degree(0,1)});
  } else{
    degree_map.at(pre)->out += 1;
  }
  if (!degree_map.contains(suf)) {
    degree_map.insert({suf, new degree(1,0)});
  } else{
    degree_map.at(suf)->in += 1;
  }
}

string flatten_vector(int *v) {
  string s;
  for (auto i = 0; i < 50; i++) {
    if (v[i] == 0) break;
    s.append(getRead(v[i]));
  }
  return s;
}

void faps(int ***paths, bool **used, std::mutex **mtx) {
  int idx = 0;
  while(true) {
    //initialize path
    int pre = -1;
    bool end = false;
    vector<int> path;

    // find node where out > in
    for (auto i = degree_map.begin(); i != degree_map.end(); i++) {
      bool expected = false;
      bool desired = true;
      if (i->second->out > i->second->in){
        pre = i->first;
        path.push_back(pre);
        break;
      }
    }

    // if no start node found return empty vector
    if (pre == -1) break;


    while(degree_map.at(pre)->out) {
      end = false;
      // iterate though 2nd map
      for (auto i = node_map.at(pre).begin(); i != node_map.at(pre).end(); i++) {
          //iterate though edges
          for (auto j = i->second.begin(); j != i->second.end(); j++) {
            bool expected = false;
            bool desired = true;
            if ((*j)->used.compare_exchange_strong(expected, desired)) {
              degree_map.at(pre)->out -= 1;
              degree_map.at(i->first)->in -= 1;
              pre = i->first;
              path.push_back(pre);
              end = true;
              break;
            }
          }
          if (end) break;
      }
    }

    bool overlap = false;
    // need to check if path overlaps with anything
    // cout << sizeof(*paths) << endl;
    for (auto i = 0; i < map_size; i++) {
      for (auto j = 0; j < 50; j++) {
        for (auto k = path.begin(); k != path.end(); k++) {
          if   ((*paths)[i][j] == *k) {
            vector<int> new_path;
            // conbine the longest part of both paths
            // (*mtx)[i].lock();
            if (j - i >= k - path.begin()) {
              new_path.insert(new_path.end(), 0, j);
            } else {
              new_path.insert(new_path.end(), path.begin(), k);
            }
            if (sizeof(paths[i]) / sizeof(paths[i][0]) - (j) >= path.size() - (k - path.begin())) {
              new_path.insert(new_path.end(), j, sizeof(paths[i]) / sizeof(paths[i][0]));
            } else {
              new_path.insert(new_path.end(), k, path.end());
            }
            // remove the old path
            (*used)[i] = true;
            //turn new_paths into an array
            for (auto l = 0; l < new_path.size(); l++) {
              (*paths)[idx][l] = new_path.at(l);
            }
            idx++;
            // (*mtx)[i].unlock();
            overlap = true;
            break;
          }
        }
        if (overlap) break;
      }
    }
    if (!overlap) {
      //add new_path to paths
      for (auto l = 0; l < path.size(); l++) {
        (*paths)[idx][l] = path.at(l);
      }
      idx++;
    }
  }

}

unordered_map<string, int> path(){
  
  // vector<vector<int>> *paths = new vector<vector<int>>(map_size);
  // malloc 2d array pahts with size 50 by map_size
  int **paths = (int **)calloc(map_size, sizeof(int *));
  for (auto i = 0; i < map_size; i++) {
    paths[i] = (int *)calloc(50, sizeof(int));
  }

  // vector<std::mutex> *locks = new vector<std::mutex>(map_size);
  std::mutex *locks = new std::mutex[map_size];
  // vector<bool> *used = new vector<bool>(map_size);
  bool *used = new bool[map_size];
  unordered_map<string, int> paths_map;

  //deploy threads that run faps
  // vector<thread> threads = vector<thread>(1);
  // for (int i = 0; i < 8; i++) {
  //   threads[i] = thread([&]() {
  //     faps(&paths, &used, &locks);
  //   });
  // }
  // for (int i = 0; i < 8; i++) {
  //   threads[i].join();
  // }

  faps(&paths, &used, &locks);
  
  for (auto i = 0; i < map_size; i++) {
    string s = flatten_vector(paths[i]);
    if (!paths_map.contains(s))
      paths_map.insert({s, 1});
  }
  return paths_map;
}


unordered_map<string, int> find_overlap(unordered_map<string, int> paths) {
  
  //initialize result and used contigs
  unordered_map<string, int> res;
  bool used = false;
  unordered_map<string, int> used_contigs;

  //iterate through all contigs
  for (auto i = paths.begin(); i != paths.end(); i++) {
    
    used = false;
    
    //if the contig has already been used skip it
    if (used_contigs.contains(i->first)) continue;
    
    //iterate though all the contigs from the current one to end
    for (auto j = i; j != paths.end(); j++) {

      if (i->first == j->first) continue;

      //if the contig has already been used skip it
      if (used_contigs.contains(j->first)) continue;

      // find the bigger and smaller contig
      string big = i->first.size() >= j->first.size() ? i->first : j->first;
      string small = i->first.size() > j->first.size() ? j->first : i->first;

      // if small is a substring of big remove small from the result
      if (big.find(small) != string::npos) {
        res.insert({big, 1});
        used_contigs.insert({j->first, 1});
        used_contigs.insert({i->first, 1});
        used = true;
        break;
      }
      
      //loop through the smaller contig and check for overlap
      for (int k = 15; k < small.size(); k++) {
        // starting at 15 overlap and moving in from the left
        if (big.substr(0, k) == small.substr(small.size() - k, k) 
            && !res.contains(small + big.substr(k))
            ) {
          res.insert({small + big.substr(k), 1});
          used_contigs.insert({j->first, 1});
          used_contigs.insert({i->first, 1});
          used = true;
          break;
        }
        // starting at 15 overlap and moving in from the right
        else if (big.substr(big.size() - k, k) == small.substr(0, k)
         && !res.contains(big + small.substr(k))
         ) {
          res.insert({big + small.substr(k), 1});
          used_contigs.insert({j->first, 1});
          used_contigs.insert({i->first, 1});
          used = true;
          break;
        }
      }
    }
    // if no overlap with anyother contig add it to the result
    if (!used) {
      res.insert({i->first, 1});
      used_contigs.insert({i->first, 1});
    }
  }
  return res;
}



void read_file(const std::string &filename) {
  FILE* fp = fopen(filename.c_str(), "r");
  if (fp == NULL)
      exit(EXIT_FAILURE);

  char* line = NULL;
  size_t len = 0;
  int i = 0;
  while ((getline(&line, &len, fp)) != -1) {
      if (line[0] == '>') continue;
      int pre = getPrefixIndex(line);
      int suf = getSufixIndex(line);

      insert(pre, suf, i++);
      map_size++;

  }
  fclose(fp);
  if (line)
      free(line);
}


void printer() {
  cout << "---- MAP ----" << endl;
  for (auto i = node_map.begin(); i != node_map.end(); i++) {
    cout << getRead(i->first) <<  " - [" << degree_map.at(i->first)->in << ", " << degree_map.at(i->first)->out << "]" <<endl;
    for (auto j = i->second.begin(); j != i->second.end(); j++){
      cout << "  " << getRead(j->first) << " - [" << j->second.size() << "]" << endl;
    }
  }
}


int main(int argc, char **argv) {
  // Parse the command-line arguments
  //
  // NB: It would be better not to put the arg_t on the heap, but then we'd need
  //     an extra level of nesting for the body of the rest of this function.
  arg_t *args;
  try {
    args = new arg_t(argc, argv);
  } catch (int i) {
    arg_t::usage(argv[0]);
    return 1;
  }


  auto start = std::chrono::high_resolution_clock::now();
  //read and insert nodes into node_map
  read_file(args->filename);
  auto end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Read File: " << elapsed.count() << std::endl;

  

  start = std::chrono::high_resolution_clock::now();
  unordered_map<string, int> paths = path();
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "Path: " << elapsed.count() << std::endl;

  int old_size = paths.size();
  int iters = 0;
  int tolerance = 1;

  start = std::chrono::high_resolution_clock::now();
  while (true) {
    paths = find_overlap(paths);
    if (paths.size() == old_size) 
      tolerance--;
    if (tolerance == 0)
      break;
    old_size = paths.size();
  }
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "Overlap: " << elapsed.count() << std::endl;
  
  
  int counter = 1;
  int total_size = 0;
  for (auto i = paths.begin(); i != paths.end(); i++) {
    cout << ">contig" << counter++ << "|size" << i->first.size() << endl;
    cout << i->first << endl;
    total_size += i->first.size();
  }
}
