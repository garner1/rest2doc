// This is my first C++ program

// g++ -std=c++11 first.cpp -o first1 
// ./first1 reads.txt ~/Work/dataset/fastq2cloud/test1K/MCmodel/

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <stdlib.h>
#include <map>
#include <vector>
#include <iterator>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <functional>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <math.h>
#include <boost/range/adaptor/sliced.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/assign.hpp>

using namespace std;

using std::cout;
using std::endl;
using std::ifstream;

using namespace boost::accumulators;

float mean_func(vector<float> v, int start, int end)
{
  double sum = accumulate(v.begin()+start,v.begin()+end,0.0);
  double mean = sum / ((end-start)*1.0);

  return mean;
}

float std_func(vector<float> v, int start, int end)
{
  double sum = accumulate(v.begin()+start,v.begin()+end,0.0);
  double mean = sum / ((end-start)*1.0);
  double sq_sum = inner_product(v.begin()+start,v.begin()+end, v.begin()+start, 0.0);
  double stdev = sqrt(sq_sum / ((end-start)*1.0) - mean * mean);

  return stdev;
}

vector<int> stop_finder( vector<float> y, int dim, int lag, float threshold, float influence )
{
  vector<int> signals(y.size());      // init the stop-words vector
  vector<float> filteredY(y.size());  // init filtering vector
  vector<float> avgFilter(y.size());  // init average vector
  vector<float> stdFilter(y.size());  // init std vector

  filteredY = y;			   // Init filteredY
  avgFilter[dim+lag] = mean_func(y,dim,dim+lag); // Init first value
  stdFilter[dim+lag] = std_func(y,dim,dim+lag); // Init first value

  for ( unsigned i = dim+lag+1; i < y.size(); i++ ) // loop over read
    {
      if ( ( y[i] - avgFilter[i-1] ) > threshold*stdFilter[i-1] && y[i-1] < y[i] )
	{
	  signals[i] = 1;
	  filteredY[i] = influence * y[i] + (1 - influence) * filteredY[i-1]; // update local data points
	}
      else                                            // if not stop word
	{
	  signals[i] = 0;                               
	  filteredY[i] = y[i];
	}
      avgFilter[i] = mean_func(filteredY,i-lag,i); // update the local mean
      stdFilter[i] = std_func(filteredY,i-lag,i);  // update the local std
    }
  return signals;
}

string revcomplement(string seq)
{
    auto lambda = [](const char c) {
        switch (c) {
        case 'A':
            return 'T';
        case 'G':
            return 'C';
        case 'C':
            return 'G';
        case 'T':
            return 'A';
        case 'N':
            return 'N';
        default:
            throw domain_error("Invalid nucleotide.");
        }
    };

    transform(seq.cbegin(), seq.cend(), seq.begin(), lambda);
    string reversed(seq.rbegin(), seq.rend());
    return reversed;
}

int main ( int argc, char* argv[] )
{ 
  string path;
  ifstream model6;
  ifstream model11;
  ifstream model12;
  map<string, float> m6,m11,m12;
  map<string, float>::iterator p6,p11,p12;
  string read,kmer;
  float value;

  path = string(argv[2]);
  model6.open( path ); // the MC model
  if (!model6.good()) 
    return 1;                         // exit if file not found
  while (!model6.eof())	// loop over model lines
    {
      model6 >> kmer >> value;
      m6[kmer] = value;
    }

  ifstream reads;
  reads.open( argv[1] );              // We assume argv[1] is the file with reads
  if (!reads.good()) 
    return 1;                         // exit if file not found

  vector<int> contextsize,signal;
  vector<float> sumvec,selection;
  float number, info, sum;
  int k,location;

  while (getline(reads, read))	// loop over reads
    {
      sumvec = {};
      vector<float> profile6(read.size(),0.0),profile11(read.size(),0.0),profile12(read.size(),0.0);

      k = 6;
      for (unsigned i = 0; i < read.length()-(k-1); i += 1)
      	{
      	  kmer = read.substr(i, k);
      	  p6 = m6.find(kmer);
      	  if (p6 != m6.end())
      	    {
      	      info = p6->second;  
      	      profile6[i+k-1] += info;
      	    }
      	}
      sum = accumulate( profile6.begin()+(k-1), profile6.end(), 0.0 );
      sumvec.push_back( sum/(1.0*(read.size()-k+1)) ); // evaluate the numb of bits per effective read length 

      if ( k == 6 )
      	selection = profile6;

      vector<float> y = selection;
      int context = k-1;
      int lag = 10;
      float threshold = 1.0;
      float influence =0.5;      
      signal = stop_finder( y, context, lag, threshold, influence );
      for( unsigned i=0; i < read.size() ; i++ )
      	{
	  if ( signal[i] == 1 )
	    cout << " " << read[i];
      	  if ( signal[i] == 0 )
      	    cout << read[i];
      	  if (i == read.size() - 1 )
      	    cout << endl;
      	}
    }// end of loop over reads
} 
