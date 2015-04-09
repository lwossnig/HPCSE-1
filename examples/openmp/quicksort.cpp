// Example codes for HPC course
// (c) 2012 Matthias Troyer, ETH Zurich

#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <iterator>


// this may not be optimal but it is simple and short
template <class It>
void quicksort(It first, It last)
{
  // empty sequence or length 1: we are done
  if (last-first <= 1)
    return;
  
  // pick a random value (here the first) and partition the sequence by it
  typedef typename std::iterator_traits<It>::value_type value_type;
  // pick a pivot
  value_type pivot = *(last-1);
  It split = std::partition(first,last,[=](value_type x) { return x < pivot;});
  // move the pivot to the center
  std::swap(*(last-1),*split);
  
  // sort the two partitions individually
  quicksort(first,split);
  quicksort(split+1,last);
}

int main()
{
  int n;
  std::cin >> n;

  // create random numbers
  std::mt19937 mt;
  std::uniform_int_distribution<int> dist(0,std::numeric_limits<int>::max());
  std::vector<int> data(n);
  std::generate(data.begin(),data.end(),std::bind(dist,mt));
  
  // check if it is sorted
  if (std::is_sorted(data.begin(), data.end()))
      std:: cout << "Initial data is sorted.\n";
  else
      std:: cout << "Initial data is not sorted.\n";

  // call quicksort
   quicksort(data.begin(),data.end());
      
  // check if it is sorted
  if (std::is_sorted(data.begin(), data.end()))
      std:: cout << "Final data is sorted.\n";
  else
      std:: cout << "Final data is not sorted.\n";

  //std::copy(data.begin(),data.end(),std::ostream_iterator<int>(std::cout,"\n"));
}