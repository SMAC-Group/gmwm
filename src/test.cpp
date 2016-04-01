#include <Rcpp.h>
#include <map>
#include <vector>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]


int mode(const std::vector<int> &v)
{
  std::map<int, size_t > nmap;

  for (const double  &elm : v)
  {
    nmap[elm]++;
  }
  
  std::map<int, size_t >::const_iterator ci = max_element(
    nmap.begin(),
    nmap.end(),
    [](const pair<double, size_t >& p1, const pair<int, size_t >& p2)
    { return p1.second < p2.second; }
  ); // find highest frequency pair
  
  return ci->first;
}