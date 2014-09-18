// // block shuffering in Cpp
// block_shuffer(std::vector<double> x, std::vector<int> id){
//   auto p0 = id.begin();
//   auto p1 = id.end();
//   std::set<int> u_id(p0, p1);
//   std::vector<double> xx;
//   for (const auto& this_id : u_id){
//     auto pp = equal_range(p0, p1, this_id);
//     copy(x.begin() + (pp - p0),
//    x.begin() + (p1 - p0),
//    back_insert(xx));
//   }
//   return xx;
// }


// // TODO: make this function generic
// // Created on 02/Jul/2014
// // Author Yi
// // functionality: create unique ids for 0s and 1s. for example
// // OT series: 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0,
// // unique id: 0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4,
// #include <Rcpp.h>
// #include <iostream>
// #include <vector>
// #include <algorithm>
// #include <iterator>
// // [[Rcpp::plugins(cpp11)]]
// using namespace std;
// // [[Rcpp::export]]
// std::vector<int> unique_id_along_01 (std::vector<int> ind){
//   int count = 1 ;
//   auto i0 = ind.begin();
//   auto i1 = ind.end();
//   auto p1 = find(i0, i1, 1); // find first 1;
//   auto p0 = find(p1, i1, 0); // find first 0 after first 1.
//   while (p1 != i1){
//     fill(p1, p0, count++);
//     p1 = find(p0, i1, 1);
//     fill(p0, p1 , count++);
//     p0 = find(p1, i1, 0);
//   }
// return ind;
// }

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
//
// [[Rcpp::export]]
double obj_cpp(NumericVector para, NumericVector yi, NumericVector yj){
  static double aLow = -1 + pow(10, -10);
  static double BigNumber = pow(10, 40);
  static double WeeNumber = pow(10, -10);
  double a = para[0];
  double b = para[1];
  NumericVector z = (yj - yi * a ) / pow(yi, b);
  double m = mean(z);
  double s = sd(z);
  NumericVector mu = yi * a + m * pow(yi, b);
  NumericVector sig = s * pow(yi, b);
  if (a < aLow | s < WeeNumber | a > 1 - WeeNumber | b > 1 - WeeNumber)
    return BigNumber;
  double res = sum(0.5 * log(2*PI) + log(sig) + 0.5 * pow((yj - mu)/sig, 2));
  if (is_infinite(NumericVector::create(res))[0] == 1)
    return(sign(NumericVector::create(res))[0] * BigNumber);
  return(res);
}
 