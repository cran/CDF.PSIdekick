#include <iostream>
#include <cfloat>
#include <fstream>
#include <sstream>
#include <string>
#include <Rcpp.h>
#include <ctime>
#include <time.h>
 
using namespace Rcpp;
using namespace std;
 
 
 
 
double sign(double d){
  if(d<0)
    return -1.0;
  else
    return 1.0;
}
 
double noise(double mu, double b){
   
  double U=as<double>(runif(1,-0.5,0.5));
  double value=mu-b*sign(U)*log(1-2*fabs(U));
   
  return value;
}
 
static double rangeH(int a, int b, std::vector<double> H, int *num) {
  double range = 0.0;
  while (a<=b){
    range += H[a];
    (*num)++;
    a++;
  }
  return range;
}
 
static double rangeT(int a, int b, int f, int height, std::vector<std::vector<double> > T, int *num) {
   
  double range = 0.0;
  //start from the leaf level
  int s = 0;
  do {
    if (s==height) {//if the range covers the whole histogram
      range += T[s][0];
      (*num)++;
       
      return range;
    }
    //add node values starting from position a and moving right
    //until you reach the last child of a's parent
    while (a%f!=0 && a<b){
      range += T[s][a];
      (*num)++;
       
      a++;
    }
    //add node values starting from position b and moving left
    //until you reach the first child of b's parent
     
    while (b%f!=f-1 && a<b){
      range += T[s][b];
      (*num)++;
       
      b--;
    }
    //move the next tree level
    if (b != a) {
      a = a/f;
      b = b/f;
       
      s++;
    } else {
      range += T[s][a];
      (*num)++;
       
      a++;
    }
  } while (a <= b);
   
  return range;
}
 
 
vector<int> dijkstra(std::vector<std::vector<double> > cost, int d) {
   std::vector<double> dist(d+1);
  //double *dist=new double[d+1];
   std::vector<int> pred(d+1);
  //int pred[d+1];  // preceeding node in path
   
  for (int i = 0; i < d+1; i++) {
    dist[i] = DBL_MAX;
     
    pred[i] = -1;
  }
   
  dist[0] = 0.0;
   
  for (int i = 0; i < d; i++) {
     
    for (int j = i+1; j < d+1; j++) {
       
      double c = dist[i] + cost[i][j];
       
      if (dist[j] > c) {
         
        dist[j] = c;
        pred[j] = i;
      }
    }
  }
   
  vector<int> path;
  vector<int>::iterator it;
  it = path.begin();
  int x = d;
  while (x != 0 && x != -1) { 
     
    path.insert(it,x);
    it = path.begin();
    x = pred[x];
  }
  path.insert(it,0);
  //delete[] dist;
  //dist=0;
  return path;
}
 
 
 
void S2(std::vector<double> H, std::vector<double> Ho, double eps1, double eps2, int d) {
  //Smoothing sub-module
  std::vector<double> Hn(d);
  //double Hn[d];//noisy histogram
  for(int i=0;i<d;i++)
    Hn[i]=0;
  //we represent each group of bins as an edge in a tournament graph
  //each bin is represented by a node
  //an edge from node 1 to node 5 represents grouping and averaging bins from 1 to 5
   std::vector<std::vector<double> > cost(d+1, std::vector<double>(d+1));
   
 // double **cost=new double *[d+1];
 // for(int i=0;i<d+1;i++)
 //   cost[i]=new double[d+1];
   
   
   
  for(int i=0;i<d+1;i++)
    for(int j=0;j<d+1;j++)
      cost[i][j]=DBL_MAX;
   
   std::vector<double> pref(d+1);
   std::vector<double> pref1(d+1);
   std::vector<double> pref2(d+1);
  //double pref1[d+1];//prefix sums array for computing the sum of histogram values
  //double pref2[d+1];//prefix sums array for computing the sum of squared histogram values
  //double pref[d+1];
  //Create Hn, pref1, and pref2
  pref1[0]=pref2[0]=pref[0]=0;
  for (int i = 0; i < d; i++) {
    Hn[i] = H[i] + noise(0, 2.0/eps1);
     
    pref[i+1]=pref[i]+H[i];
    pref1[i + 1] = pref1[i] + Hn[i];
    pref2[i + 1] = pref2[i] + pow(Hn[i], 2);
  }
   
   
   
   
  //for each edge/group calculate the error due to averaging and due to noise addition
  for (int i = 0; i < d; i++) {
    for (int j = i; j < d; j++) {
      double weight;
      //error due to averaging
      if (j != i) {
        double avg;
        avg = (pref1[j + 1] - pref1[i]);
        weight = pref2[j + 1] - pref2[i] - pow(avg, 2) / (double) (j - i + 1);
        weight -= (j - i) * pow(sqrt(2)*2.0/ eps1, 2);//fix the bias due to noise
      } else {
        weight = 0.0;
      }
      //add the error of noise addition in order to compute the total error
      //set this error as the edge's weight
      cost[i][j+1]=max(weight + 2.0*pow(1.0/eps2+1.0/(eps2*sqrt(j - i + 1)),2), 2.0*pow(1.0/eps2+1.0/(eps2*sqrt(j - i + 1)),2));//,Math.pow(Math.sqrt(2)/eps2,2)/(double)(j-i+1)));//+cost*Math.log(ts)/(Math.log(2)*ts)));//+ ((j-i+1)*cost/(double)ts)));
       
       
       
    }
  }
   
  //by computing the shortest path from the first node/bin to the last one
  //we find the grouping that minimizes the total error per bin
  vector<int> path = dijkstra(cost, d);
   
  /////////////////////////////////////////////////////////////////////////////
  //group and average
  //use the grouping/shortest path on the original histogram
  for (int k = 0; k < path.size() - 1; k++) {
    double avg = 0;
     
    avg = (pref[path[k + 1]]-pref[path[k]]) / (double) (path[k + 1] - path[k]); //group and average bin values
     
    //add reduced noise
    avg += noise(0, 1.0/eps2+1.0/(eps2 * (path[k + 1] - path[k])));
    //set the noisy smoothed value as the bin value
    for (int l = path[k]; l < path[k + 1]; l++) {
      Ho[l] = avg;
    }
  }
  //for (int j=0;j<d+1;j++)
  //  delete[] cost[j];
  //delete[] cost;
  //cost = 0;
}
 
vector<int> dijkstra2(std::vector<std::vector<double> > cost, int d, vector<int> *edges) {
   
  std::vector<double> dist(d+1);
  //double *dist=new double[d+1];
   std::vector<int> pred(d+1);
  //int pred[d+1];  // preceeding node in path
   
  for (int i = 0; i < d+1; i++) {
    dist[i] = DBL_MAX;
     
    pred[i] = -1;
  }
   
  dist[0] = 0.0;
   
   
   
  for (int i = 0; i < d; i++) {
    for (int j=0; j < edges[i].size(); j++) {
      double c = dist[i] + cost[i][edges[i][j]];
       
      if (dist[edges[i][j]] > c) {
        dist[edges[i][j]] = c;
        pred[edges[i][j]] = i;
      }
    }
  }
  vector<int> path;
  vector<int>::iterator it;
  it = path.begin();
  int x = d;
  while (x != 0 && x != -1) { 
     
    path.insert(it,x);
    it = path.begin();
    x = pred[x];
  }
  path.insert(it,0);
  //delete [] dist;
  //dist=0;
  return path;
}
 
void SUB(std::vector<double> Hin, int b, int d, std::vector<std::vector<double> > To, double eps1, double eps2) {
  //Smoothing sub-module 
   
  int height = (int) ceil(log(d) / log(b));
  vector<int> *edges=new vector<int>[(int) pow(b, height)];
  std::vector<double> H((int) pow(b, height));
  std::vector<double> Hn((int) pow(b, height));
  std::vector<std::vector<double> > cost((int) pow(b, height)+1, std::vector<double>((int) pow(b, height)+1));
  //double *H = new double[(int) pow(b, height)];
  //double *Hn = new double[(int) pow(b, height)];
  //double **cost=new double *[(int) pow(b, height)+1];
  //for(int i=0;i<(int) pow(b, height)+1;i++)
  //  cost[i]=new double[(int) pow(b, height)+1];
   
   
   
  for(int i=0;i<d;i++)
    H[i]=Hin[i];
  //in case the histogram size cannot be divided with b
  //add zero bins in order to form a complete tree
  for (int i = d; i < pow(b, height); i++) {
    H[i] = 0;
  }
   std::vector<double> pref((int) pow(b, height)+1);
   std::vector<double> pref1((int) pow(b, height)+1);
   std::vector<double> pref2((int) pow(b, height)+1);
 // double *pref1 = new double[(int) (pow(b, height)) + 1];
 // double *pref2 = new double[(int) (pow(b, height)) + 1];
 // double *pref = new double[(int) (pow(b, height)) + 1];
  pref1[0]=pref2[0]=pref[0]=0;
  for (int i = 0; i < d; i++) {
    Hn[i] = max(0.0,H[i] + noise(0,2.0/eps1));
    pref[i+1]=pref[i]+H[i];
    pref1[i + 1] = pref1[i] + Hn[i];
    pref2[i + 1] = pref2[i] + pow(Hn[i], 2);
  }
  for (int i = d; i < pow(b, height); i++) {
    pref1[i + 1] = pref1[i];
    pref2[i + 1] = pref2[i];
    pref[i+1]=pref[i];
  }
   
   
   
  //calculate the error only for the groups of bins that are leaves of complete subtrees
   
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < (pow(b, height)); i += pow(b, j)) {
      double weight;
      double avg;
      avg = (pref1[(int) (i + pow(b, j))] - pref1[i]);
      weight = pref2[(int) (i + pow(b, j))] - pref2[i] - pow(avg, 2) / (double) (pow(b, j));
      weight += 2.0*pow(ceil(log(d) / log(b)),2)*pow(1.0/eps2+1.0/(eps2*sqrt((j + 1) * pow(b, j))),2);
      weight -= (pow(b, j) - 1) * pow(sqrt(2)*2.0 / eps1, 2);
      cost[i][i + (int) pow(b, j)]=max(weight, 2.0*pow(ceil(log(d) / log(b)),2)*pow(1.0/eps2+1.0/(eps2*sqrt((j + 1) * pow(b, j))),2));
      edges[i].push_back((i+(int) pow(b, j)));
       
    }
     
  }
   
  vector<int> path = dijkstra2(cost, pow(b, height),edges);
   
  delete [] edges;
   
  ///////////////////////////////////////////////////////////////////////
  //Build tree and add noise using eps2
   
  //for the leaf level
  for (int k = 0; k < path.size() - 1; k++) {
     
    double avg = 0.0;
     
    avg = (pref[path.at(k + 1)]-pref[path.at(k)]) / (double) (path.at(k + 1) - path.at(k));
     
    double nois = 0.0;
    //calculate the noise for the root of the pruned subtree assuming 
    //we have as many noisy value estimations as the subtree height
    for (int i = 0; i < (ceil(log(path.at(k + 1) - path.at(k)) / log(b)) + 1); i++) {
      nois += noise(0, (ceil(log(d) / log(b)))/(double)eps2+(ceil(log(d) / log(b))) / (double) (eps2 * (path.at(k + 1) - path.at(k))));
    }
    //use the average of the noisy value estimations
     
    avg += nois / (ceil(log(path.at(k + 1) - path.at(k)) / log(b)) + 1);
     
    for (int l = path.at(k); l < path.at(k + 1); l++) {
      To[0][l] = avg;
       
    }
  }
   
  for (int i = d; i < pow(b, height); i++) {
    To[0][i] = 0.0;
  }
   
  //for the non-leaf levels
  for (int j = 1; j < height; j++) {
    for (int i = 0; i < pow(b, height - j); i++) {
      double sum = 0;
      for (int k = 0; k < b; k++) {
        sum += To[j - 1][b * i + k];
      }
      //check if the children are pruned/smoothed
       
      if (abs(sum - ((To[j - 1][b * i] * b))) > 0.001) {
         
        sum = pref[(int)((i+1)*pow(b,j))]-pref[(int)(i*pow(b,j))];
         
        sum += noise(0,ceil(log(d) / log(b))*2.0 / (double) (eps2));
      }
       
      To[j][i] = sum;
       
    }
  }
  //for (int j=0;j<(int) pow(b, height)+1;j++)
  //  delete[] cost[j];
  //delete[] cost;
  //cost = 0;
  //delete[] H;
  //H=0;
  //delete[] Hn;
  //Hn=0;
  //delete[] pref;
  //pref=0;
  //delete[] pref1;
  //pref1=0;
  //delete[] pref2;
  //pref2=0;
}
 
 
//' @title Creates a Tree then a CDF
//' @description This thing sure does make a fine CDF
//' @param eps An epsilon value for Differential Privacy
//' @param ds The data or something
//' @param Ks the degree of the tree
//' @param methods Either H or S2 or SUB
//' @param mins the minimum of the domain's range
//' @param maxs the maximum of the domain's range
//' @param grans The granularity
//' @param datas The data to be CDFd
//' @return A dpCDF
//' @export
// [[Rcpp::export]]
NumericVector TreeCDF(SEXP eps,SEXP ds, SEXP Ks,SEXP methods,SEXP mins, SEXP maxs, SEXP grans, SEXP datas) {
  float epsilon=as<float>(eps);
  int d=as<int>(ds);
  int f=as<int>(Ks);
  string method=as<string>(methods);
  double gran=as<double>(grans);
  int min_v=as<int>(mins);
  int max_v=as<int>(maxs);
  NumericVector data=as<NumericVector>(datas);
  int n=data.size();
   
   std::vector<double> histogram(d);
  //double histogram[d];
  for(int i=0;i<d;i++)
    histogram[i]=0;
  for(int i=0;i<n;i++){
    int myNumber=data[i];
    if(myNumber<min_v)
      myNumber=min_v;
    else if (myNumber>max_v)
      myNumber=max_v;
    histogram[(int)floor((myNumber-min_v)/(gran))]++;
  }
   
   
   
   
   
  int height=ceil(log(d)/log(f));
   
   
  double mu=0.0,b=(double)height*2.0/epsilon;
   
   std::vector<double> var(d);
   std::vector<double> cdf(d);
 // double* var = new double[d];
 // double* cdf = new double[d];
   
  if(method=="H"){
    std::vector<double> prefix((int)pow(f,height)+1);
    //std::vector<std::vector<double> > node((int)pow(f,height), std::vector<double>(height+1));
    std::vector<std::vector<double> > node(height+1, std::vector<double>((int)pow(f,height)));
    //std::array<std::array<double, (int)pow(f,height)>, height+1> node;
   // double **node=new double *[height+1];
   // for(int i=0;i<=height;i++)
   //   node[i]=new double[(int)pow(f,height-i)];
    node[height][0]=n;
     
     
     
    prefix[0]=0;
    for(int i=0;i<d;i++){
       
      prefix[i+1]=prefix[i]+histogram[i];
    }
     
    for(int i=d;i<(int)pow(f,height);i++){
       
      prefix[i+1]=prefix[i];//+histogram[i];
    }
     
     
    for(int i=height-1;i>=0;i--){
      for(int j=0;j<pow(f,height-i);j++){
        if(i==0&&j>d)
          node[i][j]=0;
        else{
          node[i][j]=(prefix[(j+1)*(int)pow(f,i)]-prefix[j*(int)pow(f,i)]+noise(mu,b));
          // node[i][j]=std::max(0.0, node[i][j]);
        }
         
      }
       
    }
     
     
     
    for(int i=0;i<d-1;i++){
      int num1=0, num2=0;
      double result1=(rangeT(0,i,f,height,node,&num1));
      double result2=((n-rangeT(i+1,d-1,f,height,node,&num2)));
       
      cdf[i]=((num2/(double)(num1+num2))*result1/(double)n+(num1/(double)(num1+num2))*result2/(double)n);
       
      var[i]= (num1*(num2/(double)(num1+num2))*(num2/(double)(num1+num2))*8.0/pow(epsilon,2)*pow(height,2)+num2*(num1/(double)(num1+num2))*(num1/(double)(num1+num2))*8.0/pow(epsilon,2)*pow(height,2))/pow(n,2);
       
    }
    cdf[d-1]=1;
    var[d-1]=0;
     
    //for (int j=0;j<=height;j++)
    //  delete[] node[j];
    //delete[] node;
     
    //node = 0;
  }
   
  else if(method=="S2"){
    std::vector<double> Ho(d);
    //double *Ho=new double[d];
    S2(histogram, Ho, epsilon/2.0, epsilon/2.0, d);
     
    for(int i=0;i<d-1;i++){
       
      int num1=0, num2=0;
      double result1=(rangeH(0,i,Ho,&num1));
      double result2=((n-rangeH(i+1,d-1,Ho,&num2)));
       
      cdf[i]=((num2/(double)(num1+num2))*result1/(double)n+(num1/(double)(num1+num2))*result2/(double)n);
      var[i]= (num1*(num2/(double)(num1+num2))*(num2/(double)(num1+num2))*8.0/pow(epsilon/2.0,2)+num2*(num1/(double)(num1+num2))*(num1/(double)(num1+num2))*8.0/pow(epsilon/2.0,2))/pow(n,2);
       
    }
    cdf[d-1]=1;
    var[d-1]=0;//sqrt(2)*8.0/epsilon;
     
     
   // delete[] Ho;
   // Ho = 0;
  }
  else if(method=="SUB"){
    std::vector<std::vector<double> > To(height+1, std::vector<double>((int)pow(f,height)));
    //std::array<std::array<double, (int)pow(f,height)>, height+1> To;
    //double **To=new double *[height+1];
     
    //for(int i=0;i<=height;i++)
    //  To[i]=new double[(int)pow(f,height-i)];
     
     
    SUB(histogram, f, d, To, epsilon/4.0, 3.0*epsilon/4.0);
    To[height][0]=n;
     
     
     
     
    for(int i=0;i<d-1;i++){
      int num1=0, num2=0;
      double result1=(rangeT(0,i,f,height,To,&num1));
      double result2=((n-rangeT(i+1,d-1,f,height,To,&num2)));
      cdf[i]=((num2/(double)(num1+num2))*result1/(double)n+(num1/(double)(num1+num2))*result2/(double)n);
      var[i]= (num1*(num2/(double)(num1+num2))*(num2/(double)(num1+num2))*8.0/pow(3.0*epsilon/4.0,2)*pow(height,2)+num2*(num1/(double)(num1+num2))*(num1/(double)(num1+num2))*8.0/pow(3.0*epsilon/4.0,2)*pow(height,2))/pow(n,2);
       
    }
    cdf[d-1]=1;
    var[d-1]=0;
     
    //for (int j=0;j<=height;j++)
    //  delete[] To[j];
    //delete[] To;
    //To = 0;
     
  }
   
   
   
  NumericVector out(2*d);
  for(int i=0;i<d;i++)
    out[i]=cdf[i];
  for(int i=d;i<2*d;i++)
    out[i]=var[i-d];
   
   
  //delete[] var;
  //delete[] cdf;
   
  return out;
}
