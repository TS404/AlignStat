
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
CharacterMatrix rcpp_prepare_alignment_matrix(CharacterMatrix ref){
//   printf("hello");
//   
  int nrow = ref.nrow();
  int ncol = ref.ncol();
  
  CharacterMatrix mat2(nrow,ncol);
  
  for (int s=0;s<ncol;s++){
    std::map<String, int> base_counts;
    for(int l=0;l<nrow;l++){
      String b = ref(l,s);
      if ( base_counts.count(b) == 0){
        base_counts[b]=0;
      }
      
      base_counts[b] +=1;
      // printf("bc %i",base_counts[b]);
      b += std::to_string(base_counts[b]);
      // printf(" %s \n",b.get_cstring());      
      mat2(l,s) = b;
    }
    
  }
  return mat2;
}

// [[Rcpp::export]]
List rcpp_align(CharacterMatrix ref,
                         CharacterMatrix aln) {

  int nrow_ref = ref.nrow();
  int nrow_aln = aln.nrow();
  int ncol_ref = ref.ncol();

  NumericMatrix results(10,nrow_ref);    
  NumericMatrix means(nrow_aln,nrow_ref);
    
  int ci = 0;
  for(int k=0;k < nrow_ref;k++){

    int maxi = -1;
    double maxval = -1;
    for(int i=0;i<nrow_aln;i++){
      int ident_count = 0;
      for(int j=0;j<ncol_ref;j++){
        String r = ref(k,j);
        String a = aln(i,j);        

        int sc = strcmp(r.get_cstring(),a.get_cstring());
        if (sc==0){
          ident_count += 1;
          ci += 1;
        }
      }
      double mn_ik = (1.0*ident_count)/(1.0*ncol_ref);
      means(i,k) = mn_ik;
      if ( mn_ik > maxval){
        maxi=i;
        maxval = mn_ik;
      }
    }
    results(0,k) = maxi+1;
  }
  
  List retlist = List::create(Rcpp::Named("results")=results,
                              Rcpp::Named("means")=means);
  return(retlist);
}
