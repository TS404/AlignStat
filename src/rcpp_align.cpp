
#include <Rcpp.h>
using namespace Rcpp;

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
