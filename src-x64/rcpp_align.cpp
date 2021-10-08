#include <sstream>

#define SSTR( x ) std::to_string(x)

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
      b += SSTR( base_counts[b] );
      // printf(" %s \n",b.get_cstring());      
      mat2(l,s) = b;
    }
    
  }
  return mat2;
}

// [[Rcpp::export]]
List rcpp_align(CharacterMatrix ref,
                         CharacterMatrix com) {

  int nrow_ref = ref.nrow();
  int nrow_com = com.nrow();
  int ncol_ref = ref.ncol();

  NumericMatrix results(9,nrow_ref);    
  NumericMatrix means(nrow_com,nrow_ref);
    
  int ci = 0;
  for(int k=0;k < nrow_ref;k++){

    int maxi = -1;
    double maxval = -1;
    for(int i=0;i<nrow_com;i++){
      int ident_count = 0;
      for(int j=0;j<ncol_ref;j++){
        String r = ref(k,j);
        String a = com(i,j);        

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
  
  CharacterMatrix cat(nrow_ref,ncol_ref);

  // Categorise differences
  for (int i=0;i<nrow_ref;i++){
    for(int j=0;j<ncol_ref;j++){
      String ref_value = ref(i,j);
      String com_value = com(results(0,i)-1,j);

      if ( ref_value == com_value){
        cat(i,j) = "M";
      } else if ( ref_value=="-"){
        if ( com_value==NA_STRING){
          cat(i,j) = "g";
        } else {
          cat(i,j) = "m";
        }
      } else {
        if (com_value==NA_STRING){
          cat(i,j) = "s";
        } else {
          cat(i,j) = "x";
        }
      }
    }
  }
  
  List retlist = List::create(Rcpp::Named("results")=results,
                              Rcpp::Named("means")=means,
                              Rcpp::Named("cat")=cat);
  return(retlist);
}
