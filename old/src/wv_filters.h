#ifndef WV_FILTERS
#define WV_FILTERS

arma::vec qmf(arma::vec g, bool inverse);

// Filters
arma::field<arma::vec> haar_filter();

// Daubechies 
arma::field<arma::vec> d4_filter() ;

arma::field<arma::vec> d6_filter() ;
arma::field<arma::vec> d8_filter() ;
arma::field<arma::vec> d16_filter() ;

arma::field<arma::vec> fk4_filter() ;
arma::field<arma::vec> fk8_filter() ;
arma::field<arma::vec> fk14_filter() ;
arma::field<arma::vec> fk22_filter() ;

arma::field<arma::vec> bl14_filter() ;
arma::field<arma::vec> bl20_filter() ;

arma::field<arma::vec> la8_filter() ;
arma::field<arma::vec> la16_filter() ;
arma::field<arma::vec> la20_filter() ;

arma::field<arma::vec> mb4_filter() ;
arma::field<arma::vec> mb8_filter() ;
arma::field<arma::vec> mb16_filter() ;
arma::field<arma::vec> mb24_filter() ;


arma::field<arma::vec> select_filter(std::string filter_name);

#endif
