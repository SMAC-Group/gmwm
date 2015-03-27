#ifndef WV_FILTERS
#define WV_FILTERS

arma::vec qmf(arma::vec g, bool inverse);

arma::field<arma::vec> haar_filter();

arma::field<arma::vec> select_filter(std::string filter_name);

#endif