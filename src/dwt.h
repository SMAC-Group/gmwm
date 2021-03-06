#ifndef DWT
#define DWT

arma::field<arma::vec> brick_wall(arma::field<arma::vec> x,  
                                  arma::field<arma::vec> wave_filter, 
                                  std::string method = "modwt") ;

arma::field<arma::vec> dwt_cpp(arma::vec x, std::string filter_name = "haar", 
                               unsigned int nlevels = 4, std::string boundary = "periodic", bool brickwall = true);
                                 
arma::field<arma::vec> modwt_cpp(arma::vec x, std::string filter_name = "haar", 
                                 unsigned int nlevels = 4, std::string boundary = "periodic", bool brickwall = true);
#endif
