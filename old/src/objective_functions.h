#ifndef OBJECTIVE_FUNCTIONS
#define OBJECTIVE_FUNCTIONS

// Used Yannick's flattening technique on guessed starting values...
double objFunStarting(const arma::vec& theta, 
                      const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type,
                      const arma::vec& wv_empir, const arma::vec& tau);

// Main objective function used by the program
double objFun(const arma::vec& theta,
              const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type,
              const arma::mat& omega,const arma::vec& wv_empir, const arma::vec& tau);
              
double getObjFunStarting(const arma::vec& theta, 
                      const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type,
                      const arma::vec& wv_empir, const arma::vec& tau);
                      
double getObjFun(const arma::vec& theta,
                    const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type,
                    const arma::mat& omega,const arma::vec& wv_empir, const arma::vec& tau);

arma::vec Rcpp_OptimStart(const arma::vec&  theta,
                          const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type,
                          const arma::vec& wv_empir, const arma::vec& tau);

arma::vec Rcpp_Optim(const arma::vec&  theta, 
                     const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type,
                     const arma::mat& omega, const arma::vec& wv_empir, const arma::vec& tau);
#endif
