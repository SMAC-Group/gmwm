#ifndef TRANSFORM_DATA
#define TRANSFORM_DATA

arma::vec pseudo_logit_inv(const arma::vec& x);

arma::vec logit_inv(const arma::vec& x);

arma::vec pseudo_logit(const arma::vec& x);

double pseudo_logit(double x);

arma::vec logit(const arma::vec& x);

double logit(double x);

arma::vec transform_values(const arma::vec& theta,
                           const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type);
                           
arma::vec untransform_values(const arma::vec& theta,
                           const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc, std::string model_type);                           

#endif
