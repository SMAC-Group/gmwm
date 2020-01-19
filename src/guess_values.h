#ifndef GUESS_VALUES
#define GUESS_VALUES
#include <map>

arma::vec ar1_draw(unsigned int draw_id, double last_phi, double sigma_tot, std::string model_type);

arma::vec arma_draws(unsigned int p, unsigned int q, double sigma2_total);

double dr_slope(const arma::vec& data);

std::string dom_process(double first_wv, double ci_low, double ci_high);

arma::vec guess_initial(const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc,
                        std::string model_type, unsigned int num_param, double expect_diff, unsigned int N,
                        const arma::mat& wv, const arma::vec& tau, double ranged, unsigned int G=1000);

arma::vec guess_initial_old(const std::vector<std::string>& desc, const arma::field<arma::vec>& objdesc,
                            std::string model_type, unsigned int num_param, double expect_diff, unsigned int N,
                            const arma::vec& wv_empir, const arma::vec& tau, unsigned int B = 1000);


#endif
