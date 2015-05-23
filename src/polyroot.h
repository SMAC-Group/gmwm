#ifndef POLYROOT_H
#define POLYROOT_H

// Function definitions
static void calct_cpp(bool &bol);
static bool fxshft_cpp(const int l2, double &zr, double &zi);
static bool vrshft_cpp(const int l3, double &zr, double &zi);
static void nexth_cpp(const bool bol);
static void noshft_cpp(const int l1);


static void polyev_cpp(const int n, const double s_r, const double s_i, 
                       std::vector<double> &p_r, std::vector<double> &p_i, 
                       std::vector<double> &q_r, std::vector<double> &q_i, 
                       double &v_r, double &v_i);
static double errev_cpp(const int n, 
                        const std::vector<double> &qr, std::vector<double> &qi,
                        const double ms, const double mp, const double a_re, const double m_re);
static double cpoly_cauchy_cpp(const int n, std::vector<double> &pot, std::vector<double> &q);
static double cpoly_scale_cpp(const int n, std::vector<double> &pot,
                              const double eps, const double BIG,
                              const double small, const double base);
static void cdivid_cpp(const double ar, const double ai, 
                       const double br, const double bi, 
                       double &cr, double &ci);

static void polyroot_cpp(const std::vector<double> &opr, const std::vector<double> &opi, int &degree,
                         std::vector<double> &zeror, std::vector<double> &zeroi, bool &fail);

static double myfmod_cpp(double x1, double x2);

double R_pow_cpp(double x, double y);

double R_pow_di_cpp(double x, int n);
#endif