sd_section("Models",
           "Supported Time Series Processes",
           c(
             "SARIMA",
             "SARMA",
             "ARIMA",
             "ARMA",
             "AR",
             "MA",
             "ARMA11",
             "AR1",
             "GM",
             "MA1",
             "WN",
             "DR",
             "QN",
             "RW"
           )
)


sd_section("GMWM Estimator",
           "Estimation Tools for the Generalized Method of Wavelet Moments estimator (GMWM)",
           c(
             "gmwm",
             "gmwm.imu",
             "rgmwm",
             "update.gmwm"
           )
)

sd_section("Model Selection",
           "Rank a group of models by using the Wavelet Information Criterion (WIC)",
           c(
             "rank.models",
             "auto.imu"
           )
)

sd_section("Computations",
           "Tools for computing information",
           c(
             "wvar",
             "wvcov",
             "avar",
             "hadam",
             "modwt",
             "dwt",
             "brickwall"
           )
)


sd_section("Generation",
           "Support for Generating Single or Composite Processes",
           c(
             "gen.gts",
             "gen.lts"
           )
)

sd_section("Inertial Measurement Units",
           "IMU specific tools",
           c("imu",
             "read.imu"
           )
)



sd_section("Visualizations",
           "Graphing Features",
           c("plot.wvar",
             "plot.wvar.imu",
             "compare.wvar",
             "plot.gmwm",
             "compare.models",
             "compare.eff",
             "compare.gmwm",
             "plot.gts",
             "plot.lts",
             "autoplot.avar",
             "autoplot.hadam",
             "plot.avar",
             "plot.hadam",
             "plot.rank.models",
             "plot.auto.imu"
           )
)


sd_section("Data",
           "Helpers to install different SMAC-Group data packages",
           c("install_datapkg",
             "install_imudata"
           )
)


sd_section("Process to Haar Wavelet Variance",
           "Functions that compute the theoretical Haar Wavelet Variance",
           c(
             "arma_to_wv",
             "arma11_to_wv",
             "ar1_to_wv",
             "ma1_to_wv",
             "wn_to_wv",
             "dr_to_wv",
             "qn_to_wv",
             "rw_to_wv"
           )
)


sd_section("Partial First Order Derivatives",
           "Functions that compute the partial first order derivatives of the theoretical Haar Wavelet Variance",
           c("derivative_first_matrix",
             "deriv_arma11",
             "deriv_ar1",
             "deriv_ma1",
             "deriv_wn",
             "deriv_dr",
             "deriv_qn",
             "deriv_rw"
           )
)

sd_section("Partial Second Order Derivatives",
           "Functions that compute the partial second order derivatives of the theoretical Haar Wavelet Variance",
           c(
             "deriv_2nd_arma11",
             "deriv_2nd_ar1",
             "deriv_2nd_ma1",
             "deriv_2nd_dr"
           )
)