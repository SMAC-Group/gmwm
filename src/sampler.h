#ifndef SAMPLE_H
#define SAMPLE_H

arma::vec rsample(const arma::vec &x, const int size, const bool replace, arma::vec prob_ );
void RSampleNoReplace(arma::vec &index, int nOrig, int size);
void RSampleReplace(arma::vec &index, int nOrig, int size);
void RProbSampleNoReplace(arma::vec &index, int nOrig, int size, arma::vec &prob);
void RProbSampleReplace(arma::vec &index, int nOrig, int size, arma::vec &prob);
void RWalkerProbSampleReplace(arma::vec &index, int nOrig, int size, arma::vec &prob);
void RFixProb(arma::vec &prob, const int size, const bool replace);

#endif