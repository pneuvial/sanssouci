// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// for_back
Rcpp::List for_back(int m, arma::mat A, arma::vec f0x, arma::vec f1x, arma::vec Pi);
RcppExport SEXP _sansSouci_for_back(SEXP mSEXP, SEXP ASEXP, SEXP f0xSEXP, SEXP f1xSEXP, SEXP PiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type f0x(f0xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type f1x(f1xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Pi(PiSEXP);
    rcpp_result_gen = Rcpp::wrap(for_back(m, A, f0x, f1x, Pi));
    return rcpp_result_gen;
END_RCPP
}
// Em_hmm
Rcpp::List Em_hmm(int m, arma::mat A, arma::vec f0x, arma::vec f1x, arma::vec Pi, double eps, int maxit);
RcppExport SEXP _sansSouci_Em_hmm(SEXP mSEXP, SEXP ASEXP, SEXP f0xSEXP, SEXP f1xSEXP, SEXP PiSEXP, SEXP epsSEXP, SEXP maxitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type f0x(f0xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type f1x(f1xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    rcpp_result_gen = Rcpp::wrap(Em_hmm(m, A, f0x, f1x, Pi, eps, maxit));
    return rcpp_result_gen;
END_RCPP
}
// Em_f1
Rcpp::List Em_f1(int m, arma::mat A, arma::vec Pi, arma::vec f0x, arma::vec f1x, Rcpp::List fw_bc_EM, arma::vec x, double eps, int maxit, double h);
RcppExport SEXP _sansSouci_Em_f1(SEXP mSEXP, SEXP ASEXP, SEXP PiSEXP, SEXP f0xSEXP, SEXP f1xSEXP, SEXP fw_bc_EMSEXP, SEXP xSEXP, SEXP epsSEXP, SEXP maxitSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type f0x(f0xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type f1x(f1xSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type fw_bc_EM(fw_bc_EMSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(Em_f1(m, A, Pi, f0x, f1x, fw_bc_EM, x, eps, maxit, h));
    return rcpp_result_gen;
END_RCPP
}
// Em_tot
Rcpp::List Em_tot(int m, arma::mat A, arma::vec Pi, arma::vec f0x, arma::vec f1x, arma::vec x, double eps, int maxit, double h);
RcppExport SEXP _sansSouci_Em_tot(SEXP mSEXP, SEXP ASEXP, SEXP PiSEXP, SEXP f0xSEXP, SEXP f1xSEXP, SEXP xSEXP, SEXP epsSEXP, SEXP maxitSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type f0x(f0xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type f1x(f1xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(Em_tot(m, A, Pi, f0x, f1x, x, eps, maxit, h));
    return rcpp_result_gen;
END_RCPP
}
// Em_tot_01
Rcpp::List Em_tot_01(int m, arma::mat A, arma::vec Pi, arma::vec f0x, arma::vec f1x, arma::vec x, double eps, int maxit, double h);
RcppExport SEXP _sansSouci_Em_tot_01(SEXP mSEXP, SEXP ASEXP, SEXP PiSEXP, SEXP f0xSEXP, SEXP f1xSEXP, SEXP xSEXP, SEXP epsSEXP, SEXP maxitSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type f0x(f0xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type f1x(f1xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(Em_tot_01(m, A, Pi, f0x, f1x, x, eps, maxit, h));
    return rcpp_result_gen;
END_RCPP
}
// colSort
arma::mat colSort(arma::mat X);
RcppExport SEXP _sansSouci_colSort(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(colSort(X));
    return rcpp_result_gen;
END_RCPP
}
// empiricalCoverageO
NumericVector empiricalCoverageO(NumericVector thr, arma::mat Z);
RcppExport SEXP _sansSouci_empiricalCoverageO(SEXP thrSEXP, SEXP ZSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type thr(thrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    rcpp_result_gen = Rcpp::wrap(empiricalCoverageO(thr, Z));
    return rcpp_result_gen;
END_RCPP
}
// get_A
arma::mat get_A(int m, arma::mat alpha, arma::mat beta, arma::mat A, arma::vec f0x, arma::vec f1x, int i);
RcppExport SEXP _sansSouci_get_A(SEXP mSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP ASEXP, SEXP f0xSEXP, SEXP f1xSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type f0x(f0xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type f1x(f1xSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(get_A(m, alpha, beta, A, f0x, f1x, i));
    return rcpp_result_gen;
END_RCPP
}
// get_L1
arma::sp_mat get_L1(arma::mat A, int m, arma::mat alpha, arma::mat beta, arma::vec f0x, arma::vec f1x);
RcppExport SEXP _sansSouci_get_L1(SEXP ASEXP, SEXP mSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP f0xSEXP, SEXP f1xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type f0x(f0xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type f1x(f1xSEXP);
    rcpp_result_gen = Rcpp::wrap(get_L1(A, m, alpha, beta, f0x, f1x));
    return rcpp_result_gen;
END_RCPP
}
// getA01
Rcpp::List getA01(int m, arma::vec li0, arma::vec f0x, arma::vec f1x, Rcpp::List Pis);
RcppExport SEXP _sansSouci_getA01(SEXP mSEXP, SEXP li0SEXP, SEXP f0xSEXP, SEXP f1xSEXP, SEXP PisSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type li0(li0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type f0x(f0xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type f1x(f1xSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Pis(PisSEXP);
    rcpp_result_gen = Rcpp::wrap(getA01(m, li0, f0x, f1x, Pis));
    return rcpp_result_gen;
END_RCPP
}
// getbound
double getbound(int m, double alpha, arma::vec li0, arma::vec f0x, arma::vec f1x, Rcpp::List Pis);
RcppExport SEXP _sansSouci_getbound(SEXP mSEXP, SEXP alphaSEXP, SEXP li0SEXP, SEXP f0xSEXP, SEXP f1xSEXP, SEXP PisSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type li0(li0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type f0x(f0xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type f1x(f1xSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Pis(PisSEXP);
    rcpp_result_gen = Rcpp::wrap(getbound(m, alpha, li0, f0x, f1x, Pis));
    return rcpp_result_gen;
END_RCPP
}
// getIC
arma::vec getIC(int m, double alpha, arma::vec li0, arma::vec f0x, arma::vec f1x, Rcpp::List Pis);
RcppExport SEXP _sansSouci_getIC(SEXP mSEXP, SEXP alphaSEXP, SEXP li0SEXP, SEXP f0xSEXP, SEXP f1xSEXP, SEXP PisSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type li0(li0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type f0x(f0xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type f1x(f1xSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Pis(PisSEXP);
    rcpp_result_gen = Rcpp::wrap(getIC(m, alpha, li0, f0x, f1x, Pis));
    return rcpp_result_gen;
END_RCPP
}
// marginalKFWER
NumericVector marginalKFWER(NumericVector thr, arma::mat Z);
RcppExport SEXP _sansSouci_marginalKFWER(SEXP thrSEXP, SEXP ZSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type thr(thrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    rcpp_result_gen = Rcpp::wrap(marginalKFWER(thr, Z));
    return rcpp_result_gen;
END_RCPP
}
// minPseudoRanks
Rcpp::NumericVector minPseudoRanks(arma::mat X, arma::mat Y);
RcppExport SEXP _sansSouci_minPseudoRanks(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(minPseudoRanks(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// partialColSortDescCpp
arma::mat partialColSortDescCpp(arma::mat X, int k);
RcppExport SEXP _sansSouci_partialColSortDescCpp(SEXP XSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(partialColSortDescCpp(X, k));
    return rcpp_result_gen;
END_RCPP
}
// posthocBySimes0Rcpp
double posthocBySimes0Rcpp(NumericVector p, NumericVector select, double alpha);
RcppExport SEXP _sansSouci_posthocBySimes0Rcpp(SEXP pSEXP, SEXP selectSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type select(selectSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(posthocBySimes0Rcpp(p, select, alpha));
    return rcpp_result_gen;
END_RCPP
}
// rowSortDesc
arma::mat rowSortDesc(arma::mat X);
RcppExport SEXP _sansSouci_rowSortDesc(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(rowSortDesc(X));
    return rcpp_result_gen;
END_RCPP
}
// sim_markov
arma::vec sim_markov(int m, arma::vec Pi, arma::mat A);
RcppExport SEXP _sansSouci_sim_markov(SEXP mSEXP, SEXP PiSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(sim_markov(m, Pi, A));
    return rcpp_result_gen;
END_RCPP
}
// sim_x_kn
NumericVector sim_x_kn(int m, NumericMatrix alpha, NumericMatrix beta, NumericMatrix A, NumericVector Pi, NumericVector f0x, NumericVector f1x);
RcppExport SEXP _sansSouci_sim_x_kn(SEXP mSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP ASEXP, SEXP PiSEXP, SEXP f0xSEXP, SEXP f1xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type f0x(f0xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type f1x(f1xSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_x_kn(m, alpha, beta, A, Pi, f0x, f1x));
    return rcpp_result_gen;
END_RCPP
}
// testBySignFlipping
arma::mat testBySignFlipping(arma::mat X, double B);
RcppExport SEXP _sansSouci_testBySignFlipping(SEXP XSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(testBySignFlipping(X, B));
    return rcpp_result_gen;
END_RCPP
}
// viterbi
arma::vec viterbi(int m, arma::mat A, arma::vec f0x, arma::vec f1x, arma::vec Pi);
RcppExport SEXP _sansSouci_viterbi(SEXP mSEXP, SEXP ASEXP, SEXP f0xSEXP, SEXP f1xSEXP, SEXP PiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type f0x(f0xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type f1x(f1xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Pi(PiSEXP);
    rcpp_result_gen = Rcpp::wrap(viterbi(m, A, f0x, f1x, Pi));
    return rcpp_result_gen;
END_RCPP
}
// viterbi_log
arma::vec viterbi_log(int m, arma::mat A_log, arma::vec f0x_log, arma::vec f1x_log, arma::vec Pi_log);
RcppExport SEXP _sansSouci_viterbi_log(SEXP mSEXP, SEXP A_logSEXP, SEXP f0x_logSEXP, SEXP f1x_logSEXP, SEXP Pi_logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A_log(A_logSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type f0x_log(f0x_logSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type f1x_log(f1x_logSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Pi_log(Pi_logSEXP);
    rcpp_result_gen = Rcpp::wrap(viterbi_log(m, A_log, f0x_log, f1x_log, Pi_log));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_sansSouci_for_back", (DL_FUNC) &_sansSouci_for_back, 5},
    {"_sansSouci_Em_hmm", (DL_FUNC) &_sansSouci_Em_hmm, 7},
    {"_sansSouci_Em_f1", (DL_FUNC) &_sansSouci_Em_f1, 10},
    {"_sansSouci_Em_tot", (DL_FUNC) &_sansSouci_Em_tot, 9},
    {"_sansSouci_Em_tot_01", (DL_FUNC) &_sansSouci_Em_tot_01, 9},
    {"_sansSouci_colSort", (DL_FUNC) &_sansSouci_colSort, 1},
    {"_sansSouci_empiricalCoverageO", (DL_FUNC) &_sansSouci_empiricalCoverageO, 2},
    {"_sansSouci_get_A", (DL_FUNC) &_sansSouci_get_A, 7},
    {"_sansSouci_get_L1", (DL_FUNC) &_sansSouci_get_L1, 6},
    {"_sansSouci_getA01", (DL_FUNC) &_sansSouci_getA01, 5},
    {"_sansSouci_getbound", (DL_FUNC) &_sansSouci_getbound, 6},
    {"_sansSouci_getIC", (DL_FUNC) &_sansSouci_getIC, 6},
    {"_sansSouci_marginalKFWER", (DL_FUNC) &_sansSouci_marginalKFWER, 2},
    {"_sansSouci_minPseudoRanks", (DL_FUNC) &_sansSouci_minPseudoRanks, 2},
    {"_sansSouci_partialColSortDescCpp", (DL_FUNC) &_sansSouci_partialColSortDescCpp, 2},
    {"_sansSouci_posthocBySimes0Rcpp", (DL_FUNC) &_sansSouci_posthocBySimes0Rcpp, 3},
    {"_sansSouci_rowSortDesc", (DL_FUNC) &_sansSouci_rowSortDesc, 1},
    {"_sansSouci_sim_markov", (DL_FUNC) &_sansSouci_sim_markov, 3},
    {"_sansSouci_sim_x_kn", (DL_FUNC) &_sansSouci_sim_x_kn, 7},
    {"_sansSouci_testBySignFlipping", (DL_FUNC) &_sansSouci_testBySignFlipping, 2},
    {"_sansSouci_viterbi", (DL_FUNC) &_sansSouci_viterbi, 5},
    {"_sansSouci_viterbi_log", (DL_FUNC) &_sansSouci_viterbi_log, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_sansSouci(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
