#include <TMB.hpp>

template <class Type>
Type ddirichlet(vector<Type> x, vector<Type> alpha, int give_log) {
  Type logres = lgamma(alpha.sum()) - lgamma(alpha).sum() + ((alpha - Type(1.0)) * log(x)).sum();
  if (give_log) return logres;
  else return exp(logres);
}

// from WHAM
template<class Type>
Type ddirmultinom(vector<Type> obs, vector<Type> alpha, int do_log) {
  int dim = obs.size();
  Type N = obs.sum();
  Type phi = sum(alpha);
  Type ll = lgamma(N + 1.0) + lgamma(phi) - lgamma(N + phi);
  for (int a = 0; a < dim; a++) ll += -lgamma(obs(a) + 1.0) + lgamma(obs(a) + alpha(a)) - lgamma(alpha(a));
  if (do_log == 1) return ll;
  else return exp(ll);
}

// from Stan
// real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
//   real sum_alpha = sum(alpha);
//   return lgamma(sum_alpha) - lgamma(sum(y) + sum_alpha)
//     // + lgamma(sum(y)+1) - sum(lgamma(to_vector(y)+1) // constant, may omit
//     + sum(lgamma(to_vector(y) + alpha)) - sum(lgamma(alpha));
// }
template <class Type>
Type ddm(vector<Type> x, vector<Type> alpha, int give_log) {
  Type sum_alpha = alpha.sum();
  Type logres = lgamma(sum_alpha) - lgamma(x.sum() + sum_alpha) + 
    lgamma(x.sum() + Type(1.0)) - lgamma(vector<Type>(x + Type(1.0))).sum() + // constant, may omit
    lgamma(vector<Type>(x + alpha)).sum() - lgamma(alpha).sum();
  if (give_log) return logres;
  else return exp(logres);
}

template <class Type>
Type ddm_wrong(vector<Type> x, vector<Type> prob, vector<Type> alpha, int give_log) {
  Type lp1 = ddirichlet(prob, alpha, give_log);
  Type lp2 = dmultinom(x, prob, give_log);
  if (give_log) return lp1 + lp2;
  else return exp(lp1 + lp2);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(x);
  DATA_MATRIX(prob);
  DATA_MATRIX(alpha);
  
  PARAMETER(par);
  
  int n = vector<Type>(x.row(0)).size();
  
  vector<Type> multinomial_dens(n);
  vector<Type> dirmult_dens(n);
  vector<Type> dirmult_dens_wham(n);
  vector<Type> dirichlet_dens(n);
  
  for (int i = 0; i < n; i++) {
    // vector<Type> y = x.col(i);
    // Type sumy = y.sum();
    // y /= sumy;
    // vector<Type> alp = prob.col(i) * sumy;
    // dirichlet_dens(i) = ddirichlet(y, alp, true);
    dirichlet_dens(i) = ddirichlet(vector<Type>(prob.col(i)), vector<Type>(alpha.col(i)), true);
    
    multinomial_dens(i) = dmultinom(vector<Type>(x.col(i)), vector<Type>(prob.col(i)), true);
    
    dirmult_dens(i) = ddm(vector<Type>(x.col(i)), vector<Type>(alpha.col(i)), true);
    dirmult_dens_wham(i) = ddirmultinom(vector<Type>(x.col(i)), vector<Type>(alpha.col(i)), true);
  }
  
  // From CCSBT ADMB code
  // lp(i) -= lf_n(i) * (x * log(pred)).sum(); // ln_like(iff) -= Nsamp(iff,iy)*((Nrobust+obs_len_freq_il(ii)(mbin,nbins)*log(Nrobust+pred_len_freq_il(ii)(mbin,nbins))));
  // x += Type(1e-6);
  // lp(i) += lf_n(i) * (x * log(x)).sum(); // mult_constant(iff) += Nsamp(iff,iy)*(1e-6+obs_len_freq_il(irec)(mbin,nbins))*log(1e-6+obs_len_freq_il(irec)(mbin,nbins));
  
  Type nll = Type(0.0);
  
  REPORT(n);
  REPORT(dirichlet_dens);
  REPORT(multinomial_dens);
  REPORT(dirmult_dens);
  REPORT(dirmult_dens_wham);
  REPORT(nll);
  
  return nll;
}
