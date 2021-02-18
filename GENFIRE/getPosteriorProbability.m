function tmp_bayes_probs = getPosteriorProbability(kSlice,calc_kSlice,noise_sigma)
diffs = abs(kSlice-calc_kSlice)./noise_sigma;
tmp_bayes_probs = 1./(2*pi*noise_sigma(:).^2).*(exp( -diffs(:).^2/2 ));
% tmp_bayes_probs = prod(tmp_bayes_probs(tmp_bayes_probs>0));

ln_tmp_bayes_probs = log(tmp_bayes_probs);
tmp_bayes_probs = sum(ln_tmp_bayes_probs);
% tmp_bayes_probs = exp(tmp_bayes_probs);
% % if (any(isnan(tmp_bayes_probs)))
% %     db =  0;
% % end
end
