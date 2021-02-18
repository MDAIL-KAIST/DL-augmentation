function GoodInds = BayesianIndSelectionFunction_topone(Probs)

GoodInds = find(Probs == max(Probs),1);

end