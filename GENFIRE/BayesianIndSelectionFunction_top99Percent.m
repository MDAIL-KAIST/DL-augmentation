function GoodInds = BayesianIndSelectionFunction_top99Percent(Probs)

[SortProb, SortInd] = sort(Probs,'descend');

currSum = 0;
for i=1:length(SortInd)
  currSum = currSum + SortProb(i);
  if currSum > 0.99
    break
  end
end

GoodInds = SortInd(1:i);

end