function err = mse(in1,in2)
% computes mean squared error between two vectors
err = mean((abs(in1-in2)).^2);
end