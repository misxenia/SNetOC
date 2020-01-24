function y = compute_ks(deg1, deg2, normalize)
	%  reweighted kolmogorov-smirnov statistic

	if nargin==2
	    normalize = true;
	end
	[f1, x1] = ecdf(deg1);
	[ff1, fx1] = fill_ecdf(f1(2:end), x1(2:end));
	[f2, x2] = ecdf(deg2);
	[ff2, fx2] = fill_ecdf(f2(2:end), x2(2:end));

	lb = max(fx1(1), fx2(1));
	ub = min(fx1(end), fx2(end));

	S = ff1(fx1>=lb & fx1 <= ub);
	P = ff2(fx2>=lb & fx2 <= ub);
	if normalize == true
	    norm = sqrt(P.*(1-P));
	    ind = norm > 0;
	    y = max(abs(S(ind)-P(ind))./norm(ind));
	else
	    y = max(abs(S-P));
	end

end