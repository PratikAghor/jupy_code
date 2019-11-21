function getH_odd(epsilon, delta, n)
	H_odd1 = zeros((n, n)); # declare Hill's matrix
	H_odd2 = zeros((n, n)); # declare Hill's matrix
	# hardcode the first two rows and the last row
	
	H_odd1[1, 1] = delta - 1*1 + 0.5*epsilon; H_odd1[1, 2] = 0.5*epsilon;

	for i = 2:n-1
		H_odd1[i, i-1]	= 0.5*epsilon;
		H_odd1[i, i]	= delta - (2*i-1)*(2*i-1);
		H_odd1[i, i+1]	= 0.5*epsilon;		
	end
	H_odd1[n, n-1]	= 0.5*epsilon;
	H_odd1[n, n]	= delta - (2*n-1)*(2*n-1);

	# H_odd2
	H_odd2 = H_odd1;
	H_odd2[1, 1] = delta - 1*1 - 0.5*epsilon;

return H_odd1, H_odd2;
end
