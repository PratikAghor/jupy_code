function getH_even(epsilon, delta, n)
	H_even1 = zeros((n, n)); # declare Hill's matrix
	H_even2 = zeros((n, n)); # declare Hill's matrix
	# hardcode the first two rows and the last row
	H_even1[1, 1] = delta; H_even1[1, 2] = 0.5*epsilon;
	H_even1[2, 1] = epsilon; H_even1[2, 2] = delta - 4; H_even1[2, 3] = 0.5*epsilon;
	for i = 3:n-1
		H_even1[i, i-1]	= 0.5*epsilon;
		H_even1[i, i]	= delta - 4.0*(i-1)*(i-1);
		H_even1[i, i+1]	= 0.5*epsilon;		
	end
	H_even1[n, n-1]	= 0.5*epsilon; H_even1[n, n]	= delta - 4.0*(n-1)*(n-1);


	H_even2 = zeros((n, n)); # declare Hill's determinant

	# hardcode the first two rows and the last row
	
	H_even2[1, 1] = delta - 4; H_even2[1, 2] = 0.5*epsilon;

	for i = 2:n-1
		H_even2[i, i-1]	= 0.5*epsilon;
		H_even2[i, i]	= delta - 4.0*(i)*(i);
		H_even2[i, i+1]	= 0.5*epsilon;		
	end
	H_even2[n, n-1]	= 0.5*epsilon;
	H_even2[n, n]	= delta - 4.0*(n)*(n);

return H_even1, H_even2;
end
