using LinearAlgebra
using DelimitedFiles
include("getH_even.jl")
include("getH_odd.jl")

# if tongue_data.nc is there, delete it

rm("tongue_data.nc")
#-------------------------------------------#
# create lists of possible values of epsilon and delta
dLength 	= 500; 	# length of delta list
eLength		= 500;	# length of epsilon list
deltaList	= range(-5,stop=20,length=dLength);
epsilonList	= range(-1,stop=60,length=eLength);
#@show epsilonList
#@show deltaList
global tongue_coords	= zeros(Float64, 1, 2); # start with [0 0], append and lose the first [0, 0] when plotting
tol 		= 1e-6; # tolerance
n 			= 5; # number of included modes, fixed. 
H_even1		= zeros(Float64, n, n);
H_even2		= zeros(Float64, n, n);
H_odd1		= zeros(Float64, n, n);
H_odd2		= zeros(Float64, n, n);
for i = 1:eLength
	for j = 2:dLength
		epsilon = epsilonList[i];
		delta	= deltaList[j];
		
		# declare det_even1_prev, det_even2_prev, det_odd1_prev, det_odd2_prev globally

		# global det_even1_prev = 1; global det_even2_prev = 1; global det_odd1_prev = 1; global det_odd2_prev = 1; 

		delta_prev	= deltaList[j-1];
		H_even1_prev, H_even2_prev	= getH_even(epsilon, delta_prev, n); # get the Hill's det
		H_odd1_prev, H_odd2_prev		= getH_odd(epsilon, delta_prev, n); # get the Hill's det  
		det_even1_prev			= det(H_even1_prev) # find det(H_even)
		det_even2_prev			= det(H_even2_prev) # find det(H_even)
		#@show det_even1, det_even2 
		det_odd1_prev			= det(H_odd1_prev) # find det(H_odd)
		det_odd2_prev			= det(H_odd2_prev) # find det(H_odd)	

		H_even1, H_even2	= getH_even(epsilon, delta, n); # get the Hill's det
		H_odd1, H_odd2		= getH_odd(epsilon, delta, n); # get the Hill's det  
		det_even1			= det(H_even1) # find det(H_even)
		det_even2			= det(H_even2) # find det(H_even)
		#@show det_even1, det_even2 
		det_odd1			= det(H_odd1) # find det(H_odd)
		det_odd2			= det(H_odd2) # find det(H_odd)
		@show det_odd1, det_odd2
		if(abs(det_even1)<= tol || abs(det_even2)<= tol || abs(det_odd1)<= tol || abs(det_odd2)<= tol)
			global tongue_coords = vcat(tongue_coords, [delta epsilon])

		# the det is a continuous function of delta, so if I have a coarse grid and there is a sign change in the determinant value, append [(delta+delta_previous)/2 epsilon] to the list. This will give the Arnold tongues!
		
		elseif( (det_even1_prev<= 0 && det_even1>= 0) || ( det_even1_prev>= 0 && det_even1<= 0) || ( det_even2_prev<= 0 && det_even2>= 0) || ( det_even2_prev>= 0 && det_even2<= 0) || (det_odd1_prev<= 0 && det_odd1>= 0) || (det_odd1_prev>= 0 && det_odd1<= 0) || (det_odd2_prev<= 0 && det_odd2>= 0) || (det_odd2_prev>= 0 && det_odd2<= 0))
			global tongue_coords = vcat(tongue_coords, [0.5*(delta + delta_prev) epsilon])

		end



		#eigs_even = eigvals(H_even);
		#eigs_odd = eigvals(H_odd);

		#@show minimum(abs.(eigs_even)), minimum(abs.(eigs_odd))
		#if(minimum(abs.(eigs_even))< tol || minimum(abs.(eigs_odd))< tol)
		#	global tongue_coords = vcat(tongue_coords, [delta epsilon])
		#end
	end	
end
#-------------------------------------------#

open("tongue_data.nc", "w") do io
	writedlm(io, tongue_coords)
end
#-------------------------------------------#
