using Plots
using DelimitedFiles
#-------------------------------------------#

tongue_coords = readdlm("tongue_data.nc", '\t', Float64, '\n');
@show tongue_coords
# now plot!
scatter(tongue_coords[2:end, 1], tongue_coords[2:end, 2], title = "Arnold Tongues")
scatter!(xlabel="\$\\delta\$", ylabel="\$\\epsilon\$")
savefig("ArnoldTongues.png")
