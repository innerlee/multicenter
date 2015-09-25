
#MM=185
MM=1
margin=.5
N=200
NN=1000
debug=false

#sites=rand(N,2)
# 161,6
# 164, 80
# 185,10
srand(MM)


include("common.jl")
include("gen_data.jl")
include("voronoi.jl")
sites,lab=GenerateTwoClasses()
#println("original labels of pts:\n $(lab)")
include("temp.jl")


println("MM=$MM")
