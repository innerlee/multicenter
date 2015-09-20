# Gnerate test data for multi-center algorisms



###
# 2D, 2 classes
###
function GenerateTwoClasses(;
    centers=[-1 0; 1 0; 0 1],
    nums=[100 100 100],
    rs=[.5 .5 .5],
    show=false)

  dim=size(centers,2)
  say("[ GenerateTwoClasses | dim=$dim,show=$show ]")

  N=sum(nums)
  pts=zeros(N,dim)
  labs=zeros(N)

  s=1
  for i=1:length(nums)
    temp=GenerateBall(center=centers[i,:],r=rs[i],num=nums[i])
    pts[s:s+nums[i]-1,:]=temp
    labs[s:s+nums[i]-1]=(ones(nums[i])*i)[:]
    s+=nums[i]
  end

  show?drawpic([pts labs]):Nothing

  return pts,labs
end


# function helpers

# generate balls
function GenerateBall(;center=[0 0],r=.5,num=100,show=false)
  say("[ GenerateBall | center=$center,r=$r,num=$num,show=$show ]")

  dim=length(center)
  pts=zeros(num,dim)
  N=1
  M=1
  while N<=num
    temp=rand(dim)-.5
    if (temp'*temp)[1]<.25
      pts[N,:]=temp
      N+=one(N)
    end
    M+=one(M)
  end

  pts=broadcast(+,center[:]',pts*(r/.5))

  #say("ran $(M-1) cycles to generate $(N-1) pts.")
  show?drawpic(pts):Nothing

  return pts
end

# generate rectangles
function GenerateRect(;dim=2,r=10,num=100,show=true)

end


### tests
#pts=GenerateBall(show=true,num=100,r=1,center=[2, 5])

pts,labs=GenerateTwoClasses(show=true)
println("hi")
