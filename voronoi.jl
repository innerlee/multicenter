say("voronoi now")

#
function composepoints(pts)
  compose(context(),circle(pts[:,1],pts[:,2],[.7mm]),fill(["black"]), stroke(["white"]))
end

#
function transformlines(lns)
  [[(lns[i,1]+lns[i,5]*lns[i,3],lns[i,2]+lns[i,5]*lns[i,4]),(lns[i,1]+lns[i,6]*lns[i,3],lns[i,2]+lns[i,6]*lns[i,4])] for i=1:size(lns,1),j=1]
end

#
function composelines(lns::AbstractArray)
  say()
  say("= in composelines() =")
  ls=transformlines(lns)
  #say(ls)
  #lines=[[(ls[i,1]+ls[i,5]*ls[i,3],ls[i,2]+ls[i,5]*ls[i,4]),(ls[i,1]+ls[i,6]*ls[i,3],ls[i,2]+ls[i,6]*ls[i,4])] for i=1:size(ls,1),j=1]
  compose(context(), line(ls),fill([colorant"black"]), stroke([colorant"black"]))
end

#
function composelines(ls)
  compose(context(), line(ls),fill([colorant"black"]), stroke([colorant"black"]))
end

# lns: n*6 (x0,y0, vx,vy, a,b)
# ls: n*4 ((x0,y0) (x1,y1))
# use borderrect
function drawvor(pts,lns,border)
  # preprocess them
  say()
  say("= drawvor =")

  # clip lines
  for i=1:size(lns,1)
    ts=[intersectlns(lns[i,:],border[k,:])[1] for k=1:size(border,1)]
    #say("*ts*")
    #say(ts)
    negs=ts[ts.<0]
    if length(negs)>0
      lns[i,5]=maximum([negs lns[i,5]])
    end
    pos=ts[ts.>0]
    if length(pos)>0
      lns[i,6]=max(lns[i,5],minimum([ts[ts.>0] lns[i,6]]))
    end
  end

  # transforms
  lns=[border;lns]
  ratio=max(b_width,b_height)
  pts=broadcast(/,broadcast(-,pts,[b_left b_bot]),[ratio ratio])
  lns=broadcast(/,broadcast(-,lns,[b_left b_bot 0 0 0 0]),[ratio ratio 1 1 ratio ratio])
  say(size(lns))

  drawcomp(compose(context(),composepoints(pts),composelines(lns)))
end

#
function intersectlns(ln1,ln2)
  A=[ln1[3] -ln2[3]
     ln1[4] -ln2[4]]
  b=(ln2[1:2]-ln1[1:2])''

  eps=0
  x=A\b

  if x[2]<=ln2[5]+eps||x[2]>=ln2[6]-eps||x[1]<=ln1[5]+eps||x[1]>=ln1[6]-eps
    return [NaN,NaN]
  else
    return x
  end
end

#
function getline(pt1,pt2,i=0,j=0)
  center=(pt1+pt2)/2
  dir=pt2-pt1
  len=sqrt(dir*(dir'))[1]
  dir=dir/len
  dir=[-dir[2],dir[1]]
  ln=[center[1] center[2] dir[1] dir[2] -Inf Inf i j]
end

# use outer vars:
# borderrect
# sites
function voron(sites)
  say()
  say("== voron() ==")
  M=3
  # sort along x axis
  p=sortperm(sites[:,1])
  sites=[sites[p[i],j] for i=1:size(sites,1),j=1:size(sites,2)]

  # init pts arrray and lns array
  pts=sites[1,:]
  lns=Array(Float64,0,8)

  # process points
  for i=2:min(M,size(sites,1))
    say("i=$i")

    # find nearest pt
    diff=broadcast(-,sites[1:i-1,:],sites[i,:])
    dist=sum(diff.*diff,2)
    ind=indmin(dist)

    # get first line
    ln=getline(sites[i,:],sites[ind,:],i,ind)
    say("$(size(lns,1)) lines to test")
    #say(lns)
    ts0=[intersectlns(ln,lns[k,:]) for k=1:size(lns,1)]
    ts=[ts0[k][1] for k=1:size(lns,1)]
    #sortts=sort(ts)
    negs=ts[ts.<0]
    if length(negs)>0
      ln[5]=max(ln[5],maximum(negs))
    end
    pos=ts[ts.>0]
    if length(pos)>0
      ln[6]=min(ln[6],minimum(pos))
    end
    lns=[lns;ln]

    # travel for new circle
    t=1
    travel=[ind]
    travelstart=[]

    if ln[5]!=-Inf
      travelstart=[5]
    elseif ln[6]!=Inf
      travelstart=[6]
    else
      travelstart=[0]
    end

    # find next intersecting line
    currentline=size(lns,1)
    currentpoint=ind
    if travelstart[t]!=0&&length(ts)>0
      nextline=find(ts.==ln[travelstart[t]])
      if length(nextline)!=1
        say("************nextline: $nextline*********************")
      end
      nextline=nextline[1]
      say("currentline: $currentline, nextline: $nextline")
      nextpoint=Int64(lns[nextline,7]==currentpoint?lns[nextline,8]:lns[nextline,7])
      say("currentpoint: $currentpoint, nextpoint: $nextpoint")
      #pts[nextpoint,:]

      # draw newline
      newln=getline(sites[i,:],sites[nextpoint,:],i,nextpoint)
      ts=[intersectlns(newln,lns[k,:])[1] for k=1:size(lns,1)]
      # disable current line and next line
      ts[currentline]=NaN
      ts[nextline]=NaN
      tnextln=intersectlns(newln,lns[nextline,:])[1]
      say(tnextln)

      #compute direction
      mid2currentp=sites[currentpoint,1:2][:]-newln[1:2][:]
      sum(mid2currentp.*newln[3:4])>0?newln[6]=tnextln:newln[5]=tnextln


      lns=[lns;newln]

      # trucate nextline
      b=ts0[nextline][2]#intersectlns(lns[currentline,:],lns[nextline,:])
      mid2p=sites[i,1:2]-lns[nextline,1:2]
      sum(mid2p.*lns[nextline,3:4])>0?lns[nextline,6]=b:lns[nextline,5]=b
      #b>0?lns[nextline,6]=b:lns[nextline,5]=b




      #say(b)

    end

    say()
    say("begin travel, first stop: $(travel[t])")

    say("ln[5]=$(ln[5]), ln[6]=$(ln[6])")

    say()


  end

  lns
end

N=3
sites=rand(N,2)
say()
say("= sites =")
say(size(sites))

# get bounds
say()
say("= bounds =")
top=maximum(sites[:,2])
bottom=minimum(sites[:,2])
left=minimum(sites[:,1])
right=maximum(sites[:,1])
say("top:$top bottom:$bottom left:$left right:$right")

width=right-left
height=top-bottom

margin=0.1
b_top=top+margin*height
b_bot=bottom-margin*height
b_left=left-margin*width
b_right=right+margin*width
b_width=width*(1+2*margin)
b_height=height*(1+2*margin)

# lines: (x0,y0, vx,vy, a,b)
# add borders
borderrect=[
b_left b_bot 1 0 0 b_width
b_right b_bot 0 1 0 b_height
b_right b_top -1 0 0 b_width
b_left b_top 0 -1 0 b_height
]

say()
say("= lines: x0,y0, vx,vy, a,b =")
say(borderrect)

lns=voron(sites)
drawvor(sites,lns[:,1:6],borderrect)
