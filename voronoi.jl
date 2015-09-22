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

#
function intersectpt(ln1,ln2)
  A=[ln1[3] -ln2[3]
     ln1[4] -ln2[4]]
  b=(ln2[1:2]-ln1[1:2])''
  x1=(A\b)[1]
  ln1[1:2]+x1*ln1[3:4]
end

# use outer vars:
# borderrect
# sites
function voron(sites)
  say()
  say("== voron() ==")
  M=size(sites,1)
  # sort along x axis
  p=sortperm(sites[:,1])
  sites=[sites[p[i],j] for i=1:size(sites,1),j=1:size(sites,2)]

  # init pts arrray and lns array
  pts=sites[1,:]
  lns=Array(Float64,0,9)
  # x,y ln1,ln2,ln3
  vorpts=Array(Float64,0,5)

  # process points
  for i=2:min(M,size(sites,1))
    say("--- new point adds in: i=$i")
say("(a)")
    # find nearest pt
    diff=broadcast(-,sites[1:i-1,:],sites[i,:])
    dist=sum(diff.*diff,2)
    ind=indmin(dist)
say("(b)")
    # get first line
    ln=[getline(sites[i,:],sites[ind,:],i,ind) [0]]
    say("$(size(lns,1)) lines to test")
    ts0=[intersectlns(ln,lns[k,:])[1] for k=1:size(lns,1)]
    negs=ts0[ts0.<0]
    if length(negs)>0
      say("will left")
      ln[5]=max(ln[5],maximum(negs))
    end
    pos=ts0[ts0.>0]
    if length(pos)>0
      say("will right")
      ln[6]=min(ln[6],minimum(pos))
    end
    say("ln[5]=$(ln[5]),ln[6]=$(ln[6])")
    say(size(lns))
    say(size(ln))
    lns=[lns;ln]

    currentline0=size(lns,1)

say("(c)")
    # travel for new circle
    MAX=10
    travel=[ind]
    t=1
    travelstart=[]
say("(d)")
    # go left and go right until unbound
    for z=5:6
      say("z=$z")
      # no existing lines, nothing to do
      length(ts0)==0?break:Nothing
say("(e)")
      say("ln[z]=$(ln[z])")
      if ln[z]!=-Inf&&ln[z]!=Inf
        travelstart=[z]
      else
        continue
      end
say("(f)")
      currentline=currentline0
      currentpoint=ind
say("(g)")
      # find next intersecting line
      nextline=find(ts0.==ln[travelstart[t]])
      if length(nextline)!=1
        say("************nextline: $nextline*********************")
      end
      nextline=nextline[1]
      nextpoint=Int64(lns[nextline,7]==currentpoint?lns[nextline,8]:lns[nextline,7])
say("(h)")
      while t<MAX
        # add to vorpts the intersection of current line and next line and new line
say("+(h1)")
        newvorpt=intersectpt(lns[currentline,:],lns[nextline,:])
        say(newvorpt)
        vorpts=[vorpts;[newvorpt[1] newvorpt[2] currentline nextline size(lns,1)+1]]
say("+(h2)")
        say(size(vorpts))
        say("[$t] currentline: $currentline, nextline: $nextline | currentpoint: $currentpoint, nextpoint: $nextpoint")
say("(i)")
        # draw newline passing between nextpoint and i-th point
        newline=[getline(sites[i,:],sites[nextpoint,:],i,nextpoint) [0]]
        # determine direction of new line
        mid2p=sites[currentpoint,1:2][:]-newline[1:2][:]
        codir=sum(mid2p.*newline[3:4])>0
say("(j)")
        # filter edges around next point
        filt=[lns[k,7]==nextpoint||lns[k,8]==nextpoint for k=1:size(lns,1)]
        ts=[filt[k]?intersectlns(newline,lns[k,:])[1]:NaN for k=1:size(lns,1)]
say("(k)")
        # trim new line
        a=ts[nextline]
        b=NaN
        say("a=$a")
        if codir
          newline[6]=a
          filterts=[ts[k]<a?ts[k]:NaN for k=1:length(ts)]
          if length(filterts)>0
            b=maximum(filterts)
            b!=Inf&&b!=-Inf&&b!=NaN?newline[5]=b:Nothing
          end
        else
          newline[5]=a
          filterts=[ts[k]>a?ts[k]:NaN for k=1:length(ts)]
          if length(filterts)>0
            b=minimum(filterts)
            b!=Inf&&b!=-Inf&&b!=NaN?newline[6]=b:Nothing
          end
        end
        lns=[lns;newline]
say("(l)")
        # trim nextline
        copycurrentline=zeros(1,size(lns,2))
        copycurrentline[:]=[lns[currentline,1:4] [-Inf Inf 0 0 0]]
        c=intersectlns(copycurrentline,lns[nextline,:])[2]
        say(c)
        mid2p=sites[i,1:2]-lns[nextline,1:2]
        #say(mid2p)
        sum(mid2p.*lns[nextline,3:4])>0?lns[nextline,6]=c:lns[nextline,5]=c
say("(m)")
        # whether ends
        isnan(b)||b==Inf||b==-Inf?break:Nothing
say("(n)")
        # now continues,
        # get the new next line which new line intersects with
        nextline2=find(filterts.==b)
        if length(nextline2)>1
          say("***********HERE************")
        end
        say(nextline2)
        say(filterts)
        say(b)
say("+(n1)")
        # delete isolated lines

        # prepare for the next round
        currentline=size(lns,1)
        currentpoint=nextpoint
        nextline=nextline2[1]
        nextpoint=Int64(lns[nextline,7]==currentpoint?lns[nextline,8]:lns[nextline,7])
say("(o)")


        # say()
        # say("begin travel, first stop: $(travel[t])")
        #
        # say("ln[5]=$(ln[5]), ln[6]=$(ln[6])")

        t+=1
say("(p)")
      end
say("(q)")
    end
say("(r)")
  end
say("(s)")
  lns
end

N=4
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

vorlns=voron(sites)
drawvor(sites,vorlns[:,1:6],borderrect)
