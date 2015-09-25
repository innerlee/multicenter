say("voronoi now")

#
function composepoints(pts,labels=[])
  if length(labels)==0
    compose(context(),circle(pts[:,1],pts[:,2],[.7mm]),fill([colorant"black"]), stroke([colorant"white"]))
  else
    dc=distinguishable_colors(10)
    #println(labels)
    colors=[dc[3+Int64(labels[k])] for k=1:length(labels)]
    compose(context(),circle(pts[:,1],pts[:,2],[.7mm]),fill(colors), stroke([colorant"white"]))
  end
end

#
function transformlines(lns)
  [[(lns[i,1]+lns[i,5]*lns[i,3],lns[i,2]+lns[i,5]*lns[i,4]),(lns[i,1]+lns[i,6]*lns[i,3],lns[i,2]+lns[i,6]*lns[i,4])] for i=1:size(lns,1),j=1]
end

#
function composelines(lns::AbstractArray,labels=[])
  say()
  say("= in composelines() =")
  ls=transformlines(lns)
  #say(ls)
  #lines=[[(ls[i,1]+ls[i,5]*ls[i,3],ls[i,2]+ls[i,5]*ls[i,4]),(ls[i,1]+ls[i,6]*ls[i,3],ls[i,2]+ls[i,6]*ls[i,4])] for i=1:size(ls,1),j=1]
  if length(labels)==0
    compose(context(), line(ls),fill([colorant"gray"]), stroke([colorant"gray"]))
  else
    colors=[lns[k,11]==0?colorant"gray":colorant"red" for k=1:size(lns,1)]
    #println("in composelines:\n $labels")
    #println("size of lns:$(size(lns,1)), size of labels:$(size(labels))")
    compose(context(), line(ls),fill(colors), stroke(colors))
  end
end

#
function composelines(ls)
  compose(context(), line(ls),fill([colorant"black"]), stroke([colorant"black"]))
end

# lns: n*6 (x0,y0, vx,vy, a,b)
# ls: n*4 ((x0,y0) (x1,y1))
# use borderrect
function drawvor(pts,lns,border,vorptlables=[])
  #println(size(lns))
  # preprocess them
  say()
  say("= drawvor =")
  say(lns[:,7:10])
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
#println("aa")
  # compute labels
  labels=[]
#  println(vorptlables)
  if length(vorptlables)>0
    #templab=[vorptlables[Int64(lns[k,7])]==vorptlables[Int64(lns[k,8])]?0:1 for
    #templab=lns[:,11]
    templab=zeros(size(lns,1))
    templab[7]=1
    #println("templab:\n  $templab")
    #k=1:size(lns,1)]
  #  println(size(templab''))
    labels=[[0 0 0 0]'; templab'']
  end
#  println("bb")
  # fileter outer lines
  lns=lns[[isout(lns[k,:]) for k=1:size(lns,1)],:]
  # transforms
  #println("size boder: $(size(tempborder[1:size(border,1),1:size(border,2)]))")
  tempborder=zeros(size(border,1),size(lns,2))
  # for k=1:size(border,1)
  #   for l=1:size(border,2)
  #     tempborder[k,l]=border[k,l]
  #   end
  # end
  tempborder[1:size(border,1),1:size(border,2)]=border[1:size(border,1),1:size(border,2)]
  lns=[tempborder;lns]
  ratio=max(b_width,b_height)
  #say(ratio)
  #say(broadcast(-,pts,[b_left b_bot])/6)
  pts=broadcast(-,pts,[b_left b_bot])/ratio
  #pts=broadcast(/,broadcast(-,pts,[b_left b_bot]),[ratio ratio])
#say("...")
  println("safe")
  lns=broadcast(/,broadcast(-,lns,[[b_left b_bot] zeros(1,size(lns,2)-2)]),[[ratio ratio 1 1 ratio ratio] ones(1,size(lns,2)-6)])
  say(size(lns))
println("cc")
  drawcomp(compose(context(),composepoints(pts,vorptlables),composelines(lns,labels)))
end

#
function isout(ln)
  #return true
  #say()
  #say("= isout =")
  #say(ln[5:6])
  twopts=[(ln[1:2]+ln[5]*ln[3:4])'; (ln[1:2]+ln[6]*ln[3:4])']
  #say(twopts)
  ep=1e-5
  for k=1:size(twopts,1)
    x=twopts[k,1]
    y=twopts[k,2]
    if x<b_left-ep||x>b_right+ep||y<b_bot-ep||y>b_top+ep
      return false
    end
  end
  return true
end

#
function intersectlns(ln1,ln2)
  A=[ln1[3] -ln2[3]
     ln1[4] -ln2[4]]
  b=(ln2[1:2]-ln1[1:2])''

  #say(A)
  #(det(A))
  if det(A)==0
    return [NaN,NaN]
  end
  eps=0
  x=A\b

  if x[2]<=ln2[5]+eps||x[2]>=ln2[6]-eps||x[1]<=ln1[5]+eps||x[1]>=ln1[6]-eps
    return [NaN,NaN]
  else
    return x
  end
end

#
function detachvor(vpts,ptno,lnno)
  if isnan(ptno)
    return
  end
  ptno=Int64(ptno)
  lnno=Int64(lnno)
  vpts[ptno,3]==lnno?vpts[ptno,3]=NaN:
  vpts[ptno,4]==lnno?vpts[ptno,4]=NaN:
  vpts[ptno,5]==lnno?vpts[ptno,5]=NaN:Nothing
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
function voron(sites,labels=[])
  say()
  #println("voron get labels:\n $labels")
  say("== voron() ==")
  M=size(sites,1)
  # perturb sites
  ep=1e-12
  #println(size(sites))
  sites=sites[:,1:2]+ep*rand(M,2)
  # sort along x axis
  p=sortperm(sites[:,1])
  #println("sortperm:\n $p")
  #println(p)
  sites=[sites[p[i],j] for i=1:size(sites,1),j=1:size(sites,2)]
  #sites[:,1:2]=aftersalt

  #####################
  sites=sites[1:min(M,NN),:]

  # init pts arrray and lns array
  pts=sites[1,:]
  # x,y vx,vy lbound,ubound pt1,pt2 vorpt1,vorpt2 isdeleted
  lns=Array(Float64,0,11)
  # x,y ln1,ln2,ln3 isdeleted
  vorpts=Array(Float64,0,6)

  # process points
  for i=2:min(M,size(sites,1))
    debug=(i>5)
    say("--- new point adds in: i=$i")
say("(a)")
    # find nearest pt
    diff=broadcast(-,sites[1:i-1,:],sites[i,:])
    dist=sum(diff.*diff,2)
    ind=indmin(dist)
    say(ind)
say("(b)")
    # get first line
    ln=[getline(sites[i,:],sites[ind,:],i,ind) [NaN NaN 0]]
    say("$(size(lns,1)) lines to test")
    ts0=[lns[k,11]==0?intersectlns(ln,lns[k,:])[1]:NaN for k=1:size(lns,1)]
    negs=ts0[ts0.<0]
    say(ts0)
    leftright=0
    if length(negs)>0
      say("will left")
      leftright+=1
      ln[5]=max(ln[5],maximum(negs))
    end
    pos=ts0[ts0.>0]
    if length(pos)>0
      say("will right")
      leftright+=1
      ln[6]=min(ln[6],minimum(pos))
    end
    say("ln[5]=$(ln[5]),ln[6]=$(ln[6])")
    say(size(lns))
    say(size(ln))
    lns=[lns;ln]

    currentline0=size(lns,1)
    say("currentline:$currentline0")
say("(c)")
    # travel for new circle
    MAX=10000
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
      #say("ln[z]=$(ln[z])")
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
      say(travelstart)
      say(t)
      nextline=find(ts0.==ln[travelstart[t]])
      if length(nextline)!=1
        say("************nextline: $nextline*********************")
      end
      nextline=nextline[1]
      nextpoint=Int64(lns[nextline,7]==currentpoint?lns[nextline,8]:lns[nextline,7])
say("+(g1)")
      # vor pt of currentline
      z==5?lns[currentline,9]=size(vorpts,1)+1:lns[currentline,10]=size(vorpts,1)+1
say("(h)")
      while t<MAX
        # add to vorpts the intersection of current line and next line and new line
say("+(h1)")
        newvorpt=intersectpt(lns[currentline,:],lns[nextline,:])
        vorpts=[vorpts;[newvorpt[1] newvorpt[2] currentline nextline size(lns,1)+1 0]]
        say(vorpts[end,3:end])
say("+(h2)")
        say(size(vorpts))
        say("[$t] currentline: $currentline, nextline: $nextline | currentpoint: $currentpoint, nextpoint: $nextpoint")
say("(i)")
        # draw newline passing between nextpoint and i-th point
        newline=[getline(sites[i,:],sites[nextpoint,:],i,nextpoint) [NaN NaN 0]]
        # determine direction of new line
        mid2p=sites[currentpoint,1:2][:]-newline[1:2][:]
        codir=sum(mid2p.*newline[3:4])>0
        say("codir:$codir")
say("(j)")
        # filter edges around next point
        filt=[lns[k,11]==0&&(lns[k,7]==nextpoint||lns[k,8]==nextpoint) for k=1:size(lns,1)]
        ts=[filt[k]?intersectlns(newline,lns[k,:])[1]:NaN for k=1:size(lns,1)]
        say(ts)
say("(k)")
        # trim new line
        a=ts[nextline]
        b=NaN
        #say("a=$a")
        if codir
          newline[6]=a
          # vor pt of new line: attach new
          newline[10]=size(vorpts,1)
          filterts=[ts[k]<a?ts[k]:NaN for k=1:length(ts)]
          if length(filterts)>0
            b=maximum(filterts)
            if b!=Inf&&b!=-Inf&&!isnan(b)
              newline[5]=b
              # vor pt 2 of new line: attach new
              newline[9]=size(vorpts,1)+1
            end
          end
        else
say("+(k0)")
          newline[5]=a
          # vor pt of new line: attach new
          newline[9]=size(vorpts,1)
          say(vorpts)
          filterts=[ts[k]>a?ts[k]:NaN for k=1:length(ts)]
          say(filterts)
          if length(filterts)>0
            b=minimum(filterts)
            if b!=Inf&&b!=-Inf&&!isnan(b)
              newline[6]=b
              # vor pt 2 of new line: attach new
              newline[10]=size(vorpts,1)+1
            end
          end
        end
say("+(k1)")
        say(newline[7:11])
        lns=[lns;newline]
say("(l)")
        # trim nextline
        copycurrentline=zeros(1,size(lns,2))
        copycurrentline[:]=[lns[currentline,1:4] [-Inf Inf 0 0 0 0 0]]
        c=intersectlns(copycurrentline,lns[nextline,:])[2]
        say(c)
        mid2p=sites[i,1:2]-lns[nextline,1:2]
        if sum(mid2p.*lns[nextline,3:4])>0
say("+(l0)")
          lns[nextline,6]=c
          # vor point of next line and next line 2
          say(nextline)
          say(lns[nextline,7:11])
          # detach old vorpt
          detachvor(vorpts,lns[nextline,10],nextline)
          #removevor=lns[nextline,10]
say("+(l1)")
          #say(lns[nextline,7:11])
          # attach new vorpt
          lns[nextline,10]=size(vorpts,1)
        else
          lns[nextline,5]=c
          # detach old vorpt
          detachvor(vorpts,lns[nextline,9],nextline)
          #removevor=lns[nextline,9]
say("+(l2)")
          say(nextline)
          say(lns[nextline,7:11])
          # attach new vorpt
          lns[nextline,9]=size(vorpts,1)
        end
say("(m)")
        # whether ends
        if isnan(b)||b==Inf||b==-Inf

          # TODO: CLEAR OUT ALL DISCONECTED VOR PTS AND LINES
          # check wheter start cleaning
          flag=z==6||(z==5&&leftright==1)
          # clean up isolated lines and vor pts
          while flag
            flag=false
            say("==^^")
            say(vorpts)
            for k=size(vorpts,1):-1:1
              say("^$k")
              sig=sum(isnan(vorpts[k,3:5]))
              if sig==1||sig==2
                flag=true
                say("here, vor pt $k has $(3-sig) lines to remove ^^^^^")
                say("vorpts $k before")
                say(vorpts[k,3:end])
                vorpts[k,6]=1
                for l=3:5
                  if !isnan(vorpts[k,l])
                    # remove line
                    lno=Int64(vorpts[k,l])
                    lns[lno,11]=1
                    say("xxxxxxxxxxxxxxxxxxxxx remove line $vorpts[k,l]")
                    # detach this line from vor pt
                    detachvor(vorpts,lns[lno,9],lno)
                    detachvor(vorpts,lns[lno,10],lno)
                    #vorpts[k,l]=NaN
                  end
                end
              end
            end
          end
          say("end left or right $z")
          break
        end
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
        say()
        say(lns[:,7:11])
        say()
        say("nextline: $nextline, nextline2: $nextline2")
#         # delete isolated lines
#         say("removevor=$removevor")
#         say(vorpts)
# say("+(n21)")
#         #say(lns[:,7:11])
#         if !isnan(removevor)&&vorpts[Int64(removevor),6]==0
#           removevor=Int64(removevor)
#           vorpts[removevor,6]=1
#           removeln=Int64(setdiff(vorpts[removevor,3:5],[nextline nextline2])[1])
#           lns[removeln,11]=1
#           say("remove line: $removeln")
#         end
        # say(vorpts)
        # say(lns[nextline,9:10])
        # say(lns[nextline2,9:10])
        # vpt=setdiff(lns[nextline,9:10],setdiff(lns[nextline,9:10],lns[nextline2,9:10]))[1]
        # say(vpt)
say("+(n22)")
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

        #t+=1
say("(p)")
      end
say("(q)")
    end
say("(r)")
  end
say("(s)")
#println("lines before:\n $(lns[:,7:11])")
say(lns[:,7:11])
  templns=lns[find(lns[:,11].==0),:]
  #println(templns[1:10,7:8])
  for k=1:size(templns,1)
    templns[k,7]=p[Int64(templns[k,7])]
    templns[k,8]=p[Int64(templns[k,8])]
    templns[k,11]=(labels[Int64(templns[k,7])]==labels[Int64(templns[k,8])]?0:1)
  end
  #println("lines after:\n $(templns[:,7:11])")
  #println(templns[1:10,7:8])
  #println([labels[Int64(templns[k,l])] for k=1:size(templns,1),l=7:8])
  #println(p)
  templns
  #lns
end
