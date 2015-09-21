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
  compose(context(), line(ls),fill(["black"]), stroke(["black"]))
end

#
function composelines(ls)
  compose(context(), line(ls),fill(["black"]), stroke(["black"]))
end

# lns: n*6 (x0,y0, vx,vy, a,b)
# ls: n*4 ((x0,y0) (x1,y1))
function drawvor(pts,lns)
  # preprocess them
  say()
  say("== drawvor =")
  ratio=max(b_width,b_height)
  pts=broadcast(/,broadcast(-,pts,[b_left b_bot]),[ratio ratio])
  lns=broadcast(/,broadcast(-,lns,[b_left b_bot 0 0 0 0]),[ratio ratio 1 1 ratio ratio])
  say(lns)

  drawcomp(compose(context(),composepoints(pts),composelines(lns)))
end


N=20
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

drawvor(sites,borderrect)
