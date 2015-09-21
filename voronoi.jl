say("voronoi now")

function drawlines(ls)
  compose(context(), line([(1,1), (0,1)]),circle(.5, .5, .3))

  lines=[[(ls[i,1]+ls[i,5]*ls[i,3],ls[i,2]+ls[i,5]*ls[i,4]),(ls[i,1]+ls[i,6]*ls[i,3],ls[i,2]+ls[i,6]*ls[i,4])] for i=1:size(ls,1),j=1]

say()
say("in drawlines()")
say(lines)
  compose(context(), line(lines),fill(["black"]), stroke(["black"])),lines
end


N=2
sites=rand(N,2)
say()
say("= sites =")
say(sites)

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


# lines: x0,y0, vx,vy, a,b)
# add borders
borderrect=[
b_left b_bot 1 0 0 width*(1+2*margin)
b_right b_bot 0 1 0 height*(1+2*margin)
b_right b_top -1 0 0 width*(1+2*margin)
b_left b_top 0 -1 0 height*(1+2*margin)
]
say()
say("= lines: x0,y0, vx,vy, a,b =")
say(borderrect)

ls,li=drawlines(borderrect)
drawcomp(ls)
