isdefined(:debug)?Nothing:debug=true

function say(something=" ")
  debug?println(something):Nothing
end

using Gadfly
using Compose

function drawpic(data::AbstractArray)
  if size(data,2)==2
    draw(SVG("temp.svg",10inch,6inch),plot(x=data[:,1],y=data[:,2],Coord.cartesian(fixed=true,aspect_ratio=1.0)))
  elseif size(data,2)==3
    draw(SVG("temp.svg",10inch,6inch),plot(x=data[:,1],y=data[:,2],color=map(x->string(round(x)),data[:,3]),Geom.point,Coord.cartesian(fixed=true,aspect_ratio=1.0)))
  end

  run(`"c:\\program files (x86)\\google\\chrome\\application\\chrome.exe" $(pwd())\\temp.svg`)
end

function showpic()
  run(`"c:\\program files (x86)\\google\\chrome\\application\\chrome.exe" $(pwd())\\temp.svg`)
end

function drawcomp(comp)
  img = SVG("temp.svg", 6inch, 6inch)
  draw(img,compose(context(),comp))
  showpic()
end

function points(pts)
  compose(context(),circle(pts[:,1],pts[:,2],[.7mm]),fill(["black"]), stroke(["white"]))
end
