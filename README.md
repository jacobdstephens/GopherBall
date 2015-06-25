# GopherBall

Screen Space Sphereical Shader 

Gopher Ball is a visual experiment to integrate Shader Toy Screen Space techniques into Unity3D. 

An anti-gravity gopher hive is floating there in front of you...just floatin'. You're not going to stand for that! Start clicking on the "gophers" to get rid of them!

##The Fun stuff:
The grass shader is a screen space technique that has been adapted to fit into Unity's World Space as Lighting Space. 
Unlike normal maps, this technique is able to take camera matricies into account to support VR.

I was able to port Inigo Quilez's Signed Distance Sphere functions (https://www.shadertoy.com/view/Xds3zN) the Rolling Hills by David Hoskins (https://www.shadertoy.com/view/Xsf3zX) to Unity to create a continuous grass texture, with with perspecive and responsive to paralax. 

![alt tag](https://raw.github.com/jacobdstephens/GopherBall/master/GopherBall.png)
