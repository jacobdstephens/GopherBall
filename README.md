# GopherBall

Screen Space Sphereical Shader 

Gopher Ball is a visual experiment to integrate Shader Toy Screen Space technigues into Unity3D. 

An anti-gravity gopher hive is floating there in front of you...just floain'. You're not going to stand for that! Start clicking on the "gophers" to get rid of them!

The Fun stuff:
The grass shader is a screen space technique that has been adapted to fit into Unity's World Space as Lighting Space. 
Unlike a normal maps, this technique is able to take camera matricies into account to support VR.

Using a Signed Distance Sphere (https://www.shadertoy.com/view/Xds3zN) as the basis for the Voronoi Plane it is possible to project spherically instead of using an XZ plane (https://www.shadertoy.com/view/Xsf3zX).

![alt tag](https://raw.github.com/jacobdstephens/GopherBall/master/GopherBall.png)
