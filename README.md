# GopherBall

###Screen Space Sphereical Shader 

Gopher Ball is a visual experiment to integrate Shader Toy Screen Space techniques into Unity3D. 

![alt tag](https://raw.github.com/jacobdstephens/GopherBall/master/GopherBall.png)

An anti-gravity gopher hive is floating there in front of you...just floatin'. You're not going to stand for that! Start clicking on the "gophers" to get rid of them!

##The Fun Stuff
The grass shader is a screen space technique that has been adapted to fit into Unity's World Space as Lighting Space. 
Unlike normal maps, this technique is able to take camera matricies into account to support VR.

I was able to port Inigo Quilez's Signed Distance Sphere functions (https://www.shadertoy.com/view/Xds3zN) and the Rolling Hills by David Hoskins (https://www.shadertoy.com/view/Xsf3zX) to Unity to create a continuous grass texture, with perspecive and responsive to parallax. 

##The Hard Stuff
Tangent Space. The real trick here is two fold: first using a signed distance function to describe a sphere in screen space instead of a 2D terrain, second (and more critically) normals of that signed distance function must be converted into tangent space (looking from the point to the camera) to bend the horizon to the outside edge of the sphere. One last cheap trick was used to compuete the the bearing of the camera to the sphere and counter animate the grass accordingly to give the illusion of orbiting around the sphere. With out this the perspective illlusion breaks and the camera becomes perpendicular to the grass blades. 
