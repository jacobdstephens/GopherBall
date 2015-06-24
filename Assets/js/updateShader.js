#pragma strict

var rend: Renderer;
var sphere: Transform;


function Start () {
	rend = GetComponent.<Renderer>();
	
	// Use the Grass shader on the material
	rend.material.shader = Shader.Find("Grass");
}


function Update () {
	// Animate the Shininess value
	var spherePos = sphere.position;
	rend.material.SetFloat("_X", spherePos.x);
	rend.material.SetFloat("_Y", spherePos.y);
	rend.material.SetFloat("_Z", spherePos.z);
}