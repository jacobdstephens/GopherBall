#pragma strict

function randomRotation () {
	transform.parent.transform.rotation = Quaternion.identity;
	transform.parent.transform.Rotate(Vector3.up, Random.Range( 0.0, 360.0 ));
	transform.parent.transform.Rotate(Vector3.right, Random.Range( 10.0, 60.0 ));
	
	
}
