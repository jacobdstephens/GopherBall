
var emitter: Transform;

private var timer : float;
private var gopher: Transform;

function Awake () 
{
	timer = Time.time + 5;
	transform.Rotate(Vector3.right, Random.Range( 10.0, 60.0 ));
	transform.Rotate(Vector3.up, Random.Range( 0.0, 360.0 ));
    	
}

function Start()
{
	
}

var randomAngle;


function Update () 
{
	
	
	//if ( timer < Time.time)
	//{
	//	gopher = Instantiate( emitter , transform.position , Quaternion.identity );
	//	gopher.Rotate(Vector3.right, Random.Range( 10.0, 60.0 ));
	//	gopher.Rotate(Vector3.up, Random.Range( 0.0, 360.0 ));
		
	//	timer = Time.time + 5;
	//}
}