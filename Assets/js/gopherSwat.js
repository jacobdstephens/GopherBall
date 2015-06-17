#pragma strict
var numHits = 0 ; 
var text: UnityEngine.UI.Text;
var dyn: Rigidbody;
var animator: Animator;
var isHitHash: int = Animator.StringToHash("isHit");

function Start () 
{
    animator = GetComponent("Animator");
    dyn = GetComponent.<Rigidbody>();
}


function Update () 
{

	var hitDetect : boolean = false ;

	if (Input.GetKey ("escape"))
	{
		Application.Quit();
	}
	
	if (Input.GetMouseButtonDown(0)) 
	{
		var hit: RaycastHit;
		var ray = Camera.main.ScreenPointToRay(Input.mousePosition) ;
		
		
		
		if (Physics.Raycast(ray,hit)) 
		{
			if ( hit.collider.tag == "Gopher") 
			{
				//hit.collider.enabled = false;
				numHits++ ;
				hitDetect = true ;
				//Destroy(gameObject);
				animator.SetTrigger(isHitHash);
				//TODO add trigger for rigidBody
				//dyn.velocity = Vector3(0,1,0);
				
				
				text.text=numHits.ToString();
				Debug.Log(hitDetect) ;
				Debug.Log(numHits);
				
			}
		}
	}
}