#pragma strict

var animator: Animator;
var isHitHash: int = Animator.StringToHash("isHit");

function gopherExit () 
{
	animator = GetComponent("Animator");
	animator.SetTrigger(isHitHash);

}

