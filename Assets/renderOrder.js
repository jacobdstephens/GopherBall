#pragma strict

function Start(){
  // get all renderers in this object and its children:
  var renders = GetComponentsInChildren(Renderer);
  for (var rendr: Renderer in renders){
    rendr.material.renderQueue = 2002; // set their renderQueue
  }
}