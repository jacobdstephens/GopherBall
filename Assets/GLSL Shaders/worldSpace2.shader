Shader "GLSL shading in world space" {
   SubShader {
      Pass {
         GLSLPROGRAM
 
         uniform mat4 _Object2World; 
            // definition of a Unity-specific uniform variable 
         uniform vec3 _WorldSpaceCameraPos;
 
         varying vec4 position_in_world_space;
 
         #ifdef VERTEX
 
         void main()
         {
            position_in_world_space = _Object2World * gl_Vertex;
               // transformation of gl_Vertex from object coordinates 
               // to world coordinates;
 
            gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
         }
 
         #endif
 
         #ifdef FRAGMENT
 
         void main()
         {
            float dist = distance(position_in_world_space, 
               vec4( _WorldSpaceCameraPos, 1.0));
               // computes the distance between the fragment position 
               // and the origin (the 4th coordinate should always be 
               // 1 for points).
 
            if (dist < 15.0)
            {
               gl_FragColor = vec4(0.0, 1.0, 0.0, 1.0); 
                  // color near origin
            }
            else
            {
               gl_FragColor = vec4(0.3, 0.3, 0.3, 1.0); 
                  // color far from origin
            }
         }
 
         #endif
 
         ENDGLSL
      }
   }
}