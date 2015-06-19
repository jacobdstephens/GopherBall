Shader "world space" {
   Properties {

   }
 
   SubShader {
      Pass {
         GLSLPROGRAM
 

 
         #include "UnityCG.glslinc" 
            // defines _Object2World and _World2Object
        uniform vec4 unity_Scale;
        uniform mat4 UNITY_MATRIX_IT_MV;
        uniform mat4 UNITY_MATRIX_V;
        uniform mat4 UNITY_MATRIX_VP;
        uniform mat4 UNITY_MATRIX_P;
        
            
		//uniform vec3 _WorldSpaceCameraPos;


 		 
         varying vec4 position_in_world_space;
         varying vec4 position;
         
 
         #ifdef VERTEX
 
         void main()
         {
            mat4 modelMatrix = _Object2World;
            mat4 worldlMatrix = _World2Object;
            
 			position_in_world_space = modelMatrix * gl_Vertex;
 			position = gl_ModelViewProjectionMatrix * gl_Vertex;
 
            gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex; 
         }
 
         #endif
 
         #ifdef FRAGMENT
         
         #define MOD2 vec2(.16632,.17369)
         
         //Convert camera and Point Position to Camera Space
         vec4 camSpacePos = UNITY_MATRIX_V * position_in_world_space;
         vec4 look = UNITY_MATRIX_V * vec4(.0,.0,.0,1.0);
         vec4 camSpace = UNITY_MATRIX_V * vec4(_WorldSpaceCameraPos,1.0);
         
         float sdSphere( vec3 p, float s )
		{
    		return length(p) -s;
		}
         
		float Hash(vec2 p)
		{
			p  = fract(p / MOD2);
   		    p += dot(p.xy, p.yx+19.19);
    		return fract(p.x * p.y);
		}
		
		float Noise( in vec2 x )
		{
    		vec2 p = floor(x);
    		vec2 f = fract(x);
    		f = f*f*(3.0-2.0*f);
    
    		float res = mix(mix( Hash(p), Hash(p+ vec2(1.0, 0.0)),f.x),
                    		mix( Hash(p+ vec2(.0, 1.0)), Hash(p+ vec2(1.0, 1.0)),f.x),f.y);
    		return res;
		}
         
		vec2 Voronoi( in vec2 x )
		{
			vec2 p = floor( x );
			vec2 f = fract( x );
			float res=100.0,id;
			for( int j=-1; j<=1; j++ )
			for( int i=-1; i<=1; i++ )
			{
				vec2 b = vec2( float(i), float(j) );
				vec2 r = b - f  + Hash( p + b );
				float d = dot(r,r);
				if( d < res )
				{
					res = d;
					id  = Hash(p+b);
				}			
   		 }
			return vec2(max(.4-sqrt(res), 0.0),id);
		}
		
		
       //Raymarching   
		vec3 DE(vec3 p)
		{
			//Signed Distance Sphere
			
			//float radius = distance(p.xyz,origin.xyz);
			float sdSphere1 = sdSphere(vec3(.0,.0,.0),1.);
			
			vec3 up = cross(p,vec3(.0,-1.0,.0));

			
			float y = .3;//sdSphere1;//distance(p.xyz,rayOrigin.xyz)*1.;// - base -height;
			y = y*y;
			
			//float z = p.z - base -height;
			//z = z*z;
			//y = z;
			
			//Pass in a 2D vector to compute the Voronoi plane
			vec2 ret = Voronoi(p.xz);
			//vec2 ret = Voronoi(vec2(0.0,1.0));
			//vec2 ret = Voronoi((p.xy*2.5+sin(y*4.0+p.yx*12.3)*.12+vec2(sin(_Time[1]*2.3+1.5*p.y),sin(_Time[1]*3.6+1.5*p.x))*y*.5));
			
			//Ret X is the height of the peaks and Y is the midline
			//float f = ret.x;
			float f = ret.x * .6 + y * .58;

			//return vec3(y,.0,.0);
			return vec3( y - f * 1.4, clamp(f * 1.5, 0.0, 1.0), ret.y);
		} 
		
		//Circle of Confusion fall off              
		float CircleOfConfusion(float t)
		{
		return max(t * .04, (2.0 / _ScreenParams.y) * (1.0+t));
		}
		
		
		// Linear Step function
		float Linstep(float a, float b, float t)
		{
			return clamp((t-a)/(b-a),0.,1.);
		}
		

		
		vec4 GrassBlades(in vec3 rO, in vec3 rD, in vec3 mat, in float dist)
		{
			float d = 0.0;
			float f;
			// Only calculate cCoC once is enough here...
			float rCoC = CircleOfConfusion(dist*.3);
			float alpha = 0.0;
	
			vec4 col = vec4(mat*0.15, 0.0);

			for (int i = 0; i < 15; i++)
			{
				if (col.w > .99) break;
				
				//Find Point on surface
				vec3 p = rO + rD * d;
				//vec3 p = position_in_world_space.xyz * d;
				
				
				//Get vector from the ray marching
				vec3 ret = DE(p);
				
				ret.x += .5 * rCoC;

				if (ret.x < rCoC)
				{
					alpha = (1.0 - col.y) * Linstep(-rCoC, rCoC, -ret.x) * 2.0;//calculate the mix like cloud density
					f = clamp(ret.y, 0.0, 1.0);
					// Mix material with white tips for grass...
					vec3 gra = mix(mat, vec3(.2, .3, min(pow(ret.z, 2.0)*3.0, .2)), pow(ret.y,20.0)*.6) * ret.y;
					col += vec4(gra * alpha, alpha);
				}
				d += max(ret.x * .7, .2);
			}
			if(col.w < .2)col.xyzw = vec4(.1, .15, .05,.0);
			return col.xyzw;
		}
		


 
         void main()
         {

            
            //Get Camera Direction
            vec3 dir = normalize(position_in_world_space.xyz-_WorldSpaceCameraPos.xyz);
            //vec3 dir = normalize(position_in_world_space.xyz-camSpace.xyz);
            vec4 rayOrigin = UNITY_MATRIX_IT_MV * vec4(position_in_world_space.xyz,1.0);
            //vec3 eye = _WorldSpaceCameraPos;
            vec4 rayDirection = UNITY_MATRIX_V * vec4(_WorldSpaceCameraPos.xyz,1.0);
            //Create default material
            vec3 mat = mix(vec3(.0,.3,.0), vec3(.2,.3,.0), Noise(position_in_world_space.xy*.025));
            
            
            //Create Distance Float
            float dist;
            
            //Render Grass 
            vec3 norm = cross(cross(rayDirection.xyz,vec3(.0,1.0,.0)), rayDirection.xyz);
            //rayDirection.xyz = reflect(rayDirection.xyz,vec3(1.0,.0,.0));
            //rayOrigin.xyz = cross(cross(rayOrigin.xyz,vec3(.0,-1.0,.0)), rayOrigin.xyz);
            //eye.xyz = reflect(eye.xyz,norm);
            //viewDir.xyz = cross(cross(viewDir.xyz,vec3(1.0,.0,.0)), viewDir.xyz);
            
            vec4 col = GrassBlades( position_in_world_space.xyz, rayOrigin.xyz , mat, dist);

			gl_FragColor = col;
            
         }
 
         #endif
 
         ENDGLSL
      }
   }
}