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
        
            
		//uniform vec3 _WorldSpaceCameraPos;


 		 
         varying vec4 position_in_world_space;
         
 
         #ifdef VERTEX
 
         void main()
         {
            mat4 modelMatrix = _Object2World;
            mat4 worldlMatrix = _World2Object;
            
 			position_in_world_space = modelMatrix * gl_Vertex;
 			
 
            gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex; 
         }
 
         #endif
 
         #ifdef FRAGMENT
         
         #define MOD2 vec2(.16632,.17369)
         
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
         
		vec3 DE(vec3 p)
		{
			vec3 up = normalize(cross(cross(p, vec3(0,1,0)), p));
			vec3 newp = up+p;
			float base = distance(position_in_world_space.xyz, vec3(.0,1.,.0))*1.5;
			float height = Noise(p.xz*2.0)*.75 + Noise(p.xz)*.35 + Noise(p.xz*.5)*.2;
			
			//p.y += height;

			float y = p.y - base -height;
			y = y*y;

			vec2 ret = Voronoi((p.xz*2.5+sin(y*4.0+p.zx*12.3)*.12+vec2(sin(_Time[1]*2.3+1.5*p.z),sin(_Time[1]*3.6+1.5*p.x))*y*.5));
			float f = ret.x * .6 + y * .58;
			return vec3( y - f*1.4, clamp(f * 1.5, 0.0, 1.0), ret.y);
		} 
		
		                
		float CircleOfConfusion(float t)
		{
		return max(t * .04, (2.0 / _ScreenParams.y) * (1.0+t));
		}
		
		
		
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
				
				vec3 p = rO + rD * d;
				//vec3 p = position_in_world_space.xyz;
				//vec3 up = normalize(cross(cross(p, vec3(0,0,1)), p));
				
				vec3 ret = DE(p);
				//ret = normalize(cross(cross(ret, vec3(0,-1,0)), ret));
				ret.x += .5 * rCoC;

				if (ret.x < rCoC)
				{
					alpha = (1.0 - col.y) * Linstep(-rCoC, rCoC, -ret.x) * 2.0;//calculate the mix like cloud density
					f = clamp(ret.y, 0.0, 1.0);
					// Mix material with white tips for grass...
					vec3 gra = mix(mat, vec3(.2, .3, min(pow(ret.z, 2.0)*3.0, .2)), pow(ret.y,20.0)*.6) * ret.y;
					col += vec4(gra * alpha, alpha);
				}
				d += max(ret.x * .7, .02);
			}
			if(col.w < .2)col.xyzw = vec4(.1, .15, .05,.0);
			return col.xyzw;
		}
		


 
         void main()
         {
         	mat4 worldlMatrix = _World2Object;
         	
         	//mat4 modelMatrixInverse = _World2Object * unity_Scale.w;
 			
      		//modelMatrixInverse[3][3] = 1.0; 
      		
      		//mat4 viewMatrix = gl_ModelViewMatrix * modelMatrixInverse;
      		
      		
            //float dir= distance(position_in_world_space, vec4(_WorldSpaceCameraPos,1.0));
            vec4 camSpacePos = UNITY_MATRIX_V * position_in_world_space;
            vec4 camSpace = UNITY_MATRIX_V * vec4(_WorldSpaceCameraPos,1.0);
            
            vec3 testDir = normalize(UNITY_MATRIX_IT_MV[2].xyz);
            
            vec4 objSpace = position_in_world_space * worldlMatrix;
            vec4 objSpaceCam = vec4(_WorldSpaceCameraPos,1.0) * worldlMatrix;
            
            //vec3 dir = normalize( camSpacePos.xyz-camSpace.xyz);
            vec3 dir = normalize( position_in_world_space.xyz+_WorldSpaceCameraPos);
            //vec3 up = cross(cross(position_in_world_space.xyz, vec3(0,1,0)), position_in_world_space.xyz);
            vec3 mat = mix(vec3(.0,.3,.0), vec3(.2,.3,.0), Noise(position_in_world_space.xy*.025));
            float dist;
            
            
            vec4 col = GrassBlades( position_in_world_space.xyz, dir , mat, dist);
            
 
			gl_FragColor = col;
            

         }
 
         #endif
 
         ENDGLSL
      }
   }
}