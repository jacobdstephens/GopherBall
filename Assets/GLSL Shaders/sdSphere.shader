Shader "Custom/sdSphere" 
{
	Properties 
	{
	}
	SubShader 
	{
		Pass 
		{
		GLSLPROGRAM
		
		#include "UnityCG.glslinc" 
		
		//----------Script Set up-----------//
		
		//Uniforms for Screen Space and WorldSpace computations
        varying vec4 worldSpacePointPosition;
        varying vec4 cameraSpacePointposition;
        
		//Signed Distance function to create a Plane
		float signedDistancePlane( vec3 p )
		{
			return p.y;
		}
		//Signed Distance function to create a Sphere
		float signedDistanceSphere( vec3 p, float s )
		{
    		return length(p)-s;
		}
		//Unite Operation
		vec2 opU( vec2 d1, vec2 d2 )
		{
			return (d1.x<d2.x) ? d1 : d2;
		}
		//Create Map
		vec2 map ( in vec3 position )
		{
			vec2 res = opU ( vec2 ( signedDistancePlane ( position ), 1.0 ),
	                vec2 ( signedDistanceSphere ( position - vec3 ( 0.0 , 0.25 , 0.0 ), 0.25 ), 46.9 ) );
	        return res;
		}
		
		float doModel ( in vec3 position )
		{
			float res = signedDistanceSphere ( position - vec3( 0.0, 0.25 , 0.0 ), 0.25 );
	        return res;
		}
		
		float calcIntersection( in vec3 ro, in vec3 rd )
		{
			const float maxd = 20.0;           // max trace distance
			const float precis = 0.001;        // precission of the intersection
    		float h = precis*2.0;
    		float t = 0.0;
			float res = -1.0;
    		for( int i=0; i<90; i++ )          // max number of raymarching iterations is 90
    		{
        		if( h<precis||t>maxd ) break;
				h = doModel( ro+rd*t );
        		t += h;
    		}

    		if ( t < maxd ) res = t;
    		return res;
		}
		
		
		//------------Vertex Shader-------------//
		#ifdef VERTEX
		
		void main()
		{
			gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
		}
		
		#endif
		
		//-------------Fragment Shader----------------//
		#ifdef FRAGMENT
		
		void main()
		{
			//Render
			//vec3 col = render( ro, rd );
			//float target = calcIntersection( rayOrigin, rayDirection );
			
			
			gl_FragColor = vec4(1.0,1.0,1.0,1.0);
		}
		
		#endif
		
		ENDGLSL
		}
	} 	
}
