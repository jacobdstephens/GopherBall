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
		#define PI 3.14159265
		
		//----------Script Set up-----------//
		
		//Uniforms for Screen Space and WorldSpace computations
		uniform mat4 UNITY_MATRIX_IT_MV;
		uniform mat4 UNITY_MATRIX_V;
		
        varying vec4 worldSpacePointPosition;
        varying vec4 cameraSpacePointPosition;
        varying vec4 modelSpacePointPosition;
        
		float orenNayarDiffuse(
		  vec3 lightDirection,
		  vec3 viewDirection,
		  vec3 surfaceNormal,
		  float roughness,
		  float albedo) {
  
		  float LdotV = dot(lightDirection, viewDirection);
		  float NdotL = dot(lightDirection, surfaceNormal);
		  float NdotV = dot(surfaceNormal, viewDirection);

		  float s = LdotV - NdotL * NdotV;
		  float t = mix(1.0, max(NdotL, NdotV), step(0.0, s));

		  float sigma2 = roughness * roughness;
		  float A = 1.0 + sigma2 * (albedo / (sigma2 + 0.13) + 0.5 / (sigma2 + 0.33));
		  float B = 0.45 * sigma2 / (sigma2 + 0.09);

		  return albedo * max(0.0, NdotL) * (A + B * s / t) / PI;
		}

		float gaussianSpecular(
		  vec3 lightDirection,
		  vec3 viewDirection,
		  vec3 surfaceNormal,
		  	float shininess) 
		  	{
		  		vec3 H = normalize(lightDirection + viewDirection);
		  		float theta = acos(dot(H, surfaceNormal));
		  		float w = theta / shininess;
		  		return exp(-w*w);
			}

        float fogFactorExp2(
		  const float dist,
		  const float density
		) {
		  const float LOG2 = -1.442695;
		  float d = density * dist;
		  return 1.0 - clamp(exp2(d * d * LOG2), 0.0, 1.0);
		}

        
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
		//Build Model in Screen Space
		float doModel ( in vec3 position )
		{
			vec3 spherePositionOffset = vec3 ( 0.0, 0.0 , 0.0 ) ;
			float scale = 3.0;
			float result = signedDistanceSphere ( position - spherePositionOffset, scale );
	        return result;
		}
		//Create Material
		vec3 doMaterial( in vec3 pos, in vec3 nor )
		{
    		return vec3 (0.125,0.1,0.2) + ( vec3 ( .6 , 0.9 , .4 ) * 3. * clamp ( length( pos ) - 0.94 , 0. , 1. ) );
		}
		
		float calcSoftshadow( in vec3 ro, in vec3 rd );

		vec3 doLighting( in vec3 pos, in vec3 nor, in vec3 rd, in float dis, in vec3 mal )
		{
    		vec3 lin = vec3(0.0);

    		// key light
    		//-----------------------------
    		vec3  view = normalize(-rd);
    		vec3  lig1 = normalize(vec3(1.0,0.7,0.9));
    		vec3  lig2 = normalize(vec3(1.0,0.9,0.9)*-1.);
    
    		float spc1 = gaussianSpecular(lig1, view, nor, 0.95)*0.5;
    		float dif1 = max(0., orenNayarDiffuse(lig1, view, nor, -20.1, 1.0));
    		float sha1 = 0.0; if( dif1>0.01 ) sha1=calcSoftshadow( pos+0.01*nor, lig1 );
    		vec3  col1 = vec3(2.,4.2,4.);
    		lin += col1*spc1+dif1*col1*sha1;
    
    		float spc2 = gaussianSpecular(lig2, view, nor, 0.95);
    		float dif2 = max(0., orenNayarDiffuse(lig2, view, nor, -20.1, 1.0));
    		float sha2 = 0.0; if( dif2>0.01 ) sha2=calcSoftshadow( pos+0.01*nor, lig2 );
    		vec3  col2 = vec3(2.00,0.05,0.15);
    		lin += col2*spc2+dif2*col2*sha1;

    		// ambient light
    		//-----------------------------
    		lin += vec3(0.05);

    
    		// surface-light interacion
    		//-----------------------------
    		vec3 col = mal*lin;

    		return col;
		}
		
		//Calculate Ray Intersection in Screen Space
		float calcIntersection( in vec3 rayOrigin, in vec3 rayDirection )
		{
			const float maximumTraceDiscance = 20.0;       
			const float intersectionPrecision = 0.001;        
    		float heightInScreenSpace = intersectionPrecision * 2.0;
    		float traceDistance = 0.0;
			float cameraTargetDistance = -1.0;
    		for( int i = 0 ; i < 90; i++ )          // max number of raymarching iterations is 90
    		{
        		if( heightInScreenSpace < intersectionPrecision || traceDistance > maximumTraceDiscance ) break;
				heightInScreenSpace = doModel ( rayOrigin + rayDirection * traceDistance );
        		traceDistance += heightInScreenSpace;
    		}

    		if ( traceDistance < maximumTraceDiscance ) cameraTargetDistance = traceDistance;
    		return cameraTargetDistance;
		}
		
		vec2 castRay( in vec3 ro, in vec3 rd )
		{
    		float tmin = 1.0;
    		float tmax = 20.0;
    
			#if 0
   			float tp1 = (0.0-ro.y)/rd.y; if( tp1>0.0 ) tmax = min( tmax, tp1 );
    		float tp2 = (1.6-ro.y)/rd.y; if( tp2>0.0 ) { if( ro.y>1.6 ) tmin = max( tmin, tp2 );
                                                 else           tmax = min( tmax, tp2 ); }
			#endif
    
			float precis = 0.002;
    		float t = tmin;
    		float m = -1.0;
    		for( int i=0; i<50; i++ )
    		{
	    		vec2 res = map( ro+rd*t );
        		if( res.x<precis || t>tmax ) break;
        		t += res.x;
	    		m = res.y;
    		}

    		if( t>tmax ) m=-1.0;
    		return vec2( t, m );
		}
		
		vec3 calcNormal( in vec3 pos )
		{
    		const float eps = 0.002;             // precision of the normal computation

    		const vec3 v1 = vec3( 1.0,-1.0,-1.0);
    		const vec3 v2 = vec3(-1.0,-1.0, 1.0);
    		const vec3 v3 = vec3(-1.0, 1.0,-1.0);
    		const vec3 v4 = vec3( 1.0, 1.0, 1.0);

			return normalize( v1*doModel( pos + v1*eps ) + 
					  		v2*doModel( pos + v2*eps ) + 
					  		v3*doModel( pos + v3*eps ) + 
					  		v4*doModel( pos + v4*eps ) );
		}
		
		float calcSoftshadow( in vec3 ro, in vec3 rd )
		{
		    float res = 1.0;
		    float t = 0.0001;                 // selfintersection avoidance distance
			float h = 1.0;
		    for( int i=0; i<5; i++ )         // 40 is the max numnber of raymarching steps
		    {
		        h = doModel(ro + rd*t);
		        res = min( res, 4.0*h/t );   // 64 is the hardness of the shadows
				t += clamp( h, 0.02, 2.0 );   // limit the max and min stepping distances
		    }
		    return clamp(res,0.0,1.0);
		}
		
		mat3 calcLookAtMatrix( in vec3 ro, in vec3 ta, in float roll )
		{
    		vec3 ww = normalize( ta - ro );
    		vec3 uu = normalize( cross(ww,vec3(sin(roll),cos(roll),0.0) ) );
   	 		vec3 vv = normalize( cross(uu,ww));
    		return mat3( uu, vv, ww );
		}
		
		
		
		
		//------------Vertex Shader-------------//
		#ifdef VERTEX
		
		void main()
		{
            mat4 modelMatrix = _Object2World;
            
 			worldSpacePointPosition = modelMatrix * gl_Vertex;
 			cameraSpacePointPosition = gl_ModelViewProjectionMatrix * gl_Vertex;
 			modelSpacePointPosition = gl_Vertex;
 			
			gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
		}
		
		#endif
		
		//-------------Fragment Shader----------------//
		#ifdef FRAGMENT
		
		void main()
		{
			vec2 p = -1.0 + 2.0 * gl_FragCoord.xy / _ScreenParams.xy;
			p.x *= _ScreenParams.x/_ScreenParams.y;
			//Render
			vec3 color = vec3(.3,.3,.3);
			
			vec4 rayOrigin = worldSpacePointPosition;//_WorldSpaceCameraPos;//UNITY_MATRIX_IT_MV[0].xyz;
			vec4 rO = UNITY_MATRIX_IT_MV * worldSpacePointPosition;//_WorldSpaceCameraPos;
			vec3 camTarget = vec3 (.0,.0,.0);
			mat3 camMat = calcLookAtMatrix( rO.xyz, camTarget, 0.0 );  // 0.0 is the camera roll
			vec3 rd = normalize( camMat * vec3(p.xy,2.0) );
			
			
			vec3 rayDirection = rd;//UNITY_MATRIX_IT_MV[2].xyz;
			float target = calcIntersection( rayOrigin.xyz, rayDirection.xyz );
			if ( target > -0.5 )
			{
				vec3 targetPosition = rO.xyz + target * rayDirection.xyz;
				vec3 targetNormal = calcNormal(targetPosition);
				vec3 materialContribution = doMaterial( targetPosition, targetNormal );
				vec3 lightContribution = doLighting( targetPosition, targetNormal, rayDirection.xyz, target, materialContribution );
				color = mix ( lightContribution, color, fogFactorExp2 (target, .01));
			}
			
			
			gl_FragColor = vec4(color,1.0);
		}
		
		#endif
		
		ENDGLSL
		}
	} 	
}
