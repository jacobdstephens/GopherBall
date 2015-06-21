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
		uniform mat4 UNITY_MATRIX_MV;
		uniform mat4 UNITY_MATRIX_V;
		uniform mat4 UNITY_MATRIX_VP;
        uniform mat4 UNITY_MATRIX_P;
        
        uniform mat4 _Projector;
		
        varying vec4 worldSpacePointPosition;
        varying vec4 worldSpaceOrigin;
        varying vec4 cameraSpacePointPosition;
        varying vec4 modelSpacePointPosition;
        
        varying vec4 positionInProjSpace; 
        
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
			vec2 res = vec2 ( signedDistanceSphere ( position - vec3(0.0,0.0,0.0), 5.0 ), 46.9 );
	        return res;
		}
		//Build Model in Screen Space
		float doModel ( in vec3 position )
		{
			vec3 spherePositionOffset = vec3 (  0.0 ,  0.0 ,  0.0 ) * 1.;
			//vec4 relPos = ;
			//vec4 ( _WorldSpaceCameraPos, 1.0))
			vec4 relPos = (vec4 ( 0.0,0.0,0.0, 1.0)-vec4 ( _WorldSpaceCameraPos, 1.0))*1.;
			float scale = 7.25;
			float result = signedDistanceSphere ( position - spherePositionOffset.xyz, scale );
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
    		vec3  lig1 = (vec3(100.0,70.7,90.9));
    		vec3  lig2 = (vec3(100.0,90.9,90.9)*-1.);
    
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
		
		float softshadow( in vec3 ro, in vec3 rd, in float mint, in float tmax )
		{
			float res = 1.0;
    		float t = mint;
    		for( int i=0; i<16; i++ )
    		{
				float h = map( ro + rd*t ).x;
        		res = min( res, 8.0*h/t );
        		t += clamp( h, 0.02, 0.10 );
        		if( h<0.001 || t>tmax ) break;
    		}
    		return clamp( res, 0.0, 1.0 );

		}
		
		mat3 calcLookAtMatrix( in vec3 ro, in vec3 ta, in float roll )
		{
    		vec3 ww = normalize( ta - ro );
    		vec3 uu = normalize( cross(ww,vec3(sin(roll),cos(roll),0.0) ) );
   	 		vec3 vv = normalize( cross(uu,ww));
    		return mat3( uu, vv, ww );
		}
		
		mat3 setCamera( in vec3 ro, in vec3 ta, float cr )
		{
			vec3 cw = normalize(ta-ro);
			vec3 cp = vec3(sin(cr), cos(cr),0.0);
			vec3 cu = normalize( cross(cw,cp) );
			vec3 cv = normalize( cross(cu,cw) );
		    return mat3( cu, cv, cw );
		}
		
		float calcAO( in vec3 pos, in vec3 nor )
		{
			float occ = 0.0;
    		float sca = 1.0;
    		for( int i=0; i<5; i++ )
    		{
        		float hr = 0.01 + 0.12*float(i)/4.0;
        		vec3 aopos =  nor * hr + pos;
        		float dd = map( aopos ).x;
        		occ += -(dd-hr)*sca;
        		sca *= 0.95;
    		}
    		return clamp( 1.0 - 3.0*occ, 0.0, 1.0 );    
		}

		
		vec3 render( in vec3 ro, in vec3 rd )
		{ 
    		vec3 col = vec3(0.8, 0.9, 1.0);
    		vec2 res = castRay(ro,rd);
    		float t = res.x;
			float m = res.y;
    		if( m>-0.5 )
    		{
        		vec3 pos = ro + t*rd;
        		vec3 nor = calcNormal( pos );
        		vec3 ref = reflect( rd, nor );
        
        		// material        
				col = 0.45 + 0.3*sin( vec3(0.05,0.08,0.10)*(m-1.0) );
		
        		if( m<1.5 )
        		{
            
            		float f = mod( floor(5.0*pos.z) + floor(5.0*pos.x), 2.0);
            		col = 0.4 + 0.1*f*vec3(1.0);
        		}

        		// lighitng        
        		float occ = calcAO( pos, nor );
				vec3  lig = normalize( vec3(-0.6, 0.7, -0.5) );
				float amb = clamp( 0.5+0.5*nor.y, 0.0, 1.0 );
        		float dif = clamp( dot( nor, lig ), 0.0, 1.0 );
        		float bac = clamp( dot( nor, normalize(vec3(-lig.x,0.0,-lig.z))), 0.0, 1.0 )*clamp( 1.0-pos.y,0.0,1.0);
        		float dom = smoothstep( -0.1, 0.1, ref.y );
        		float fre = pow( clamp(1.0+dot(nor,rd),0.0,1.0), 2.0 );
				float spe = pow(clamp( dot( ref, lig ), 0.0, 1.0 ),16.0);
        
        		dif *= softshadow( pos, lig, 0.02, 2.5 );
        		dom *= softshadow( pos, ref, 0.02, 2.5 );

				vec3 brdf = vec3(0.0);
        		brdf += 1.20*dif*vec3(1.00,0.90,0.60);
				brdf += 1.20*spe*vec3(1.00,0.90,0.60)*dif;
        		brdf += 0.30*amb*vec3(0.50,0.70,1.00)*occ;
        		brdf += 0.40*dom*vec3(0.50,0.70,1.00)*occ;
        		brdf += 0.30*bac*vec3(0.25,0.25,0.25)*occ;
        		brdf += 0.40*fre*vec3(1.00,1.00,1.00)*occ;
				brdf += 0.02;
				col = col*brdf;

    			col = mix( col, vec3(0.8,0.9,1.0), 1.0-exp( -0.0005*t*t ) );

    		}

			return vec3( clamp(col,0.0,1.0) );
		}
		
		
		
		
		//------------Vertex Shader-------------//
		#ifdef VERTEX
		
		void main()
		{
            mat4 modelMatrix = _Object2World;
            positionInProjSpace = _Projector * gl_Vertex;
            
 			worldSpacePointPosition = modelMatrix * gl_Vertex;
 			worldSpaceOrigin = modelMatrix * vec4(.0,.0,.0,1.0);
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
			
			float lensDistance = gl_FragCoord.y * .5/tan(radians(60.0) * .5 );
			
			//Render
			vec3 color = vec3(.3,.3,.3);
			
			vec4 rayOrigin = worldSpacePointPosition - vec4(0.0,0.0,0.0,1.0) ;//normalize( worldSpacePointPosition );//_WorldSpaceCameraPos;//UNITY_MATRIX_IT_MV[0].xyz;
			vec3 rayNorm = calcNormal(rayOrigin.xyz);
			vec4 rO = rayOrigin;//normalize(vec4 ( _WorldSpaceCameraPos, 1.0));//_WorldSpaceCameraPos;
			vec3 camTarget = _WorldSpaceCameraPos; //UNITY_MATRIX_V[2];// + vec4 (_WorldSpaceCameraPos, 1.0);//vec4 ( 0.0,.0,0.0, 1.0);
			mat3 ca = setCamera( rO.xyz, camTarget.xyz, .0 );
			//mat3 camMat = calcLookAtMatrix( rO.xyz, camTarget.xyz, 0.0 );  // 0.0 is the camera roll
			
			vec3 rd = ca * normalize(vec3(p.xy,-lensDistance) );//Flip Z Axis to project onto sphere
			vec3 rayDirection = rd.xyz;//normalize((vec3(p.xy,2.0)));
			
			vec3 col = render( rO.xyz, rayDirection.xyz );

			col = pow( col, vec3(0.4545) );
			//vec3 rd = ca * normalize(worldSpacePointPosition.xyz );
			
			//vec4 rd = UNITY_MATRIX_P * vec4(.0,.0,1.0, 1.0);
			
			
			float target = calcIntersection( rO.xyz, rayDirection.xyz );
			if ( target > -0.5 )
			{
				vec3 targetPosition = rO.xyz + target * rayDirection.xyz;
				vec3 targetNormal = calcNormal(targetPosition);
				vec3 materialContribution = doMaterial( targetPosition, targetNormal );
				vec3 lightContribution = doLighting( targetPosition, targetNormal, rayDirection.xyz, target, materialContribution );
				color = mix ( lightContribution, color, fogFactorExp2 (target, .01));
			}
			
			
			gl_FragColor = vec4(col,0.0);
		}
		
		#endif
		
		ENDGLSL
		}
	} 	
}
