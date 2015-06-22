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
		#define MOD2 vec2(3.07965, 7.4235)
		
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
        		
		//float PI  = 4.0*atan(1.0);
		vec3 sunLight  = normalize( vec3(  0.35, 0.2,  0.3 ) );
		vec3 cameraPos;
		vec3 sunColour = vec3(1.0, .75, .6);
		const mat2 rotate2D = mat2(1.932, 1.623, -1.623, 1.952);
		float gTime = 0.0;
        
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
			vec2 res = vec2 ( signedDistanceSphere ( position - vec3(0.0,0.0,0.0), 6.0 ), 46.9 );
	        return res;
		}
		
		


		//--------------------------------------------------------------------------
		// Noise functions...
		float Hash( float p )
		{
			vec2 p2 = fract(vec2(p) / MOD2);
		    p2 += dot(p2.yx, p2.xy+19.19);
			return fract(p2.x * p2.y);
		}

		//--------------------------------------------------------------------------
		float Hash(vec2 p)
		{
			p  = fract(p / MOD2);
		    p += dot(p.xy, p.yx+19.19);
		    return fract(p.x * p.y);
		}


		//--------------------------------------------------------------------------
		float Noise( in vec2 x )
		{
		    vec2 p = floor(x);
		    vec2 f = fract(x);
		    f = f*f*(3.0-2.0*f);
		    float n = p.x + p.y*57.0;
		    float res = mix(mix( Hash(n+  0.0), Hash(n+  1.0),f.x),
		                    mix( Hash(n+ 57.0), Hash(n+ 58.0),f.x),f.y);
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
				vec2 r = vec2( b ) - f  + Hash( p + b );
				float d = dot(r,r);
				if( d < res )
				{
					res = d;
					id  = Hash(p+b);
				}			
		    }
			return vec2(max(.4-sqrt(res), 0.0),id);
		}


		//--------------------------------------------------------------------------
		vec2 Terrain( in vec2 p)
		{
			float type = 0.0;
			vec2 pos = p*0.003;
			float w = 50.0;
			float f = .0;
			for (int i = 0; i < 3; i++)
			{
				f += Noise(pos) * w;
				w = w * 0.62;
				pos *= 2.5;
			}
	        

			return vec2(f, type);
		}

		//--------------------------------------------------------------------------
		vec2 Map(in vec3 p)
		{
			vec2 h = Terrain(p.xz);
		    return vec2(p.y - h.x, h.y);
		}

		//--------------------------------------------------------------------------
		float FractalNoise(in vec2 xy)
		{
			float w = .7;
			float f = 0.0;

			for (int i = 0; i < 3; i++)
			{
				f += Noise(xy) * w;
				w = w*0.6;
				xy = 2.0 * xy;
			}
			return f;
		}

		//--------------------------------------------------------------------------
		// Grab all sky information for a given ray from camera
		vec3 GetSky(in vec3 rd)
		{
			float sunAmount = max( dot( rd, sunLight), 0.0 );
			float v = pow(1.0-max(rd.y,0.0),6.);
			vec3  sky = mix(vec3(.1, .2, .3), vec3(.32, .32, .32), v);
			sky = sky + sunColour * sunAmount * sunAmount * .25;
			sky = sky + sunColour * min(pow(sunAmount, 800.0)*1.5, .3);
			return clamp(sky, 0.0, 1.0);
		}

		//--------------------------------------------------------------------------
		// Merge grass into the sky background for correct fog colouring...
		vec3 ApplyFog( in vec3  rgb, in float dis, in vec3 dir)
		{
			float fogAmount = clamp(dis*dis* 0.0000012, 0.0, 1.0);
			return mix( rgb, GetSky(dir), fogAmount );
		}

		//--------------------------------------------------------------------------
		vec3 DE(in vec3 p, in vec3 normal, in vec2 camPlane)
		{
			vec4 camNorm = UNITY_MATRIX_IT_MV * vec4(0.0,0.0,0.0,1.0);
			p = normal - p;
			float base = signedDistanceSphere ( p - vec3(0.0,0.0,0.0), 3.0 )*1.5;
			
			float height = Noise(p.xz*2.0)*.75 + Noise(p.xz)*.35 + Noise(p.xz*.5)*.2;
			//p.y += height;
			float y = base;
			y = y*y;
			//TODO create vec2 that represents the camera plane
			
			
			vec2 ret = Voronoi((p.xy*2.5+sin(y*4.0+p.yx*12.3)*.12+vec2(sin(_Time[0]*2.3+1.5*p.y),sin(_Time[0]*3.6+1.5*p.x))*y*.5));
			float f = ret.x * .6 + y * .58;
			return vec3( y - f*1.4, ret.y,clamp(f * 1.5, 0.0, 1.0));
		}

		//--------------------------------------------------------------------------
		// eiffie's code for calculating the aperture size for a given distance...
		float CircleOfConfusion(float t)
		{
			return max(t * .04, (2.0 / _ScreenParams.y) * (1.0+t));
		}

		//--------------------------------------------------------------------------
		float Linstep(float a, float b, float t)
		{
			return clamp((t-a)/(b-a),0.,1.);
		}

		//--------------------------------------------------------------------------
		vec3 GrassBlades(in vec3 rO, in vec3 rD, in vec3 mat, in float dist , in vec3 normal, in vec2 camPlane)
		{
			float d = 0.0;
			// Only calculate cCoC once is enough here...
			float rCoC = CircleOfConfusion(dist*.3);
			float alpha = 0.0;
			
			vec4 col = vec4(mat*0.15, 0.0);

			for (int i = 0; i < 20; i++)
			{
				if (col.w > .99) break;
				vec3 p = rO + rD * d;
				
				vec3 ret = DE(p, normal, camPlane);
				ret.x += .5 * rCoC;

				if (ret.x < rCoC)
				{
					alpha = (1.0 - col.y) * Linstep(-rCoC, rCoC, -ret.x);//calculate the mix like cloud density
					// Mix material with white tips for grass...
					vec3 gra = mix(mat, vec3(.35, .35, min(pow(ret.z, 4.0)*35.0, .35)), pow(ret.y, 9.0)*.7) * ret.y;
					col += vec4(gra * alpha, alpha);
				}
				d += max(ret.x * .7, .1);
			}
			if(col.w < .2)
				col.xyz = vec3(0.1, .15, 0.05);
			return col.xyz;
		}

		//--------------------------------------------------------------------------
		// Calculate sun light...
		void DoLighting(inout vec3 mat, in vec3 pos, in vec3 normal, in vec3 eyeDir, in float dis)
		{
			float h = dot(sunLight,normal);
			mat = mat * sunColour*(max(h, 0.0)+.2);
		}

		//--------------------------------------------------------------------------
		vec3 TerrainColour(vec3 pos, vec3 dir,  vec3 normal, float dis, in vec2 camPlane,float type)
		{
			vec3 mat;
			if (type == 0.0)
			{
				// Random colour...
				mat = mix(vec3(.0,.3,.0), vec3(.2,.3,.0), Noise(pos.xz*.025));
				// Random shadows...
				float t = FractalNoise(pos.xz * .1)+.5;
				// Do grass blade tracing...
				mat = GrassBlades(pos, dir, mat, dis, normal, camPlane) * t;
				DoLighting(mat, pos, normal,dir, dis);
			}
			mat = ApplyFog(mat, dis, dir);
			return mat;
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

		
		vec3 render( in vec3 ro, in vec3 rd, in vec2 camPlane )
		{ 
    		vec3 col = vec3(0.8, 0.9, 1.0);
    		vec3 grass;
    		vec2 res = castRay(ro,rd);
    		float t = res.x;
			float m = res.y;
			float distance;
			float type;
    		if( m>-0.5 )
    		{
        		vec3 pos = ro + t*rd;
        		vec3 nor = calcNormal( pos );
        		
        		vec3 ref = reflect( rd, nor );
        
        		grass = TerrainColour(ro, rd, nor, distance, camPlane, type);
        		
        		// material        
				//col = 0.45 + 0.3*sin( vec3(0.05,0.08,0.10)*(m-1.0) );
		
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

			return vec3( clamp(grass,0.0,1.0) );
		}
		
		
		
		
		//------------Vertex Shader-------------//
		#ifdef VERTEX
		
		void main()
		{
            mat4 modelMatrix = _Object2World;
            
 			worldSpacePointPosition = modelMatrix * gl_Vertex;
			gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
		}
		
		#endif
		
		//-------------Fragment Shader----------------//
		#ifdef FRAGMENT
		
		void main()
		{
			vec2 camPlane = -1.0 + 2.0 * gl_FragCoord.xy / _ScreenParams.xy;
			camPlane.x *= _ScreenParams.x/_ScreenParams.y;
			
			float lensDistance = gl_FragCoord.y * .5/tan(radians(60.0) * .5 );//60 is the FOV in Unity
			
			//Render
			vec3 color = vec3(.3,.3,.3);
			
			vec4 rayOrigin = worldSpacePointPosition ;
			vec3 camTarget = _WorldSpaceCameraPos;
			mat3 ca = setCamera( rayOrigin.xyz, camTarget.xyz, .0 );
			vec3 rayDirection = ca * normalize(vec3(camPlane.xy,-lensDistance) );//Flip Z Axis to project onto sphere
			//vec3 surfaceNormal = calcNormal();
			
			vec3 col = render( rayOrigin.xyz, rayDirection.xyz, camPlane );

			col = pow( col, vec3(0.4545) );

			gl_FragColor = vec4(col,0.0);
		}
		
		#endif
		
		ENDGLSL
		}
	} 	
}
