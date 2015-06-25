Shader "Custom/Grass" 
{
	Properties 
	{
		_Color ("Color", Color) = (1,1,1,1)
		_X ("X", Range(-10,10)) = 0.0
		_Y ("Y", Range(-10,10)) = 0.0
		_Z ("Z", Range(-10,10)) = 0.0
		
		_Xa ("Xa", Range(-10,10)) = 0.0
		_Ya ("Ya", Range(-10,10)) = 0.0
		_Za ("Za", Range(-10,10)) = 0.0
		_GrassBall ("Grass_Ball", Range(0,10)) = 7.25
	}
	SubShader 
	{
		Tags { "RenderType"="Transparent" }
		Pass 
		{
		GLSLPROGRAM
		
		
		
		//----------Script Set up-----------//
		
		#include "UnityCG.glslinc" 
		#define PI 3.14159265
		#define MOD2 vec2(3.07965, 7.4235)
   
        //UI Slider inputs from Shader
        
        uniform float _X;
        uniform float _Y;
        uniform float _Z;
        
        uniform float _Xa;
        uniform float _Ya;
        uniform float _Za;
        
        uniform float _GrassBall;

		//Position vectors for world splace computation
        varying vec4 worldSpacePointPosition;
        varying vec4 cameraSpacePointPosition;
        varying vec4 tempCol;

		//Constants and Defaults
		vec3 sunLight  = normalize(_WorldSpaceLightPos0.xyz);
		vec3 cameraPos;
		vec3 sunColour = vec3(1.0, .75, .6);
		const mat2 rotate2D = mat2(1.932, 1.623, -1.623, 1.952);
		float gTime = 0.0;

		// Noise functions
		//--------------------------------------------------------------------------
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
			return vec2(max(.4-sqrt(res), 0.01),id);
		}

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

		//Build Model in Screen Space
		float doModel ( in vec3 position )
		{
			//Todo hook up values to UI
			vec3 spherePositionOffset = vec3 (  0.0 , 0.0, 0.0 ) * 1.;
			//float scale = 7.25;
			float result = signedDistanceSphere ( position - spherePositionOffset.xyz, _GrassBall );
	        return result;
		}

		
		//Calculate Ray Intersection in Screen Space
		vec2 castRay( in vec3 ro, in vec3 rd )
		{
    		float tmin = 1.0;
    		float tmax = 80.0;//Todo add control for distance
    
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
	    		vec2 res = vec2(doModel( ro+rd*t ),1.0);
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
				float h = doModel( ro + rd*t );
        		res = min( res, 8.0*h/t );
        		t += clamp( h, 0.02, 0.10 );
        		if( h<0.001 || t>tmax ) break;
    		}
    		return clamp( res, 0.0, 1.0 );

		}
		
		//Todo rename SetCamera to create Look At Matrix or something better
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
        		float dd = doModel( aopos );
        		occ += -(dd-hr)*sca;
        		sca *= 0.95;
    		}
    		return clamp( 1.0 - 3.0*occ, 0.0, 1.0 );    
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
			//return vec4(0.0,0.0,0.0,1.0);
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
		vec3 DE(in vec3 p, in vec3 normal, in vec3 camSpaceP)
		{
			//Convert to Tangent Space
			mat3 XcameraMatrix = setCamera( vec3(0.0,1.0,.0) , _WorldSpaceCameraPos, 0.0);
			p = p * XcameraMatrix ;
			 
			//Set distance to bottom of grass
			float y = signedDistanceSphere ( p - vec3(_Xa,_Ya,_Za), 5.5 );
			//y = y*y;
			//Convert _X and _Y to elevation azmuth. Todo Fix popping when atan flips

			p.xy = p.xy + vec2(atan(_WorldSpaceCameraPos.z , _WorldSpaceCameraPos.x )*-5.0, normalize(_WorldSpaceCameraPos).y*-5.0).xy;
			//Pass Tangent Space Points into Voronoi function to create grass in perspective
			//vec2 ret = Voronoi(p.xy);
			vec2 ret = Voronoi((p.xy*2.5+sin(y*4.0+p.yx*12.3)*.12+vec2(sin(_Time[1]*1.3+1.5*p.y),sin(_Time[1]*2.6+1.5*p.x))*y*.5));

			float f = ret.x * .6 + y * .58;
			return vec3( y - f*1.4, clamp(f * 1.5, 0.0, 1.0), ret.y);
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
		vec4 GrassBlades(in vec3 rO, in vec3 rD, in vec3 mat, in float dist , in vec3 normal, in vec3 camSpaceP)
		{
			float d = 0.0;
			// Only calculate cCoC once is enough here...
			float rCoC = CircleOfConfusion(dist*.3);
			float alpha = 0.0;
			
			vec4 col = vec4(mat*0.15, 0.0);

			for (int i = 0; i < 99; i++)
			{
				if (col.w > .99) break;
				
				//p is the point of intersection on the grass plane
				vec3 p = rO + rD * d;
				
				vec3 nor = calcNormal( p );
				
				//Todo need to pass in tangent and binormal for the RayDirection to get the Voronoi Plane to orient correctly
				vec3 ret = DE(p, nor, camSpaceP);
				ret.x += .5 * rCoC;

				if (ret.x < rCoC)
				{
					alpha = (1.0 - col.y) * Linstep(-rCoC, rCoC, -ret.x);//calculate the mix like cloud density
					// Mix material with white tips for grass...
					vec3 gra = mix(mat, vec3(.35, .35, min(pow(ret.z, 4.0)*35.0, .35)), pow(ret.y, 9.0)*.7) * ret.y;
					col += vec4(gra * alpha, alpha);
				}
				d += max(ret.x *2., .01);
			}
			if(col.w < .2)
				col = vec4(GetSky(rD), 0.0);
				
			return col;
		}

		//--------------------------------------------------------------------------
		// Calculate sun light...
		void DoLighting(inout vec3 mat, in vec3 pos, in vec3 normal, in vec3 eyeDir, in float dis)
		{
			float h = dot(sunLight,normal);
			mat = mat * sunColour*(max(h, 0.0)+.2);
		}


		
		vec4 render( in vec3 ro, in vec3 rd, in vec3 camSpaceP )
		{ 
    		vec4 col;
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
        		col.xyz = mix(vec3(.0,.3,.0), vec3(.2,.3,.0), Noise(pos.xz*.025));
        		
        		col = GrassBlades(ro, rd, col.xyz, distance, nor, camSpaceP);
        		DoLighting(col.xyz, pos, nor,rd, distance);
        		col.xyz = ApplyFog(col.xyz, distance, rd);
        		    
    		}

			return col;
		}
		
		//--------------------------------------------------------------------------
		// Home in on the surface by dividing by two and split...
		float BinarySubdivision(in vec3 rO, in vec3 rD, float t, float oldT)
		{
			float halfwayT = 0.0;
			for (int n = 0; n < 5; n++)
			{
				halfwayT = (oldT + t ) * .5;
				if (doModel(rO + halfwayT*rD) < .05)
				{
					t = halfwayT;
				}else
				{
					oldT = halfwayT;
				}
			}
			return t;
		}
		
		
		bool Scene(in vec3 rO, in vec3 rD, out float resT, out float type )
		{
		    float t = 5.;
			float oldT = 0.0;
			float delta = 0.;
			vec2 h = vec2(1.0, 1.0);
			bool hit = false;
			for( int j=0; j < 70; j++ )
			{
			    vec3 p = rO + t*rD;
			    
				h = vec2(doModel(p), 1.0); // ...Get this position's height mapping.

				// Are we inside, and close enough to fudge a hit?...
				if( h.x < 0.05)
				{
					hit = true;
		            break;
				}
			        
				delta = h.x + (t*0.03);
				oldT = t;
				t += delta;
			}
		    type = h.y;
		    resT = BinarySubdivision(rO, rD, t, oldT);
			return hit;
		}
		
		vec3 PostEffects(vec3 rgb, vec2 xy)
		{
			// Gamma first...
			rgb = pow(rgb, vec3(0.45));
			
			// Then...
			#define CONTRAST 1.1
			#define SATURATION 1.3
			#define BRIGHTNESS 1.3
			rgb = mix(vec3(.5), mix(vec3(dot(vec3(.2125, .7154, .0721), rgb*BRIGHTNESS)), rgb*BRIGHTNESS, SATURATION), CONTRAST);
			// Vignette...
			rgb *= .4+0.5*pow(40.0*xy.x*xy.y*(1.0-xy.x)*(1.0-xy.y), 0.2 );	
			return rgb;
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
			vec2 xy = gl_FragCoord.xy / _ScreenParams.xy;
			

			
			float lensDistance = gl_FragCoord.y * .5/tan(radians(60.0) * .5 );//60 is the FOV in Unity
			
			//Render

			vec3 rayOrigin = worldSpacePointPosition.xyz ;
			vec3 camTarget = _WorldSpaceCameraPos;
			mat3 cameraMatrix = setCamera( camTarget.xyz,rayOrigin.xyz, .0 );
			vec3 rayDirection = cameraMatrix * normalize(vec3(camPlane.xy,lensDistance) );
						
			vec3 cw = normalize(_WorldSpaceCameraPos);
			//vec3 cp = vec3(sin(roll), cos(roll),0.0);
			vec3 cu = cross(camTarget,rayOrigin);
			vec3 cv = cross(cu,camTarget);
			vec3 dir = normalize(camPlane.x*cu + camPlane.y*cv + 1.3*cw);
			mat3 camMat = mat3(cu, cv, cw);
			
			
			vec4 col;

			float distance;
			float type;
			if( !Scene(camTarget, rayDirection, distance, type) )
			{
				// Missed scene, now just get the sky...
				col = vec4(GetSky(rayDirection), 1.0);
			}
				else
			{
				col = render( camTarget, rayDirection, vec3(camPlane.xy,0.0) );
			
			}

			col.xyz = PostEffects(col.xyz, xy);	
			
			gl_FragColor = col;
		
		}
		#endif
		ENDGLSL
		}
	} 	
}
