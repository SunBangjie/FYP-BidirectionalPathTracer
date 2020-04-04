// This is a modification of the default Falcor routine in Shading.slang
ShadingData simplePrepareShadingData(VertexOut v, MaterialData m, float3 camPosW)
{
	ShadingData sd = initShadingData();

    ExplicitLodTextureSampler lodSampler = { 0 };  // Specify the tex lod/mip to use here

	// Sample the diffuse texture and apply the alpha test
	float4 baseColor = sampleTexture(m.resources.baseColor, m.resources.samplerState, v.texC, m.baseColor, EXTRACT_DIFFUSE_TYPE(m.flags), lodSampler);
	sd.opacity = m.baseColor.a;

	sd.posW = v.posW;
	sd.uv = v.texC;
	sd.V = normalize(camPosW - v.posW);
	sd.N = normalize(v.normalW);
	sd.B = normalize(v.bitangentW - sd.N * (dot(v.bitangentW, sd.N)));
	sd.T = normalize(cross(sd.B, sd.N));

	// Sample the spec texture
	float4 spec = sampleTexture(m.resources.specular, m.resources.samplerState, v.texC, m.specular, EXTRACT_SPECULAR_TYPE(m.flags), lodSampler);
	if (EXTRACT_SHADING_MODEL(m.flags) == ShadingModelMetalRough)
	{
		sd.diffuse = lerp(baseColor.rgb, float3(0), spec.b);
		sd.specular = lerp(float3(0.04f), baseColor.rgb, spec.b);
		sd.linearRoughness = spec.g;
	}
	else // if (EXTRACT_SHADING_MODEL(m.flags) == ShadingModelSpecGloss)
	{
		sd.diffuse = baseColor.rgb;
		sd.specular = spec.rgb;
		sd.linearRoughness = 1 - spec.a;
	}

	sd.linearRoughness = max(0.08, sd.linearRoughness); // Clamp the roughness so that the BRDF won't explode
	sd.roughness = sd.linearRoughness * sd.linearRoughness;
    sd.emissive = sampleTexture(m.resources.emissive, m.resources.samplerState, v.texC, float4(m.emissive, 1), EXTRACT_EMISSIVE_TYPE(m.flags), lodSampler).rgb;
    sd.IoR = m.IoR;
	sd.doubleSidedMaterial = EXTRACT_DOUBLE_SIDED(m.flags);

	// We're not applying normal mapping on 2ndary surfaces; saves quite a bit of cost
	//applyNormalMap(m, sd, lodSampler);
	sd.NdotV = dot(sd.N, sd.V);

	// Flip the normal if it's backfacing
	if (sd.NdotV <= 0 && sd.doubleSidedMaterial)
	{
		sd.N = -sd.N;
		sd.NdotV = -sd.NdotV;
	}

	return sd;
}

// Encapsulates a bunch of Falcor stuff into one simpler function. 
//    -> This can only be called within a closest hit or any hit shader
ShadingData getHitShadingData(BuiltInTriangleIntersectionAttributes attribs, float3 cameraPos )
{
	// Run a pair of Falcor helper functions to compute important data at the current hit point
	VertexOut  vsOut = getVertexAttributes(PrimitiveIndex(), attribs);
	return simplePrepareShadingData(vsOut, gMaterial, cameraPos);
}

// A work-around function because some DXR drivers seem to have broken atan2() implementations
float atan2_WAR(float y, float x)
{
	if (x > 0.f)
		return atan(y / x);
	else if (x < 0.f && y >= 0.f)
		return atan(y / x) + M_PI;
	else if (x < 0.f && y < 0.f)
		return atan(y / x) - M_PI;
	else if (x == 0.f && y > 0.f)
		return M_PI / 2.f;
	else if (x == 0.f && y < 0.f)
		return -M_PI / 2.f;
	return 0.f; // x==0 && y==0 (undefined)
}

// Convert our world space direction to a (u,v) coord in a latitude-longitude spherical map
float2 wsVectorToLatLong(float3 dir)
{
	float3 p = normalize(dir);

	// atan2_WAR is a work-around due to an apparent compiler bug in atan2
	float u = (1.f + atan2_WAR(p.x, -p.z) * M_1_PI) * 0.5f;
	float v = acos(p.y) * M_1_PI;
	return float2(u, v);
}

// Generates a seed for a random number generator from 2 inputs plus a backoff
uint initRand(uint val0, uint val1, uint backoff = 16)
{
	uint v0 = val0, v1 = val1, s0 = 0;

	[unroll]
	for (uint n = 0; n < backoff; n++)
	{
		s0 += 0x9e3779b9;
		v0 += ((v1 << 4) + 0xa341316c) ^ (v1 + s0) ^ ((v1 >> 5) + 0xc8013ea4);
		v1 += ((v0 << 4) + 0xad90777d) ^ (v0 + s0) ^ ((v0 >> 5) + 0x7e95761e);
	}
	return v0;
}

// Takes our seed, updates it, and returns a pseudorandom float in [0..1]
float nextRand(inout uint s)
{
	s = (1664525u * s + 1013904223u);
	return float(s & 0x00FFFFFF) / float(0x01000000);
}

// This function tests if the alpha test fails, given the attributes of the current hit. 
//   -> Can legally be called in a DXR any-hit shader or a DXR closest-hit shader, and 
//      accesses Falcor helpers and data structures to extract and perform the alpha test.
bool alphaTestFails(BuiltInTriangleIntersectionAttributes attribs)
{
	// Run a Falcor helper to extract the current hit point's geometric data
	VertexOut  vsOut = getVertexAttributes(PrimitiveIndex(), attribs);

    // Extracts the diffuse color from the material (the alpha component is opacity)
    ExplicitLodTextureSampler lodSampler = { 0 };  // Specify the tex lod/mip to use here
    float4 baseColor = sampleTexture(gMaterial.resources.baseColor, gMaterial.resources.samplerState,
        vsOut.texC, gMaterial.baseColor, EXTRACT_DIFFUSE_TYPE(gMaterial.flags), lodSampler);

	// Test if this hit point fails a standard alpha test.  
	return (baseColor.a < gMaterial.alphaThreshold);
}

uint2 getLaunchIndexFromDirection(float3 dir, uint2 dim, float2 jitter)
{
    float d1 = dot(dir, gCamera.cameraU) / dot(gCamera.cameraU, gCamera.cameraU);
    float d2 = dot(dir, gCamera.cameraV) / dot(gCamera.cameraV, gCamera.cameraV);
    float d3 = dot(dir, gCamera.cameraW) / dot(gCamera.cameraW, gCamera.cameraW);
    float2 ndc = float2(d1 / d3, -d2 / d3);
    float2 pixelCenter = ndc * (0.5, 0.5) + float2(0.5, 0.5);
    uint2 id = uint2(round(pixelCenter * dim - jitter));
    return id;
}

void sampleLight(inout uint rndSeed, out float3 origin, out float3 dir, out float3 lightIntensity)
{
    int index = min(int(nextRand(rndSeed) * gLightsCount), gLightsCount - 1);
    
    origin = gLights[index].posW;
    lightIntensity = gLights[index].intensity;

    if (gLights[index].type == LightDirectional)
        dir = gLights[index].dirW;
    else
        dir = sampleUnitSphere(rndSeed);
    dir = getCosHemisphereSample(rndSeed, dir);
}

float evalG(PathVertex a, PathVertex b)
{
    float3 posA = a.posW;
    float3 posB = b.posW;
    float3 vecAB = posB - posA;
    float lengthAB = length(vecAB);
    float invLengthAB = 1.0f / lengthAB;
    float3 dirAB = vecAB * invLengthAB;

    float cosA = abs(dot(a.N, dirAB));
    float cosB = abs(dot(b.N, dirAB));
	
    float visibilityTestDistance = distance(posA, posB);
    bool V = shadowRayVisibility(posA, dirAB, gMinT, lengthAB);
    
    return V ? cosA * cosB * invLengthAB * invLengthAB : 0.f;
}

float evalGWithoutV(PathVertex a, PathVertex b)
{
    float3 posA = a.posW;
    float3 posB = b.posW;
    float3 vecAB = posB - posA;
    float invLengthAB = 1.0f / length(vecAB);
    float3 dirAB = vecAB * invLengthAB;

    float cosA = abs(dot(a.N, dirAB));
    float cosB = abs(dot(b.N, dirAB));
    
    return cosA * cosB * invLengthAB * invLengthAB;
}

float3 getUnweightedContribution(PathVertex cameraPath[9], PathVertex lightPath[9], uint cameraIndex, uint lightIndex, float G)
{
	// this function only computes unweighted contribution for path length > 1
    if (cameraIndex <= 1 || lightIndex <= 1)
        return float3(0, 0, 0);
	
	// get end vertex from each path
    PathVertex cameraEndV = cameraPath[cameraIndex - 1];
    PathVertex lightEndV = lightPath[lightIndex - 1];
	
	// get throughput from the last vertex of each path
    float3 aE = cameraPath[cameraIndex - 2].color;
    float3 aL = lightPath[cameraIndex - 2].color;
	
	// connect two end vertices
    float3 connectDir = normalize(cameraEndV.posW - lightEndV.posW); // from light to camera

    float3 wi, wo;
    
	// compute fsL
    wi = connectDir; // towards camera
    wo = normalize(lightPath[lightIndex - 2].posW - lightEndV.posW); // towards light
    float3 fsL = evalBRDF(wi, wo, lightEndV.posW, lightEndV.N, lightEndV.N, lightEndV.dif, lightEndV.spec, lightEndV.rough, lightEndV.isSpecular);

    if (all(fsL == 0))
        return fsL;
	
	// compute fsE
    wi = -connectDir; // towards light
    wo = normalize(cameraPath[cameraIndex - 2].posW - cameraEndV.posW); // towards camera
    float3 fsE = evalBRDF(wi, wo, cameraEndV.posW, cameraEndV.N, cameraEndV.N, cameraEndV.dif, cameraEndV.spec, cameraEndV.rough, cameraEndV.isSpecular);
    
    if (all(fsE == 0))
        return fsE;

    float3 cst = fsL * G * fsE;
	
    return aL * cst * aE;
}
