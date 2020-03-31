// Define 1/pi
#define M_1_PI  0.318309886183790671538
#define M_1_4_PI  0.07957747154

struct PathVertex
{
    float3 color;
	
    float3 posW;
    float3 N;
    float3 V;
	
    float3 dif;
    float3 spec;
    float  rough;
    bool   isSpecular;
    
    float  pdfForward;
	
    static PathVertex init()
    {
        PathVertex v;
        v.color = float3(0, 0, 0);
        v.posW = float3(0, 0, 0);
        v.N = float3(0, 0, 0);
        v.V = float3(0, 0, 0);
        v.dif = float3(0, 0, 0);
        v.spec = float3(0, 0, 0);
        v.rough = 0.0f;
        v.isSpecular = false;
        v.pdfForward = 0.0f;
        return v;
    }
	
    static PathVertex create(float3 color, float3 posW, float3 N, float3 V, float3 dif, float3 spec, float rough, bool isSpecular, float pdfForward)
    {
        PathVertex v;
        v.color = color;
        v.posW = posW;
        v.N = N;
        v.V = V;
        v.dif = dif;
        v.spec = spec;
        v.rough = rough;
        v.isSpecular = isSpecular;
        v.pdfForward = pdfForward;
        return v;
    }
};

float3 ggxDirectWrapper(PathVertex v, inout uint rndSeed)
{
    return ggxDirect(rndSeed, v.posW, v.N, v.V, v.dif, v.spec, v.rough);
}

float3 connectToCamera(PathVertex v)
{
    return evalGGXBRDF(v.V, normalize(gCamera.posW - v.posW), v.posW, v.N, v.N, v.dif, v.spec, v.rough, v.isSpecular);
}

float3 clampVec(float3 v)
{
    return float3(clamp(v.x, 0, 0.5), clamp(v.y, 0, 0.5), clamp(v.z, 0, 0.5));
}

// Our material has have both a diffuse and a specular lobe.  
//     With what probability should we sample the diffuse one?
float probabilityToSampleDiffuse(float3 difColor, float3 specColor)
{
    float lumDiffuse = max(0.01f, luminance(difColor.rgb));
    float lumSpecular = max(0.01f, luminance(specColor.rgb));
    return lumDiffuse / (lumDiffuse + lumSpecular);
}

// A helper to extract important light data from internal Falcor data structures.  What's going on isn't particularly
//     important -- any framework you use will expose internal scene data in some way.  Use your framework's utilities.
void getLightData(in int index, in float3 hitPos, out float3 toLight, out float3 lightIntensity, out float distToLight)
{
	// Use built-in Falcor functions and data structures to fill in a LightSample data structure
	//   -> See "Lights.slang" for it's definition
    LightSample ls;

	// Is it a directional light?
    if (gLights[index].type == LightDirectional)
        ls = evalDirectionalLight(gLights[index], hitPos);

	// No?  Must be a point light.
    else
        ls = evalPointLight(gLights[index], hitPos);

	// Convert the LightSample structure into simpler data
    toLight = normalize(ls.L);
    lightIntensity = ls.diffuse;
    distToLight = length(ls.posW - hitPos);
}

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

// Utility function to get a vector perpendicular to an input vector 
//    (from "Efficient Construction of Perpendicular Vectors Without Branching")
float3 getPerpendicularVector(float3 u)
{
	float3 a = abs(u);
	uint xm = ((a.x - a.y)<0 && (a.x - a.z)<0) ? 1 : 0;
	uint ym = (a.y - a.z)<0 ? (1 ^ xm) : 0;
	uint zm = 1 ^ (xm | ym);
	return cross(u, float3(xm, ym, zm));
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

// Get a cosine-weighted random vector centered around a specified normal direction.
float3 getCosHemisphereSample(inout uint randSeed, float3 hitNorm)
{
	// Get 2 random numbers to select our sample with
    float2 randVal = float2(nextRand(randSeed), nextRand(randSeed));

	// Cosine weighted hemisphere sample from RNG
    float3 bitangent = getPerpendicularVector(hitNorm);
    float3 tangent = cross(bitangent, hitNorm);
    float r = sqrt(randVal.x);
    float phi = 2.0f * 3.14159265f * randVal.y;

	// Get our cosine-weighted hemisphere lobe sample direction
    return tangent * (r * cos(phi).x) + bitangent * (r * sin(phi)) + hitNorm.xyz * sqrt(max(0.0, 1.0f - randVal.x));
}

float3 sampleUnitSphere(inout uint rndSeed)
{
	// rejection sampling
    float3 p = float3(2.0f, 2.0f, 2.0f);
    while (length(p) > 1.0f)
        p = float3(nextRand(rndSeed) * 2.0f - 1.0f, nextRand(rndSeed) * 2.0f - 1.0f, nextRand(rndSeed) * 2.0f - 1.0f);
    return p;
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

float3 evalGGXBRDF(float3 V, float3 L, float3 hit, float3 N, float3 noNormalN, float3 dif, float3 spec, float rough, bool isSpecular)
{
    // check wi and wo, if wo is the reflected vector of wi, then we take specular
    if (!isSpecular)
    {
        // Check to make sure our randomly selected, normal mapped diffuse ray didn't go below the surface.
        if (dot(noNormalN, L) <= 0.0f)
            return float3(0.0f);
    
        return dif * M_1_PI;
    }
    else
    {
        float3 H = normalize(L + V);
        if (dot(noNormalN, L) <= 0.0f)
            return float3(0, 0, 0);
        float NdotL = saturate(dot(N, L));
        float NdotV = saturate(dot(N, V));
        float ggxProb;
        return ggxLighting(H, L, N, NdotL, NdotV, rough, spec, ggxProb);
    }
}

float evalGGXPdf(float3 V, float3 L, float3 hit, float3 N, float3 noNormalN, float3 dif, float3 spec, float rough, bool isSpecular)
{
    float probDiffuse = probabilityToSampleDiffuse(dif, spec);
    
    if (!isSpecular)
    {
        // Check to make sure our randomly selected, normal mapped diffuse ray didn't go below the surface.
        if (dot(noNormalN, L) <= 0.0f)
            return 0.0f;
        
        float NdotL = saturate(dot(N, L));
        return (NdotL * M_1_PI) * probDiffuse;
    }
    else
    {
        float3 H = normalize(L + V);
        if (dot(noNormalN, L) <= 0.0f)
            return 0.f;
        float NdotL = saturate(dot(N, L));
        float NdotV = saturate(dot(N, V));
        float ggxProb;
        ggxLighting(H, L, N, NdotL, NdotV, rough, spec, ggxProb);
        
        return ggxProb * (1.0f - probDiffuse);
    }
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

float3 getUnweightedContribution(PathVertex cameraPath[4], PathVertex lightPath[4], uint cameraIndex, uint lightIndex, float G)
{
	// this function only computes unweighted contribution for path length > 1
    if (cameraIndex <= 1 || lightIndex <= 1)
        return float3(0, 0, 0);
	
	// get end vertex from each path
    PathVertex cameraEndV = cameraPath[cameraIndex - 1];
    PathVertex lightEndV = lightPath[lightIndex - 1];
	
	// get throughput from the last vertex of each path
    float3 aE = cameraEndV.color;
    float3 aL = lightEndV.color;
	
	// connect two end vertices
    float3 connectDir = normalize(cameraEndV.posW - lightEndV.posW); // from light to camera

    float3 wi, wo;
    
	// compute fsL
    wi = connectDir; // towards camera
    wo = normalize(lightPath[lightIndex - 2].posW - lightEndV.posW); // towards light
    float3 fsL = evalGGXBRDF(wi, wo, lightEndV.posW, lightEndV.N, lightEndV.N, lightEndV.dif, lightEndV.spec, lightEndV.rough, lightEndV.isSpecular);

    if (all(fsL == 0))
        return fsL;
	
	// compute fsE
    wi = -connectDir; // towards light
    wo = normalize(cameraPath[cameraIndex - 2].posW - cameraEndV.posW); // towards camera
    float3 fsE = evalGGXBRDF(wi, wo, cameraEndV.posW, cameraEndV.N, cameraEndV.N, cameraEndV.dif, cameraEndV.spec, cameraEndV.rough, cameraEndV.isSpecular);
    
    if (all(fsE == 0))
        return fsE;

    float3 cst = fsL * G * fsE;
	
    return aL * cst * aE;
}
