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
    float rough;
    
    float pdfForward;
    float pdfBackward;
	
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
        v.pdfForward = 0.0f;
        v.pdfBackward = 0.0f;
        return v;
    }
	
    static PathVertex create(float3 color, float3 posW, float3 N, float3 V, float3 dif, float3 spec, float rough, float pdfForward, float pdfBackward)
    {
        PathVertex v;
        v.color = color;
        v.posW = posW;
        v.N = N;
        v.V = V;
        v.dif = dif;
        v.spec = spec;
        v.rough = rough;
        v.pdfForward = pdfForward;
        v.pdfBackward = pdfBackward;
        return v;
    }
};

struct MisNode
{
    float pTowardCamera;
    float pTowardLight;
    bool  isSpecular;
    
    static MisNode init()
    {
        MisNode n;
        n.pTowardCamera = 0.0;
        n.pTowardLight = 0.0;
        n.isSpecular = false;
        return n;
    }
};

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

float3 evalGGXBRDF(float3 wi, float3 wo, float3 hit, float3 N, float3 noNormalN, float3 dif, float3 spec, float rough)
{
    float3 L = wo;
    float3 V = wi;
    
    if (!all(spec > 0.f))
    {    
        // Check to make sure our randomly selected, normal mapped diffuse ray didn't go below the surface.
        if (dot(noNormalN, L) <= 0.0f)
            return float3(0.0f);
    
        return dif * M_1_PI;
    }
    else
    {
        return float3(0, 0, 0);
    }
}

float evalGGXPdf(float3 wi, float3 wo, float3 hit, float3 N, float3 noNormalN, float3 dif, float3 spec, float rough)
{
    float3 L = wo;
    float3 V = wi;
    float probDiffuse = probabilityToSampleDiffuse(dif, spec);
    
    if (!all(spec > 0.f))
    {
        // Check to make sure our randomly selected, normal mapped diffuse ray didn't go below the surface.
        if (dot(noNormalN, L) <= 0.0f)
            return 0.0f;
        
        float NdotL = saturate(dot(N, L));
        return (NdotL * M_1_PI) * probDiffuse;
    }
    else
    {
        return 0.f;
    }
}

float evalG(PathVertex a, PathVertex b)
{
    float3 posA = a.posW;
    float3 posB = b.posW;
    float3 vecAB = posB - posA;
    float invLengthAB = 1.0f / length(vecAB);
    float3 dirAB = vecAB * invLengthAB;

    float cosA = abs(dot(a.N, dirAB));
    float cosB = abs(dot(b.N, dirAB));
	
    float3 visibilityTestDir = normalize(posA - posB);
    float visibilityTestDistance = distance(posA, posB);
    float V = shadowRayVisibility(posB, visibilityTestDir, gMinT, visibilityTestDistance);
    return V * cosA * cosB * invLengthAB * invLengthAB;
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

float3 getUnweightedContribution(PathVertex cameraPath[4], PathVertex lightPath[4], uint cameraIndex, uint lightIndex)
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
    wi = connectDir;
    wo = normalize(lightPath[lightIndex - 2].posW - lightEndV.posW); // towards light
    float3 fsL = evalGGXBRDF(wi, wo, lightEndV.posW, lightEndV.N, lightEndV.N, lightEndV.dif, lightEndV.spec, lightEndV.rough);

    if (all(fsL == 0))
        return fsL;
	
	// compute fsE
    wi = -connectDir;
    wo = normalize(cameraPath[cameraIndex - 2].posW - cameraEndV.posW); // towards camera
    float3 fsE = evalGGXBRDF(wi, wo, cameraEndV.posW, cameraEndV.N, cameraEndV.N, cameraEndV.dif, cameraEndV.spec, cameraEndV.rough);
    
    if (all(fsE == 0))
        return fsE;
	
	// compute cst by fsL * G * fsE
    float G = evalG(lightEndV, cameraEndV);

    float3 cst = fsL * G * fsE;
	
    return aL * cst * aE;
}

float getMISWeight(PathVertex cameraPath[4], PathVertex lightPath[4], uint cameraIndex, uint lightIndex, MisNode misNodes[8])
{
    if (cameraIndex < 1 || lightIndex < 1)
        return 0;
    
    PathVertex lightEndV = lightPath[lightIndex - 1];
    PathVertex cameraEndV = cameraPath[cameraIndex - 1];
    
    float pdfLightEndForward, pdfLightEndBackward, pdfCameraEndForward, pdfCameraEndBackward;
    
    float3 pLightEnd = lightEndV.posW;
    float3 nLightEnd = lightEndV.N;
    float3 pCameraEnd = cameraEndV.posW;
    float3 nCameraEnd = cameraEndV.N;
    
    float3 dLightToCamera = normalize(pCameraEnd - pLightEnd);
    if (lightIndex == 1)
    {
        // TODO: change the probability to different light source type
        float pdfW = M_1_4_PI;
        pdfLightEndBackward = pdfW;
        pdfLightEndForward = pdfW;
    }
    else
    {
        float3 dLightEndTowardsLight = normalize(lightPath[lightIndex - 2].posW - pLightEnd);
        pdfLightEndForward = evalGGXPdf(dLightToCamera, dLightEndTowardsLight, lightEndV.posW, 
                lightEndV.N, lightEndV.N, lightEndV.dif, lightEndV.spec, lightEndV.rough) / dot(dLightToCamera, nLightEnd);
        
        pdfLightEndBackward = evalGGXPdf(dLightEndTowardsLight, dLightToCamera, lightEndV.posW,
                lightEndV.N, lightEndV.N, lightEndV.dif, lightEndV.spec, lightEndV.rough) / dot(dLightEndTowardsLight, nLightEnd);
    }
    
    float3 dCameraToLight = -dLightToCamera;
    if (cameraIndex == 1)
    {
        // TODO: change the probability to different camera model
        float pdfW = 1.0f;
        pdfCameraEndForward = pdfW / dot(nCameraEnd, dCameraToLight);
        pdfCameraEndBackward = cameraEndV.pdfBackward;
    }
    else
    {
        float3 dCameraEndTowardsCamera = normalize(cameraPath[cameraIndex - 2].posW - pCameraEnd);
        pdfCameraEndForward = evalGGXPdf(dCameraToLight, dCameraEndTowardsCamera, cameraEndV.posW,
                cameraEndV.N, cameraEndV.N, cameraEndV.dif, cameraEndV.spec, cameraEndV.rough) / dot(dCameraToLight, nCameraEnd);
        pdfCameraEndBackward = evalGGXPdf(dCameraEndTowardsCamera, dCameraToLight, cameraEndV.posW,
                cameraEndV.N, cameraEndV.N, cameraEndV.dif, cameraEndV.spec, cameraEndV.rough) / dot(dCameraEndTowardsCamera, nCameraEnd);
    }
    
    // construct misNodes
    
    int total = lightIndex + cameraIndex - 1;
    for (uint i = 0; i < lightIndex - 1; ++i)
    {
        misNodes[i].pTowardLight = (i == 0) ? lightPath[0].pdfBackward : lightPath[i].pdfBackward * evalGWithoutV(lightPath[i], lightPath[i - 1]);
        misNodes[i].pTowardCamera = lightPath[i].pdfForward * evalGWithoutV(lightPath[i + 1], lightPath[i]);
        misNodes[i].isSpecular = all(lightPath[i].spec > 0.f);
    }
    if (lightIndex > 0)
    {
        misNodes[lightIndex - 1].pTowardLight = (lightIndex == 1) ? pdfLightEndBackward : pdfLightEndBackward * evalGWithoutV(lightPath[lightIndex - 1], lightPath[lightIndex - 2]);
        misNodes[lightIndex - 1].pTowardCamera = (lightIndex - 1 == total) ? pdfLightEndForward : pdfLightEndForward * evalGWithoutV(lightEndV, cameraEndV);
        misNodes[lightIndex - 1].isSpecular = all(lightPath[lightIndex - 1].spec > 0.f);
    }
    
    for (uint i = 0; i < cameraIndex - 1; ++i)
    {
        misNodes[total - i].pTowardCamera = (i == 0) ? cameraPath[0].pdfBackward : cameraPath[i].pdfBackward * evalGWithoutV(cameraPath[i], cameraPath[i - 1]);
        misNodes[total - i].pTowardLight = cameraPath[i].pdfForward * evalGWithoutV(cameraPath[i + 1], cameraPath[i]);
        misNodes[total - i].isSpecular = all(cameraPath[i].spec > 0.f);
    }
    if (cameraIndex > 0)
    {
        misNodes[total - lightIndex + 1].pTowardCamera = (cameraIndex == 1) ? pdfCameraEndBackward : pdfCameraEndBackward * evalGWithoutV(cameraPath[lightIndex - 1], cameraPath[lightIndex - 2]);
        misNodes[total - lightIndex + 1].pTowardLight = (cameraIndex - 1 == total) ? pdfCameraEndForward : pdfCameraEndForward * evalGWithoutV(lightEndV, cameraEndV);
        misNodes[total - lightIndex + 1].isSpecular = all(cameraPath[cameraIndex - 1].spec > 0.f);
    }
    
    // iterate from the connection end point to calculate relative pdfA and add it to misWeightSum (power heuristic)
    float pK = 1.0f;
    float misWeightSum = 1.0f;
    for (int i = lightIndex; i < total; ++i)
    {
        if (i == 0)
        {
            pK *= misNodes[0].pTowardLight / misNodes[1].pTowardLight;
            if (misNodes[1].isSpecular)
                continue;
        }
        else
        {
            pK *= misNodes[i - 1].pTowardCamera / misNodes[i + 1].pTowardLight;
            if (misNodes[i].isSpecular || misNodes[i+1].isSpecular)
                continue;
        }
    }
    misWeightSum += pK * pK;
    
    pK = 1.0f;
    for (int i = lightIndex; i > 0; --i)
    {
        if (i == (total + 1))
        {
            pK *= misNodes[total].pTowardCamera / misNodes[total - 1].pTowardCamera;
            if (misNodes[total - 1].isSpecular)
                continue;
        }
        else
        {
            pK *= misNodes[i].pTowardLight / misNodes[i - 2].pTowardCamera;
            if (misNodes[i - 1].isSpecular || misNodes[i - 2].isSpecular)
                continue;
        }
        misWeightSum += pK * pK;
    }
    
    return 1.0f / misWeightSum;
}
