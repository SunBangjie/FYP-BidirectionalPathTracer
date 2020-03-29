// Define pi
#define M_1_PI  0.318309886183790671538

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

    ExplicitLodTextureSampler lodSampler = { 0 }; // Specify the tex lod/mip to use here

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
ShadingData getHitShadingData(BuiltInTriangleIntersectionAttributes attribs, float3 cameraPos)
{
	// Run a pair of Falcor helper functions to compute important data at the current hit point
    VertexOut vsOut = getVertexAttributes(PrimitiveIndex(), attribs);
    return simplePrepareShadingData(vsOut, gMaterial, cameraPos);
}

// Utility function to get a vector perpendicular to an input vector 
//    (from "Efficient Construction of Perpendicular Vectors Without Branching")
float3 getPerpendicularVector(float3 u)
{
    float3 a = abs(u);
    uint xm = ((a.x - a.y) < 0 && (a.x - a.z) < 0) ? 1 : 0;
    uint ym = (a.y - a.z) < 0 ? (1 ^ xm) : 0;
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

// This function tests if the alpha test fails, given the attributes of the current hit. 
//   -> Can legally be called in a DXR any-hit shader or a DXR closest-hit shader, and 
//      accesses Falcor helpers and data structures to extract and perform the alpha test.
bool alphaTestFails(BuiltInTriangleIntersectionAttributes attribs)
{
	// Run a Falcor helper to extract the current hit point's geometric data
    VertexOut vsOut = getVertexAttributes(PrimitiveIndex(), attribs);

    // Extracts the diffuse color from the material (the alpha component is opacity)
    ExplicitLodTextureSampler lodSampler = { 0 }; // Specify the tex lod/mip to use here
    float4 baseColor = sampleTexture(gMaterial.resources.baseColor, gMaterial.resources.samplerState,
        vsOut.texC, gMaterial.baseColor, EXTRACT_DIFFUSE_TYPE(gMaterial.flags), lodSampler);

	// Test if this hit point fails a standard alpha test.  
    return (baseColor.a < gMaterial.alphaThreshold);
}