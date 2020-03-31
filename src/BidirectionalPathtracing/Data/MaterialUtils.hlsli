// Define 1/pi
#define M_1_PI  0.318309886183790671538
#define M_1_4_PI  0.07957747154

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

float3 ggxDirect(inout uint rndSeed, float3 hit, float3 N, float3 V, float3 dif, float3 spec, float rough)
{
	// Pick a random light from our scene to shoot a shadow ray towards
    int lightToSample = min(int(nextRand(rndSeed) * gLightsCount), gLightsCount - 1);

	// Query the scene to find info about the randomly selected light
    float distToLight;
    float3 lightIntensity;
    float3 L;
    getLightData(lightToSample, hit, L, lightIntensity, distToLight);

	// Compute our lambertion term (N dot L)
    float NdotL = saturate(dot(N, L));

	// Shoot our shadow ray to our randomly selected light
    bool vis = shadowRayVisibility(hit, L, gMinT, distToLight);
    float shadowMult = vis ? float(gLightsCount) : 0.f;

	// Compute half vectors and additional dot products for GGX
    float3 H = normalize(V + L);
    float NdotH = saturate(dot(N, H));
    float LdotH = saturate(dot(L, H));
    float NdotV = saturate(dot(N, V));

	// Evaluate terms for our GGX BRDF model
    float D = ggxNormalDistribution(NdotH, rough);
    float G = ggxSchlickMaskingTerm(NdotL, NdotV, rough);
    float3 F = schlickFresnel(spec, LdotH);

	// Evaluate the Cook-Torrance Microfacet BRDF model
	//     Cancel out NdotL here & the next eq. to avoid catastrophic numerical precision issues.
    float3 ggxTerm = D * G * F / (4 * NdotV /* * NdotL */);

	// Compute our final color (combining diffuse lobe plus specular GGX lobe)
    return shadowMult * lightIntensity * ( /* NdotL * */ggxTerm + NdotL * dif / M_PI);
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

float3 lambertianDirect(inout uint rndSeed, float3 hit, float3 norm, float3 difColor)
{
	// Pick a random light from our scene to shoot a shadow ray towards
    int lightToSample = min(int(nextRand(rndSeed) * gLightsCount), gLightsCount - 1);

	// Query the scene to find info about the randomly selected light
    float distToLight;
    float3 lightIntensity;
    float3 toLight;
    getLightData(lightToSample, hit, toLight, lightIntensity, distToLight);

	// Compute our lambertion term (L dot N)
    float LdotN = saturate(dot(norm, toLight));

	// Shoot our shadow ray to our randomly selected light
    float shadowMult = float(gLightsCount) * shadowRayVisibility(hit, toLight, gMinT, distToLight);

	// Return the Lambertian shading color using the physically based Lambertian term (albedo / pi)
    return shadowMult * LdotN * lightIntensity * difColor / M_PI;
}

float3 lambertianBRDF(inout uint rndSeed, float3 hit, float3 norm, float3 difColor, out float3 L)
{
	// Shoot a randomly selected cosine-sampled diffuse ray.
    L = getCosHemisphereSample(rndSeed, norm);

	// Accumulate the color: (NdotL * incomingLight * difColor / pi) 
	// Probability of sampling:  (NdotL / pi)
    return difColor;
}

float lambertianPdf(float3 N, float3 L)
{
    return saturate(dot(N, L) * M_1_PI);
}

