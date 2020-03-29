// The payload structure for our indirect rays
struct LightRayPayload
{
    float3 accumulated;
    float3 color;    // throughput at each path vertex

    float3 posW;
    float3 N;

    uint   rndSeed; // random seed, so we pick uncorrelated RNGs along our ray
    float3 rayOrigin;
    float3 rayDir;
    
    bool   terminated;
    
};

void shootLightRay(inout LightRayPayload payload)
{
	// Setup our indirect ray
	RayDesc rayColor;
	rayColor.Origin = payload.rayOrigin;  // Where does it start?
	rayColor.Direction = payload.rayDir;  // What direction do we shoot it?
	rayColor.TMin = gMinT;         // The closest distance we'll count as a hit
	rayColor.TMax = 1.0e38f;      // The farthest distance we'll count as a hit

	// Trace our ray to get a color in the indirect direction.  Use hit group #1 and miss shader #1
	TraceRay(gRtScene, 0, 0xFF, 1, hitProgramCount, 1, rayColor, payload);
}

[shader("miss")]
void LightRayMiss(inout LightRayPayload rayData)
{	
    rayData.terminated = true;
}

[shader("anyhit")]
void LightRayAnyHit(inout LightRayPayload rayData, BuiltInTriangleIntersectionAttributes attribs)
{
	// Is this a transparent part of the surface?  If so, ignore this hit
	if (alphaTestFails(attribs))
		IgnoreHit();
}

void handleLightRayHit(float3 hit, float3 N, float3 noNormalN, float3 V, float3 dif, float3 spec, float rough, inout LightRayPayload rayData)
{
	// We have to decide whether we sample our diffuse or specular/ggx lobe.
	float probDiffuse = probabilityToSampleDiffuse(dif, spec);
	float chooseDiffuse = (nextRand(rayData.rndSeed) < probDiffuse);

	// We'll need NdotV for both diffuse and specular...
	float NdotV = saturate(dot(N, V));

	// If we randomly selected to sample our diffuse lobe...
	if (chooseDiffuse)
	{
        // Shoot a randomly selected cosine-sampled diffuse ray.
        float3 L = getCosHemisphereSample(rayData.rndSeed, N);
		
		// Check to make sure our randomly selected, normal mapped diffuse ray didn't go below the surface.
        if (dot(noNormalN, L) <= 0.0f)
            return;
        
        rayData.accumulated = rayData.color * dif;
        rayData.color *= dif / probDiffuse;
        rayData.rayOrigin = hit;
        rayData.rayDir = L;
        rayData.posW = hit;
        rayData.N = N;
    }
	// Otherwise we randomly selected to sample our GGX lobe
	else
	{
        float3 H = getGGXMicrofacet(rayData.rndSeed, rough, N);
        
        float3 L = normalize(gCamera.posW - hit);

        if (dot(noNormalN, L) <= 0.0f)
            return;
        
        float ggxProb;
        float NdotL = saturate(dot(N, L));
        float3 ggxTerm = ggxLighting(H, L, N, NdotL, NdotV, rough, spec, ggxProb);
        
        rayData.accumulated = rayData.color * NdotL * ggxTerm;
        

        L = normalize(2.f * dot(V, H) * H - V);
        if (dot(noNormalN, L) <= 0.0f)
            return;
        NdotL = saturate(dot(N, L));
        ggxTerm = ggxLighting(H, L, N, NdotL, NdotV, rough, spec, ggxProb);

        rayData.color *= NdotL * ggxTerm / (ggxProb * (1.0f - probDiffuse));
        rayData.posW = hit;
        rayData.N = N;
        rayData.rayOrigin = hit;
        rayData.rayDir = L;
    }
}

[shader("closesthit")]
void LightRayClosestHit(inout LightRayPayload rayData, BuiltInTriangleIntersectionAttributes attribs)
{
	// Run a helper functions to extract Falcor scene data for shading
	ShadingData shadeData = getHitShadingData( attribs, WorldRayOrigin() );
    handleLightRayHit(shadeData.posW, shadeData.N, shadeData.N, shadeData.V, 
				shadeData.diffuse, shadeData.specular, shadeData.roughness, rayData);
}
