void shootRay(inout RayPayload payload)
{
	// Setup our indirect ray
    RayDesc rayColor;
    rayColor.Origin = payload.rayOrigin; // Where does it start?
    rayColor.Direction = payload.rayDir; // What direction do we shoot it?
    rayColor.TMin = gMinT; // The closest distance we'll count as a hit
    rayColor.TMax = 1.0e38f; // The farthest distance we'll count as a hit

	// Trace our ray to get a color in the indirect direction.  Use hit group #2 and miss shader #2
    TraceRay(gRtScene, 0, 0xFF, 1, hitProgramCount, 1, rayColor, payload);
}

[shader("miss")]
void RayMiss(inout RayPayload rayData)
{
    rayData.color = float3(0, 0, 0);
    rayData.terminated = true;
}

[shader("anyhit")]
void RayAnyHit(inout RayPayload rayData, BuiltInTriangleIntersectionAttributes attribs)
{
	// Is this a transparent part of the surface?  If so, ignore this hit
    if (alphaTestFails(attribs))
        IgnoreHit();
}

void handleIndirectRayHit(float3 hit, float3 N, float3 noNormalN, float3 V, float3 dif, float3 spec, float rough, inout RayPayload rayData)
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
        
        float3 color = dif / probDiffuse;
        float NdotL = saturate(dot(N, L));
        float pdfForward = (NdotL * M_1_PI) * probDiffuse;
        updateRayData(rayData, color, hit, L, N, V, dif, spec, rough, false, pdfForward);
    }
	// Otherwise we randomly selected to sample our GGX lobe
    else
    {
        float3 H = getGGXMicrofacet(rayData.rndSeed, rough, N);
        float3 L = normalize(2.f * dot(V, H) * H - V);
        if (dot(noNormalN, L) <= 0.0f)
            return;
        float NdotL = saturate(dot(N, L));
        float ggxProb;
        float3 ggxTerm = ggxLighting(H, L, N, NdotL, NdotV, rough, spec, ggxProb);

        float3 color = NdotL * ggxTerm / (ggxProb * (1.0f - probDiffuse));
        float pdfForward = ggxProb * (1.0f - probDiffuse);
        updateRayData(rayData, color, hit, L, N, V, dif, spec, rough, true, pdfForward);
    }
}

[shader("closesthit")]
void RayClosestHit(inout RayPayload rayData, BuiltInTriangleIntersectionAttributes attribs)
{
	// Run a helper functions to extract Falcor scene data for shading
    ShadingData shadeData = getHitShadingData(attribs, WorldRayOrigin());
    handleIndirectRayHit(shadeData.posW, shadeData.N, shadeData.N, shadeData.V,
			shadeData.diffuse, shadeData.specular, shadeData.roughness, rayData);
}
