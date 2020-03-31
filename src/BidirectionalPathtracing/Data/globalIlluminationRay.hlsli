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
    float3 L;
    float pdf;
    bool isSpecular;
    float3 color = sampleBRDF(rayData.rndSeed, hit, N, noNormalN, V, dif, spec, rough, L, pdf, isSpecular);
    updateRayData(rayData, color, hit, L, N, V, dif, spec, rough, isSpecular, pdf);
}

[shader("closesthit")]
void RayClosestHit(inout RayPayload rayData, BuiltInTriangleIntersectionAttributes attribs)
{
	// Run a helper functions to extract Falcor scene data for shading
    ShadingData shadeData = getHitShadingData(attribs, WorldRayOrigin());
    handleIndirectRayHit(shadeData.posW, shadeData.N, shadeData.N, shadeData.V,
			shadeData.diffuse, shadeData.specular, shadeData.roughness, rayData);
}
