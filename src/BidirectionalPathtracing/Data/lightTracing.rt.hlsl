// Some shared Falcor stuff for talking between CPU and GPU code
#include "HostDeviceSharedMacros.h"
#include "HostDeviceData.h"           

// Include and import common Falcor utilities and data structures
import Raytracing;                   // Shared ray tracing specific functions & data
import ShaderCommon;                 // Shared shading data structures
import Shading;                      // Shading functions, etc     
import Lights;                       // Light structures for our current scene

// A constant buffer we'll populate from our C++ code  (used for our ray generation shader)
shared cbuffer GlobalCB
{
	float gMinT;           // Min distance to start a ray to avoid self-occlusion
	uint  gFrameCount;     // An integer changing every frame to update the random number
    bool  gDoDirectGI;
	bool  gDoIndirectGI;   // A boolean determining if we should shoot indirect GI rays
	uint  gMaxDepth;       // Maximum number of recursive bounces to allow
    float gEmitMult;       // Multiply emissive amount by this factor (set to 1, usually)
    float2 gPixelJitter;   // Camera jitter
}

// Input and out textures that need to be set by the C++ code (for the ray gen shader)
shared Texture2D<float4> gPos;
shared Texture2D<float4> gNorm;
shared Texture2D<float4> gDiffuseMatl;
shared Texture2D<float4> gSpecMatl;
shared Texture2D<float4> gExtraMatl;
shared Texture2D<float4> gEnvMap;
shared Texture2D<float4> gEmissive;
shared RWTexture2D<float4> gOutput;

// A separate file with some simple utility functions: getPerpendicularVector(), initRand(), nextRand()
#include "lightTracingUtils.hlsli"

// Include implementations of GGX normal distribution function, Fresnel approx,
//     masking term and function to sampl NDF 
#include "microfacetBRDFUtils.hlsli"

// Include shader entries, data structures, and utility functions to spawn rays
#include "standardShadowRay.hlsli"
#include "pathTracingRay.hlsli"

IndirectRayPayload initPayload(float3 rayOrigin, float3 rayDir, float3 color, uint randSeed)
{
    IndirectRayPayload payload;
    payload.rayOrigin = rayOrigin;
    payload.rayDir = rayDir;
    payload.rndSeed = randSeed;
    payload.color = color;
    payload.posW = rayOrigin;
    payload.N = float3(0.0f);
    payload.V = float3(0.0f);
    payload.dif = float3(0.0f);
    payload.spec = float3(0.0f);
    payload.rough = 0.0;
    payload.pdfForward = 0.0;
    payload.pdfBackward = 0.0;
    payload.terminated = false;
    return payload;
}

float3 ggxDirectWrapper(PathVertex v, inout uint rndSeed)
{
    return ggxDirect(rndSeed, v.posW, v.N, v.V, v.dif, v.spec, v.rough);
}

float3 connectToCamera(PathVertex v)
{
    return evalGGXBRDF(v.V, normalize(gCamera.posW - v.posW), v.posW, v.N, v.N, v.dif, v.spec, v.rough);
}

// How do we shade our g-buffer and spawn indirect and shadow rays?
[shader("raygeneration")]
void SimpleDiffuseGIRayGen()
{
	// Where is this ray on screen?
	uint2 launchIndex    = DispatchRaysIndex().xy;
	uint2 launchDim      = DispatchRaysDimensions().xy;
    
    // Load g-buffer data
	// g-buffer is used to optimize camera path construction
    float4 worldPos = gPos[launchIndex];
    float4 worldNorm = gNorm[launchIndex];
    float4 difMatlColor = gDiffuseMatl[launchIndex];
    float4 specMatlColor = gSpecMatl[launchIndex];
    float4 extraData = gExtraMatl[launchIndex];
    float4 pixelEmissive = gEmissive[launchIndex];
    
    // Does this g-buffer pixel contain a valid piece of geometry?  (0 in pos.w for invalid)
    bool isGeometryValid = (worldPos.w != 0.0f);

	// Extract and compute some material and geometric parameters
    float roughness = specMatlColor.a * specMatlColor.a;
    float3 V = normalize(gCamera.posW - worldPos.xyz);

	// Make sure our normal is pointed the right direction
    if (dot(worldNorm.xyz, V) <= 0.0f)
        worldNorm.xyz = -worldNorm.xyz;
    float NdotV = dot(worldNorm.xyz, V);
    float3 noMapN = normalize(extraData.yzw);
    if (dot(noMapN, V) <= 0.0f)
        noMapN = -noMapN;

	// If we don't hit any geometry, our difuse material contains our background color.
    float3 shadeColor = isGeometryValid ? float3(0, 0, 0) : difMatlColor.rgb;
    
    uint randSeed = initRand(launchIndex.x + launchIndex.y * launchDim.x, gFrameCount, 16);
    
    /**
	 ** Construct camera path by shooting rays from the camera
	 ** TODO: Now assume pinhole camera is used, but we can extend it to thin lens camera
	 */
	
    PathVertex cameraPath[4];
    PathVertex lightPath[4];
    for (uint i = 0; i < 4; i++)
    {
        cameraPath[i] = PathVertex.init();
        lightPath[i] = PathVertex.init();
    }
    
	// The first vertex is the camera
    cameraPath[0] = PathVertex.init();
    cameraPath[0].posW = gCamera.posW;
    cameraPath[0].color = float3(1.0f);
    
	// Do shading, if we have geoemtry here (otherwise, output the background color)
    if (isGeometryValid)
    {	
        IndirectRayPayload payload = initPayload(gCamera.posW, normalize(worldPos.xyz - gCamera.posW), float3(1.0f), randSeed);
        
        for (uint depth = 0; depth < gMaxDepth && !payload.terminated; depth++)
        {
            shootIndirectRay(payload);
            cameraPath[depth + 1] = PathVertex.create(payload.color, payload.posW, payload.N, payload.V, payload.dif, payload.spec, payload.rough, payload.pdfForward, payload.pdfBackward);
        }
        
        randSeed = payload.rndSeed;
    }
    
    /**
	 ** Construct light path by shooting rays from the light source
	 ** TODO: Now assume point light source is used, but we can extend it to area and directional light
	 */
    
    float3 lightOrigin, lightDir, lightIntensity;
    sampleLight(randSeed, lightOrigin, lightDir, lightIntensity);
    bool takeContribution[4] = { true, true, true, true };
    
    // first vertex is the light sample
    lightPath[0].posW = lightOrigin;
    lightPath[0].color = lightIntensity;
    
    IndirectRayPayload payload = initPayload(lightOrigin, lightDir, lightIntensity, randSeed);
        
    for (uint depth = 0; depth < gMaxDepth && !payload.terminated; depth++)
    {
        shootIndirectRay(payload);
        lightPath[depth + 1] = PathVertex.create(payload.color, payload.posW, payload.N, payload.V, payload.dif, payload.spec, payload.rough, payload.pdfForward, payload.pdfBackward);
        takeContribution[depth + 1] = !payload.terminated;
    }
        
    randSeed = payload.rndSeed;
    
    // add contributions
    
    MisNode misNodes[8];
    for (uint i = 0; i < 8; i++)
    {
        misNodes[i] = MisNode.init();
    }
    
    for (uint i = 0; i < gMaxDepth; i++)
    {
        shadeColor = cameraPath[i].color * ggxDirectWrapper(cameraPath[i + 1], randSeed);
        float weight = getMISWeight(cameraPath, lightPath, i + 2, 1, misNodes);
        shadeColor = saturate(shadeColor);
        bool colorsNan = any(isnan(shadeColor));
        gOutput[launchIndex] = saturate(gOutput[launchIndex] + float4(colorsNan ? float3(0, 0, 0) : shadeColor, 1.0f));
    }
    
    for (uint i = 0; i < gMaxDepth && takeContribution[i+1]; i++)
    {
        shadeColor = lightPath[i].color * connectToCamera(lightPath[i + 1]);
        float3 lastHitPos = lightPath[i + 1].posW;
        float3 lastHitN = lightPath[i + 1].N;

        float3 dirToCamera = normalize(gCamera.posW - lastHitPos);
        float disToCamera = length(gCamera.posW - lastHitPos);
        float3 cameraN = normalize(gCamera.cameraW);

        if (dot(cameraN, dirToCamera) < 0 && takeContribution[i+1])
        {
            // test visibility towards camera
            float vis = shadowRayVisibility(lastHitPos, dirToCamera, gMinT, disToCamera);
            if (vis > 0)
            {
                uint2 id = getLaunchIndexFromDirection(dirToCamera, launchDim, gPixelJitter);
                if (gMaxDepth > 0)
                {
                    float theta1 = saturate(abs(dot(dirToCamera, cameraN)));
                    float theta2 = saturate(abs(dot(dirToCamera, lastHitN)));
                    shadeColor *= theta1 * theta2;
                }
                else
                {
                    gOutput[id] = float4(shadeColor, 1.0f);
                }
                float weight = getMISWeight(cameraPath, lightPath, 1, i + 2, misNodes);
                bool colorsNan = any(isnan(shadeColor));
                gOutput[id] = saturate(gOutput[id] + float4(colorsNan ? float3(0, 0, 0) : shadeColor, 1.0f));
            }
        }
    }
    
    for (uint lightLength = 2; lightLength <= gMaxDepth + 1; lightLength++)
    {
        for (uint cameraLength = 2; cameraLength <= gMaxDepth + 1; cameraLength++)
        {
            shadeColor = getUnweightedContribution(cameraPath, lightPath, lightLength, cameraLength);
            float weight = getMISWeight(cameraPath, lightPath, lightLength, cameraLength, misNodes);
            shadeColor *= weight;
            bool colorsNan = any(isnan(shadeColor));
            gOutput[launchIndex] = saturate(gOutput[launchIndex] + float4(colorsNan ? float3(0, 0, 0) : shadeColor, 1.0f));
        }
    }
}
