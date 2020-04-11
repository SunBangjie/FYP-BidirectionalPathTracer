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
    uint  gMatIndex;       // The index of material brdf to use
    float gRefractiveIndex; // The refractive index of dielectric materials
	uint  gMaxDepth;       // Maximum number of recursive bounces to allow
    float gEmitMult;       // Multiply emissive amount by this factor (set to 1, usually)
    float gClampUpper;     // Clamping upper bound
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

#include "RayPathData.hlsli"
#include "BRDFUtils.hlsli"
#include "MaterialUtils.hlsli"
#include "BDPTUtils.hlsli"
#include "standardShadowRay.hlsli"
#include "globalIlluminationRay.hlsli"

// How do we shade our g-buffer and spawn indirect and shadow rays?
[shader("raygeneration")]
void SimpleDiffuseGIRayGen()
{
    // get pixel index
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
    
    // Optimization: skip if hit background
    if (!isGeometryValid)
    {
        gOutput[launchIndex] = float4(difMatlColor.rgb, 1.0f);
        return;
    }

	// Extract and compute some material and geometric parameters
    float roughness = specMatlColor.a * specMatlColor.a;
    float3 V = normalize(gCamera.posW - worldPos.xyz);
    
    // Initialize random seed
    uint randSeed = initRand(launchIndex.x + launchIndex.y * launchDim.x, gFrameCount, 16);
	
    // Initialize camera and light paths
    PathVertex cameraPath[9];
    PathVertex lightPath[9];
    for (uint i = 0; i < 9; i++)
    {
        cameraPath[i] = PathVertex.init();
        lightPath[i] = PathVertex.init();
    }
    
     /**
	 ** Construct camera path by shooting rays from the camera
	 ** TODO: Now assume pinhole camera is used, but we can extend it to thin lens camera
	 */
    
	// The first vertex is the camera
    cameraPath[0].posW = gCamera.posW;
    cameraPath[0].N = normalize(gCamera.cameraW);
    cameraPath[0].color = float3(1.0f); // We = 1 for pinhole camera
    cameraPath[0].pdfForward = 1.0f; // the probability of shooting in the direction
    
    // We can make use of the G-buffer to update the second vertex
    float3 outDir;
    float pdf;
    bool isHitSpecular;
    float3 hitThrouhput = sampleBRDF(randSeed, worldPos.xyz, worldNorm.xyz, worldNorm.xyz, V, difMatlColor.xyz, specMatlColor.xyz, roughness, outDir, pdf, isHitSpecular);
    cameraPath[1] = PathVertex.create(hitThrouhput, worldPos.xyz, worldNorm.xyz, V, difMatlColor.xyz, specMatlColor.xyz, roughness, isHitSpecular, pdf);
    
    // Initialize ray payload
    RayPayload payload = initPayload(worldPos.xyz, outDir, hitThrouhput, randSeed);
    
    // Start tracing from camera
    for (uint depth = 1; depth < gMaxDepth && !payload.terminated; depth++)
    {
        shootRay(payload);
        cameraPath[depth + 1] = PathVertex.create(payload.color, payload.posW, payload.N, payload.V, 
                                                  payload.dif, payload.spec, payload.rough, 
                                                  payload.isSpecular, payload.pdfForward);
    }
    
    // Carry forward the random seed
    randSeed = payload.rndSeed;
    
    /**
	 ** Construct light path by shooting rays from the light source
	 ** TODO: Now assume point light source is used, but we can extend it to area and directional light
	 */
    
    // Sample light position and light direction
    float3 lightOrigin, lightDir, lightIntensity;
    sampleLight(randSeed, lightOrigin, lightDir, lightIntensity);
    bool takeContribution[9];
    for (uint i = 0; i < 9; i++)
        takeContribution[i] = true;
    
    // first vertex is the light sample
    lightPath[0].posW = lightOrigin;
    lightPath[0].color = lightIntensity;
    lightPath[0].pdfForward = 1.0f / gLightsCount;
    
    // initialize light ray payload
    RayPayload lightRayPayload = initPayload(lightOrigin, lightDir, lightIntensity, randSeed);

    // Start tracing from light sample
    for (uint depth = 0; depth < gMaxDepth && !lightRayPayload.terminated; depth++)
    {
        shootRay(lightRayPayload);
        lightPath[depth + 1] = PathVertex.create(lightRayPayload.color, lightRayPayload.posW, lightRayPayload.N, lightRayPayload.V, 
                                                 lightRayPayload.dif, lightRayPayload.spec, lightRayPayload.rough, 
                                                 lightRayPayload.isSpecular, lightRayPayload.pdfForward);
        takeContribution[depth + 1] = !lightRayPayload.terminated;
    }
    
    // Carry forward the random seed
    randSeed = lightRayPayload.rndSeed;
    
    // Add weighted contributions to the frame
    
    // Initialize shadeColor
    float3 shadeColor = float3(0, 0, 0);
    
    // Add path-tracing weighted contributions
    for (uint i = 0; i < gMaxDepth; i++)
    {
        shadeColor = cameraPath[i].color * evalDirectWrapper(cameraPath[i + 1], randSeed);
        //shadeColor = clampVec(shadeColor / (i + 2));
        bool colorsNan = any(isnan(shadeColor));
        gOutput[launchIndex] = gOutput[launchIndex] + float4(colorsNan ? float3(0, 0, 0) : shadeColor, 1.0f);
    }
    
    // add light-tracing weighted contributions
    shadeColor = float3(0, 0, 0);
    for (uint i = 0; i < gMaxDepth && takeContribution[i+1]; i++)
    {
        float3 lastHitPos = lightPath[i + 1].posW;
        float3 lastHitN = lightPath[i + 1].N;
        
        float3 cameraN = normalize(gCamera.cameraW);
        float3 dirToCamera = normalize(gCamera.posW - lastHitPos);
        float disToCamera = length(gCamera.posW - lastHitPos);

        if (dot(cameraN, dirToCamera) < 0 && takeContribution[i+1])
        {
            // test visibility towards camera
            bool vis = shadowRayVisibility(lastHitPos, dirToCamera, gMinT, disToCamera);
            if (vis)
            {
                uint2 id = getLaunchIndexFromDirection(dirToCamera, launchDim, gPixelJitter);
                if (gMaxDepth > 0)
                {
                    // calculate G
                    float theta1 = saturate(abs(dot(dirToCamera, cameraN)));
                    float theta2 = saturate(abs(dot(dirToCamera, lastHitN)));
                    float invDisToCamera = 1.0 / disToCamera;
                    float G = theta1 * theta2 * invDisToCamera * invDisToCamera;
                    
                    // color = thp * brdf * G
                    shadeColor = (lightPath[i].color * connectToCamera(lightPath[i + 1])) * G;
                    shadeColor = clampVec(shadeColor / (i + 2));
                    bool colorsNan = any(isnan(shadeColor));
                    //gOutput[id] = saturate(gOutput[id] + float4(colorsNan ? float3(0, 0, 0) : shadeColor, 1.0f));
                }
                else
                {
                    // point light source occupies at most one pixel
                    gOutput[id] = float4(shadeColor, 1.0f);
                }
            }
        }
    }

    // add connections of camera and light path vertices
    shadeColor = float3(0, 0, 0);
    for (uint totalLength = 2; totalLength <= gMaxDepth; totalLength++)
    {
        for (uint cameraLength = 1; cameraLength <= gMaxDepth - 1; cameraLength++)
        {
            uint lightLength = totalLength - cameraLength;
            float G = evalGWithoutV(cameraPath[cameraLength], lightPath[lightLength]);
            
            float3 posA = cameraPath[cameraLength].posW;
            float3 posB = lightPath[lightLength].posW;
            float lengthAB = length(posB - posA);
            float3 dirAB = (posB - posA) / lengthAB;
            bool V = shadowRayVisibility(posA, dirAB, gMinT, lengthAB);
            
            if (V)
            {
                shadeColor = getUnweightedContribution(cameraPath, lightPath, cameraLength, lightLength, G);
                //shadeColor *= getWeight(cameraPath, lightPath, cameraLength, lightLength);
                bool colorsNan = any(isnan(shadeColor));
                //gOutput[launchIndex] = saturate(gOutput[launchIndex] + float4(colorsNan ? float3(0, 0, 0) : shadeColor, 1.0f));
            }
        }
    }
}

