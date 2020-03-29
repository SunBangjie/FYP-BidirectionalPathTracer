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
	bool  gDoIndirectGI;   // A boolean determining if we should shoot indirect GI rays
	bool  gDoDirectGI;     // A boolean determining if we should compute direct lighting
	uint  gMaxDepth;       // Maximum number of recursive bounces to allow
    float gEmitMult;       // Multiply emissive amount by this factor (set to 1, usually)
}

// Input and out textures that need to be set by the C++ code (for the ray gen shader)
shared Texture2D<float4>   gPos;
shared Texture2D<float4>   gNorm;
shared Texture2D<float4>   gDiffuseMatl;
shared Texture2D<float4>   gSpecMatl;
shared Texture2D<float4>   gExtraMatl;
shared Texture2D<float4>   gEnvMap;
shared Texture2D<float4>   gEmissive;
shared RWTexture2D<float4> gOutput;

// A separate file with some simple utility functions: getPerpendicularVector(), initRand(), nextRand()
#include "pathTracingUtils.hlsli"

// Include implementations of GGX normal distribution function, Fresnel approx,
//     masking term and function to sampl NDF 
#include "microfacetBRDFUtils.hlsli"

// Our material has have both a diffuse and a specular lobe.  
//     With what probability should we sample the diffuse one?
float probabilityToSampleDiffuse(float3 difColor, float3 specColor)
{
	float lumDiffuse = max(0.01f, luminance(difColor.rgb));
	float lumSpecular = max(0.01f, luminance(specColor.rgb));
	return lumDiffuse / (lumDiffuse + lumSpecular);
}

// Include shader entries, data structures, and utility functions to spawn rays
#include "standardShadowRay.hlsli"
#include "pathTracingRay.hlsli"

struct PathVertex
{
    float3 color;
    float3 posW;
	
	static PathVertex init()
	{
        PathVertex v;
        v.color = float3(0, 0, 0);
        v.posW = float3(0, 0, 0);
        return v;
    }
	
    static PathVertex create(float3 color, float3 posW)
    {
        PathVertex v;
        v.color = color;
        v.posW = posW;
        return v;
    }
};

// How do we shade our g-buffer and spawn indirect and shadow rays?
[shader("raygeneration")]
void SimpleDiffuseGIRayGen()
{
	// Where is this ray on screen?
	uint2 launchIndex    = DispatchRaysIndex().xy;
	uint2 launchDim      = DispatchRaysDimensions().xy;

	// Load g-buffer data
	// g-buffer is used to optimize camera path construction
	float4 worldPos      = gPos[launchIndex];
	float4 worldNorm     = gNorm[launchIndex];
	float4 difMatlColor  = gDiffuseMatl[launchIndex];
	float4 specMatlColor = gSpecMatl[launchIndex];
	float4 extraData     = gExtraMatl[launchIndex];
    float4 pixelEmissive = gEmissive[launchIndex];
	
	// Does this g-buffer pixel contain a valid piece of geometry?  (0 in pos.w for invalid)
	bool isGeometryValid = (worldPos.w != 0.0f);

	// Extract and compute some material and geometric parameters
	float roughness = specMatlColor.a * specMatlColor.a;
	float3 V = normalize(gCamera.posW - worldPos.xyz);

	// Make sure our normal is pointed the right direction
	if (dot(worldNorm.xyz, V) <= 0.0f) worldNorm.xyz = -worldNorm.xyz;
	float NdotV = dot(worldNorm.xyz, V);

	// Grab our geometric normal.  Also make sure this points the right direction.
	//     This is badly hacked into our G-buffer for now.  We need this because 
	//     sometimes, when normal mapping, our randomly selected indirect ray will 
	//     be *below* the surface (due to the normal map perturbations), which will 
	//     cause light leaking.  We solve by ignoring the ray's contribution if it
	//     is below the horizon.  
	float3 noMapN = normalize(extraData.yzw);
	if (dot(noMapN, V) <= 0.0f) noMapN = -noMapN;

	// If we don't hit any geometry, our difuse material contains our background color.
	float3 shadeColor    = isGeometryValid ? float3(1,0,0) : difMatlColor.rgb;

	// Initialize our random number generator
	uint randSeed        = initRand(launchIndex.x + launchIndex.y * launchDim.x, gFrameCount, 16);

	/**
	 ** Construct camera path by shooting rays from the camera
	 ** TODO: Now assume pinhole camera is used, but we can extend it to thin lens camera
	 */
	
    PathVertex cameraPath[5];
	for (uint i = 0; i < 5; i++) cameraPath[i] = PathVertex.init();
	
	// The first vertex is the camera
    cameraPath[0] = PathVertex.create(float3(1, 1, 1), gCamera.posW);
	
	// The second vertex is from g-buffer
    cameraPath[1] = PathVertex.create(float3(1, 1, 1) /* We / P = 1 */, worldPos.xyz);
    uint cameraPathLength = 2;
	
	// Do shading, if we have geoemtry here (otherwise, output the background color)
	if (isGeometryValid)
	{
        // Add any emissive color from primary rays
        shadeColor = gEmitMult * pixelEmissive.rgb;
		
		// (Optionally) do explicit direct lighting to a random light in the scene
		// This is not needed and it is already accounted for when light path length = 1

		if (gDoDirectGI)
			shadeColor += ggxDirect(randSeed, worldPos.xyz, worldNorm.xyz, V,
				                   difMatlColor.rgb, specMatlColor.rgb, roughness);

        float3 posW;
		// (Optionally) do indirect lighting for global illumination
		if (gDoIndirectGI && (gMaxDepth > 0))
			shadeColor += ggxIndirect(randSeed, worldPos.xyz, worldNorm.xyz, noMapN,
										V, difMatlColor.rgb, specMatlColor.rgb, roughness, 0, posW);
    }
	
    bool colorsNan = any(isnan(shadeColor));
	// Store out the color of this shaded pixel
	gOutput[launchIndex] = float4(colorsNan?float3(0,0,0):shadeColor, 1.0f);
}