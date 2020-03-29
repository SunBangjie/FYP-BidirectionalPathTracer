#include "LightTracing.h"

// Some global vars, used to simplify changing shader location & entry points
namespace {
	// Where is our shaders located?
	const char* kFileRayTrace = "Tutorial14\\lightTracing.rt.hlsl";

	// What are the entry points in that shader for various ray tracing shaders?
	const char* kEntryPointRayGen        = "SimpleDiffuseGIRayGen";

	const char* kEntryPointMiss0         = "ShadowMiss";
	const char* kEntryShadowAnyHit       = "ShadowAnyHit";
	const char* kEntryShadowClosestHit   = "ShadowClosestHit";

	const char* kEntryPointMiss1         = "IndirectMiss";
	const char* kEntryIndirectAnyHit     = "IndirectAnyHit";
	const char* kEntryIndirectClosestHit = "IndirectClosestHit";

	const float kMSAA[8][2] = { { 1,-3 },{ -1,3 },{ 5,1 },{ -3,-5 },{ -5,5 },{ -7,-1 },{ 3,7 },{ 7,-7 } };
};

bool LightTracingPass::initialize(RenderContext* pRenderContext, ResourceManager::SharedPtr pResManager)
{
	// Stash a copy of our resource manager so we can get rendering resources
	mpResManager = pResManager;
	mpResManager->requestTextureResources({ "WorldPosition", "WorldNormal", "MaterialDiffuse", "MaterialSpecRough", "MaterialExtraParams", "Emissive" });
	mpResManager->requestTextureResource(mOutputTextureName);
	mpResManager->requestTextureResource(ResourceManager::kEnvironmentMap);

	// Set the default scene to load
	mpResManager->setDefaultSceneName("Data/pink_room/pink_room.fscene");

	// Create our wrapper around a ray tracing pass.  Tell it where our ray generation shader and ray-specific shaders are
	mpRays = RayLaunch::create(kFileRayTrace, kEntryPointRayGen);

	// Add ray type #0 (shadow rays)
	mpRays->addMissShader(kFileRayTrace, kEntryPointMiss0);
	mpRays->addHitShader(kFileRayTrace, kEntryShadowClosestHit, kEntryShadowAnyHit);

	// Add ray type #1 (GI rays)
	mpRays->addMissShader(kFileRayTrace, kEntryPointMiss1);
	mpRays->addHitShader(kFileRayTrace, kEntryIndirectClosestHit, kEntryIndirectAnyHit);

	// Now that we've passed all our shaders in, compile and (if available) setup the scene
	mpRays->compileRayProgram(120u);
	mpRays->setMaxRecursionDepth(uint32_t(mMaxPossibleRayDepth));
	if (mpScene) mpRays->setScene(mpScene);
    return true;
}

void LightTracingPass::initScene(RenderContext* pRenderContext, Scene::SharedPtr pScene)
{
	// Stash a copy of the scene and pass it to our ray tracer (if initialized)
    mpScene = std::dynamic_pointer_cast<RtScene>(pScene);
	if (mpRays) mpRays->setScene(mpScene);
}

void LightTracingPass::renderGui(Gui* pGui)
{
	int dirty = 0;
	dirty |= (int)pGui->addIntVar("Max RayDepth", mUserSpecifiedRayDepth, 0, mMaxPossibleRayDepth);
	dirty |= (int)pGui->addCheckBox(mDoDirectGI ? "Compute direct illumination" : "Skipping direct illumination",
									mDoDirectGI);
	dirty |= (int)pGui->addCheckBox(mDoIndirectGI ? "Shooting global illumination rays" : "Skipping global illumination", 
		                            mDoIndirectGI);
	if (dirty) setRefreshFlag();
}


void LightTracingPass::execute(RenderContext* pRenderContext)
{
	// Get the output buffer we're writing into
	Texture::SharedPtr pDstTex = mpResManager->getClearedTexture(mOutputTextureName, vec4(0.0f, 0.0f, 0.0f, 0.0f));

	// Do we have all the resources we need to render?  If not, return
	if (!pDstTex || !mpRays || !mpRays->readyToRender()) return;

	// Set our variables into the global HLSL namespace
	auto globalVars = mpRays->getGlobalVars();
	globalVars["GlobalCB"]["gMinT"]         = mpResManager->getMinTDist();
	globalVars["GlobalCB"]["gFrameCount"]   = mFrameCount++;
	globalVars["GlobalCB"]["gDoIndirectGI"] = mDoIndirectGI;
	globalVars["GlobalCB"]["gDirectGI"]     = mDoDirectGI;
	globalVars["GlobalCB"]["gMaxDepth"]     = mUserSpecifiedRayDepth;
    globalVars["GlobalCB"]["gEmitMult"]     = 1.0f;
	globalVars["gPos"] = mpResManager->getTexture("WorldPosition");
	globalVars["gNorm"] = mpResManager->getTexture("WorldNormal");
	globalVars["gDiffuseMatl"] = mpResManager->getTexture("MaterialDiffuse");
	globalVars["gSpecMatl"] = mpResManager->getTexture("MaterialSpecRough");
	globalVars["gExtraMatl"] = mpResManager->getTexture("MaterialExtraParams");
	globalVars["gEmissive"] = mpResManager->getTexture("Emissive");
	globalVars["gOutput"]      = pDstTex;
	globalVars["gEnvMap"] = mpResManager->getTexture(ResourceManager::kEnvironmentMap);

	// Determine our offset in the pixel
	float xOff = kMSAA[mFrameCount % 8][0] * 0.0625f;
	float yOff = kMSAA[mFrameCount % 8][1] * 0.0625f;

	// Set our shader and the scene camera to use the computed jitter
	globalVars["GlobalCB"]["gPixelJitter"] = vec2(xOff + 0.5f, yOff + 0.5f);
	mpScene->getActiveCamera()->setJitter(xOff / float(pDstTex->getWidth()), yOff / float(pDstTex->getHeight()));

	// Shoot our rays and shade our primary hit points
	mpRays->execute( pRenderContext, mpResManager->getScreenSize() );

}
