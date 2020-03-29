#pragma once
#include "../SharedUtils/RenderPass.h"
#include "../SharedUtils/RayLaunch.h"

class LightTracingPass : public ::RenderPass, inherit_shared_from_this<::RenderPass, LightTracingPass>
{
public:
    using SharedPtr = std::shared_ptr<LightTracingPass>;
    using SharedConstPtr = std::shared_ptr<const LightTracingPass>;

    static SharedPtr create(const std::string &outChannel) { return SharedPtr(new LightTracingPass(outChannel)); }
    virtual ~LightTracingPass() = default;

protected:
	LightTracingPass(const std::string &outChannel) : mOutputTextureName(outChannel),
		::RenderPass("Light Tracing", "Light Tracing Options") {}

    // Implementation of RenderPass interface
    bool initialize(RenderContext* pRenderContext, ResourceManager::SharedPtr pResManager) override;
    void initScene(RenderContext* pRenderContext, Scene::SharedPtr pScene) override;
    void execute(RenderContext* pRenderContext) override;
	void renderGui(Gui* pGui) override;

	// Override some functions that provide information to the RenderPipeline class
	bool requiresScene() override { return true; }
	bool usesRayTracing() override { return true; }

    // Rendering state
	RayLaunch::SharedPtr    mpRays;                       ///< Our wrapper around a DX Raytracing pass
    RtScene::SharedPtr      mpScene;                      ///< Our scene file (passed in from app)  

	// Recursive ray tracing can be slow.  Add a toggle to disable, to allow you to manipulate the scene
	bool                    mDoIndirectGI = true;
	bool                    mDoDirectGI = true;

	int32_t                 mUserSpecifiedRayDepth = 1;   ///<  What is the current maximum ray depth
	const int32_t           mMaxPossibleRayDepth = 3;     ///<  The largest ray depth we support (without recompile)

	// What texture should was ask the resource manager to store our result in?
	std::string             mOutputTextureName;
    
	// Various internal parameters
	uint32_t                mFrameCount = 0x1337u;        ///< A frame counter to vary random numbers over time
};