#include "Falcor.h"
#include "../SharedUtils/RenderingPipeline.h"
#include "../CommonPasses/LightProbeGBufferPass.h"
#include "Passes/BDPTPass.h"
#include "Passes/DenoisePass.h"

int WINAPI WinMain(_In_ HINSTANCE hInstance, _In_opt_ HINSTANCE hPrevInstance, _In_ LPSTR lpCmdLine, _In_ int nShowCmd)
{
	// Create our rendering pipeline
	RenderingPipeline *pipeline = new RenderingPipeline();

	// Add passes into our pipeline
	pipeline->setPass(0, LightProbeGBufferPass::create());
	pipeline->setPass(1, BDPTPass::create(ResourceManager::kOutputChannel));
	pipeline->setPass(2, BlockwiseMultiOrderFeatureRegression::create(ResourceManager::kOutputChannel));
   
	// Define a set of config / window parameters for our program
    SampleConfig config;
	config.windowDesc.resizableWindow = true;
	config.windowDesc.width = 1280; 
	config.windowDesc.height = 720;
	config.windowDesc.title = "Bidirectional Path Tracing (HD resolution 1280 x 720)";

	// Start our program!
	RenderingPipeline::run(pipeline, config);
}
