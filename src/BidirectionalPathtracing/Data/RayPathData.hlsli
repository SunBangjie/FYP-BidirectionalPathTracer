struct PathVertex
{
    float3 color;
	
    float3 posW;
    float3 N;
    float3 V;
	
    float3 dif;
    float3 spec;
    float rough;
    bool isSpecular;
    
    float pdfForward;
	
    static PathVertex init()
    {
        PathVertex v;
        v.color = float3(0, 0, 0);
        v.posW = float3(0, 0, 0);
        v.N = float3(0, 0, 0);
        v.V = float3(0, 0, 0);
        v.dif = float3(0, 0, 0);
        v.spec = float3(0, 0, 0);
        v.rough = 0.0f;
        v.isSpecular = false;
        v.pdfForward = 0.0f;
        return v;
    }
	
    static PathVertex create(float3 color, float3 posW, float3 N, float3 V, float3 dif, float3 spec, float rough, bool isSpecular, float pdfForward)
    {
        PathVertex v;
        v.color = color;
        v.posW = posW;
        v.N = N;
        v.V = V;
        v.dif = dif;
        v.spec = spec;
        v.rough = rough;
        v.isSpecular = isSpecular;
        v.pdfForward = pdfForward;
        return v;
    }
};

// The payload structure for our indirect rays
struct RayPayload
{
    float3 color;
    uint rndSeed;

    float3 posW;
    float3 N;
    float3 V;
    float3 dif;
    float3 spec;
    float rough;
    bool isSpecular;
    
    float pdfForward;
    
    float3 rayOrigin;
    float3 rayDir;
    
    bool terminated;
};

RayPayload initPayload(float3 rayOrigin, float3 rayDir, float3 color, uint randSeed)
{
    RayPayload payload;
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
    payload.isSpecular = false;
    payload.pdfForward = 0.0;
    payload.terminated = false;
    return payload;
}

void updateRayData(inout RayPayload rayData, float3 color, float3 hit, float3 L, float3 N, float3 V, float3 dif, float3 spec, float rough, bool isSpecular, float pdfF)
{
    // update brdf
    rayData.color *= color;
        
    // update next ray segment
    rayData.rayOrigin = hit;
    rayData.rayDir = L;
        
    // update geometry data
    rayData.posW = hit;
    rayData.N = N;
    rayData.V = V;
        
    // update material data
    rayData.dif = dif;
    rayData.spec = spec;
    rayData.rough = rough;
    rayData.isSpecular = isSpecular;
    
    // update pdf
    rayData.pdfForward = pdfF;
}