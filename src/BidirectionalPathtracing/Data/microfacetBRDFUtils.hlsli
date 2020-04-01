// The NDF for GGX, see Eqn 19 from 
//    http://blog.selfshadow.com/publications/s2012-shading-course/hoffman/s2012_pbs_physics_math_notes.pdf
//
// This function can be used for "D" in the Cook-Torrance model:  D*G*F / (4*NdotL*NdotV)
float ggxNormalDistribution(float NdotH, float roughness)
{
	float a2 = roughness * roughness;
	float d = ((NdotH * a2 - NdotH) * NdotH + 1);
	return a2 / max(0.001f, (d * d * M_PI));
}

// This from Schlick 1994, modified as per Karas in SIGGRAPH 2013 "Physically Based Shading" course
//
// This function can be used for "G" in the Cook-Torrance model:  D*G*F / (4*NdotL*NdotV)
float ggxSchlickMaskingTerm(float NdotL, float NdotV, float roughness)
{
	// Karis notes they use alpha / 2 (or roughness^2 / 2)
	float k = roughness*roughness / 2;

	// Karis also notes they can use the following equation, but only for analytical lights
	//float k = (roughness + 1)*(roughness + 1) / 8; 

	// Compute G(v) and G(l).  These equations directly from Schlick 1994
	//     (Though note, Schlick's notation is cryptic and confusing.)
	float g_v = NdotV / (NdotV*(1 - k) + k);
	float g_l = NdotL / (NdotL*(1 - k) + k);

	// Return G(v) * G(l)
	return g_v * g_l;
}

// Traditional Schlick approximation to the Fresnel term (also from Schlick 1994)
//
// This function can be used for "F" in the Cook-Torrance model:  D*G*F / (4*NdotL*NdotV)
float3 schlickFresnel(float3 f0, float u)
{
	return f0 + (float3(1.0f, 1.0f, 1.0f) - f0) * pow(1.0f - u, 5.0f);
}

// Get a GGX half vector / microfacet normal, sampled according to the distribution computed by
//     the function ggxNormalDistribution() above.  
//
// When using this function to sample, the probability density is pdf = D * NdotH / (4 * HdotV)
float3 getGGXMicrofacet(inout uint randSeed, float roughness, float3 hitNorm)
{
	// Get our uniform random numbers
	float2 randVal = float2(nextRand(randSeed), nextRand(randSeed));

	// Get an orthonormal basis from the normal
	float3 B = getPerpendicularVector(hitNorm);
	float3 T = cross(B, hitNorm);

	// GGX NDF sampling
	float a2 = roughness * roughness;
	float cosThetaH = sqrt(max(0.0f, (1.0 - randVal.x) / ((a2 - 1.0) * randVal.x + 1)));
	float sinThetaH = sqrt(max(0.0f, 1.0f - cosThetaH * cosThetaH));
	float phiH = randVal.y * M_PI * 2.0f;

	// Get our GGX NDF sample (i.e., the half vector)
	return T * (sinThetaH * cos(phiH)) + B * (sinThetaH * sin(phiH)) + hitNorm * cosThetaH;
}

float3 ggxLighting(float3 H, float3 L, float3 N, float NdotL, float NdotV, float rough, float3 spec, out float ggxProb)
{
    float NdotH = saturate(dot(N, H));
    float LdotH = saturate(dot(L, H));

    float D = ggxNormalDistribution(NdotH, rough); // The GGX normal distribution
    float G = ggxSchlickMaskingTerm(NdotL, NdotV, rough); // Use Schlick's masking term approx
    float3 F = schlickFresnel(spec, LdotH); // Use Schlick's approx to Fresnel
    ggxProb = D * NdotH / (4 * LdotH);
    return D * G * F / (4 * NdotL * NdotV); // The Cook-Torrance microfacet BRDF
}

float specularReflect(float3 N, float3 V, float etai, float etat, out float3 L)
{
    float cosi = dot(N, V);
    float ei = etai;
    float et = etat;
    bool entering = cosi > 0.0f;
    if (!entering)
    {
        float temp = et;
        et = ei;
        ei = temp;
        N = -N;
        cosi = -cosi;
    }
    float f = fresnelDielectric(cosi, ei, et);
    L = 2.0f * cosi * N - V;
    float cosr = cosi;
    return f / cosr;
}

float fresnelDielectric(float cosi, float etai, float etat)
{
    cosi = clamp(cosi, -1.0f, 1.0f);
    float sint = (etai / etat) * sqrt(max(0.0f, 1.0f - cosi * cosi));
    if (sint >= 1.0f)
        return 1.0f;
    float cost = sqrt(max(0.0f, 1.0f - sint * sint));
    cosi = abs(cosi);
    float rParl = ((etat * cosi) - (etai * cost)) /
        ((etat * cosi) + (etai * cost));
    float rPerp = ((etai * cosi) - (etat * cost)) /
        ((etai * cosi) + (etat * cost));
    return (rParl * rParl + rPerp * rPerp) / 2.0f;
}

float specularRefract(float3 N, float3 V, float etao, float etai, out float3 L)
{
    float coso = dot(N, V);
    float et = etao;
    float ei = etai;
    bool entering = coso > 0.0f;
    if (!entering)
    {
        float temp = et;
        et = ei;
        ei = temp;
        N = -N;
        coso = -coso;
    }
    float f = fresnelDielectric(coso, et, ei);
    // total reflection
    if (f == 1.0f)
    {
        return 0.0f;
    }
    /*
        * Wi = -N * cosi - WoPerpN * sini / sino =
        * -N * cosi - (sini / sino) * (Wo - dot(N, Wo) * N) =
        * -N * sqrt(1 - (etao / etai)^2 * (1 - (dot(N, Wo))^2))) +
        * etao / etai * (Wo - dot(N, Wo) * N) =
        * N * (etao / etai * dot(N, Wo) -
        * sqrt(1 - (etao / etai)^2(1 - (dot(N, Wo))^2)) -
        * etao / etai * Wo
        */
    float eta = et / ei;
    L = normalize(N * (eta * coso -
        sqrt(max(0.0f, 1.0f - eta * eta * (1.0f - coso * coso)))) -
        eta * V);
    /*
        * see Veach 97 Chapter 5: The Sources of Non-Symmetric Scattering
        * for detail derivation, solid angle would be "squeezed" from
        * smaller IOR to bigger IOR and cause radiance increase, while
        * it doesn't apply to importance
        */
    return eta * eta * (1.0f - f) / abs(dot(V, N));
}