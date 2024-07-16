
//- Built In Unity Shader for Metallic Roughness
//- ===============================================
//-
//- Import from libraries.


import lib-sampler.glsl
import lib-pbr.glsl
import lib-sparse.glsl
import lib-alpha.glsl
import lib-defines.glsl
import lib-env.glsl
import lib-random.glsl
import lib-defines.glsl
import lib-normal.glsl
import lib-utils.glsl
import lib-vectors.glsl
import lib-sss.glsl



//: param auto environment_rotation
uniform float environment_rotation;
//- Show back faces as there may be holes in front faces.
//: state cull_face off

//- Enable alpha blending
//: state blend over

//- Channels needed for metal/rough workflow are bound here.
//: param auto channel_basecolor
uniform SamplerSparse basecolor_tex;
//: param auto channel_roughness
uniform SamplerSparse roughness_tex;
//: param auto channel_metallic
uniform SamplerSparse metallic_tex;

//: param auto is_2d_view
uniform bool uniform_2d_view;


//--------------------------------------------------------------------------------------------------
vec3 pow3(vec3 a,float b)
{
    return vec3(pow(a.x,b),pow(a.y,b),pow(a.z,b));
}
//--------------------------------------------------------------------------------------------------
vec3 mix3(vec3 a,vec3 b,float c)
{
    return a*(1-c)+c*b;
}
//--------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------
// Appoximation of joint Smith term for GGX
// [Heitz 2014, "Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs"]
float Vis_SmithJointApprox(float Roughness,float NoV,float NoL)
{
    float a=Roughness*Roughness;
    float Vis_SmithV=NoL*(NoV*(1-a)+a);
    float Vis_SmithL=NoV*(NoL*(1-a)+a);
    return.5/(Vis_SmithV+Vis_SmithL+1e-4);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// [Schlick 1994, "An Inexpensive BRDF Model for Physically-Based Rendering"]
vec3 F_Schlick(vec3 SpecularColor,float VoH)
{
    float Fc=pow(1.-VoH,5.);
    // Anything less than 2% is physically impossible and is instead considered to be shadowing
    return min(50.*SpecularColor.g,1.)*Fc+(1.-Fc)*SpecularColor;
}

//--------------------------------------------------------------------------------------------------
float Vis_Schlick(float Roughness,float NoV,float NoL)
{
    float k=Roughness*Roughness*.5;
    float Vis_SchlickV=NoV*(1-k)+k;
    float Vis_SchlickL=NoL*(1-k)+k;
    return.25/(Vis_SchlickV*Vis_SchlickL);
}
//--------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------
//- Compute the outgoing light radiance to the viewer's eye
vec4 ComputeUnityBRDF(V2F inputs,vec3 diffColor,float metallic,float roughness,float occlusion,int SamplesNum,LocalVectors vectors)
{
    //- Get world space normal
    vec3 normal_vec=computeWSNormal(inputs.sparse_coord,inputs.tangent,inputs.bitangent,inputs.normal);
    
    vec3 fresnelDielect=vec3(pow(.04,1));
    vec3 specColor=mix3(fresnelDielect,diffColor,metallic);
    //--invers gamma, intented to fake unity's default gamma space lighting
    float gamma=1;
    float gammaNorm=.87604;
    diffColor=pow3(diffColor,1/gamma);
    metallic=pow(metallic,gamma);
    
    vec3 eye_vec=normalize(camera_pos-inputs.position);
    
    //- Diffuse contribution
    vec3 diffuseIrradiance=pow3(envIrradiance(normal_vec),1/gamma);
    vec3 contribE=pow(occlusion,1/gamma)*diffuseIrradiance*diffColor*(1.-metallic)*gammaNorm;
    
    float lerp2uniform=max(2*(roughness-.5),.0001);
    
    //- Create a local basis for BRDF work
    vec3 Tp=normalize(inputs.tangent-normal_vec*dot(inputs.tangent,normal_vec));// local tangent
    vec3 Bp=normalize(inputs.bitangent-normal_vec*dot(inputs.bitangent,normal_vec)-Tp*dot(inputs.bitangent,Tp));// local bitangent
    vec3 contribS=vec3(0.);
    
    vec2 Random=fibonacci2D(alg_random_seed,SamplesNum);
    
    vec3 FullIS=vec3(0.);
    
    float ndv=dot(vectors.eye,vectors.normal);
    
    for(int i=0;i<SamplesNum;++i)
    {
        vec2 Xi=fibonacci2D(i,nbSamples);
        vec3 Hn=importanceSampleGGX(
        Xi,
        vectors.tangent,
        vectors.bitangent,
        vectors.normal,
        roughness);

        vec3 Ln=-reflect(vectors.eye,Hn);
        
        float fade=horizonFading(dot(vectors.vertexNormal,Ln),horizonFade);
        
        float ndl=dot(vectors.normal,Ln);
        ndl=max(1e-8,ndl);
        float vdh=max(1e-8,dot(vectors.eye,Hn));
        float ndh=max(1e-8,dot(vectors.normal,Hn));
        float PDF=(4.*Vis_SmithJointApprox(roughness,ndv,ndl)*ndl*vdh/ndh);
        //4.0f * Vis * NdotL * VdotH / NdotH;
        
        float lodS=roughness<.01?0.:computeLOD(Ln,probabilityGGX(ndh,vdh,roughness));
        
        vec3 irr=envSampleLOD(Ln,lodS);
        //miplevel probablityggx
        
        // FullIS+=irr*F_Schlick(specColor,vdh)*clamp(visibility(ndl,ndv,roughness)*PDF*ndl,0,1);
        //FullIS+=fade*envSampleLOD(Ln,lodS)*cook_torrance_contrib(vdh,ndh,ndl,ndv,specColor,roughness);
        FullIS+=irr*clamp(F_Schlick(specColor,vdh)*Vis_Schlick(roughness,ndv,ndl)*PDF,0.,.1);
    }
    
    FullIS/=SamplesNum;
    contribS=pow3(FullIS,1/gamma);
    
    //- Emissive
    vec3 contribEm=emissive_intensity*textureSparse(emissive_tex,inputs.sparse_coord).rgb;
    
    //- Sum diffuse + spec + emissive
    return vec4(pow3(contribS+contribE+contribEm,gamma),1);
    
}

//- Shader entry point.
vec4 shade(V2F inputs)
{
    //float glossiness=getRoughness(roughness_tex,inputs.sparse_coord);
    float roughness=sRGB2linear(getRoughness(roughness_tex,inputs.sparse_coord));
    vec4 baseColor=textureSparse(basecolor_tex,inputs.sparse_coord);
    alphaKill(inputs.sparse_coord);
    vec3 metallic=textureSparse(metallic_tex,inputs.sparse_coord).xyz;
    float occlusion=getAO(inputs.sparse_coord)*getShadowFactor();
    LocalVectors vectors=computeLocalFrame(inputs);
    sssCoefficientsOutput(getSSSCoefficients(inputs.sparse_coord));
    // Feed parameters for a physically based BRDF integration
    return vec4(
    ComputeUnityBRDF(
    inputs,
    baseColor.xyz,
    metallic.x,
    roughness,
    occlusion,
    nbSamples,
    vectors).rgb,1.);
}
