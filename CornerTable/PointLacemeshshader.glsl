#version 450 core
#extension GL_NV_mesh_shader : require
#define GROUP_SIZE 32

#define WIREEND 0xFF
#define EMPTYWIRE 0xFFFFFFFF
#define INGEOMASK 0xFC000000


layout(local_size_x=GROUP_SIZE) in;
layout(max_vertices=64, max_primitives=126) out;//error may occur
layout(points) out;

layout(std430,binding = 1) readonly buffer layoutDesLoc{
	uint DesLoc[];
};

layout(std430,binding = 2) readonly buffer layoutDesInfo{
	uint DesInfo[];
};

layout(std430,binding = 3) readonly buffer layoutInterCon{
	uint InterCon[];
};

layout(std430,binding = 4) readonly buffer layoutExterCon{
	uint ExterCon[];
};

layout(std430,binding = 5) readonly buffer layoutInterGeo{
	float InterGeo[];
};

layout(std430,binding = 6) readonly buffer layoutExterGeo{
	float ExterGeo[];
};

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

// mesh output
layout (location = 0) out myPerVertexData
{
  vec3 color;
} v_out[]; //[max_vertices]


uint AnaUint(uint value,uint seq){
    value = (value>>(8*(3-seq)))&0x000000FF;
    return value;
}

vec3 extractColorFromUint(uint value) {
    uint r = AnaUint(value,1);
    uint g = AnaUint(value,2);
    uint b = AnaUint(value,3);

    vec3 color = vec3(float(r), float(g), float(b)) / 255.0;
    return color;
}

uint GetInterConInfo(uint constart,uint seq){
    constart = constart + seq/4;
    seq = seq%4;
    uint value = InterCon[constart];
    return AnaUint(value,seq);
}

uint GetExterConInfo(uint constart,uint seq){
    constart = constart + seq/4;
    seq = seq%4;
    uint value = ExterCon[constart];
    return AnaUint(value,seq);
}


uint presum[GROUP_SIZE];
bool reverse[GROUP_SIZE];//true ½âÎöright false ½âÎöleft
uint exvernum;

uint FindWireId(uint idx,uint ewirenum){
	for(int i=0;i<ewirenum;i++)
		if(idx<presum[i])
			return i;
	return ewirenum;
}

void main(){


	uint mi = gl_WorkGroupID.x;


    if(mi == 1)
    {
	uint threadid = gl_LocalInvocationID.x;

	uint start = DesLoc[mi];

	uint ewirenum = AnaUint(DesInfo[start],0);
	vec3 meshletcolor = extractColorFromUint(DesInfo[start]);
	uint irrnum = AnaUint(DesInfo[start+1],0);
	uint numvertex = AnaUint(DesInfo[start+1],1);

	uint intergeolocation = DesInfo[start+2];
    uint interconlocation = DesInfo[start+3];
    uint exterstartgeolocation = DesInfo[start+4];
    uint exterstartconlocation = DesInfo[start+4+ewirenum];

	uint intergeonum = (intergeolocation & INGEOMASK)>>26;
	uint intergeorealloc = intergeolocation & 0x03FFFFFF;
    exvernum = numvertex - intergeonum;

	uint pretemp = 0;
	for(int i=0;i<ewirenum;i++){
		uint exloc = DesInfo[start+4 + i];
		uint exnum = (exloc>>26) & (0x001F);

		presum[i] = exnum+pretemp-1;
		pretemp = presum[i];
		reverse[i] = (exloc & 0x80000000)!=0;
	}


//vertex part
	for(int i = 0; i+threadid<intergeonum;i+=GROUP_SIZE){
			uint ingeostart = intergeorealloc + (i+threadid)*3;
			vec4 vergeo = vec4(InterGeo[ingeostart],InterGeo[ingeostart+1],InterGeo[ingeostart+2],1.0f);
			gl_MeshVerticesNV[i+threadid].gl_Position = projection*view*model*vergeo;
			v_out[i+threadid].color = meshletcolor;
	}

	for(int i = 0; i+threadid < numvertex - intergeonum;i+=GROUP_SIZE){
		uint wireid = FindWireId(i+threadid,ewirenum);
        uint temp = 0;
        if(wireid!=0)
            temp = presum[wireid-1];
        uint vertexid;
        if(reverse[wireid]==false)
            vertexid = i+threadid - temp;
        else
            vertexid = presum[wireid] - (i+threadid);
        uint geoloc = (DesInfo[start+4+wireid])&0x03FFFFFF;
        geoloc = geoloc + vertexid*3;
        vec4 vergeo = vec4(ExterGeo[geoloc],ExterGeo[geoloc+1],ExterGeo[geoloc+2],1.0f);
        gl_MeshVerticesNV[i+threadid+intergeonum].gl_Position = projection*view*model*vergeo;
        v_out[i+threadid+intergeonum].color = meshletcolor;
	}

	for(int i = 0; i+threadid < numvertex;i+=GROUP_SIZE){
		gl_PrimitiveIndicesNV[i+threadid] = i+threadid;
		}
    
	if(threadid==0)
		gl_PrimitiveCountNV = numvertex;
	
	}
}
