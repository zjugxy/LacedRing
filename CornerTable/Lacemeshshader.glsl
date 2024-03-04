#version 450 core
#extension GL_NV_mesh_shader : require
#define GROUP_SIZE 32

#define WIREEND 0xFF
#define EMPTYWIRE 0xFFFFFFFF



layout(local_size_x=GROUP_SIZE) in;
layout(max_vertices=64, max_primitives=126) out;//error may occur
layout(triangles) out;

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
uint ewirevernum[GROUP_SIZE];
uint exvernum;

uint FindWireId(uint idx,uint ewirenum){
	for(int i=0;i<ewirenum;i++)
		if(idx<presum[i])
			return i;
	return ewirenum;
}

void main(){


	uint mi = gl_WorkGroupID.x;

	uint threadid = gl_LocalInvocationID.x;

	uint start = DesLoc[mi];

	uint ewirenum = AnaUint(DesInfo[start],0);
	vec3 meshletcolor = extractColorFromUint(DesInfo[start]);
	uint irrnum = AnaUint(DesInfo[start+1],0);
	uint numvertex = AnaUint(DesInfo[start+1],1);

    uint intergeonum = AnaUint(DesInfo[start+1],2);
    uint ingeostart = AnaUint(DesInfo[start+1],3);


	uint intergeolocation = DesInfo[start+ingeostart];
    uint interconlocation = DesInfo[start+ingeostart+1];
    uint exterstartgeolocation = DesInfo[start+ingeostart+2];
    uint exterstartconlocation = DesInfo[start+2+ingeostart+ewirenum];

    exvernum = numvertex - intergeonum;

	uint pretemp = 0;
	for(int i=0;i<ewirenum;i++){

		uint exnum = AnaUint(DesInfo[start+2+i/4],i%4);

        ewirevernum[i] = exnum&0x7F;
		presum[i] = (exnum&0x7F) +pretemp-1;
		pretemp = presum[i];
		reverse[i] = (exnum & 0x80)!=0;
	}

    if(mi>=0){

//vertex part
	for(int i = 0; i+threadid<intergeonum;i+=GROUP_SIZE){
			uint ingeostart = intergeolocation + (i+threadid)*3;
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
        uint geoloc = (DesInfo[start+2+ingeostart+wireid]);
        geoloc = geoloc + vertexid*3;
        vec4 vergeo = vec4(ExterGeo[geoloc],ExterGeo[geoloc+1],ExterGeo[geoloc+2],1.0f);
        gl_MeshVerticesNV[i+threadid+intergeonum].gl_Position = projection*view*model*vergeo;
        v_out[i+threadid+intergeonum].color = meshletcolor;
	}


//con
    for(int i=0;i+threadid<intergeonum*2;i+=GROUP_SIZE){
            uint id = i+threadid;
            uint vertexid = id;
            if(vertexid>=intergeonum)
                vertexid -= intergeonum;
            
            uint idx = GetInterConInfo(interconlocation,i+threadid);
            if(idx==0xFF){
                gl_PrimitiveIndicesNV[id*3] = 0xFFFFFFFF;
                gl_PrimitiveIndicesNV[id*3+1] = 0xFFFFFFFF;
                gl_PrimitiveIndicesNV[id*3+2] = 0xFFFFFFFF;
            }else{
                gl_PrimitiveIndicesNV[id*3] = vertexid;
                gl_PrimitiveIndicesNV[id*3+1] = (vertexid+1)%intergeonum;
                gl_PrimitiveIndicesNV[id*3+2] = idx;
            }
    }

    for(int i = 0; i+threadid < numvertex - intergeonum;i+=GROUP_SIZE){
        uint wireid = FindWireId(i+threadid,ewirenum);
        uint constart = DesInfo[start+2+ingeostart+ewirenum+wireid];
        uint temp = 0;
        if(wireid!=0)
            temp = presum[wireid-1];
        uint vertexid = i+threadid - temp;
        uint idx;

        if(reverse[wireid]==false)
            idx = GetExterConInfo(constart,vertexid);
        else
            idx = GetExterConInfo(constart,vertexid+ewirevernum[wireid]);

        uint triid = 2*intergeonum+i+threadid;
    
        if(idx==0xFF){
            gl_PrimitiveIndicesNV[triid*3] = 0;
            gl_PrimitiveIndicesNV[triid*3+1] = 0;
            gl_PrimitiveIndicesNV[triid*3+2] = 0;
        }else{
            gl_PrimitiveIndicesNV[triid*3] = intergeonum+vertexid+temp;
            gl_PrimitiveIndicesNV[triid*3+1] = intergeonum+(vertexid+temp+1)%exvernum;
            gl_PrimitiveIndicesNV[triid*3+2] = idx;
        }
    }
    
    
    for(int i = 0; i+threadid < irrnum;i+=GROUP_SIZE){
        uint id = 2*intergeonum + (i+threadid)*3;
        uint idx0 = GetInterConInfo(interconlocation,id);
        uint idx1 = GetInterConInfo(interconlocation,id+1);
        uint idx2 = GetInterConInfo(interconlocation,id+2);
        uint triid = 2*intergeonum + exvernum + i+threadid;
        gl_PrimitiveIndicesNV[triid*3] = idx0;
        gl_PrimitiveIndicesNV[triid*3+1] = idx1;
        gl_PrimitiveIndicesNV[triid*3+2] = idx2;
    }
    }

	if(threadid==0)
			gl_PrimitiveCountNV = intergeonum*2+exvernum+irrnum;
    
}
