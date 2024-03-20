#version 450 core
#extension GL_NV_mesh_shader : require
#define GROUP_SIZE 32
#define MEM_LOCAT 16

#define WIREEND 0xFF
#define EMPTYWIRE 0xFFFFFFFF

layout(local_size_x=GROUP_SIZE) in;
layout(max_vertices=64, max_primitives=126) out;//error may occur
layout(points) out;

layout(std430,binding = 1) readonly buffer layoutDesLoc{
uint DesLoc[];
};

layout(std430,binding = 2) readonly buffer layoutDesInfo{
uint DesInfo[];
};

layout(std430,binding = 3) readonly buffer layoutInter{
uint InterWireData[];
};

layout(std430,binding = 4) readonly buffer layoutExter{
uint ExterWireData[];
};

layout(std430,binding = 5)readonly buffer layoutCorner{
uint CornerVerData[];
};

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform float MeshGloData[12];

uniform float CornerGlodata[12];
uniform uint Cornerxyznum[3];

// mesh output
layout (location = 0) out myPerVertexData
{
vec3 color;
} v_out[]; //[max_vertices]


//*****common function
mat3 matrotateX(float angle) {
    float c = cos(angle);
    float s = sin(angle);
    mat3 result;
    result[0] = vec3(1.0,0.0,0.0);
    result[1] = vec3(0.0,c,s);
    result[2] = vec3(0.0,-s,c);

    return result;
}

mat3 matrotateY(float angle) {
    float c = cos(angle);
    float s = sin(angle);

    mat3 result;
    result[0] = vec3(c,0.0,-s);
    result[1] = vec3(0.0,1.0,0.0);
    result[2] = vec3(s,0.0,c);

    return result;
}

mat3 matrotateZ(float angle) {
    float c = cos(angle);
    float s = sin(angle);

    mat3 result;

    result[0] = vec3(c,s,0.0);
    result[1] = vec3(-s,c,0.0);
    result[2] = vec3(0.0,0.0,1.0);
    return result;
}

//
uint LowAnaUint(uint value,uint seq){
  value = (value>>(8*seq))&0x000000FF;
return value;
}

uint GetInterWireByuint8(uint location,uint start){
    //start 记录了 uint8_t的个数
    uint offset = start/4;
    uint seq = start%4;
    return LowAnaUint(InterWireData[location+offset],seq);
}
uint GetExterWireByuint8(uint location,uint start){
    //start 记录了 uint8_t的个数
    uint offset = start/4;
    uint seq = start%4;
    return LowAnaUint(ExterWireData[location+offset],seq);
}
vec3 extractColorFromUint(uint value) {
	uint r = LowAnaUint(value,1);
	uint g = LowAnaUint(value,2);
	uint b = LowAnaUint(value,3);

	vec3 color = vec3(float(r), float(g), float(b)) / 255.0;
	return color;
}

shared mat3 DequanMatrixs[MEM_LOCAT];
shared vec3 TransVec[MEM_LOCAT];
shared uint xlength[MEM_LOCAT];
shared uint ylength[MEM_LOCAT];
shared uint zlength[MEM_LOCAT];
shared uint pntlength[MEM_LOCAT];

uint presum[MEM_LOCAT];
bool reverse[MEM_LOCAT];//true 解析right false 解析left
uint ewirevernum[MEM_LOCAT];
uint exvernum;


void ParseInterTempGeo(uint intervnum,uint irrnum,uint location,
            inout mat3 dequanmatrix,inout vec3 transvec,out uint xl,out uint yl,out uint zl,out uint pntl){
    //
    uint intergeostart = 2*intervnum+3*irrnum;

    uint rotatx = GetInterWireByuint8(location,intergeostart+0);
    uint rotaty = GetInterWireByuint8(location,intergeostart+1);
    uint rotatz = GetInterWireByuint8(location,intergeostart+2);

    uint translatex = GetInterWireByuint8(location,intergeostart+3);
    uint translatey = GetInterWireByuint8(location,intergeostart+4);
    uint translatez = GetInterWireByuint8(location,intergeostart+5);

    uint scalex = GetInterWireByuint8(location,intergeostart+6);
    uint scaley = GetInterWireByuint8(location,intergeostart+7);
    uint scalez = GetInterWireByuint8(location,intergeostart+8);

    uint xnum = GetInterWireByuint8(location,intergeostart+9);
    uint ynum = GetInterWireByuint8(location,intergeostart+10);
    uint znum = GetInterWireByuint8(location,intergeostart+11);
    //
    xl = xnum;
    yl = ynum;
    zl = znum;
    pntl = xnum+ynum+znum;

    // get rotation matrix
    mat3 Meulerx = matrotateX(float(rotatx)/255.0* 6.28318530718);
    mat3 Meulery = matrotateY(float(rotaty)/255.0* 6.28318530718);
    mat3 Meulerz = matrotateZ(float(rotatz)/255.0* 6.28318530718);

    //***********

    float descalex = MeshGloData[6]*pow(MeshGloData[7],float(scalex)/255.0);
    float descaley = MeshGloData[8]*pow(MeshGloData[9],float(scaley)/255.0);
    float descalez = MeshGloData[10]*pow(MeshGloData[11],float(scalez)/255.0);

    mat3 scaleMatrix = mat3(
          descalex, 0.0, 0.0,
          0.0, descaley, 0.0,
          0.0, 0.0, descalez
        );

    dequanmatrix = Meulerx*Meulery*Meulerz*scaleMatrix;
    // get transvec
    float transx = MeshGloData[0]+(float(translatex)/255.0)*MeshGloData[1];
    float transy = MeshGloData[2]+(float(translatey)/255.0)*MeshGloData[3];
    float transz = MeshGloData[4]+(float(translatez)/255.0)*MeshGloData[5];
    transvec = vec3(transx,transy,transz);
}

void ParseExterTempGeo(uint extervnum,uint location ,
    inout mat3 dequanmatrix,inout vec3 transvec,out uint xl,out uint yl,out uint zl,out uint pntl){
    //
    uint exgeostart = 2*extervnum + 8;

    uint rotatx = GetExterWireByuint8(location,exgeostart+0);
    uint rotaty = GetExterWireByuint8(location,exgeostart+1);
    uint rotatz = GetExterWireByuint8(location,exgeostart+2);

    uint translatex = GetExterWireByuint8(location,exgeostart+3);
    uint translatey = GetExterWireByuint8(location,exgeostart+4);
    uint translatez = GetExterWireByuint8(location,exgeostart+5);

    uint scalex = GetExterWireByuint8(location,exgeostart+6);
    uint scaley = GetExterWireByuint8(location,exgeostart+7);
    uint scalez = GetExterWireByuint8(location,exgeostart+8);

    uint xnum = GetExterWireByuint8(location,exgeostart+9);
    uint ynum = GetExterWireByuint8(location,exgeostart+10);
    uint znum = GetExterWireByuint8(location,exgeostart+11);
    //
    xl = xnum;
    yl = ynum;
    zl = znum;
    pntl = xnum+ynum+znum;

    // get rotation matrix
    mat3 Meulerx = matrotateX(float(rotatx)/255.0* 6.28318530718);
    mat3 Meulery = matrotateY(float(rotaty)/255.0* 6.28318530718);
    mat3 Meulerz = matrotateZ(float(rotatz)/255.0* 6.28318530718);

    //***********

    float descalex = MeshGloData[6]*pow(MeshGloData[7],float(scalex)/255.0);
    float descaley = MeshGloData[8]*pow(MeshGloData[9],float(scaley)/255.0);
    float descalez = MeshGloData[10]*pow(MeshGloData[11],float(scalez)/255.0);

    mat3 scaleMatrix = mat3(
          descalex, 0.0, 0.0,
          0.0, descaley, 0.0,
          0.0, 0.0, descalez
        );

    dequanmatrix = Meulerx*Meulery*Meulerz*scaleMatrix;
    // get transvec
    float transx = MeshGloData[0]+(float(translatex)/255.0)*MeshGloData[1];
    float transy = MeshGloData[2]+(float(translatey)/255.0)*MeshGloData[3];
    float transz = MeshGloData[4]+(float(translatez)/255.0)*MeshGloData[5];
    transvec = vec3(transx,transy,transz);
}
// 16 * intergeonum + 3*8*irnum + pntlength* (i+threadid)
vec4 ParseInterPnt(uint geoidx,uint startloc,uint bitoffset){
    uint temppntlength = pntlength[geoidx];
    uint tempxnum = xlength[geoidx];
    uint tempynum = ylength[geoidx];
    uint tempznum = zlength[geoidx];




    uint start = bitoffset/32;
    uint end = (bitoffset+temppntlength - 1)/32;
    uint offset = bitoffset%32;

    uint value;
    if(start==end){
        uint mask = (1<<temppntlength)-1;
        value = (InterWireData[start+startloc]>>offset)&mask;
        }
    else{
        uint lowvalue = InterWireData[start+startloc]>>offset;
        uint highnum = temppntlength +offset - 32;
        uint highmask = (1<<highnum)-1;
        uint highvalue = InterWireData[start+1+startloc]&highmask;
        value = lowvalue+(highvalue<<(temppntlength - highnum));
    }

    uint xmask = (1<<tempxnum)-1;
    uint ymask = (1<<tempynum)-1;
    uint zmask = (1<<tempznum)-1;

    uint xvalue = value&xmask;
    uint yvalue = (value>>tempxnum)&ymask;
    uint zvalue = (value>>(tempxnum+tempynum))&zmask;

    vec3 rawdata =  vec3(
        float(xvalue)/float(xmask)*2.0-1.0,
        float(yvalue)/float(ymask)*2.0-1.0,
        float(zvalue)/float(zmask)*2.0-1.0
    );
    rawdata = DequanMatrixs[geoidx]*rawdata+TransVec[geoidx];
    return vec4(rawdata,1.0);
}

vec4 ParseExterPnt(uint geoidx,uint startloc,uint bitoffset){
    uint temppntlength = pntlength[geoidx];
    uint tempxnum = xlength[geoidx];
    uint tempynum = ylength[geoidx];
    uint tempznum = zlength[geoidx];




    uint start = bitoffset/32;
    uint end = (bitoffset+temppntlength - 1)/32;
    uint offset = bitoffset%32;

    uint value;
    if(start==end){
        uint mask = (1<<temppntlength)-1;
        value = (ExterWireData[start+startloc]>>offset)&mask;
        }
    else{
        uint lowvalue = ExterWireData[start+startloc]>>offset;
        uint highnum = temppntlength +offset - 32;
        uint highmask = (1<<highnum)-1;
        uint highvalue = ExterWireData[start+1+startloc]&highmask;
        value = lowvalue+(highvalue<<(temppntlength - highnum));
    }

    uint xmask = (1<<tempxnum)-1;
    uint ymask = (1<<tempynum)-1;
    uint zmask = (1<<tempznum)-1;

    uint xvalue = value&xmask;
    uint yvalue = (value>>tempxnum)&ymask;
    uint zvalue = (value>>(tempxnum+tempynum))&zmask;

    vec3 rawdata =  vec3(
        float(xvalue)/float(xmask)*2.0-1.0,
        float(yvalue)/float(ymask)*2.0-1.0,
        float(zvalue)/float(zmask)*2.0-1.0
    );
    rawdata = DequanMatrixs[geoidx]*rawdata+TransVec[geoidx];
    return vec4(rawdata,1.0);
}

float getCornerDataElem(uint startloc,uint elemlength){
    uint mask = (1<<elemlength)-1;
    uint start = startloc/32;
    uint offset = startloc%32;
    uint end = (startloc+elemlength-1)/32;
    uint value;

    if(start==end){
        value = (CornerVerData[start]>>offset)&mask;
        }
    else{
        uint lowvalue = CornerVerData[start]>>offset;
        uint highnum = elemlength +offset - 32;
        uint highmask = (1<<highnum)-1;
        uint highvalue = CornerVerData[start+1]&highmask;
        value = lowvalue+(highvalue<<(elemlength - highnum));
    }
    return float(value)/float(mask)*2.0 - 1.0;
}

uint FindWireId(uint idx,uint ewirenum){
	for(int i=0;i<ewirenum;i++)
		if(idx<presum[i])
			return i;
            //最大返回ewirenum - 1
	return ewirenum;
}



void main(){


	uint mi = gl_WorkGroupID.x;
	uint threadid = gl_LocalInvocationID.x;

	uint start = DesLoc[mi];

	uint ewirenum = LowAnaUint(DesInfo[start],0);
	vec3 meshletcolor = extractColorFromUint(DesInfo[start]);
	uint irrnum = LowAnaUint(DesInfo[start+1],0);
	uint numvertex = LowAnaUint(DesInfo[start+1],1);
    uint intergeonum = LowAnaUint(DesInfo[start+1],2);
    uint ingeostart = LowAnaUint(DesInfo[start+1],3);

	uint interwirelocation = DesInfo[start+ingeostart];
    uint exterwirestartlocation = DesInfo[start+ingeostart+1];

    exvernum = numvertex - intergeonum;

	uint pretemp = 0;
	for(int i=0;i<ewirenum;i++){
		uint exnum = LowAnaUint(DesInfo[start+2+i/4],i%4);
        ewirevernum[i] = exnum&0x7F;
		presum[i] = (exnum&0x7F) +pretemp-1;
		pretemp = presum[i];
		reverse[i] = (exnum & 0x80)!=0;
	}


//vertex part

    if(intergeonum!=0){
        if(threadid==0){
            ParseInterTempGeo(intergeonum,irrnum,interwirelocation,
            DequanMatrixs[0],TransVec[0],xlength[0],ylength[0],zlength[0],pntlength[0]);
        }
    }

    for(int i = 0;i+threadid<ewirenum;i+=GROUP_SIZE){
        if(ewirevernum[i+threadid]==2)continue;
        uint exterwirelocation = DesInfo[start+ingeostart+1+i+threadid];
        ParseExterTempGeo(
            ewirevernum[i+threadid],exterwirelocation,
            DequanMatrixs[i+threadid+1],TransVec[i+threadid+1],xlength[i+threadid+1],ylength[i+threadid+1],
            zlength[i+threadid+1],pntlength[i+threadid+1]);
        
    }


	for(int i = 0; i+threadid<intergeonum;i+=GROUP_SIZE){
                uint bitoffset = 2*8*intergeonum + 3*8*irrnum + 96 + (i+threadid)*pntlength[0];
                vec4 vergeo = ParseInterPnt(0,interwirelocation,bitoffset);
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
    
        if((vertexid == 0) || (vertexid == ewirevernum[wireid] - 1))//这是corner
            continue;
    
        uint bitoffset = 64+2*8*ewirevernum[wireid]+96+(vertexid-1)*pntlength[wireid+1];
        vec4 vergeo = ParseExterPnt(wireid+1,DesInfo[start+ingeostart+1+wireid],bitoffset);

    
        gl_MeshVerticesNV[i+threadid+intergeonum].gl_Position = projection*view*model*vergeo;
        v_out[i+threadid+intergeonum].color =meshletcolor;
	}
    //
//corner data
    //
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
    
        if((vertexid == 0) || (vertexid == ewirevernum[wireid] - 1)){
            //
            uint idx = (vertexid==0?0:1);
            uint loc = DesInfo[start+ingeostart+1+wireid];
            uint corneridx = ExterWireData[loc+idx];
            uint cpntlength = Cornerxyznum[0]+Cornerxyznum[1]+Cornerxyznum[2];
            uint startcornerloc = corneridx*cpntlength;
        
            float x = getCornerDataElem(startcornerloc,Cornerxyznum[0]);
            float y = getCornerDataElem(startcornerloc+Cornerxyznum[0],Cornerxyznum[1]);
            float z = getCornerDataElem(startcornerloc+Cornerxyznum[0]+Cornerxyznum[1],
            Cornerxyznum[2]);
    
            mat3 demat = mat3(
            CornerGlodata[0],CornerGlodata[1],CornerGlodata[2],
            CornerGlodata[3],CornerGlodata[4],CornerGlodata[5],
            CornerGlodata[6],CornerGlodata[7],CornerGlodata[8]);
            vec3 newcentroid = vec3(CornerGlodata[9],CornerGlodata[10],CornerGlodata[11]);
    
            vec4 vergeo = vec4(demat*vec3(x,y,z)+newcentroid,1.0);
            gl_MeshVerticesNV[i+threadid+intergeonum].gl_Position = projection*view*model*vergeo;
            v_out[i+threadid+intergeonum].color =meshletcolor;
        }
	}


//con data










    for(int i=0;i+threadid<numvertex;i+=GROUP_SIZE)
        gl_PrimitiveIndicesNV[i+threadid] = i+threadid;

//con
    	if(threadid==0)
			gl_PrimitiveCountNV = numvertex;
}

