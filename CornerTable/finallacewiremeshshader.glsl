#version 450 core
#extension GL_NV_mesh_shader : require
#define GROUP_SIZE 32
#define MEM_LOCAT 16

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
	uint InterGeo[];
};

layout(std430,binding = 6) readonly buffer layoutExterGeo{
	uint ExterGeo[];
};

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform float MeshGloData[12];

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


//*************
// 8bit��connect info��uint32�ǰ��ո�λ����λ�ķ�ʽ�洢��
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
//**************

shared mat3 DequanMatrixs[MEM_LOCAT];
shared vec3 TransVec[MEM_LOCAT];
shared uint xlength[MEM_LOCAT];
shared uint ylength[MEM_LOCAT];
shared uint zlength[MEM_LOCAT];
shared uint pntlength[MEM_LOCAT];


uint presum[GROUP_SIZE];
bool reverse[GROUP_SIZE];//true 解析right false 解析left
uint ewirevernum[GROUP_SIZE];
uint exvernum;


void ParseTempGeo(uint parsevalue1,uint parsevalue2,uint parsevalue3,
            inout mat3 dequanmatrix,inout vec3 transvec,out uint xl,out uint yl,out uint zl,out uint pntl){
    uint rotatx = AnaUint(parsevalue1,0);
    uint rotaty = AnaUint(parsevalue1,1);
    uint rotatz = AnaUint(parsevalue1,2);
    uint translatex = AnaUint(parsevalue1,3);

    uint translatey = AnaUint(parsevalue2,0);
    uint translatez = AnaUint(parsevalue2,1);
    uint scalex = AnaUint(parsevalue2,2);
    uint scaley = AnaUint(parsevalue2,3);

    uint scalez = AnaUint(parsevalue3,0);
    uint xnum = AnaUint(parsevalue3,1);
    uint ynum = AnaUint(parsevalue3,2);
    uint znum = AnaUint(parsevalue3,3);
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
    //�������е��������Ӱ��
    mat3 scaleMatrix = mat3(
          descalex, 0.0, 0.0,
          0.0, descaley, 0.0,
          0.0, 0.0, descalez
        );

    dequanmatrix = ((Meulerx*Meulery)*Meulerz)*scaleMatrix;
    // get transvec
    float transx = MeshGloData[0]+(float(translatex)/255.0)*MeshGloData[1];
    float transy = MeshGloData[2]+(float(translatey)/255.0)*MeshGloData[3];
    float transz = MeshGloData[4]+(float(translatez)/255.0)*MeshGloData[5];
    transvec = vec3(transx,transy,transz);
}

vec4 ParseInterPnt(uint geoidx,uint startidx,uint startloc){
    uint temppntlength = pntlength[geoidx];
    uint tempxnum = xlength[geoidx];
    uint tempynum = ylength[geoidx];
    uint tempznum = zlength[geoidx];

    uint start = startidx/32;
    uint end = (startidx+temppntlength - 1)/32;
    uint offset = startidx%32;

    uint value;
    if(start==end){
        uint mask = (1<<temppntlength)-1;
        value = (InterGeo[start+startloc]>>offset)&mask;
        }
    else{
        uint lowvalue = InterGeo[start+startloc]>>offset;
        uint highnum = temppntlength +offset - 32;
        uint highmask = (1<<highnum)-1;
        uint highvalue = InterGeo[start+1+startloc]&highmask;
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

vec4 ParseExterPnt(uint geoidx,uint startidx,uint startloc){
    uint temppntlength = pntlength[geoidx];
    uint tempxnum = xlength[geoidx];
    uint tempynum = ylength[geoidx];
    uint tempznum = zlength[geoidx];

    uint start = startidx/32;
    uint end = (startidx+temppntlength - 1)/32;
    uint offset = startidx%32;

    uint value;
    if(start==end){
        uint mask = (1<<temppntlength)-1;
        value = (ExterGeo[start+startloc]>>offset)&mask;
        }
    else{
        uint lowvalue = ExterGeo[start+startloc]>>offset;
        uint highnum = temppntlength +offset - 32;
        uint highmask = (1<<highnum)-1;
        uint highvalue = ExterGeo[start+1+startloc]&highmask;
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


//vertex part
    if(intergeonum!=0){
        if(threadid==0)
            ParseTempGeo(InterGeo[intergeolocation],InterGeo[intergeolocation+1],InterGeo[intergeolocation+2],
            DequanMatrixs[0],TransVec[0],xlength[0],ylength[0],zlength[0],pntlength[0]);
    }

    for(int i=0;i+threadid<ewirenum;i+=GROUP_SIZE){
        uint loc = DesInfo[start+ingeostart+2+i+threadid];
        ParseTempGeo(ExterGeo[loc],ExterGeo[loc+1],ExterGeo[loc+2],
        DequanMatrixs[i+threadid+1],TransVec[i+threadid+1],xlength[i+threadid+1],ylength[i+threadid+1],
        zlength[i+threadid+1],pntlength[i+threadid+1]);
    }

	for(int i = 0; i+threadid<intergeonum;i+=GROUP_SIZE){
                vec4 vergeo = ParseInterPnt(0,96+(i+threadid)*pntlength[0],intergeolocation);
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
        vec4 vergeo = ParseExterPnt(wireid+1,96+(vertexid*pntlength[wireid+1]),geoloc);
        gl_MeshVerticesNV[i+threadid+intergeonum].gl_Position = projection*view*model*vergeo;
        v_out[i+threadid+intergeonum].color =meshletcolor;
	}

    //for(int i=0;i+threadid<numvertex;i+=GROUP_SIZE)
    //    gl_PrimitiveIndicesNV[i+threadid] = i+threadid;

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
             gl_PrimitiveIndicesNV[triid*3] = intergeonum+i+threadid;
             gl_PrimitiveIndicesNV[triid*3+1] = intergeonum+(i+threadid+1)%exvernum;
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


	if(threadid==0)
			gl_PrimitiveCountNV = intergeonum*2+exvernum+irrnum;
    
}
