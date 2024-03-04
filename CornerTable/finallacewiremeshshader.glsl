#version 450 core
#extension GL_NV_mesh_shader : require
#define GROUP_SIZE 32

#define WIREEND 0xFF
#define EMPTYWIRE 0xFFFFFFFF
//后面可以改小
#define MEM_LOCAT 16
#define TWOPI 6.28318530718
//需要这里解析geo数据


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
mat3 rotateX(float angle) {
  float c = cos(angle);
  float s = sin(angle);
  return mat3(
    1.0, 0.0, 0.0,
    0.0, c, -s,
    0.0, s, c
  );
}

mat3 rotateY(float angle) {
  float c = cos(angle);
  float s = sin(angle);
  return mat3(
    c, 0.0, s,
    0.0, 1.0, 0.0,
    -s, 0.0, c
  );
}

mat3 rotateZ(float angle) {
  float c = cos(angle);
  float s = sin(angle);
  return mat3(
    c, -s, 0.0,
    s, c, 0.0,
    0.0, 0.0, 1.0
  );
}

//


//*************
// 8bit的connect info在uint32是按照高位到低位的方式存储的
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

shared mat3 INDeqmat;
shared vec3 INtranvec;
shared uint inx,iny,inz,inpnt;

shared mat3 DequanMatrixs[MEM_LOCAT];
shared vec3 TransVec[MEM_LOCAT];
shared uint xlength[MEM_LOCAT];
shared uint ylength[MEM_LOCAT];
shared uint zlength[MEM_LOCAT];
shared uint pntlength[MEM_LOCAT];


uint presum[MEM_LOCAT];
bool reverse[MEM_LOCAT];
uint ewirevernum[MEM_LOCAT];
uint exvernum;

//***************

void ParseTempGeo(uint parsevalue1,uint parsevalue2,uint parsevalue3,
            out mat3 dequanmatrix,out vec3 transvec,out uint xl,out uint yl,out uint zl,out uint pntl){
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
    mat3 Meulerx = rotateX(float(rotatx)/255.0*TWOPI);
    mat3 Meulery = rotateY(float(rotaty)/255.0*TWOPI);
    mat3 Meulerz = rotateZ(float(rotatz)/255.0*TWOPI);

    //get scale matrix
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

//geo数据是先填充低位
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


vec3 newParse(uint startidx,uint startloc,uint temppntlength,uint tempxnum,uint tempynum,uint tempznum){


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
    return rawdata;
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
//************




uint FindWireId(uint idx,uint ewirenum){
	for(int i=0;i<ewirenum;i++)
		if(idx<presum[i])
			return i;
	return ewirenum;
}


void main(){

    uint mid = gl_WorkGroupID.x;
    uint threadid = gl_LocalInvocationID.x;

    uint start = DesLoc[mid];
    uint ewirenum = AnaUint(DesInfo[start],0);
	vec3 meshletcolor = extractColorFromUint(DesInfo[start]);

    uint irrnum = AnaUint(DesInfo[start+1],0);
	uint numvertex = AnaUint(DesInfo[start+1],1);
    uint intergeonum = AnaUint(DesInfo[start+1],2);
    uint ingeostart = AnaUint(DesInfo[start+1],3);


    uint intergeolocation = DesInfo[start+ingeostart];
    uint interconlocation = DesInfo[start+ingeostart+1];
    uint exterstartgeolocation = DesInfo[start+ingeostart+2];
    uint extergeostartindes = start+2+ingeostart;
    uint exterstartconlocation = DesInfo[start+2+ingeostart+ewirenum];

    //可以优化的部分
    exvernum = numvertex - intergeonum;
    uint pretemp = 0;
	for(int i=0;i<ewirenum;i++){
		uint exnum = AnaUint(DesInfo[start+2+i/4],i%4);
        ewirevernum[i] = exnum&0x7F;
		presum[i] = (exnum&0x7F) +pretemp-1;
		pretemp = presum[i];
		reverse[i] = (exnum & 0x80)!=0;
	}

    //vertex part matrix
    if(threadid == 0){
        if(intergeonum!=0){
            ParseTempGeo(InterGeo[intergeolocation],InterGeo[intergeolocation+1],InterGeo[intergeolocation+2],INDeqmat,INtranvec,inx,iny,inz,inpnt);

            for(uint i = 0;i<intergeonum;++i){
                vec3 rawdata = newParse(96+i*inpnt,intergeolocation,inpnt,inx,iny,inz);
                rawdata = INDeqmat*rawdata+INtranvec;
                vec4 vergeo = vec4(rawdata,1.0);
                
                gl_MeshVerticesNV[i].gl_Position = projection*view*model*vergeo;
			    v_out[i+threadid].color = meshletcolor;
            }
        }
    }

    if(threadid==0){
        for(int i=0;i<ewirenum;++i){
            uint loc = DesInfo[start+2+ingeostart+i];
            ParseTempGeo(ExterGeo[loc],ExterGeo[loc+1],ExterGeo[loc+2],DequanMatrixs[i],TransVec[i],xlength[i],ylength[i],zlength[1],pntlength[i]);
        }

        for(int i = 0; i < numvertex - intergeonum;i++){
		    uint wireid = FindWireId(i,ewirenum);
            uint temp = 0;
            if(wireid!=0)
                temp = presum[wireid-1];
            uint vertexid;
            if(reverse[wireid]==false)
                vertexid = i - temp;
            else
                vertexid = presum[wireid] - (i);
            uint geoloc = (DesInfo[start+2+ingeostart+wireid]);




            vec4 vergeo = ParseExterPnt(wireid,96+pntlength[wireid]*vertexid,geoloc);
            gl_MeshVerticesNV[i+intergeonum].gl_Position = projection*view*model*vergeo;
            v_out[i+intergeonum].color = meshletcolor;
	    }

    }



    // con part
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


	if(threadid==0)
			gl_PrimitiveCountNV = intergeonum*2+exvernum+irrnum;
    
}