#include"MyLR.h"


MyLR::MyLR(const CornerTable& ct)
{


	m = ct.Ctable.size();
	mring = ct.ringvertex.size();
	mi = m - mring;

	//填充vertices 并且生成 ctvertexidx -> lrvertexidx的映射
	std::vector<int> ctv2lrv(m, 0);
	int tempid = 0;
	for ( tempid = 0; tempid < ct.ringvertex.size(); tempid++) {
		int ctvid = ct.ringvertex[tempid];
		ctv2lrv[ctvid] = tempid;
		vertices.push_back(ct.vertices[ctvid]);
	}
	//此时tempid = mring,处理孤立顶点
	for(int i=0;i<ct.Mvertex.size();i++)
		if (ct.Mvertex[i] == 0) {
			vertices.push_back(ct.vertices[i]);
			ctv2lrv[i] = tempid;
			tempid++;
		}
	//simple check
	if (tempid != m) {
		std::cout << "error in fillin vertices" << std::endl;
		return;
	}
	if (vertices.size() != m) {
		std::cout << "error in fillin vertices" << std::endl;
		return;
	}
	
	//处理corner 完成 ctcorner -> lrcorner的映射
	std::vector<int> ctc2lrc = std::vector<int>(ct.Vtable.size(), 0);
	//T1 T2三角形corner映射
	for (tempid = 0; tempid < ct.ringvertex.size(); tempid++) {
		//处理left triangle corner
		int leftcorner = ct.RingLeftCorners[tempid];
		int rightcorner = ct.RingRightCorners[tempid];
		int start = tempid * 8;
		//由于三角形的遍历顺序，后一个T2会被记载如ctc2lrc中
		//left tri corner 
		ctc2lrc[ct.Prev(leftcorner)] = start;
		ctc2lrc[leftcorner] = start + 1;
		ctc2lrc[ct.Next(leftcorner)] = start + 2;
		//right tri corner
		ctc2lrc[ct.Prev(rightcorner)] = start + 4;
		ctc2lrc[rightcorner] = start + 5;
		ctc2lrc[ct.Next(rightcorner)] = start + 6;
	}
	//T0三角形corner映射
	int T0cornerstart = ct.ringvertex.size() * 8;
	for (int i = 0; i < ct.T0trilist.size(); i++) {
		int triidx = ct.T0trilist[i];
		ctc2lrc[3 * triidx]   = T0cornerstart;
		ctc2lrc[3 * triidx+1] = T0cornerstart+1;
		ctc2lrc[3 * triidx+2] = T0cornerstart+2;
		T0cornerstart += 4;
	}
	// mr*2 + totrilist.size()个逻辑上的三角形
	 

	//LR中T0对应的ct表示
	for (int i = 0; i < ct.T0trilist.size(); i++) {
		int triidx = ct.T0trilist[i];
		int v0 = ct.Vtable[triidx * 3];
		int v1 = ct.Vtable[triidx * 3+1];
		int v2 = ct.Vtable[triidx * 3+2];
		T0Vtable.push_back(ctv2lrv[v0]);
		T0Vtable.push_back(ctv2lrv[v1]);
		T0Vtable.push_back(ctv2lrv[v2]);

		int corner0 = ct.Otable[triidx * 3];
		int corner1 = ct.Otable[triidx * 3+1];
		int corner2 = ct.Otable[triidx * 3+2];
		T0Otable.push_back(ctc2lrc[corner0]);
		T0Otable.push_back(ctc2lrc[corner1]);
		T0Otable.push_back(ctc2lrc[corner2]);
	}

	//处理孤立点的Ctable
	if (mring != m) {
		//遍历ct.Mvertex
		//插入的顺序是与isolated points插入顺序一致的
		for (int i = 0; i < ct.Mvertex.size(); i++)
			if (ct.Mvertex[i] == 0) {
				Ctable.push_back(ctc2lrc[ct.Ctable[i]]);
			}
	}


	//处理LR中的T1,T2
	int leftidx = 0;
	int rightidx = 0;
	for (tempid = 0; tempid < ct.ringvertex.size(); tempid++) {
		int cornerright = ct.RingRightCorners[tempid];
		int cornerleft = ct.RingLeftCorners[tempid];

		TypeTri typeright = ct.triangletypes[ct.Tri(cornerright)];
		TypeTri typeleft = ct.triangletypes[ct.Tri(cornerleft)];

		//先处理left
		if ((typeleft == TypeTri::T1)||(typeleft == TypeTri::T2)) {
			int vertexidxLR = ctv2lrv[ct.Vtable[cornerleft]];
			vertexidxLR = vertexidxLR << 1;
			vertexidxLR = vertexidxLR | SETNORMAL;//设置末尾1bit为0
			LeftTable.push_back(vertexidxLR);
		}
		else if ((typeleft == TypeTri::T1_irregular) || (typeleft == TypeTri::T2_irregular)) {
			VOleft temp;
			temp.vleft = ctv2lrv[ct.Vtable[cornerleft]];
			temp.v0opp = ctc2lrc[ct.Opp(ct.Prev(cornerleft))];
			temp.v2opp = ctc2lrc[ct.Opp(ct.Next(cornerleft))];
			int idx = leftidx; leftidx++;
			idx = idx << 1;
			idx = idx | SETIRREGULAR;
			LeftTable.push_back(idx);
			VOLlist.push_back(temp);
		}
		else {
			std::cout << "address unexpected triangle T0" << std::endl;
			return;
		}
		
		//处理right
		if ((typeright == TypeTri::T1) || (typeright == TypeTri::T2)) {
			int vertexidxLR = ctv2lrv[ct.Vtable[cornerright]];
			vertexidxLR = vertexidxLR << 1;
			vertexidxLR = vertexidxLR | SETNORMAL;//设置末尾1bit为0
			RightTable.push_back(vertexidxLR);
		}
		//看到这里
		else if ((typeright == TypeTri::T1_irregular) || (typeright == TypeTri::T2_irregular)) {
			VOright temp;
			temp.vright = ctv2lrv[ct.Vtable[cornerright]];
			temp.v4opp = ctc2lrc[ct.Opp(ct.Prev(cornerright))];
			temp.v6opp = ctc2lrc[ct.Opp(ct.Next(cornerright))];
			int idx = rightidx; rightidx++;
			idx = idx << 1;
			idx = idx | SETIRREGULAR;
			RightTable.push_back(idx);
			VORlist.push_back(temp);
		} else{
			std::cout << "address unexpected triangle T0" << std::endl;
			return;
		}
	}




	//check
	checkLR();

	GenerateEBO();

}
//not fast design for check by rendering
void MyLR::RenderByLRCHECK()
{
	//vertices不变生成一个rennderlist

}

void MyLR::Decompress()
{
}

//根据原论文先处理left
void MyLR::Addresstri(int idx,TypeTri leftt, TypeTri rightt) {

//unused
}


bool MyLR::ISirregular(int idx) {
	return idx & ISIRREGULAR;//取末尾
}

bool MyLR::IsNormaltri(int idx) {
	return ((idx & ISIRREGULAR) == 0);
}

bool MyLR::ParseLeft(int leftcode, int& leftvertexidx)
{
	if (leftcode & ISIRREGULAR == 1) {
		int index = leftcode >> 1;
		leftvertexidx = VOLlist[index].vleft;
		return 1;
	}
	else
		leftvertexidx = leftcode >> 1;

	return 0;
}

bool MyLR::ParseRight(int rightcode, int& rightvertexidx)
{
	if (rightcode & ISIRREGULAR == 1) {
		int index = rightcode >> 1;
		rightvertexidx = VORlist[index].vright;
		return 1;
	}
	else
		rightvertexidx = rightcode >> 1;
	return 0;
}

void MyLR::checkLR()
{
	std::cout << "lr whether is right needed to be checked" << std::endl;

}

int MyLR::CornertoVertexidx(int corner)
{
	if (corner >= 8 * mring) {
		int i = corner - (corner / 4) - 6 * mring;
		return T0Vtable[i];
	}
	else {
		int i = corner / 8;
		int judgeidx = corner % 8;

		if ((judgeidx==0)||(judgeidx==6))
			return i;
		else if ((judgeidx == 2) || (judgeidx == 4))
			return i + 1;
		else if (judgeidx == 1) {
				int returnid;
				ParseLeft(LeftTable[i], returnid);
				return returnid;
		}
		else if (judgeidx == 5) {
			int returnid;
			ParseRight(RightTable[i], returnid);
			return returnid;
		}
		else {
			std::cout << "error in LR c.V operator" << std::endl;
			return -1;
		}
	}

	return 0;
}

vec3 MyLR::CornertoVertexvalue(int corner)
{
	int vertexidx = CornertoVertexidx(corner);
	return vertices[vertexidx];
}

void MyLR::GenerateEBO(){
	for (int i = 0; i < mring; i++) {
		//left first and not first T2
		//有问题
		int leftv;
		int type;
		type = ParseLeft(LeftTable[i], leftv);

		//type 0 normal 1 i
		if (leftv != (i + 2)%mring) {
			//填充循序不对就不渲染我日
			EBO_check.push_back(i);
			EBO_check.push_back(leftv);
			EBO_check.push_back((i + 1) % mring);
			facetypeid.push_back(type);
		}
		int rightv;
		type = ParseRight(RightTable[i], rightv);
		if (rightv != (i + 2)%mring) {
			EBO_check.push_back(i);
			EBO_check.push_back((i + 1) % mring);
			EBO_check.push_back(rightv);
			facetypeid.push_back(type);
		}
		//right tri
	}

	for (int i = 0; i < T0Vtable.size(); i += 3) {
		EBO_check.push_back(T0Vtable[i+2]);
		EBO_check.push_back(T0Vtable[i+1]);
		EBO_check.push_back(T0Vtable[i]);
		facetypeid.push_back(2);
	}
	std::cout << "lr ebo finished" << std::endl;
}
