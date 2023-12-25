#pragma once
#include<vector>
#include<string>
#include<iostream>
#include<fstream>
using uint = unsigned int;
class MeshFile
{
public:
	std::vector<uint> DesLoc;
	std::vector<uint> Desinfo;
	std::vector<uint> newintercon;
	std::vector<uint> newextercon;
	std::vector<float> intergeo;
	std::vector<float> extergeo;


public:
	bool ImportFile(std::string& filename) {
		size_t lastDotpos = filename.find_last_of('.');
		std::string partfilename = filename.substr(0, lastDotpos+1);

		const std::string& DesLocfilename = partfilename + std::string{ "DesLoc" };
		std::string Desinfofilename = partfilename + std::string{ "Desinfo" };
		std::string newinterconfilename = partfilename + std::string{ "newintercon" };
		std::string newexterconfilename = partfilename + std::string{ "newextercon" };
		std::string intergeofilename = partfilename + std::string{ "intergeo" };
		std::string extergeofilename = partfilename + std::string{ "extergeo" };


		std::cout <<Desinfofilename << std::endl;
		std::cout << Desinfofilename << std::endl;

		{
			std::ifstream file;
			file.open(DesLocfilename, std::ios::binary);
			if (!file.is_open()) {
				std::cerr << "�޷����ļ���" << DesLocfilename << std::endl;
				return false;
			}

			// ��ȡ������С
			size_t size;
			file.read(reinterpret_cast<char*>(&size), sizeof(size));

			// ��ȡ��������
			DesLoc.clear();
			DesLoc.resize(size);

			file.read(reinterpret_cast<char*>(DesLoc.data()), size * sizeof(uint));
			file.close();
		}
		{
			std::ifstream file;
			file.open(Desinfofilename, std::ios::binary);
			if (!file.is_open()) {
				std::cerr << "�޷����ļ���" << Desinfofilename << std::endl;
				return false;
			}

			// ��ȡ������С
			size_t size;
			file.read(reinterpret_cast<char*>(&size), sizeof(size));

			// ��ȡ��������
			Desinfo.clear();
			Desinfo.resize(size);

			file.read(reinterpret_cast<char*>(Desinfo.data()), size * sizeof(uint));
			file.close();
		}
		{
			std::ifstream file;
			file.open(newinterconfilename, std::ios::binary);
			if (!file.is_open()) {
				std::cerr << "�޷����ļ���" << newinterconfilename << std::endl;
				return false;
			}

			// ��ȡ������С
			size_t size;
			file.read(reinterpret_cast<char*>(&size), sizeof(size));

			// ��ȡ��������
			newintercon.clear();
			newintercon.resize(size);

			file.read(reinterpret_cast<char*>(newintercon.data()), size * sizeof(uint));
			file.close();
		}
		{
			std::ifstream file;
			file.open(newexterconfilename, std::ios::binary);
			if (!file.is_open()) {
				std::cerr << "�޷����ļ���" << newexterconfilename << std::endl;
				return false;
			}

			// ��ȡ������С
			size_t size;
			file.read(reinterpret_cast<char*>(&size), sizeof(size));

			// ��ȡ��������
			newextercon.clear();
			newextercon.resize(size);

			file.read(reinterpret_cast<char*>(newextercon.data()), size * sizeof(uint));
			file.close();
		}
		{
			std::ifstream file;
			file.open(extergeofilename, std::ios::binary);
			if (!file.is_open()) {
				std::cerr << "�޷����ļ���" << extergeofilename << std::endl;
				return false;
			}

			// ��ȡ������С
			size_t size;
			file.read(reinterpret_cast<char*>(&size), sizeof(size));

			// ��ȡ��������
			extergeo.clear();
			extergeo.resize(size);

			file.read(reinterpret_cast<char*>(extergeo.data()), size * sizeof(float));
			file.close();
		}
		{
			std::ifstream file;
			file.open(intergeofilename, std::ios::binary);
			if (!file.is_open()) {
				std::cerr << "�޷����ļ���" << intergeofilename << std::endl;
				return false;
			}

			// ��ȡ������С
			size_t size;
			file.read(reinterpret_cast<char*>(&size), sizeof(size));

			// ��ȡ��������
			intergeo.clear();
			intergeo.resize(size);

			file.read(reinterpret_cast<char*>(intergeo.data()), size * sizeof(float));
			file.close();
		}

		return true;
	}


};

