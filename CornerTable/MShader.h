#pragma once


#include <glad/glad.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include<glm/glm.hpp>

class MShader
{
public:
    unsigned int ID;
    MShader() {}
    // constructor generates the shader on the fly
    // ------------------------------------------------------------------------
    MShader(const char* meshPath, const char* fragmentPath)
    {
        // 1. retrieve the vertex/fragment source code from filePath
        std::string meshCode;
        std::string fragmentCode;
        std::ifstream mShaderFile;
        std::ifstream fShaderFile;
        // ensure ifstream objects can throw exceptions:
        mShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        fShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        try
        {
            // open files
            mShaderFile.open(meshPath);
            fShaderFile.open(fragmentPath);
            std::stringstream vShaderStream, fShaderStream;
            // read file's buffer contents into streams
            vShaderStream << mShaderFile.rdbuf();
            fShaderStream << fShaderFile.rdbuf();
            // close file handlers
            mShaderFile.close();
            fShaderFile.close();
            // convert stream into string
            meshCode = vShaderStream.str();
            fragmentCode = fShaderStream.str();
        }
        catch (std::ifstream::failure& e)
        {
            std::cout << "ERROR::SHADER::FILE_NOT_SUCCESFULLY_READ: " << e.what() << std::endl;
        }
        const char* mShaderCode = meshCode.c_str();
        const char* fShaderCode = fragmentCode.c_str();
        // 2. compile shaders
        unsigned int mesh, fragment;
        // vertex shader
        mesh = glCreateShader(GL_MESH_SHADER_NV);
        glShaderSource(mesh, 1, &mShaderCode, NULL);
        glCompileShader(mesh);
        checkCompileErrors(mesh, "MESH");
        // fragment Shader
        fragment = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fragment, 1, &fShaderCode, NULL);
        glCompileShader(fragment);
        checkCompileErrors(fragment, "FRAGMENT");
        // shader Program
        ID = glCreateProgram();
        glAttachShader(ID, mesh);
        glAttachShader(ID, fragment);
        glLinkProgram(ID);
        checkCompileErrors(ID, "PROGRAM");
        // delete the shaders as they're linked into our program now and no longer necessary
        glDeleteShader(mesh);
        glDeleteShader(fragment);
    }
    // activate the shader
    // ------------------------------------------------------------------------
    void use()
    {
        glUseProgram(ID);
    }
    // utility uniform functions
    // ------------------------------------------------------------------------
    void setBool(const std::string& name, bool value) const
    {
        glUniform1i(glGetUniformLocation(ID, name.c_str()), (int)value);
    }
    // ------------------------------------------------------------------------
    void setInt(const std::string& name, int value) const
    {
        glUniform1i(glGetUniformLocation(ID, name.c_str()), value);
    }
    // ------------------------------------------------------------------------
    void setFloat(const std::string& name, float value) const
    {
        glUniform1f(glGetUniformLocation(ID, name.c_str()), value);
    }
    void setVec3(const std::string& name, const glm::vec3& value) const
    {
        glUniform3fv(glGetUniformLocation(ID, name.c_str()), 1, &value[0]);
    }
    void setVec3(const std::string& name, float x, float y, float z) const
    {
        glUniform3f(glGetUniformLocation(ID, name.c_str()), x, y, z);
    }

    void setMat4(const std::string& name, glm::mat4 matrix) const {
        unsigned int matLoc = glGetUniformLocation(ID, name.c_str());
        glUniformMatrix4fv(matLoc, 1, GL_FALSE, glm::value_ptr(matrix));
    }

    void setVec12(const std::string& name,const std::vector<float> data) {
        GLint uniformLocation = glGetUniformLocation(ID, name.c_str());
        glUniform1fv(uniformLocation, 12, data.data());
    }

    void setUintArray3(const std::string& name, const std::vector<uint32_t>& data) {
        GLint uniformLocation = glGetUniformLocation(ID, name.c_str());
        glUniform1uiv(uniformLocation, 3, reinterpret_cast<const GLuint*>(data.data()));
    }

private:
    // utility function for checking shader compilation/linking errors.
    // ------------------------------------------------------------------------
    void checkCompileErrors(unsigned int shader, std::string type)
    {
        int success;
        char infoLog[1024];
        if (type != "PROGRAM")
        {
            glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
            if (!success)
            {
                glGetShaderInfoLog(shader, 1024, NULL, infoLog);
                std::cout << "ERROR::SHADER_COMPILATION_ERROR of type: " << type << "\n" << infoLog << "\n -- --------------------------------------------------- -- " << std::endl;
            }
        }
        else
        {
            glGetProgramiv(shader, GL_LINK_STATUS, &success);
            if (!success)
            {
                glGetProgramInfoLog(shader, 1024, NULL, infoLog);
                std::cout << "ERROR::PROGRAM_LINKING_ERROR of type: " << type << "\n" << infoLog << "\n -- --------------------------------------------------- -- " << std::endl;
            }
        }
    }
};


