#ifndef SHADER_H
#define SHADER_H

#include <glad/glad.h>
#include <glm/glm.hpp>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

class Shader
{
public:
    unsigned int ID;
    // constructor generates the shader on the fly
    // ------------------------------------------------------------------------
    Shader() {};
    Shader(const char* vertexPath, const char* fragmentPath, const char* geometryPath = nullptr)
    {
        // 1. retrieve the vertex/fragment source code from filePath
        std::string vertexCode;
        std::string fragmentCode;
        std::string geometryCode;
        std::ifstream vShaderFile;
        std::ifstream fShaderFile;
        std::ifstream gShaderFile;
        // ensure ifstream objects can throw exceptions:
        vShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        fShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        gShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        try
        {
            // open files
            vShaderFile.open(vertexPath);
            fShaderFile.open(fragmentPath);
            std::stringstream vShaderStream, fShaderStream;
            // read file's buffer contents into streams
            vShaderStream << vShaderFile.rdbuf();
            fShaderStream << fShaderFile.rdbuf();
            // close file handlers
            vShaderFile.close();
            fShaderFile.close();
            // convert stream into string
            vertexCode = vShaderStream.str();
            fragmentCode = fShaderStream.str();
            // if geometry shader path is present, also load a geometry shader
            if (geometryPath != nullptr)
            {
                gShaderFile.open(geometryPath);
                std::stringstream gShaderStream;
                gShaderStream << gShaderFile.rdbuf();
                gShaderFile.close();
                geometryCode = gShaderStream.str();
            }
        }
        catch (std::ifstream::failure& e)
        {
            std::cout << "ERROR::SHADER::FILE_NOT_SUCCESFULLY_READ: " << e.what() << std::endl;
        }
        const char* vShaderCode = vertexCode.c_str();
        const char* fShaderCode = fragmentCode.c_str();
        // 2. compile shaders
        unsigned int vertex, fragment;
        // vertex shader
        vertex = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vertex, 1, &vShaderCode, NULL);
        glCompileShader(vertex);
        checkCompileErrors(vertex, "VERTEX");
        // fragment Shader
        fragment = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fragment, 1, &fShaderCode, NULL);
        glCompileShader(fragment);
        checkCompileErrors(fragment, "FRAGMENT");
        // if geometry shader is given, compile geometry shader
        unsigned int geometry;
        if (geometryPath != nullptr)
        {
            const char* gShaderCode = geometryCode.c_str();
            geometry = glCreateShader(GL_GEOMETRY_SHADER);
            glShaderSource(geometry, 1, &gShaderCode, NULL);
            glCompileShader(geometry);
            checkCompileErrors(geometry, "GEOMETRY");
        }
        // shader Program
        ID = glCreateProgram();
        glAttachShader(ID, vertex);
        glAttachShader(ID, fragment);
        if (geometryPath != nullptr)
            glAttachShader(ID, geometry);
        glLinkProgram(ID);
        checkCompileErrors(ID, "PROGRAM");
        // delete the shaders as they're linked into our program now and no longer necessery
        glDeleteShader(vertex);
        glDeleteShader(fragment);
        if (geometryPath != nullptr)
            glDeleteShader(geometry);

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
    // ------------------------------------------------------------------------
    void setVec2(const std::string& name, const glm::vec2& value) const
    {
        glUniform2fv(glGetUniformLocation(ID, name.c_str()), 1, &value[0]);
    }
    void setVec2(const std::string& name, float x, float y) const
    {
        glUniform2f(glGetUniformLocation(ID, name.c_str()), x, y);
    }
    // ------------------------------------------------------------------------
    void setVec3(const std::string& name, const glm::vec3& value) const
    {
        glUniform3fv(glGetUniformLocation(ID, name.c_str()), 1, &value[0]);
    }
    void setVec3(const std::string& name, float x, float y, float z) const
    {
        glUniform3f(glGetUniformLocation(ID, name.c_str()), x, y, z);
    }
    // ------------------------------------------------------------------------
    void setVec4(const std::string& name, const glm::vec4& value) const
    {
        glUniform4fv(glGetUniformLocation(ID, name.c_str()), 1, &value[0]);
    }
    void setVec4(const std::string& name, float x, float y, float z, float w)
    {
        glUniform4f(glGetUniformLocation(ID, name.c_str()), x, y, z, w);
    }
    // ------------------------------------------------------------------------
    void setMat2(const std::string& name, const glm::mat2& mat) const
    {
        glUniformMatrix2fv(glGetUniformLocation(ID, name.c_str()), 1, GL_FALSE, &mat[0][0]);
    }
    // ------------------------------------------------------------------------
    void setMat3(const std::string& name, const glm::mat3& mat) const
    {
        glUniformMatrix3fv(glGetUniformLocation(ID, name.c_str()), 1, GL_FALSE, &mat[0][0]);
    }
    // ------------------------------------------------------------------------
    void setMat4(const std::string& name, const glm::mat4& mat) const
    {
        glUniformMatrix4fv(glGetUniformLocation(ID, name.c_str()), 1, GL_FALSE, &mat[0][0]);
    }

    void setArrayInt(const std::string& name, const std::vector<int>& vec) const
    {
        GLint location = glGetUniformLocation(ID, name.c_str());
        glUniform1iv(location, vec.size(), vec.data());
    }

    void setDirlight(std::string name, const glm::vec3 dir, const glm::vec3 ambient, const glm::vec3 diffuse, const glm::vec3 specular) {
        std::string lightdir = name + ".dir";
        std::string lightambient = name + ".ambient";
        std::string lightdiffuse = name + ".diffuse";
        std::string lightspecular = name + ".specular";

        glUniform3fv(glGetUniformLocation(ID, lightdir.c_str()), 1, &dir[0]);
        glUniform3fv(glGetUniformLocation(ID, lightambient.c_str()), 1, &ambient[0]);
        glUniform3fv(glGetUniformLocation(ID, lightdiffuse.c_str()), 1, &diffuse[0]);
        glUniform3fv(glGetUniformLocation(ID, lightspecular.c_str()), 1, &specular[0]);
    }

    void setPointlight(std::string name, const glm::vec3 pos, const glm::vec3 ambient, const glm::vec3 diffuse, const glm::vec3 specular, const glm::vec3 diselem = glm::vec3(1.0, 0.9, 0.13)) {
        std::string lightpos = name + ".pos";
        std::string lightambient = name + ".ambient";
        std::string lightdiffuse = name + ".diffuse";
        std::string lightspecular = name + ".specular";
        std::string constant = name + ".constant";
        std::string linear = name + ".linear";
        std::string quadratic = name + ".quadratic";

        glUniform3fv(glGetUniformLocation(ID, lightpos.c_str()), 1, &pos[0]);
        glUniform3fv(glGetUniformLocation(ID, lightambient.c_str()), 1, &ambient[0]);
        glUniform3fv(glGetUniformLocation(ID, lightdiffuse.c_str()), 1, &diffuse[0]);
        glUniform3fv(glGetUniformLocation(ID, lightspecular.c_str()), 1, &specular[0]);
        setFloat(constant, diselem[0]);
        setFloat(linear, diselem[1]);
        setFloat(quadratic, diselem[2]);

    }

    //struct Pointlight {
    //    vec3 pos;
    //    vec3 ambient;
    //    vec3 diffuse;
    //    vec3 specular;
    //    float constant, linear, quadratic;
    //};
    //
    //struct Spotlight {
    //    vec3 pos;
    //    vec3 dir;
    //    vec3 ambient;
    //    vec3 diffuse;
    //    vec3 specular;
    //    float incutoff, outcutoff;
    //    float constant, linear, quadratic;
    //};

    void setSpotlight(std::string name, const glm::vec3 pos, const glm::vec3 dir, const float incutoff, const float outcutoff, const glm::vec3 ambient, const glm::vec3 diffuse, const glm::vec3 specular, const glm::vec3 diselem = glm::vec3(1.0, 0.09, 0.013)) {
        std::string lightpos = name + ".pos";
        std::string  lightdir = name + ".dir";
        std::string  incf = name + ".incutoff";
        std::string  outcf = name + ".outcutoff";

        std::string lightambient = name + ".ambient";
        std::string lightdiffuse = name + ".diffuse";
        std::string lightspecular = name + ".specular";
        std::string constant = name + ".constant";
        std::string linear = name + ".linear";
        std::string quadratic = name + ".quadratic";

        glUniform3fv(glGetUniformLocation(ID, lightpos.c_str()), 1, &pos[0]);
        glUniform3fv(glGetUniformLocation(ID, lightdir.c_str()), 1, &dir[0]);
        glUniform3fv(glGetUniformLocation(ID, lightambient.c_str()), 1, &ambient[0]);
        glUniform3fv(glGetUniformLocation(ID, lightdiffuse.c_str()), 1, &diffuse[0]);
        glUniform3fv(glGetUniformLocation(ID, lightspecular.c_str()), 1, &specular[0]);

        setFloat(constant, diselem[0]);
        setFloat(linear, diselem[1]);
        setFloat(quadratic, diselem[2]);
        setFloat(incf, incutoff);
        setFloat(outcf, outcutoff);
    }
private:
    // utility function for checking shader compilation/linking errors.
    // ------------------------------------------------------------------------
    void checkCompileErrors(GLuint shader, std::string type)
    {
        GLint success;
        GLchar infoLog[1024];
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
#endif