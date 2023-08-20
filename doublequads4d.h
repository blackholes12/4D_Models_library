#pragma once
//#include<GLFW/glfw3.h>
//#include<vec4.hpp>
namespace doublequads4d {
	static glm::vec4 point4D[] =
	{
	glm::vec4(-0.5f, -0.f, 0.5f,-0.5f),
	glm::vec4(0.5f, -0.f, 0.5f,-0.5f),
	glm::vec4(0.5f, -0.f, -0.5f,-0.5f),
	glm::vec4(-0.5f, -0.f, -0.5f,-0.5f),
	glm::vec4(-0.5f, 0.f, 0.5f,-0.5f),
	glm::vec4(0.5f, 0.f, 0.5f,-0.5f),
	glm::vec4(0.5f, 0.f, -0.5f,-0.5f),
	glm::vec4(-0.5f, 0.f, -0.5f,-0.5f),
	glm::vec4(-0.5f, -0.f, 0.5f,0.5f),
	glm::vec4(0.5f, -0.f, 0.5f,0.5f),
	glm::vec4(0.5f, -0.f, -0.5f,0.5f),
	glm::vec4(-0.5f, -0.f, -0.5f,0.5f),
	glm::vec4(-0.5f, 0.f, 0.5f,0.5f),
	glm::vec4(0.5f, 0.f, 0.5f,0.5f),
	glm::vec4(0.5f, 0.f, -0.5f,0.5f),
	glm::vec4(-0.5f, 0.f, -0.5f,0.5f)
	};
	static GLuint indices4D0[] =
	{
		3,8,10,11,
		2,3,1,10,
		0,8,1,3,
		1,8,9,10,
		1,3,8,10,

		15,4,6,7,
		14,15,13,6,
		12,4,13,15,
		13,4,5,6,
		13,15,4,6,
	};
}
