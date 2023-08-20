#pragma once

#include<vector>
#include"Vertex.h"
#include"triPrism.h"
#include"hypercube.h"
#include"quad4d.h"
#include"doublequads4d.h"
#include"box4d.h"
#include"halfbox4d.h"
#include"pentachoron.h"
#include"pyramid4d.h"
#include"hexadecachoron.h"
#include"Geometry.h"

class Primitive4D
{
private:
	std::vector<Vertex4D> vertices4D;
	std::vector<GLuint> indices4D;
	glm::vec4 normal4D;
public:
	glm::vec4* vertexData4D;
	Primitive4D() {}
	virtual ~Primitive4D() {}
	//Functions
	void set(const Vertex4D* vertices4D,
		const unsigned int sizeOfVertices4D,
		const GLuint* indices4D,
		const unsigned int sizeOfIndices4D)
	{
		this->vertices4D.clear();
		this->indices4D.clear();
		for (size_t i = 0; i < sizeOfVertices4D; i++)
		{this->vertices4D.push_back(vertices4D[i]);}
		for (size_t i = 0; i < sizeOfIndices4D; i++)
		{this->indices4D.push_back(indices4D[i]);}
	}
	inline Vertex4D* getVertices4D() { return this->vertices4D.data(); }
	inline GLuint* getIndices4D() { return this->indices4D.data(); }
	inline size_t getSizeOfVertices4D() { return this->vertices4D.size(); }
	inline size_t getSizeOfIndices4D() { return this->indices4D.size(); }
	
	void transform_vertex4D(glm::vec4* point4d, const unsigned int* indices4d0, const unsigned int telax_size, Vertex4D* vertices4d, GLuint* indices4d) {
		for (GLuint j = 0;j < telax_size;j++) {
			this->normal4D = normalize(cross4d(point4d[indices4d0[4 * j]] - point4d[indices4d0[4 * j + 1]], point4d[indices4d0[4 * j]] - point4d[indices4d0[4 * j + 2]], point4d[indices4d0[4 * j]] - point4d[indices4d0[4 * j + 3]]));
			glm::vec4 point_A = point4d[indices4d0[4 * j]];
			glm::vec4 point_B = point4d[indices4d0[4 * j + 1]];
			glm::vec4 point_C = point4d[indices4d0[4 * j + 2]];
			glm::vec4 point_D = point4d[indices4d0[4 * j + 3]];
			glm::vec3 texcoord_A = texcoord_position(this->normal4D, point4d[indices4d0[4 * j]]);
			glm::vec3 texcoord_B = texcoord_position(this->normal4D, point4d[indices4d0[4 * j + 1]]);
			glm::vec3 texcoord_C = texcoord_position(this->normal4D, point4d[indices4d0[4 * j + 2]]);
			glm::vec3 texcoord_D = texcoord_position(this->normal4D, point4d[indices4d0[4 * j + 3]]);
			glm::vec4 normal_A = this->normal4D;
			glm::vec4 normal_B = this->normal4D;
			glm::vec4 normal_C = this->normal4D;
			glm::vec4 normal_D = this->normal4D;
			for (GLuint i = 0;i < 4;i++) {

				vertices4d[i + 4 * j].point_A = point_A;
				vertices4d[i + 4 * j].point_B = point_B;
				vertices4d[i + 4 * j].point_C = point_C;
				vertices4d[i + 4 * j].point_D = point_D;

				vertices4d[i + 4 * j].texcoord_A = texcoord_A;
				vertices4d[i + 4 * j].texcoord_B = texcoord_B;
				vertices4d[i + 4 * j].texcoord_C = texcoord_C;
				vertices4d[i + 4 * j].texcoord_D = texcoord_D;

				vertices4d[i + 4 * j].normal_A = normal_A;
				vertices4d[i + 4 * j].normal_B = normal_B;
				vertices4d[i + 4 * j].normal_C = normal_C;
				vertices4d[i + 4 * j].normal_D = normal_D;

				vertices4d[i + 4 * j].id = i;
			}
			indices4d[6 * j] = 4 * j;
			indices4d[6 * j + 1] = 4 * j + 1;
			indices4d[6 * j + 2] = 4 * j + 2;
			indices4d[6 * j + 3] = 4 * j;
			indices4d[6 * j + 4] = 4 * j + 2;
			indices4d[6 * j + 5] = 4 * j + 3;
		}
		this->set(vertices4d, 4 * telax_size, indices4d, 6 * telax_size);
	}

	void transform_vertex4D(glm::vec4* point4d, glm::vec4* normal4d, const unsigned int* indices4d0, const unsigned int telax_size, Vertex4D* vertices4d, GLuint* indices4d) {
		for (GLuint j = 0; j < telax_size; j++) {
			glm::vec4 point_A = point4d[indices4d0[4 * j]];
			glm::vec4 point_B = point4d[indices4d0[4 * j + 1]];
			glm::vec4 point_C = point4d[indices4d0[4 * j + 2]];
			glm::vec4 point_D = point4d[indices4d0[4 * j + 3]];
			glm::vec3 texcoord_A = texcoord_position(glm::vec4(0.f,0.f,1.f,0.f), point4d[indices4d0[4 * j]]);
			glm::vec3 texcoord_B = texcoord_position(glm::vec4(0.f, 0.f, 1.f, 0.f), point4d[indices4d0[4 * j + 1]]);
			glm::vec3 texcoord_C = texcoord_position(glm::vec4(0.f, 0.f, 1.f, 0.f), point4d[indices4d0[4 * j + 2]]);
			glm::vec3 texcoord_D = texcoord_position(glm::vec4(0.f, 0.f, 1.f, 0.f), point4d[indices4d0[4 * j + 3]]);
			glm::vec4 normal_A = normal4d[indices4d0[4 * j]];
			glm::vec4 normal_B = normal4d[indices4d0[4 * j + 1]];
			glm::vec4 normal_C = normal4d[indices4d0[4 * j + 2]];
			glm::vec4 normal_D = normal4d[indices4d0[4 * j + 3]];
			for (GLuint i = 0; i < 4; i++) {

				vertices4d[i + 4 * j].point_A = point_A;
				vertices4d[i + 4 * j].point_B = point_B;
				vertices4d[i + 4 * j].point_C = point_C;
				vertices4d[i + 4 * j].point_D = point_D;

				vertices4d[i + 4 * j].texcoord_A = texcoord_A;
				vertices4d[i + 4 * j].texcoord_B = texcoord_B;
				vertices4d[i + 4 * j].texcoord_C = texcoord_C;
				vertices4d[i + 4 * j].texcoord_D = texcoord_D;

				vertices4d[i + 4 * j].normal_A = normal_A;
				vertices4d[i + 4 * j].normal_B = normal_B;
				vertices4d[i + 4 * j].normal_C = normal_C;
				vertices4d[i + 4 * j].normal_D = normal_D;

				vertices4d[i + 4 * j].id = i;
			}
			indices4d[6 * j] = 4 * j;
			indices4d[6 * j + 1] = 4 * j + 1;
			indices4d[6 * j + 2] = 4 * j + 2;
			indices4d[6 * j + 3] = 4 * j;
			indices4d[6 * j + 4] = 4 * j + 2;
			indices4d[6 * j + 5] = 4 * j + 3;
		}
		this->set(vertices4d, 4 * telax_size, indices4d, 6 * telax_size);
	}

	void transform_terrain4D(glm::vec4* point4d, glm::vec4* normal4d, glm::ivec3* cubeIndex_XZW,const unsigned int* indices4d0,glm::vec4 offset4d, const unsigned int telax_size, Vertex4D* vertices4d, GLuint* indices4d) {
		for (GLuint j = 0;j < telax_size;j++) {
			glm::vec4 point_A = point4d[indices4d0[4 * j]];
			glm::vec4 point_B = point4d[indices4d0[4 * j + 1]];
			glm::vec4 point_C = point4d[indices4d0[4 * j + 2]];
			glm::vec4 point_D = point4d[indices4d0[4 * j + 3]];
			glm::vec3 texcoord_A = texcoord_position2(point4d[indices4d0[4 * j]]+ offset4d);
			glm::vec3 texcoord_B = texcoord_position2(point4d[indices4d0[4 * j + 1]] + offset4d);
			glm::vec3 texcoord_C = texcoord_position2(point4d[indices4d0[4 * j + 2]] + offset4d);
			glm::vec3 texcoord_D = texcoord_position2(point4d[indices4d0[4 * j + 3]] + offset4d);
			glm::vec4 normal_A = normal4d[indices4d0[4 * j]];
			glm::vec4 normal_B = normal4d[indices4d0[4 * j+1]];
			glm::vec4 normal_C = normal4d[indices4d0[4 * j+2]];
			glm::vec4 normal_D = normal4d[indices4d0[4 * j+3]];
			glm::ivec3 cubeIndex_Xzw = cubeIndex_XZW[4 * j];
			for (GLuint i = 0;i < 4;i++) {

				vertices4d[i + 4 * j].point_A = point_A;
				vertices4d[i + 4 * j].point_B = point_B;
				vertices4d[i + 4 * j].point_C = point_C;
				vertices4d[i + 4 * j].point_D = point_D;

				vertices4d[i + 4 * j].texcoord_A = texcoord_A;
				vertices4d[i + 4 * j].texcoord_B = texcoord_B;
				vertices4d[i + 4 * j].texcoord_C = texcoord_C;
				vertices4d[i + 4 * j].texcoord_D = texcoord_D;

				vertices4d[i + 4 * j].normal_A = normal_A;
				vertices4d[i + 4 * j].normal_B = normal_B;
				vertices4d[i + 4 * j].normal_C = normal_C;
				vertices4d[i + 4 * j].normal_D = normal_D;

				vertices4d[i + 4 * j].cubeIndex_XZW = cubeIndex_Xzw;

				vertices4d[i + 4 * j].id = i;
			}
			indices4d[6 * j] = 4 * j;
			indices4d[6 * j + 1] = 4 * j + 1;
			indices4d[6 * j + 2] = 4 * j + 2;
			indices4d[6 * j + 3] = 4 * j;
			indices4d[6 * j + 4] = 4 * j + 2;
			indices4d[6 * j + 5] = 4 * j + 3;
		}
		this->set(vertices4d, 4 * telax_size, indices4d, 6 * telax_size);
	}
};


class TriPrism : public Primitive4D
{
public:
	TriPrism()
		: Primitive4D()
	{
		const unsigned int telax_size = sizeof(triPrism::indices4D0) / sizeof(GLuint) / 4;
		this->transform_vertex4D(triPrism::point4D, triPrism::indices4D0, telax_size, new Vertex4D[4 * telax_size], new GLuint[6 * telax_size]);
	}
	~TriPrism() {
	}
};
#define PI glm::pi<float>()
class Hypersphere : public Primitive4D
{
	glm::vec4* point4D;
	glm::vec4* normal4D;
	GLuint* indices4D0;
public:
	Hypersphere()
		: Primitive4D()
	{
		const int length = 8;
		const int width = 16;
		this->point4D =new glm::vec4[width*length* length];
		this->normal4D = new glm::vec4[width * length * length];
		for (int i = 0; i < length; i++)
		{
			for (int j = 0; j < length; j++)
			{
				for (int k = 0; k < width; k++)
				{
					glm::vec2 xz = glm::vec2(cos(2.f*PI*(float)(k)/ (float)(width)), sin(2.f * PI * (float)(k) / (float)(width)));
					float factor = (float)length / (float)(length - 1);
					float y = sin(PI  * (((float)j*factor - (float)length / 2.f) / (float)length));
					float xzRy = cos(PI * (((float)j* factor - (float)length / 2.f) / (float)length));
					float w = sin(PI * (((float)i*factor - (float)length / 2.f) / (float)length));
					float xyzRw = cos(PI * (((float)i* factor - (float)length / 2.f) / (float)length));
					this->point4D[k+width*j+width*length*i] = glm::vec4(0.5*xyzRw*xzRy*xz.x, 0.5 * xyzRw*y, 0.5 * xyzRw*xzRy*xz.y, 0.5 * w);
				}
			}
		}
		for (int i = 0; i < width * length * length; i++)
		{
			this->normal4D[i] = normalize(this->point4D[i]);
		}
		this->indices4D0 = new GLuint[width * (length-1) * (length-1) * 4 * 6];
		int index[8];
		int coordinate[24];
		for (int i = 0; i < length - 1; i++)
		{
			for (int j = 0; j < length - 1; j++)
			{
				for (int k = 0; k < width; k++)
				{
					if (k == width - 1)
					{
						index[0] = width - 1, index[1] = 0, index[2] = width - 1+width, index[3] = 0+width,
						index[4] = width - 1+width*length, index[5] = 0 + width * length, 
						index[6] = width - 1 + width + width * length, index[7] = 0 + width + width * length;
						for (int l = 0; l < 8; l++)
						{
							index[l] += width * j + width * length * i;
						}
					}
					if (k != width - 1)
					{
						index[0] = 0, index[1] = 1, index[2] = 0 + width, index[3] = 1 + width,
						index[4] = 0 + width * length, index[5] = 1 + width * length,
						index[6] = 0 + width + width * length, index[7] = 1 + width + width * length;
						for (int l = 0; l < 8; l++)
						{
							index[l] += k+width * j + width * length * i;
						}
					}
					coordinate[0] = index[5], coordinate[1] = index[3], coordinate[2] = index[1], coordinate[3] = index[0],
					coordinate[4] = index[4], coordinate[5] = index[2], coordinate[6] = index[3], coordinate[7] = index[0],
					coordinate[8] = index[5], coordinate[9] = index[3], coordinate[10] = index[0], coordinate[11] = index[4],
					coordinate[12] = index[7], coordinate[13] = index[6], coordinate[14] = index[2], coordinate[15] = index[4],
					coordinate[16] = index[3], coordinate[17] = index[7], coordinate[18] = index[2], coordinate[19] = index[4],
					coordinate[20] = index[7], coordinate[21] = index[3], coordinate[22] = index[5], coordinate[23] = index[4];
					for (int l = 0; l < 24; l++)
					{
						this->indices4D0[24 * (k + width * j + width * (length-1) * i)+l] = coordinate[l];
					}
				}
			}
		}
		const unsigned int telax_size = 6 * width * (length - 1) * (length - 1);//sizeof(this->indices4D0) / sizeof(GLuint) / 4;
		this->transform_vertex4D(this->point4D, this->normal4D, this->indices4D0, telax_size, new Vertex4D[4 * telax_size], new GLuint[6 * telax_size]);
	}
	~Hypersphere() {
		delete[] this->point4D;
		delete[] this->normal4D;
		delete[] this->indices4D0;
	}
};

class Hypercylinder : public Primitive4D
{
	glm::vec4* point4D;
	glm::vec4* normal4D;
	GLuint* indices4D0;
public:
	Hypercylinder()
		: Primitive4D()
	{
		const int length = 8;
		const int width = 16;
		this->point4D = new glm::vec4[width*length * 2+2];
		this->normal4D = new glm::vec4[width * length * 2 + 2];
		this->indices4D0 = new GLuint[width * (length - 1) * 4 * 6 + width * (length - 1) * 4 * 2 * 2];
		for (int i = 0; i < length; i++)
		{
			for (int j = 0; j < width; j++)
			{
				glm::vec2 xz = glm::vec2(cos(2.f * PI * (float)(j) / (float)(width)), sin(2.f * PI * (float)(j) / (float)(width)));
				float factor = (float)length / (float)(length - 1);
				float w = sin(PI * (((float)i * factor - (float)length / 2.f) / (float)length));
				float xzRw = cos(PI * (((float)i * factor - (float)length / 2.f) / (float)length));
				this->point4D[j + width * i] = glm::vec4(xzRw * xz.x, xzRw * xz.y, 1.f, w)*0.5f;
				this->point4D[width * length+j + width * i] = glm::vec4(xzRw * xz.x, xzRw * xz.y, -1.f,  w)*0.5f;
			}
		}
		this->point4D[width * length * 2 + 0]= glm::vec4(0, 0, 1.f, 0) *0.5f;
		this->point4D[width * length * 2 + 1]= glm::vec4(0, 0, -1.f, 0) *0.5f;
		for (int i = 0; i < width * length * 2; i++)
		{
			this->normal4D[i] = normalize(glm::vec4(this->point4D[i].x, this->point4D[i].y,0.f, this->point4D[i].w));
		}
		this->normal4D[width * length * 2 + 0] = glm::vec4(0, 0,1.f, 0);
		this->normal4D[width * length * 2 + 1] = glm::vec4(0, 0,-1.f, 0);
		
		{
			int index[8];
			int coordinate[24];
			for (int i = 0; i < length - 1; i++)
			{
				for (int j = 0; j < width; j++)
				{
					if (j == width - 1)
					{
						index[0] = width - 1, index[1] = 0, index[2] = width - 1 + width, index[3] = 0 + width,
							index[4] = width - 1 + width * length, index[5] = 0 + width * length,
							index[6] = width - 1 + width + width * length, index[7] = 0 + width + width * length;
						for (int l = 0; l < 8; l++)
						{
							index[l] += width * i;
						}
					}
					if (j != width - 1)
					{
						index[0] = 0, index[1] = 1, index[2] = 0 + width, index[3] = 1 + width,
							index[4] = 0 + width * length, index[5] = 1 + width * length,
							index[6] = 0 + width + width * length, index[7] = 1 + width + width * length;
						for (int l = 0; l < 8; l++)
						{
							index[l] += j + width * i;
						}
					}
					coordinate[0] = index[3], coordinate[1] = index[5], coordinate[2] = index[1], coordinate[3] = index[0],
						coordinate[4] = index[2], coordinate[5] = index[4], coordinate[6] = index[3], coordinate[7] = index[0],
						coordinate[8] = index[3], coordinate[9] = index[5], coordinate[10] = index[0], coordinate[11] = index[4],
						coordinate[12] = index[6], coordinate[13] = index[7], coordinate[14] = index[2], coordinate[15] = index[4],
						coordinate[16] = index[7], coordinate[17] = index[3], coordinate[18] = index[2], coordinate[19] = index[4],
						coordinate[20] = index[3], coordinate[21] = index[7], coordinate[22] = index[5], coordinate[23] = index[4];
					for (int l = 0; l < 24; l++)
					{
						this->indices4D0[24 * (j + width * i) + l+ width * (length - 1) * 4 * 2 * 2] = coordinate[l];
					}
				}
			}
		}
		
		{
			int index[4];
			int coordinate[8];
			for (int i = 0; i < length - 1; i++)
			{
				for (int j = 0; j < width; j++)
				{
					if (j == width - 1)
					{
						index[0] = width - 1, index[1] = 0, index[2] = width - 1 + width, index[3] = 0 + width;
						for (int l = 0; l < 4; l++)
						{
							index[l] += width * i;
						}
					}
					if (j != width - 1)
					{
						index[0] = 0, index[1] = 1, index[2] = 0 + width, index[3] = 1 + width;
						for (int l = 0; l < 4; l++)
						{
							index[l] += j + width * i;
						}
					}
					coordinate[0] = index[0], coordinate[1] = index[1], coordinate[2] = index[2], coordinate[3] = width * length * 2 + 0,
					coordinate[4] = index[2], coordinate[5] = index[1], coordinate[6] = index[3], coordinate[7] = width * length * 2 + 0;
					for (int l = 0; l < 8; l++)
					{
						this->indices4D0[8*(j + width * i) + l]= coordinate[l];
					}
					coordinate[0] = index[1]+ width * length, coordinate[1] = index[0] + width * length, coordinate[2] = index[2] + width * length, coordinate[3] = width * length * 2 + 1,
					coordinate[4] = index[1] + width * length, coordinate[5] = index[2] + width * length, coordinate[6] = index[3] + width * length, coordinate[7] = width * length * 2 + 1;
					for (int l = 0; l < 8; l++)
					{
						this->indices4D0[8*width * (length - 1) + 8*(j + width * i) + l] = coordinate[l];
					}
				}
			}
		}
		

		const unsigned int telax_size = (width * (length - 1) * 4 * 6 + width * (length - 1) * 4 * 2 * 2 )/4;
		this->transform_vertex4D(this->point4D, this->normal4D,this->indices4D0, telax_size, new Vertex4D[4 * telax_size], new GLuint[6 * telax_size]);
	}
	~Hypercylinder() {
		delete[] this->point4D;
		delete[] this->normal4D;
		delete[] this->indices4D0;
	}
};

class Hypercube : public Primitive4D
{
public:
	Hypercube()
		: Primitive4D()
	{
		const unsigned int telax_size = sizeof(hypercube::indices4D0) / sizeof(GLuint) / 4;
		this->transform_vertex4D(hypercube::point4D, hypercube::indices4D0, telax_size, new Vertex4D[4 * telax_size], new GLuint[6 * telax_size]);
	}
	~Hypercube() {
	}
};

class Quad4d : public Primitive4D
{
public:
	Quad4d()
		: Primitive4D()
	{
		const unsigned int telax_size = sizeof(quad4d::indices4D0) / sizeof(GLuint) / 4;
		this->transform_vertex4D(quad4d::point4D, quad4d::indices4D0, telax_size, new Vertex4D[4 * telax_size], new GLuint[6 * telax_size]);
	}
	~Quad4d() {
	}
};

class Doublequads4d : public Primitive4D
{
public:
	Doublequads4d()
		: Primitive4D()
	{
		const unsigned int telax_size = sizeof(doublequads4d::indices4D0) / sizeof(GLuint) / 4;
		this->transform_vertex4D(doublequads4d::point4D, doublequads4d::indices4D0, telax_size, new Vertex4D[4 * telax_size], new GLuint[6 * telax_size]);
	}
	~Doublequads4d() {
	}
};

class Box4d : public Primitive4D
{
public:
	Box4d()
		: Primitive4D()
	{
		const unsigned int telax_size = sizeof(box4d::indices4D0) / sizeof(GLuint) / 4;
		this->transform_vertex4D(box4d::point4D, box4d::indices4D0, telax_size, new Vertex4D[4 * telax_size], new GLuint[6 * telax_size]);
	}
	~Box4d() {
	}
};

class Halfbox4d : public Primitive4D
{
public:
	Halfbox4d()
		: Primitive4D()
	{
		const unsigned int telax_size = sizeof(halfbox4d::indices4D0) / sizeof(GLuint) / 4;
		this->transform_vertex4D(halfbox4d::point4D, halfbox4d::indices4D0, telax_size, new Vertex4D[4 * telax_size], new GLuint[6 * telax_size]);
	}
	~Halfbox4d() {
	}
};

class Pentachoron : public Primitive4D
{

public:
	Pentachoron()
		: Primitive4D()
	{
		const unsigned int telax_size = sizeof(pentachoron::indices4D0) / sizeof(GLuint) / 4;
		this->transform_vertex4D(pentachoron::point4D, pentachoron::indices4D0, telax_size, new Vertex4D[4 * telax_size], new GLuint[6 * telax_size]);
	}
	~Pentachoron() {
	}
};

class Pyramid4D : public Primitive4D
{
public:
	Pyramid4D()
		: Primitive4D()
	{
		const unsigned int telax_size = sizeof(pyramid4d::indices4D0) / sizeof(GLuint) / 4;
		this->transform_vertex4D(pyramid4d::point4D, pyramid4d::indices4D0, telax_size, new Vertex4D[4 * telax_size], new GLuint[6 * telax_size]);
	}
	~Pyramid4D() {
	}
};

class Hexadecachoron : public Primitive4D
{
public:
	Hexadecachoron()
		: Primitive4D()
	{
		const unsigned int telax_size = sizeof(hexadecachoron::indices4D0) / sizeof(GLuint) / 4;
		this->transform_vertex4D(hexadecachoron::point4D, hexadecachoron::indices4D0, telax_size, new Vertex4D[4 * telax_size], new GLuint[6 * telax_size]);
	}
	~Hexadecachoron() {
	}
};

class Terrain4d : public Primitive4D
{
	float* terrainHeight2;
	glm::vec4* terrainNormal2;
	glm::ivec3* coordinate_XZW;
	glm::ivec3* cubeIndex_XZW;
	glm::vec4* point4D;
	glm::vec4* normal4D;
	GLuint* indices4D0;
public:
	Terrain4d(glm::ivec4 pos, const int size)
		: Primitive4D()
	{
		int half_size = int((float)size/ 2.f);
		this->terrainHeight2 = new float[(size+1)* (size + 1)* (size + 1)];
		this->terrainNormal2 = new glm::vec4[(size + 1) * (size + 1) * (size + 1)];
		for (int w = 0;w < size+1;w++) {
			for (int z = 0;z < size+1;z++) {
				for (int x = 0;x < size+1;x++) {
					this->terrainHeight2[x + (size + 1) * z + (size + 1) * (size + 1) * w] =
						terrain_Height2(glm::vec4(x + pos.x - half_size, 0, z + pos.z - half_size, w + pos.w - half_size));
					
					this->terrainNormal2[x + (size + 1) * z + (size + 1) * (size + 1) * w] =
						getTerrainNormal(glm::vec4(x + pos.x - half_size, 0, z + pos.z - half_size, w + pos.w - half_size));
				}
			}
		}

		this->coordinate_XZW = new glm::ivec3[24];
		this->point4D = new glm::vec4[4 * 6 * size * size * size];
		this->normal4D = new glm::vec4[4 * 6 * size * size * size];
		this->cubeIndex_XZW = new glm::ivec3[4 * 6 * size * size * size];
		for (size_t w = 0;w < size;w++) {
			for (size_t z = 0;z < size;z++) {
				for (size_t x = 0;x < size;x++) {					
					unsigned int cubeIndex = 4 * 6 * (size * size * w + size * z + x);
					    this->coordinate_XZW[0] = glm::ivec3(1, 0, 1), this->coordinate_XZW[3] = glm::ivec3(0, 0, 0),
						this->coordinate_XZW[2] = glm::ivec3(1, 0, 0), this->coordinate_XZW[1] = glm::ivec3(1, 1, 0),

						this->coordinate_XZW[4] = glm::ivec3(0, 0, 1), this->coordinate_XZW[7] = glm::ivec3(0, 0, 0),
						this->coordinate_XZW[6] = glm::ivec3(1, 0, 1), this->coordinate_XZW[5] = glm::ivec3(0, 1, 0),

						this->coordinate_XZW[8] = glm::ivec3(1, 0, 1), this->coordinate_XZW[11] = glm::ivec3(0, 1, 0),
						this->coordinate_XZW[10] = glm::ivec3(0, 0, 0), this->coordinate_XZW[9] = glm::ivec3(1, 1, 0),

						this->coordinate_XZW[12] = glm::ivec3(0, 1, 1), this->coordinate_XZW[15] = glm::ivec3(0, 1, 0),
						this->coordinate_XZW[14] = glm::ivec3(0, 0, 1), this->coordinate_XZW[13] = glm::ivec3(1, 1, 1),

						this->coordinate_XZW[16] = glm::ivec3(1, 1, 1), this->coordinate_XZW[19] = glm::ivec3(0, 1, 0),
						this->coordinate_XZW[18] = glm::ivec3(0, 0, 1), this->coordinate_XZW[17] = glm::ivec3(1, 0, 1),

						this->coordinate_XZW[20] = glm::ivec3(1, 0, 1), this->coordinate_XZW[23] = glm::ivec3(0, 1, 0),
						this->coordinate_XZW[22] = glm::ivec3(1, 1, 0), this->coordinate_XZW[21] = glm::ivec3(1, 1, 1);
					for (size_t i = 0;i < 24;i++) {
						int x0 = x + this->coordinate_XZW[i].x;
						int z0 = z + this->coordinate_XZW[i].y;
						int w0 = w + this->coordinate_XZW[i].z;
						this->cubeIndex_XZW[cubeIndex + i] = glm::ivec3(x - half_size + pos.x, z - half_size + pos.z, w - half_size + pos.w);
						this->point4D[cubeIndex + i] = glm::vec4(float(x0 - half_size), this->terrainHeight2[x0 + (size + 1) * z0 + (size + 1) * (size + 1) * w0], float(z0 - half_size), float(w0 - half_size));
						this->normal4D[cubeIndex + i] = this->terrainNormal2[x0 + (size + 1) * z0 + (size + 1) * (size + 1) * w0];
					}
				}
			}
		}

		this->indices4D0 = new GLuint[4 * 6 * size * size * size];

		for (int i = 0;i < 4 * 6 * size * size * size;i++)
		{
			this->indices4D0[i] = i;
		}
		
		const unsigned int telax_size = 6 * size * size * size;
		this->transform_terrain4D(this->point4D, this->normal4D, this->cubeIndex_XZW, this->indices4D0, pos, telax_size, new Vertex4D[4 * telax_size], new GLuint[6 * telax_size]);

		this->vertexData4D = new glm::vec4[4 * telax_size];
		for (size_t i = 0;i < 4 * telax_size;i++) {
			this->vertexData4D[i] = this->point4D[this->indices4D0[i]];
		}
	}
	~Terrain4d() {
		delete[] this->terrainHeight2;
		delete[] this->terrainNormal2;

		delete[] this->coordinate_XZW;
		delete[] this->cubeIndex_XZW;
		delete[] this->point4D;
		delete[] this->normal4D;
		delete[] this->indices4D0;
	}
};