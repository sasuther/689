#pragma once

//#include "face.h"
#include <vector>
#include "Vector.h"
#include <glm/vec3.hpp> // glm::vec3
#include <glm/vec2.hpp> // glm::vec2

namespace vmr {

	class Face {
	public:

		void setface(glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, glm::vec2 uuvv0, glm::vec2 uuvv1, glm::vec2 uuvv2, glm::vec3 n0, glm::vec3 n1, glm::vec3 n2) {

			vertex0.set(v0.x, v0.y, v0.z);
			vertex1.set(v1.x, v1.y, v1.z);
			vertex2.set(v2.x, v2.y, v2.z);

			uv0.set(uuvv0.x, uuvv0.y, 0);
			uv1.set(uuvv1.x, uuvv1.y, 0);
			uv2.set(uuvv2.x, uuvv2.y, 0);

			norm0.set(n0.x, n0.y, n0.z);
			norm1.set(n1.x, n1.y, n1.z);
			norm2.set(n2.x, n2.y, n2.z);

		}
		const Vector& getVertex0() const { return vertex0; }

		const Vector& getVertex1() const { return vertex1; }

		const Vector& getVertex2() const { return vertex2; }

		const Vector& getNorm0() const { return norm0; }

		const Vector& getNorm1() const { return norm1; }

		const Vector& getNorm2() const { return norm2; }


		Vector getUV0() {

			return uv0;
		}

		Vector getUV1() {

			return uv1;
		}

		Vector getUV2() {

			return uv2;
		}



		Vector getFaceNorm() {

			Vector faceNorm = (vertex1 - vertex0) ^ (vertex2 - vertex0);
			faceNorm = faceNorm.unitvector();
			return faceNorm;
		}

	private:
		Vector vertex0, vertex1, vertex2, norm0, norm1, norm2, uv0, uv1, uv2;

	};



	class Polygon {

	public:

		Polygon() {}

		bool loadOBJ(
			const char * path,
			std::vector < glm::vec3 > & out_vertices,
			std::vector < glm::vec2 > & out_uvs,
			std::vector < glm::vec3 > & out_normals
		)
		{

			int faceCount = 0;
			FILE * file = fopen(path, "r");
			if (file == NULL) {
				printf("Impossible to open the file !\n");
				return false;
			}

			Face newFace;

			while (1) {

				char lineHeader[128]; //basically a made up length
									  // read the first word of the line
				int res = fscanf(file, "%s", lineHeader);
				if (res == EOF)
					break; // EOF = End Of File. Quit the loop.

						   // else : parse lineHeader

				if (strcmp(lineHeader, "v") == 0) {
					glm::vec3 vertex;
					fscanf(file, "%f %f %f\n", &vertex.x, &vertex.y, &vertex.z);
					temp_vertices.push_back(vertex);

				}
				else if (strcmp(lineHeader, "vt") == 0) {
					glm::vec2 uv;
					fscanf(file, "%f %f\n", &uv.x, &uv.y);
					temp_uvs.push_back(uv);
				}
				else if (strcmp(lineHeader, "vn") == 0) {
					glm::vec3 normal;
					fscanf(file, "%f %f %f\n", &normal.x, &normal.y, &normal.z);
					temp_normals.push_back(normal);

				}
				else if (strcmp(lineHeader, "f") == 0) {
					std::string vertex1, vertex2, vertex3;
					unsigned int vertexIndex[3], uvIndex[3], normalIndex[3];
					int matches = fscanf(file, "%d/%d/%d %d/%d/%d %d/%d/%d\n", &vertexIndex[0], &uvIndex[0], &normalIndex[0], &vertexIndex[1], &uvIndex[1], &normalIndex[1], &vertexIndex[2], &uvIndex[2], &normalIndex[2]);
					if (matches != 9) {
						printf("File can't be read by our simple parser : ( Try exporting with other options\n");
						return false;
					}

					vertexIndices.push_back(vertexIndex[0]);
					vertexIndices.push_back(vertexIndex[1]);
					vertexIndices.push_back(vertexIndex[2]);
					uvIndices.push_back(uvIndex[0]);
					uvIndices.push_back(uvIndex[1]);
					uvIndices.push_back(uvIndex[2]);
					normalIndices.push_back(normalIndex[0]);
					normalIndices.push_back(normalIndex[1]);
					normalIndices.push_back(normalIndex[2]);

					faceCount++;
				}

			}
			// For each vertex of each triangle
			for (unsigned int i = 0; i < vertexIndices.size(); i++) {
				unsigned int vertexIndex = vertexIndices[i];
				glm::vec3 vertex = temp_vertices[vertexIndex - 1];
				out_vertices.push_back(vertex);


			}

			for (unsigned int i = 0; i < uvIndices.size(); i++) {
				unsigned int uvIndex = uvIndices[i];
				glm::vec2 uv = temp_uvs[uvIndex - 1];
				out_uvs.push_back(uv);
			}

			for (unsigned int i = 0; i < uvIndices.size(); i++) {
				unsigned int normalIndex = normalIndices[i];
				glm::vec3 normal = temp_normals[normalIndex - 1];
				out_normals.push_back(normal);
			}

			int k = 0;
			for (int j = 0; j < faceCount; j++) {

				newFace.setface(out_vertices.at(k), out_vertices.at(k + 1), out_vertices.at(k + 2),
					out_uvs.at(k), out_uvs.at(k + 1), out_uvs.at(k + 2),
					out_normals.at(k), out_normals.at(k + 1), out_normals.at(k + 2));

				myFaces.push_back(newFace);

				k = k + 3;
			}
		}


		const std::vector<Face>& getFaces() const { return myFaces; }

	private:
		std::vector< unsigned int > vertexIndices, uvIndices, normalIndices;
		std::vector< glm::vec3 > temp_vertices;
		std::vector< glm::vec2 > temp_uvs;
		std::vector< glm::vec3 > temp_normals;

		std::vector<Face> myFaces;
	};
}