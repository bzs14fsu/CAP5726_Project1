// C++ include
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cfloat>
#include <Eigen/Geometry> //Had to include this to do cross products

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"
#include "utils.h"

// Shortcut to avoid Eigen:: and std:: everywhere, DO NOT USE IN .h
using namespace std;
using namespace Eigen;

void allparts()
{
	//Part 1 Declarations
	cout << "Part 1: Simple ray tracer, multiple spheres with orthographic projection\n";

	const string filename1("part1.png");
	MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
	MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask
	
	//Create x spheres
	int number;
	cout << "Enter number of spheres to make (max. 6), or enter -1\n";
	cout << "to use the same inputs as in the README: ";
	cin >> number;
	
	while(number > 6)
	{
		cout << "Max. 6: ";
		cin >> number;
	}
	
	double sph[6][4];
	double posmax = 0;
	
	//Input for spheres
	if(number != -1)
	{
		cout << "Create the spheres to draw. Enter values as x y z and radius, separated only by the enter key -\n";
		for(int i = 0; i < number; i++)
		{
			cout << "Circle " << i + 1 << ": ";
			cin >> sph[i][0] >> sph[i][1] >> sph[i][2] >> sph[i][3];
		}
	}
	else
	{
		number = 5;
		sph[0][0] = 0; sph[0][1] = 0; sph[0][2] = 0; sph[0][3] = 2;
		sph[1][0] = -3; sph[1][1] = 3; sph[1][2] = 3; sph[1][3] = 0.25;
		sph[2][0] = 2; sph[2][1] = -2; sph[2][2] = 2; sph[2][3] = 0.5;
		sph[3][0] = 4; sph[3][1] = -4; sph[3][2] = 0; sph[3][3] = 0.5;
		sph[4][0] = -2; sph[4][1] = -3.5; sph[4][2] = -3.5; sph[4][3] = 2;
	}
	
	
	//Part 2 Declarations
	//We reuse the same spheres but recolor/retexture
	//User input would be a bit tedious for this part, so I am hard-coding these values
	cout << "Part 2: Add colors and more advanced shading\n";
	const string filename2("part2.png");
	
	MatrixXd R = MatrixXd::Zero(800, 800);
	MatrixXd G = MatrixXd::Zero(800, 800);
	MatrixXd B = MatrixXd::Zero(800, 800);
	
	double sphcolors[6][3];
	double sphshading[6];
	
	for(int i = 0; i < number; i++)
	{
		sphcolors[i][0] = max(255 * ((i % 3) - 1), 0); //Red
		sphcolors[i][1] = max(255 * (((i + 1) % 3) - 1), 0); //Green
		sphcolors[i][2] = max(255 * (((i + 2) % 3) - 1), 0); //Blue
		
		sphshading[i] = pow(1.5, i); //Phong exponent
	}
	
	//Part 3 declarations
	cout << "Part 3: Perspective viewing\n";
	const string filename3("part3.png");
	
	MatrixXd R3 = MatrixXd::Zero(800, 800);
	MatrixXd G3 = MatrixXd::Zero(800, 800);
	MatrixXd B3 = MatrixXd::Zero(800, 800);
	MatrixXd A3 = MatrixXd::Zero(800, 800);
	
	//Part 5 declarations
	cout << "Part 5: Shadows\n";
	const string filename5("part5.png");
	
	MatrixXd R5 = MatrixXd::Zero(800, 800);
	MatrixXd G5 = MatrixXd::Zero(800, 800);
	MatrixXd B5 = MatrixXd::Zero(800, 800);
	MatrixXd A5 = MatrixXd::Zero(800, 800);
	
	//Part 6 declarations
	cout << "Part 6: Mirror\n";
	const string filename6("part6.png");
	
	MatrixXd R6 = MatrixXd::Zero(800, 800);
	MatrixXd G6 = MatrixXd::Zero(800, 800);
	MatrixXd B6 = MatrixXd::Zero(800, 800);
	MatrixXd A6 = MatrixXd::Zero(800, 800);
	
	
	//Find maxes and mins for camera view range
	for(int i = 0; i < number; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(sph[i][j] + sph[i][3] > posmax)
				posmax = sph[i][j] + sph[i][3];
			if(sph[i][j] - sph[i][3] < -posmax)
				posmax = -(sph[i][j] - sph[i][3]);
		}
	}
	
	//The camera is orthographic, pointing in the direction -z and covering the needed distance to see the spheres
	Vector3d origin(-posmax - 0.5, posmax + 0.5, posmax + 0.5);
	Vector3d x_displacement((2 * posmax + 1)/C.cols(),0,0);
	Vector3d y_displacement(0,(2 * -posmax - 1)/C.rows(),0);

	//Light sources
	Vector3d light_position(-posmax - 0.5, posmax + 0.5, posmax + 0.5);
	Vector3d second_light(posmax + 0.5, -posmax - 0.5, 0); //This light is only used after part 2
	
	//For perspective view in part 3
	Vector3d perspective_view(0, 0, 2 * (posmax + 0.5));
	
	for (int i = 0; i < C.cols(); i++)
	{
		for (int j = 0; j < C.rows(); j++)
		{
			//Prepare the ray
			Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
			Vector3d ray_direction = RowVector3d(0, 0, -1);
			
			//Does it intersect?
			double true_int = DBL_MAX;
			int true_sphere = -1;
			for(int k = 0; k < number; k++)
			{
				Vector3d ray_center = RowVector3d(sph[k][0] - ray_origin(0), sph[k][1] - ray_origin(1), sph[k][2] - ray_origin(2));
				double projection_length = ray_center.dot(ray_direction);
				double dist_apart = sqrt(ray_center.dot(ray_center) - projection_length * projection_length);
				
				//Which intersection is closest?
				if (dist_apart <= sph[k][3])
				{
					double intersect_sec = sqrt(sph[k][3] * sph[k][3] - dist_apart * dist_apart);
					double intersect = projection_length - intersect_sec;
					
					if(intersect < true_int)
					{
						true_int = intersect;
						true_sphere = k;
					}
				}
			}
			
			//Part 1 shading
			if (true_int != DBL_MAX)
			{	
				// The ray hit the sphere, compute the exact intersection point
				Vector3d ray_intersection(ray_origin(0), ray_origin(1), ray_origin(2) - true_int);
				
				Vector3d ray_normal(ray_intersection(0) - sph[true_sphere][0], ray_intersection(1) - sph[true_sphere][1],
									ray_intersection(2) - sph[true_sphere][2]);
				
				// Simple diffuse model
				C(i, j) = (light_position - ray_intersection).normalized().dot(ray_normal);

				// Clamp to zero
				C(i, j) = max(C(i, j), 0.);

				// Disable the alpha mask for this pixel
				A(i, j) = 1;
			}
			//Part 1 end
			
			
			//Part 2 shading, with lambertian, specular, and ambient coeff. as .5, .2, and .25, repectively
			//Light 2 has an intensity of .5
			if(true_int != DBL_MAX)
			{
				Vector3d ray_intersection(ray_origin(0), ray_origin(1), ray_origin(2) - true_int);
				Vector3d light_ray1 = (light_position - ray_intersection).normalized();
				Vector3d light_ray2 = (second_light - ray_intersection).normalized();
				
				Vector3d ray_normal(ray_intersection(0) - sph[true_sphere][0], ray_intersection(1) - sph[true_sphere][1],
									ray_intersection(2) - sph[true_sphere][2]);
				
				//Shading source 1
				double lambert = (light_ray1).dot(ray_normal);
				lambert = max(lambert, 0.0);
				double specular = (light_ray1 - ray_direction).normalized().dot(ray_normal);
				specular = pow(max(specular, 0.0), sphshading[true_sphere]);
				double ambient = 1;
				
				double total_light = lambert * 0.5 + specular * 0.2 + ambient * 0.25;
				
				//Shading source 2
				lambert = (light_ray2).dot(ray_normal);
				lambert = max(lambert, 0.0);
				specular = (light_ray2 - ray_direction).normalized().dot(ray_normal);
				specular = pow(max(specular, 0.0), sphshading[true_sphere]);
				
				total_light = min(total_light + (lambert * 0.5 + specular * 0.2) * 0.5, 1.0);
				
				R(i, j) = sphcolors[true_sphere][0]/255.0 * total_light;
				G(i, j) = sphcolors[true_sphere][1]/255.0 * total_light;
				B(i, j) = sphcolors[true_sphere][2]/255.0 * total_light;
				A(i, j) = 1;
			}
			//Part 2 end
			
			
			//Part 3 shading, which requires that intersection points are recalculated
			//Still uses many of the same functions and calculations
			
			//Prepare the ray
			Vector3d screen_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
			Vector3d perspective_ray = screen_origin - perspective_view;
			ray_direction = perspective_ray.normalized();
			
			//Does it intersect?
			true_int = DBL_MAX;
			true_sphere = -1;
			for(int k = 0; k < number; k++)
			{
				Vector3d ray_center = RowVector3d(sph[k][0] - perspective_view(0), sph[k][1] - perspective_view(1),
												sph[k][2] - perspective_view(2));
				double projection_length = ray_center.dot(ray_direction);
				double dist_apart = sqrt(ray_center.dot(ray_center) - projection_length * projection_length);
				
				//Which intersection is closest?
				if (dist_apart <= sph[k][3])
				{
					double intersect_sec = sqrt(sph[k][3] * sph[k][3] - dist_apart * dist_apart);
					double intersect = projection_length - intersect_sec;
					
					if(intersect < true_int)
					{
						true_int = intersect;
						true_sphere = k;
					}
				}
			}
			
			//Shading for parts 3, 5, and 6
			double total_light_shadow = 0;
			bool light_intersect = false;
			
			if(true_int != DBL_MAX)
			{
				Vector3d ray_intersection = perspective_view + ray_direction * true_int;
				Vector3d light_ray1 = (light_position - ray_intersection).normalized();
				Vector3d light_ray2 = (second_light - ray_intersection).normalized();
				
				Vector3d ray_normal(ray_intersection(0) - sph[true_sphere][0], ray_intersection(1) - sph[true_sphere][1],
									ray_intersection(2) - sph[true_sphere][2]);
				
				//Shading source 1
				double lambert = (light_ray1).dot(ray_normal);
				lambert = max(lambert, 0.0);
				double specular = (light_ray1 - ray_direction).normalized().dot(ray_normal);
				specular = pow(max(specular, 0.0), sphshading[true_sphere]);
				double ambient = 1;
				
				//For part 3
				double total_light = min(lambert * 0.5 + specular * 0.2 + ambient * 0.25, 1.0);
				
				
				//Part 5 determines if I should count this light source
				//Find if there is an intersection between the point and the light
				for(int k = 0; k < number; k++)
				{
					Vector3d ray_center = RowVector3d(sph[k][0] - ray_intersection(0), sph[k][1] - ray_intersection(1),
													sph[k][2] - ray_intersection(2));
					double projection_length = ray_center.dot(light_ray1);
					double dist_apart = DBL_MAX;
					if(projection_length > 0)
						dist_apart = sqrt(ray_center.dot(ray_center) - projection_length * projection_length);
					
					//if there is an intersection, break
					if (dist_apart <= sph[k][3] && k != true_sphere)
					{
						light_intersect = true;
						break;
					}
				}
				
				//Don't add if there was an intersection (except ambient)
				if(!light_intersect)
					total_light_shadow = min(lambert * 0.5 + specular * 0.2 + ambient * 0.25, 1.0);
				else
					total_light_shadow = min(ambient * 0.25, 1.0);

				
				//Shading source 2
				lambert = (light_ray2).dot(ray_normal);
				lambert = max(lambert, 0.0);
				specular = (light_ray2 - ray_direction).normalized().dot(ray_normal);
				specular = pow(max(specular, 0.0), sphshading[true_sphere]);
				
				//For part 3
				total_light = min(total_light + (lambert * 0.5 + specular * 0.2) * 0.5, 1.0);
				
				
				//Do part 5 for this as well
				light_intersect = false;
				for(int k = 0; k < number; k++)
				{
					Vector3d ray_center = RowVector3d(sph[k][0] - ray_intersection(0), sph[k][1] - ray_intersection(1),
													sph[k][2] - ray_intersection(2));
					double projection_length = ray_center.dot(light_ray2);
					double dist_apart = DBL_MAX;
					if(projection_length > 0)
						dist_apart = sqrt(ray_center.dot(ray_center) - projection_length * projection_length);
					
					//if there is an intersection, break
					if (dist_apart <= sph[k][3] && k != true_sphere)
					{
						light_intersect = true;
						break;
					}
				}
				
				if(!light_intersect)
					total_light_shadow = min(total_light_shadow + (lambert * 0.5 + specular * 0.2) * 0.5, 1.0);
				
				
				R3(i, j) = sphcolors[true_sphere][0]/255.0 * total_light;
				G3(i, j) = sphcolors[true_sphere][1]/255.0 * total_light;
				B3(i, j) = sphcolors[true_sphere][2]/255.0 * total_light;
				A3(i, j) = 1;
				
				R5(i, j) = sphcolors[true_sphere][0]/255.0 * total_light_shadow;
				G5(i, j) = sphcolors[true_sphere][1]/255.0 * total_light_shadow;
				B5(i, j) = sphcolors[true_sphere][2]/255.0 * total_light_shadow;
				A5(i, j) = 1;
				
				R6(i, j) = sphcolors[true_sphere][0]/255.0 * total_light_shadow;
				G6(i, j) = sphcolors[true_sphere][1]/255.0 * total_light_shadow;
				B6(i, j) = sphcolors[true_sphere][2]/255.0 * total_light_shadow;
				A6(i, j) = 1;
			}
			//Part 3 end
			//Part 5 end
			
			//Part 6 creates a mirror at the bottom of the scene
			//This part was not done generally, and only mirrors the bottom
			else if(ray_direction(1) < 0)
			{
				//Where it intersects the bottom
				Vector3d new_perspective = perspective_view + ray_direction * (-posmax - 0.5)/ray_direction(1);
				Vector3d new_direction = RowVector3d(ray_direction(0), -ray_direction(1), ray_direction(2));
				
				//If it's in the scene bounds, look for intersections again
				if(new_perspective(2) > (-posmax - 0.5) && new_perspective(0) > (-posmax - 0.5) && new_perspective(0) < (posmax + 0.5))
				{
					true_int = DBL_MAX;
					true_sphere = -1;
					for(int k = 0; k < number; k++)
					{
						Vector3d ray_center = RowVector3d(sph[k][0] - new_perspective(0), sph[k][1] - new_perspective(1),
														sph[k][2] - new_perspective(2));
						double projection_length = ray_center.dot(new_direction);
						double dist_apart = DBL_MAX;
						if(projection_length > 0)
							dist_apart = sqrt(ray_center.dot(ray_center) - projection_length * projection_length);
						
						//Which intersection is closest?
						if (dist_apart <= sph[k][3])
						{
							double intersect_sec = sqrt(sph[k][3] * sph[k][3] - dist_apart * dist_apart);
							double intersect = projection_length - intersect_sec;
							
							if(intersect < true_int)
							{
								true_int = intersect;
								true_sphere = k;
							}
						}
					}
				}
				
				if(true_int != DBL_MAX)
				{
					Vector3d ray_intersection = new_perspective + new_direction * true_int;
					Vector3d light_ray1 = (light_position - ray_intersection).normalized();
					Vector3d light_ray2 = (second_light - ray_intersection).normalized();
					
					Vector3d ray_normal(ray_intersection(0) - sph[true_sphere][0], ray_intersection(1) - sph[true_sphere][1],
										ray_intersection(2) - sph[true_sphere][2]);
					
					//Shading source 1
					double lambert = (light_ray1).dot(ray_normal);
					lambert = max(lambert, 0.0);
					double specular = (light_ray1 - new_direction).normalized().dot(ray_normal);
					specular = pow(max(specular, 0.0), sphshading[true_sphere]);
					double ambient = 1;
					
					for(int k = 0; k < number; k++)
					{
						Vector3d ray_center = RowVector3d(sph[k][0] - ray_intersection(0), sph[k][1] - ray_intersection(1),
														sph[k][2] - ray_intersection(2));
						double projection_length = ray_center.dot(light_ray1);
						double dist_apart = DBL_MAX;
						if(projection_length > 0)
							dist_apart = sqrt(ray_center.dot(ray_center) - projection_length * projection_length);
						
						//if there is an intersection, break
						if (dist_apart <= sph[k][3] && k != true_sphere)
						{
							light_intersect = true;
							break;
						}
					}
					
					//Don't add if there was an intersection (except ambient)
					if(!light_intersect)
						total_light_shadow = min(lambert * 0.5 + specular * 0.2 + ambient * 0.25, 1.0);
					else
						total_light_shadow = min(ambient * 0.25, 1.0);

					
					//Shading source 2
					lambert = (light_ray2).dot(ray_normal);
					lambert = max(lambert, 0.0);
					specular = (light_ray2 - new_direction).normalized().dot(ray_normal);
					specular = pow(max(specular, 0.0), sphshading[true_sphere]);
					
					
					light_intersect = false;
					for(int k = 0; k < number; k++)
					{
						Vector3d ray_center = RowVector3d(sph[k][0] - ray_intersection(0), sph[k][1] - ray_intersection(1),
														sph[k][2] - ray_intersection(2));
						double projection_length = ray_center.dot(light_ray2);
						double dist_apart = DBL_MAX;
						if(projection_length > 0)
							dist_apart = sqrt(ray_center.dot(ray_center) - projection_length * projection_length);
						
						//if there is an intersection, break
						if (dist_apart <= sph[k][3] && k != true_sphere)
						{
							light_intersect = true;
							break;
						}
					}
					
					if(!light_intersect)
						total_light_shadow = min(total_light_shadow + (lambert * 0.5 + specular * 0.2) * 0.5, 1.0);
					
					
					R6(i, j) = sphcolors[true_sphere][0]/255.0 * total_light_shadow;
					G6(i, j) = sphcolors[true_sphere][1]/255.0 * total_light_shadow;
					B6(i, j) = sphcolors[true_sphere][2]/255.0 * total_light_shadow;
					A6(i, j) = 1;
				}
			}
			//Part 6 end
		}
	}
	
	//Part 4 declarations
	cout << "Part 4: Bumpy cube and bunny in OFF format\n";
	const string filename4("part4.png");
	
	MatrixXd R4 = MatrixXd::Zero(800, 800);
	MatrixXd G4 = MatrixXd::Zero(800, 800);
	MatrixXd B4 = MatrixXd::Zero(800, 800);
	MatrixXd A4 = MatrixXd::Zero(800, 800);
	
	string trashstring;
	int cubevertices = 0, cubefaces = 0, bunnyvertices = 0, bunnyfaces = 0, trashint;
	posmax = 0;
	
	int trianglecolors[4][2]; //RGB and shading (4) x cube and bunny (2)
	trianglecolors[0][0] = 255;
	trianglecolors[1][0] = 0;
	trianglecolors[2][0] = 0;
	trianglecolors[3][0] = 10;
	trianglecolors[0][1] = 0;
	trianglecolors[1][1] = 255;
	trianglecolors[2][1] = 0;
	trianglecolors[3][1] = 50;
	
	ifstream read;
	read.open("../data/bumpy_cube.off");
	
	//Part 4 new scene with shapes loaded in OFF format (triangles only)
	//Read in the data for the cube, offset by 4 in the x direction
	read >> trashstring;
	read >> cubevertices >> cubefaces >> trashint;
	Matrix<double, Dynamic, Dynamic> cubevertexarray = Matrix<double, Dynamic, Dynamic>::Zero(cubevertices, 3);
	Matrix<int, Dynamic, Dynamic> cubefacearray = Matrix<int, Dynamic, Dynamic>::Zero(cubefaces, 3);
	
	for(int i = 0; i < cubevertices; i++)
	{
		read >> cubevertexarray(i, 0) >> cubevertexarray(i, 1) >> cubevertexarray(i, 2);
		cubevertexarray(i, 0) = cubevertexarray(i, 0) + 4.0;
	}
	
	for(int i = 0; i < cubefaces; i++)
		read >> trashint >> cubefacearray(i, 0) >> cubefacearray(i, 1) >> cubefacearray(i, 2);
	
	//Find maxes and mins for camera view range
	for(int i = 0; i < cubevertices; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(cubevertexarray(i, j) > posmax)
				posmax = cubevertexarray(i, j);
			if(cubevertexarray(i, j) < -posmax)
				posmax = -(cubevertexarray(i, j));
		}
	}
	
	read.close();
	read.open("../data/bunny.off");
	
	//Repeat for the bunny, but offset by 4 in the -x direction and also scaled by 20 times
	read >> trashstring;
	read >> bunnyvertices >> bunnyfaces >> trashint;
	Matrix<double, Dynamic, Dynamic> bunnyvertexarray = Matrix<double, Dynamic, Dynamic>::Zero(bunnyvertices, 3);
	Matrix<int, Dynamic, Dynamic> bunnyfacearray = Matrix<int, Dynamic, Dynamic>::Zero(bunnyfaces, 3);
	
	for(int i = 0; i < bunnyvertices; i++)
	{
		read >> bunnyvertexarray(i, 0) >> bunnyvertexarray(i, 1) >> bunnyvertexarray(i, 2);
		bunnyvertexarray(i, 0) = bunnyvertexarray(i, 0) * 20 - 4.0;
		bunnyvertexarray(i, 1) = bunnyvertexarray(i, 1) * 20;
		bunnyvertexarray(i, 2) = bunnyvertexarray(i, 2) * 20;
	}
	
	for(int i = 0; i < bunnyfaces; i++)
		read >> trashint >> bunnyfacearray(i, 0) >> bunnyfacearray(i, 1) >> bunnyfacearray(i, 2);
	
	for(int i = 0; i < bunnyvertices; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(bunnyvertexarray(i, j) > posmax)
				posmax = bunnyvertexarray(i, j);
			if(bunnyvertexarray(i, j) < -posmax)
				posmax = -(bunnyvertexarray(i, j));
		}
	}
	
	read.close();
	
	origin = RowVector3d(-posmax - 0.5, posmax + 0.5, posmax + 0.5);
	x_displacement = RowVector3d((2 * posmax + 1)/R4.cols(), 0, 0);
	y_displacement = RowVector3d(0, (2 * -posmax - 1)/R4.rows(), 0);

	//Light sources
	light_position = RowVector3d(-posmax - 0.5, posmax + 0.5, posmax + 0.5);
	second_light = RowVector3d(posmax + 0.5, -posmax - 0.5, 0);
	
	//For perspective view
	perspective_view = RowVector3d(0, 0, 2 * (posmax + 0.5));
	
	Vector3d screen_origin, perspective_ray, ray_direction;
	Vector3d v0, v1, v2, e0, e1, e2, trinormal;
	Vector3d to_intersect, vertexint0, vertexint1, vertexint2;
	int v0index, v1index, v2index;
	
	for (int i = 0; i < R4.cols(); i++)
	{
		for (int j = 0; j < R4.rows(); j++)
		{
			//Prepare the ray
			screen_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
			perspective_ray = screen_origin - perspective_view;
			ray_direction = perspective_ray.normalized();
			
			//Does it intersect?
			double true_int = DBL_MAX;
			int true_triangle = -1;
			Vector3d true_normal(0, 0, 0);
			
			//Cube first
			for(int k = 0; k < cubefaces; k++)
			{
				v0index = cubefacearray(k, 0);
				v1index = cubefacearray(k, 1);
				v2index = cubefacearray(k, 2);
				v0 = RowVector3d(cubevertexarray(v0index, 0), cubevertexarray(v0index, 1), cubevertexarray(v0index, 2));
				v1 = RowVector3d(cubevertexarray(v1index, 0), cubevertexarray(v1index, 1), cubevertexarray(v1index, 2));
				v2 = RowVector3d(cubevertexarray(v2index, 0), cubevertexarray(v2index, 1), cubevertexarray(v2index, 2));
				
				//Find the triangle normal
				e1 = v1 - v0;
				e2 = v2 - v0;
				trinormal = e1.cross(e2);
				trinormal = trinormal.normalized();
				
				//Not parallel
				if(trinormal.dot(ray_direction) == 0)
					continue;
				
				//Intersection with triangle?
				double dist_apart = (trinormal.dot(v0) - trinormal.dot(perspective_view)) / trinormal.dot(ray_direction);
				to_intersect = perspective_view + ray_direction * dist_apart;
				vertexint0 = to_intersect - v0;
				vertexint1 = to_intersect - v1;
				vertexint2 = to_intersect - v2;
				
				e0 = v1 - v0;
				e1 = v2 - v1;
				e2 = v0 - v2;
				
				//Find closest intersection
				if(trinormal.dot(e0.cross(vertexint0)) >= 0.0 && trinormal.dot(e1.cross(vertexint1)) >= 0.0 && trinormal.dot(e2.cross(vertexint2)) >= 0.0)
				{
					if(dist_apart < true_int)
					{
						//cout << "Here " << dist_apart << '\n';
						true_int = dist_apart;
						true_normal = trinormal;
						true_triangle = 0;
					}
				}
			}
			
			//Then bunny
			for(int k = 0; k < bunnyfaces; k++)
			{
				v0index = bunnyfacearray(k, 0);
				v1index = bunnyfacearray(k, 1);
				v2index = bunnyfacearray(k, 2);
				v0 = Vector3d(bunnyvertexarray(v0index, 0), bunnyvertexarray(v0index, 1), bunnyvertexarray(v0index, 2));
				v1 = Vector3d(bunnyvertexarray(v1index, 0), bunnyvertexarray(v1index, 1), bunnyvertexarray(v1index, 2));
				v2 = Vector3d(bunnyvertexarray(v2index, 0), bunnyvertexarray(v2index, 1), bunnyvertexarray(v2index, 2));
				
				//Find the triangle normal
				e1 = v1 - v0;
				e2 = v2 - v0;
				trinormal = e1.cross(e2);
				trinormal = trinormal.normalized();
				
				//Not parallel
				if(trinormal.dot(ray_direction) == 0)
					continue;
				
				//Intersection with triangle?
				double dist_apart = (trinormal.dot(v0) - trinormal.dot(perspective_view)) / trinormal.dot(ray_direction);
				to_intersect = perspective_view + ray_direction * dist_apart;
				vertexint0 = to_intersect - v0;
				vertexint1 = to_intersect - v1;
				vertexint2 = to_intersect - v2;
				
				e0 = v1 - v0;
				e1 = v2 - v1;
				e2 = v0 - v2;
				
				//Find closest intersection
				if(trinormal.dot(e0.cross(vertexint0)) >= 0.0 && trinormal.dot(e1.cross(vertexint1)) >= 0.0 && trinormal.dot(e2.cross(vertexint2)) >= 0.0)
				{
					if(dist_apart < true_int)
					{
						//cout << "bunny " << dist_apart << '\n';
						true_int = dist_apart;
						true_normal = trinormal;
						true_triangle = 1;
					}
				}
			}
			
			if(true_int < DBL_MAX)
			{
				Vector3d ray_intersection = perspective_view + ray_direction * true_int;
				Vector3d light_ray1 = (light_position - ray_intersection).normalized();
				Vector3d light_ray2 = (second_light - ray_intersection).normalized();
				
				Vector3d ray_normal = true_normal;
				
				//Shading source 1
				double lambert = (light_ray1).dot(ray_normal);
				lambert = max(lambert, 0.0);
				double specular = (light_ray1 - ray_direction).normalized().dot(ray_normal);
				specular = pow(max(specular, 0.0), trianglecolors[3][true_triangle]);
				double ambient = 1;
				
				double total_light = lambert * 0.5 + specular * 0.2 + ambient * 0.25;
				
				//Shading source 2
				lambert = (light_ray2).dot(ray_normal);
				lambert = max(lambert, 0.0);
				specular = (light_ray2 - ray_direction).normalized().dot(ray_normal);
				specular = pow(max(specular, 0.0), trianglecolors[3][true_triangle]);
				
				total_light = min(total_light + (lambert * 0.5 + specular * 0.2) * 0.5, 1.0);
				
				R4(i, j) = trianglecolors[0][true_triangle]/255.0 * total_light;
				G4(i, j) = trianglecolors[1][true_triangle]/255.0 * total_light;
				B4(i, j) = trianglecolors[2][true_triangle]/255.0 * total_light;
				A4(i, j) = 1;
			}
		}
		cout << double(i) * 100/R4.cols() << "%\n";
	}
	
	
	
	
	// Save parts 1-6 to png
	write_matrix_to_png(C, C, C, A, filename1);
	write_matrix_to_png(R, G, B, A, filename2);
	write_matrix_to_png(R3, G3, B3, A3, filename3);
	write_matrix_to_png(R4, G4, B4, A4, filename4);
	write_matrix_to_png(R5, G5, B5, A5, filename5);
	write_matrix_to_png(R6, G6, B6, A6, filename6);
}


int main()
{
	allparts();
	
	return 0;
}