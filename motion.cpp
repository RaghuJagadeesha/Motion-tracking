// Motion.cpp : Defines the entry point for the console application.
//

#include <afxwin.h>  // necessary for MFC to work properly
#include "motion.h"
#include "../../src/blepo.h"
#include <math.h>
#include <stack>
#ifdef _DEBUG
#define new DEBUG_NEW
#endif

using namespace blepo;

class moved_pos{
public:
	float x_pos;
	float y_pos;
	moved_pos()
	{
		x_pos = 0;
		y_pos = 0;
	}
};

class feature_points{
public:
	moved_pos position;
	int frame;
	feature_points()
	{
		frame = 0;
		position = moved_pos();
	}
};
int gethalfwidth(float sigma);
void computegaussian(float kern[], int width, float sigma);
void computedervgaussian(float dervkern[], float kern[], int width, float sigma);
void compute_gradient(ImgGray img, int width, int halfwidth, ImgFloat *outx, ImgFloat *outy, float dervkern[], float kern[]);
void compute_covariance(ImgFloat outx, ImgFloat outy, int window_size, int x_pos, int y_pos, float covar_mat[2][2]);
void compute_covariance(ImgFloat outx, ImgFloat outy, int window_size, float x_pos, float y_pos, float covar_mat[2][2]);
moved_pos lucaskanade(ImgGray prev_img, ImgGray current_img, feature_points good_features, int window_size, ImgFloat gradx, ImgFloat grady, float covar_mat[2][2]);
void compute_err_vect(ImgGray prev_img, ImgGray current_img, ImgFloat gradx, ImgFloat grady, float x_pos, float y_pos, int window_size, moved_pos *u, float err[2][1]);
moved_pos compute_disp(float covar_mat[2][2], float err[2][1]);
float interpolate(ImgGray img, float x, float y);
float interpolate(ImgFloat img, float x, float y);


int main(int argc, const char* argv[], const char* envp[])
{
	// Initialize MFC and return if failure
	HMODULE hModule = ::GetModuleHandle(NULL);
	if (hModule == NULL || !AfxWinInit(hModule, NULL, ::GetCommandLine(), 0))
	{
		printf("Fatal Error: MFC initialization failed (hModule = %x)\n", hModule);
		return 1;
	}
	if (argc != 6)
	{
		printf("Usage : filename-format first-frame last-frame sigma window-size\n");
		exit(0);
	}
	try
	{
		CString file;
		//const char* format = "img%04d.bmp";  // format string – could be CString instead
	//	const char* format = "img%03d.pgm";  // format string – could be CString instead
		const char* format = argv[1];
		file.Format(format, atoi(argv[2]));     // ‘filename’ is now “img0400.bmp”
		const char *filename;
		//CString path = "../../../statue_sequence/statue_seq/";
		CString path = "../../images/";
		string fn;
		fn = path + file;
		filename = fn.c_str();
		float kern[100];
		float dervkern[100];
		int index = 0;
		//printf("%s", filename);
		for (int i = 0; i < 100; i++)
		{
			kern[i] = 0;
			dervkern[i] = 0;
		}

		int width;
		/*Sigma hardcoded to 1 and window size to 3 */
		width = 2 * gethalfwidth(1) + 1;

		computegaussian(kern, width, 1);
		computedervgaussian(dervkern, kern, width, 1);

		ImgGray input_img;
		Load(filename, &input_img);

		ImgFloat outx;
		ImgFloat outy;
		ImgFloat cornerness;

		outx.Reset(input_img.Width(), input_img.Height());
		Set(&outx, 0);
		outy.Reset(input_img.Width(), input_img.Height());
		Set(&outy, 0);
		cornerness.Reset(input_img.Width(), input_img.Height());
		Set(&cornerness, 0);
		
		for (int i = 0; i < width / 2; i++)
		{
			float temp = kern[i];
			kern[i] = kern[width - i-1];
			kern[width - i-1] = temp;
		}
		for (int i = 0; i < width / 2; i++)
		{
			float temp = dervkern[i];
			dervkern[i] = dervkern[width - i-1];
			dervkern[width - i-1] = temp;
		}
		//printf("kernel is \n");
		//for (int i = 0; i < width; i++)
		//{
		//	printf("%f\n", kern[i]);
		//}
		//printf("Derived kernel is \n");
		//for (int i = 0; i < width; i++)
		//{
		//	printf("%f\n", dervkern[i]);
		//}
		compute_gradient(input_img, width, (width - 1) / 2, &outx, &outy, dervkern, kern);

		int i, j;
		float covar_mat[2][2] = { 0, 0, 0, 0 };
		float eigen_1, eigen_2;
		float threshold = 1350;

		for (j = 1; j < input_img.Height() - 1; j++)
		{
			for (i = 1; i < input_img.Width() - 1; i++)
			{
				compute_covariance(outx, outy, 3, i, j, covar_mat);
				eigen_1 = (0.5*(covar_mat[0][0] + covar_mat[1][1] + sqrt((pow(covar_mat[0][0] - covar_mat[1][1], 2)) + (4 * (pow(covar_mat[0][1], 2))))));
				eigen_2 = (0.5*(covar_mat[0][0] + covar_mat[1][1] - sqrt((pow(covar_mat[0][0] - covar_mat[1][1], 2)) + (4 * (pow(covar_mat[0][1], 2))))));
				if (eigen_2 >= threshold)
				{
					cornerness(i, j) = eigen_2;
				}
			}
		}

		int x, y;

		for (j = 1; j < input_img.Height() - 1; j++)
		{
			for (i = 1; i < input_img.Width() - 1; i++)
			{
				for (y = -1; y <= 1; y++)
				{
					for (x = -1; x <= 1; x++)
					{
						if (x == 0 && y == 0)
							continue;
						if (cornerness(i, j) < cornerness(i + x, j + y))
							cornerness(i, j) = 0;
					}
				}
			}
		}
		for (j = 3; j < input_img.Height() - 3; j++)
		{
			for (i = 3; i < input_img.Width() - 3; i++)
			{
				for (y = -3; y <= 3; y++)
				{
					for (x = -3; x <= 3; x++)
					{
						if (cornerness(i, j) != 0 && x != 0 && y != 0)
						{
							cornerness(i + x, j + y) = 0;
						}
					}
				}
			}
		}

		std::vector<feature_points> good_features;
		for (j = 0; j < input_img.Height(); j++)
		{
			for (i = 0; i < input_img.Width(); i++)
			{
				if (cornerness(i, j) != 0)
				{
					feature_points temp;
					temp.position.x_pos = i;
					temp.position.y_pos = j;
					temp.frame = atoi(argv[2]);
					good_features.push_back(temp);
				}
			}
		}

		moved_pos u;
		int window_size;
		float sigma = atof(argv[4]);
		window_size = atoi(argv[5]);
		ImgGray previous_img;
		ImgGray current_img;
		ImgBgr initial_img;
		Load(filename, &initial_img);
		for (j = 0; j < good_features.size(); j++)
		{
			printf("%f,%f\n", good_features[j].position.x_pos, good_features[j].position.y_pos);
			for (y = -1; y <= 1; y++)
			{
				for (x = -1; x <= 1; x++)
				{
					if (x == 0 && y == 0)
						continue;
					if (good_features[j].position.x_pos + x < 0 || good_features[j].position.y_pos + y < 0)
						continue;
					if (good_features[j].position.x_pos + x >= initial_img.Width() || good_features[j].position.y_pos + y >= initial_img.Height())
						continue;
					initial_img(good_features[j].position.x_pos + x, good_features[j].position.y_pos + y).r = 255;
					initial_img(good_features[j].position.x_pos + x, good_features[j].position.y_pos + y).g = 0;
					initial_img(good_features[j].position.x_pos + x, good_features[j].position.y_pos + y).b = 0;
				}
			}
		}
		Figure fig1;
		fig1.SetTitle("Initial Image with features");
		fig1.Draw(initial_img);
		ImgBgr display_img;

		Load(filename, &previous_img);

		width = 2 * gethalfwidth(sigma) + 1;

		for (int i = 0; i < 100; i++)
		{
			kern[i] = 0;
			dervkern[i] = 0;
		}
		computegaussian(kern, width, sigma);
		computedervgaussian(dervkern, kern, width, sigma);
		for (int i = 0; i < width / 2; i++)
		{
			float temp = kern[i];
			kern[i] = kern[width - i-1];
			kern[width - i-1] = temp;
		}
		for (int i = 0; i < width / 2; i++)
		{
			float temp = dervkern[i];
			dervkern[i] = dervkern[width - i-1];
			dervkern[width - i-1] = temp;
		}
		//printf("kernel is \n");
		//for (int i = 0; i < width; i++)
		//{
		//	printf("%f\n", kern[i]);
		//}
		//printf("Derived kernel is \n");
		//for (int i = 0; i < width; i++)
		//{
		//	printf("%f\n", dervkern[i]);
		//}
		Figure fig2;
		fig2.SetTitle("Final Image with tracking");
		for (int frames = atoi(argv[2]) + 1; frames <= atoi(argv[3]); frames++)
		{
			file.Format(format, frames);
			fn = path + file;
			filename = fn.c_str();
			Load(filename, &current_img);
			Load(filename, &display_img);
			outx.Reset(input_img.Width(), input_img.Height());
			Set(&outx, 0);
			outy.Reset(input_img.Width(), input_img.Height());
			Set(&outy, 0);
			compute_gradient(previous_img, width, ((width - 1) / 2), &outx, &outy, dervkern, kern);
			for (j = 0; j < good_features.size(); j++)
			{
				u = moved_pos();
				compute_covariance(outx, outy, window_size, good_features[j].position.x_pos, good_features[j].position.y_pos, covar_mat);
				u = lucaskanade(previous_img, current_img, good_features[j], window_size, outx, outy, covar_mat);
				good_features[j].position.x_pos += u.x_pos;
				good_features[j].position.y_pos += u.y_pos;
				good_features[j].frame = frames;
				//if (j == 38)
				//	printf("%f,%f,%d\n", good_features[j].position.x_pos, good_features[j].position.y_pos, frames);
			}
			previous_img = current_img;
			//for (j = 0; j < good_features.size(); j++)
			//{
			//	if (round(good_features[j].position.x_pos) < 0 || round(good_features[j].position.y_pos) < 0)
			//	{
			//		good_features.erase(good_features.begin() + j);
			//		j = 0;
			//		continue;
			//	}
			//	if (round(good_features[j].position.x_pos) > initial_img.Width() || round(good_features[j].position.y_pos) > initial_img.Height())
			//	{
			//		good_features.erase(good_features.begin() + j);
			//		j = 0;
			//		continue;
			//	}
			//}
			for (j = 0; j < good_features.size(); j++)
			{
			//	printf("%f,%f, %d\n", good_features[j].position.x_pos, good_features[j].position.y_pos, good_features[j].frame);
				for (y = -1; y <= 1; y++)
				{
					for (x = -1; x <= 1; x++)
					{
						if (x == 0 && y == 0)
							continue;
						if (round(good_features[j].position.x_pos) + x < 0 || round(good_features[j].position.y_pos) + y < 0)
							continue;
						if (round(good_features[j].position.x_pos) + x >= initial_img.Width() || round(good_features[j].position.y_pos) + y >= initial_img.Height())
							continue;
						display_img(round(good_features[j].position.x_pos) + x, round(good_features[j].position.y_pos) + y).r = 255;
						display_img(round(good_features[j].position.x_pos) + x, round(good_features[j].position.y_pos) + y).g = 0;
						display_img(round(good_features[j].position.x_pos) + x, round(good_features[j].position.y_pos) + y).b = 0;
					}
				}
			}
			fig2.Draw(display_img);
		}
		EventLoop();
	}
	catch (const Exception& e)
	{
		e.Display();
	}
}

int gethalfwidth(float sigma)
{
	float val = round(2.5*sigma - 0.5);
	return (int)(round(val));
}

void computegaussian(float kern[], int width, float sigma)
{
	int i = 0;
	float sum = 0;
	for (i = 0; i < width; i++)
	{
		kern[i] = exp(-pow((i - (width - 1) / 2), 2) / (2 * pow(sigma, 2)));
		sum += kern[i];
	}
	for (i = 0; i < width; i++)
		kern[i] /= sum;
}

void computedervgaussian(float dervkern[], float kern[], int width, float sigma)
{
	int i = 0;
	float sum = 0;
	for (i = 0; i < width; i++)
	{
		dervkern[i] = (i - (width - 1) / 2) * kern[i];
		sum += (i*dervkern[i]);
	}
	for (i = 0; i < width; i++)
		dervkern[i] /= sum;
}
void compute_gradient(ImgGray img, int width, int halfwidth, ImgFloat *outx, ImgFloat *outy, float dervkern[], float kern[])
{
	float sum = 0;
	int x, y, i, j;
	int d;
	float *temp;
	temp = (float *)calloc(img.Height()*img.Width(), sizeof(float));
	for (y = 0; y < img.Height(); y++)
	{
		for (x = 0; x < img.Width(); x++)
		{
			sum = 0;
			for (i = 0; i < width; i++)
			{
				if (x + halfwidth - i < 0)
					d = (x + halfwidth - i)*-1;
				else if (x + halfwidth - i >= img.Width())
					d = img.Width() - (x + halfwidth - i) - 1;
				else
					d = 0;
				sum += ((img(x + halfwidth - i + d, y))*(dervkern[i]));
			}
			temp[y*img.Width() + x] = sum;
		}
	}
	for (y = halfwidth; y < img.Height() - halfwidth - 1; y++)
	{
		for (x = 0; x < img.Width(); x++)
		{
			sum = 0;
			for (i = 0; i < width; i++)
			{
				if (y + halfwidth - i < 0)
					d = (y + halfwidth - i)*-1;
				else if (y + halfwidth - i >= img.Width())
					d = img.Height() - (y + halfwidth - i) - 1;
				else
					d = 0;
				sum += ((temp[(y + halfwidth - i + d)* img.Width() + x])*(kern[i]));
			}
			(*outx)(x, y) = sum;
		}
	}

	for (y = 0; y < img.Height(); y++)
	{
		for (x = halfwidth; x < img.Width() - halfwidth - 1; x++)
		{
			sum = 0;
			for (i = 0; i < width; i++)
			{
				if (x + halfwidth - i < 0)
					d = (x + halfwidth - i)*-1;
				else if (x + halfwidth - i >= img.Width())
					d = img.Width() - (x + halfwidth - i) - 1;
				else
					d = 0;
				sum += ((img(x + halfwidth - i + d, y))*(kern[i]));
			}
			temp[y*img.Width() + x] = sum;
		}
	}
	for (y = halfwidth; y < img.Height() - halfwidth - 1; y++)
	{
		for (x = 0; x < img.Width(); x++)
		{
			sum = 0;
			for (i = 0; i < width; i++)
			{
				if (y + halfwidth - i < 0)
					d = (y + halfwidth - i)*-1;
				else if (y + halfwidth - i >= img.Width())
					d = img.Height() - (y + halfwidth - i) - 1;
				else
					d = 0;
				sum += ((temp[(y + halfwidth - i + d)*img.Width() + x])*(dervkern[i]));
			}
			(*outy)(x, y) = sum;
		}
	}

}
void compute_covariance(ImgFloat outx, ImgFloat outy, int window_size, int x_pos, int y_pos, float covar_mat[2][2])
{
	covar_mat[0][0] = 0;
	covar_mat[0][1] = 0;
	covar_mat[1][0] = 0;
	covar_mat[1][1] = 0;
	for (int y = -(window_size - 1) / 2; y <= (window_size - 1) / 2; y++)
	{
		for (int x = -(window_size - 1) / 2; x <= (window_size - 1) / 2; x++)
		{
			covar_mat[0][0] += pow(outx(x_pos + x, y_pos + y), 2);
			covar_mat[1][1] += pow(outy(x_pos + x, y_pos + y), 2);
			covar_mat[1][0] += outx(x_pos + x, y_pos + y)*outy(x_pos + x, y_pos + y);
			covar_mat[0][1] = covar_mat[1][0];
		}
	}
}

void compute_covariance(ImgFloat outx, ImgFloat outy, int window_size, float x_pos, float y_pos, float covar_mat[2][2])
{
	covar_mat[0][0] = 0;
	covar_mat[0][1] = 0;
	covar_mat[1][0] = 0;
	covar_mat[1][1] = 0;
	for (int y = -(window_size - 1) / 2; y <= (window_size - 1) / 2; y++)
	{
		for (int x = -(window_size - 1) / 2; x <= (window_size - 1) / 2; x++)
		{
			covar_mat[0][0] += pow(interpolate(outx, x_pos + x, y_pos + y), 2);
			covar_mat[1][1] += pow(interpolate(outy, x_pos + x, y_pos + y), 2);
			covar_mat[1][0] += interpolate(outx, x_pos + x, y_pos + y)*interpolate(outy, x_pos + x, y_pos + y);
			covar_mat[0][1] = covar_mat[1][0];
		}
	}
}

moved_pos lucaskanade(ImgGray prev_img, ImgGray current_img, feature_points good_features, int window_size, ImgFloat gradx, ImgFloat grady, float covar_mat[2][2])
{
	float err[2][1] = { 0, 0 };
	moved_pos u = moved_pos();
	u.x_pos = 0;
	u.y_pos = 0;
	float threshold = 0.001;
	for (int iter = 0; iter < 3; iter++)
	{
		moved_pos temp = moved_pos();
		compute_err_vect(prev_img, current_img, gradx, grady, good_features.position.x_pos, good_features.position.y_pos, window_size, &u, err);
		temp = compute_disp(covar_mat, err);
		u.x_pos += temp.x_pos;
		u.y_pos += temp.y_pos;
		if (u.x_pos < threshold && u.y_pos < threshold)
			break;
	}
	return u;
}
void compute_err_vect(ImgGray prev_img, ImgGray current_img, ImgFloat gradx, ImgFloat grady, float x_pos, float y_pos, int window_size, moved_pos *u, float err[2][1])
{
	err[0][0] = 0;
	err[1][0] = 0;
	for (int y = -(window_size - 1) / 2; y <= (window_size - 1) / 2; y++)
	{
		for (int x = -(window_size - 1) / 2; x <= (window_size - 1) / 2; x++)
		{
			err[0][0] += interpolate(gradx, x_pos + x, y_pos + y) * ((interpolate(prev_img, x_pos + x, y_pos + y)) - interpolate(current_img, x_pos + x + u->x_pos, y_pos + y + u->y_pos));
			err[1][0] += interpolate(grady, x_pos + x, y_pos + y) * ((interpolate(prev_img, x_pos + x, y_pos + y)) - interpolate(current_img, x_pos + x + u->x_pos, y_pos + y + u->y_pos));
		}
	}
}
moved_pos compute_disp(float covar_mat[2][2], float err[2][1])
{
	moved_pos u = moved_pos();
	float det = (covar_mat[0][0] * covar_mat[1][1]) - (covar_mat[0][1] * covar_mat[0][1]);
	if (det == 0)
		det = 1;
	u.x_pos = (1 / det)*(covar_mat[1][1] * err[0][0] - covar_mat[0][1] * err[1][0]);
	u.y_pos = (1 / det)*(covar_mat[0][0] * err[1][0] - covar_mat[1][0] * err[0][0]);
	return u;
}

float interpolate(ImgGray img, float x, float y)
{
	int x0 = x;
	int y0 = y;
	int x1 = x0 + 1;
	int y1 = y0 + 1;
	float r1 = x - x0;
	float r2 = y - y0;
	if (x0 < 0)
	{
		x0 = 0;
		x1 = 0;
	}
	if (y0 < 0)
	{
		y0 = 0;
		y1 = 0;
	}
	if (x0 >= img.Width() - 1)
	{
		x0 = img.Width() - 1;
		x1 = img.Width() - 1;
	}
	if (y0 >= img.Height() - 1)
	{
		y0 = img.Height() - 1;
		y1 = img.Height() - 1;
	}
//	printf(" ingray img %d,%d\n", x0, y0);
	float val = (((1 - r1)*(1 - r2)*(img(x0, y0))) + (r1 * (1 - r2)*img(x1, y0)) + ((1 - r1)*r2*img(x0, y1)) + (r1*r2*img(x1, y1)));
	return val;
}

float interpolate(ImgFloat img, float x, float y)
{
	int x0 = x;
	int y0 = y;
	int x1 = x0 + 1;
	int y1 = y0 + 1;
	float r1 = x - x0;
	float r2 = y - y0; 
	if (x0 < 0)
	{
		x0 = 0;
		x1 = 0;
	}
	if (y0 < 0)
	{
		y0 = 0;
		y1 = 0;
	}	
	if (x0 >= img.Width() - 1)
	{
		x0 = img.Width() - 1;
		x1 = img.Width() - 1;
	}
	if (y0 >= img.Height() - 1)
	{
		y0 = img.Height() - 1;
		y1 = img.Height() - 1;
	}
//	printf("in float img %d,%d\n", x0, y0);
	float val = (((1 - r1)*(1 - r2)*(img(x0, y0))) + (r1 * (1 - r2)*img(x1, y0)) + ((1 - r1)*r2*img(x0, y1)) + (r1*r2*img(x1, y1)));
	return val;
}
