#include <qstring.h>
#include "myWidget.h"

myWidget* mywidget;

int TCSL = 0/*layer_num*/;
int sssflag = 0;

QList<Mdl_Face> FACELIST;
QList<QList<double>> XingBianChangList;

struct Mdl_Point
{
	int pNum;
	float x, y, z;
	int CF;
	Mdl_Point() {}
	Mdl_Point(int num, double X, double Y, double Z)
	{
		pNum = num;
		x = X;
		y = Y;
		z = Z;
	}
	int yflag;
	double sx1;
	~Mdl_Point() {}
};

struct Mdl_Triangle
{
	Mdl_Triangle() {}
	~Mdl_Triangle() {}

	int triNum;
	Mdl_Point p1, p2, p3;
	bool isOverlap;
	int AdjacentTriID_1, AdjacentTriID_2, AdjacentTriID_3;
	int displayOrder;
	int Level;
};

struct Mdl_Tetrahedra
{
	Mdl_Tetrahedra() {}
	Mdl_Tetrahedra(Mdl_Point P1, Mdl_Point P2, Mdl_Point P3, Mdl_Point P4)
	{
		p1 = P1;
		p2 = P2;
		p3 = P3;
		p4 = P4;
	}
	~Mdl_Tetrahedra() {}

	int tetraNum;
	Mdl_Point p1, p2, p3, p4;

	// center
	double Tet_Centerx = (p1.x + p2.x + p3.x + p4.x) / 4;
	double Tet_Centery = (p1.y + p2.y + p3.y + p4.y) / 4;
	double Tet_Centerz = (p1.z + p2.z + p3.z + p4.z) / 4;

	int Sum_yflag;
	double TET_HFTJ;

	double Top_Tet_Porosity = 1;
	double Tet_Porosity = (p1.sx1 + p2.sx1 + p3.sx1 + p4.sx1) / 4;

	double Get_Tet_Volume()
	{
		double A1, A2, A3, A4, value;
		A1 = p2.x * (p3.y * p4.z - p4.y * p3.z) - p3.x * (p2.y * p4.z - p4.y * p2.z) + p4.x * (p2.y * p3.z - p3.y * p2.z);
		A2 = p1.x * (p3.y * p4.z - p4.y * p3.z) - p3.x * (p1.y * p4.z - p4.y * p1.z) + p4.x * (p1.y * p3.z - p3.y * p1.z);
		A3 = p1.x * (p2.y * p4.z - p4.y * p2.z) - p2.x * (p1.y * p4.z - p4.y * p1.z) + p4.x * (p1.y * p2.z - p2.y * p1.z);
		A4 = p1.x * (p2.y * p3.z - p3.y * p2.z) - p2.x * (p1.y * p3.z - p3.y * p1.z) + p3.x * (p1.y * p2.z - p2.y * p1.z);
		value = (A1 - A2 + A3 - A4) / 6;
		return value;
	}
};

struct Mdl_Face
{
	int faceNum; 
	QList<Mdl_Point> triPointList; 
	QList<Mdl_Triangle> triList; 
	double Face_Volume; 

};

void on_timeslider_clicked()
{
	mywidget->editShow(QString::fromLocal8Bit("= = = = = = = = = = ="));
	mywidget->editShow(QString::fromLocal8Bit("select number£º") + QString::number(mywidget->timeslider->value()));
	mywidget->editShow(QString::fromLocal8Bit("= = = = = = = = = = ="));

	if (mywidget->timeslider->value() != sssflag)
	{
		mywidget->editShow(QString::fromLocal8Bit("start..."));

		XingBianChangList.clear();

		for (int tnum = 0; tnum < TCSL + 2; tnum++)
		{
			long double SUM = 0;
			int colornum = tnum;
			for (int i = 0; i < FACELIST[tnum].triPointList.size(); i++)
			{
				SUM += FACELIST[tnum].triPointList[i].z;
			}
			double avz = SUM / FACELIST[tnum].triPointList.size();
			QList<Mdl_Point> pointlist;
			for (int i = 0; i < FACELIST[tnum].triPointList.size(); i++)
			{
				Mdl_Point p;
				p.pNum = FACELIST[tnum].triPointList[i].pNum;
				p.CF = FACELIST[tnum].triPointList[i].CF;
				p.x = FACELIST[tnum].triPointList[i].x;
				p.y = FACELIST[tnum].triPointList[i].y;
				p.z = avz;
				pointlist.push_back(p);
			}
			QList<double> detlist;
			detlist.clear();
			for (int i = 0; i < FACELIST[tnum].triPointList.size(); i++)
			{
				double det = pointlist[i].z - FACELIST[tnum].triPointList[i].z;
				detlist.push_back(det);
			}
			XingBianChangList.push_back(detlist);
		}
		QString trioutpath1 = "file.txt"; 
		
		for (int tnum = mywidget->timeslider->value(); tnum <= mywidget->timeslider->value(); tnum++)
		{
			colorflag = tnum;

			mywidget->setProgressValue(0); 
			myview->rootSeparator->removeAllChildren();
			if (mytree->model()->rowCount() > 0) mytree->model()->removeRows(0, mytree->model()->rowCount(), mytree->currentIndex().parent());
			if (myview->rootSeparator->findChild(Horizons) == -1)
			{
				Horizons = new SoSeparator();
				Horizons->setName(SbName(QString("layer").toLocal8Bit().data()));
				myview->rootSeparator->addChild(Horizons);
			}
			if (myview->rootSeparator->findChild(Tetrahedron) == -1)
			{
				Tetrahedron = new SoSeparator();
				Tetrahedron->setName(SbName(QString("tetrahedron").toLocal8Bit().data()));
				myview->rootSeparator->addChild(Tetrahedron);
			}

			QList<Mdl_Face> FACELIST1;
			for (int i = tnum; i < TCSL + 2; i++)
			{
				mywidget->setProgressValue((float(i - tnum) / float(mywidget->timeslider->value() - tnum)) * 100);
				Mdl_Face FACE;
				
				if (tnum == 0)
				{
					for (int j = 0; j < FACELIST[i].triPointList.size(); j++)
					{
						FACELIST[i].triPointList[j].z += XingBianChangList[tnum][j];
					}
				}
				else
				{
					for (int j = 0; j < FACELIST[i].triPointList.size(); j++)
					{
						FACELIST[i].triPointList[j].z += XingBianChangList[tnum][j];
						FACELIST[i].triPointList[j].z -= XingBianChangList[tnum - 1][j];
					}
				}

				double det_por = res_volume(FACELIST[i], FACELIST[i + 1], Tetrahedron);

				FACE.triPointList.clear();
				for (int j = 0; j < FACELIST[i].triPointList.size(); j++)
				{
					Mdl_Point p;
					p.CF = FACELIST[i].triPointList[j].CF;
					p.pNum = FACELIST[i].triPointList[j].pNum;
					p.x = FACELIST[i].triPointList[j].x;
					p.y = FACELIST[i].triPointList[j].y;
					p.z = FACELIST[i].triPointList[j].z - det_por;
					FACE.triPointList.push_back(p);
				}
				
				if (FACE.triPointList.size() == 0) mywidget->editShow("FACE.triPointList is empty");
				VisualTriList.clear();
				Triangulation(FACE.triPointList, mPoly, trioutpath1);
				ReadTriangulation(VisualTriList, trioutpath1);

				// save triList 
				FACE.triList.clear();
				for (int j = 0; j < VisualTriList.size(); j++) {
					int have = 0;
					for (int k = 0; k < mPoly.polylist.size(); k++)
					{
						if (IsEqualPoint(VisualTriList[j].p1, mPoly.polylist[k]) || IsEqualPoint(VisualTriList[j].p2, mPoly.polylist[k]) || IsEqualPoint(VisualTriList[j].p3, mPoly.polylist[k])) have = 1;
					}
					if (have == 0)
					{
						FACE.triList.push_back(VisualTriList[j]);
					}
				}
				FACELIST1.append(FACE);
			}

			// clear window
			if (tnum != mywidget->timeslider->value())
			{
				myview->rootSeparator->removeAllChildren();
				if (mytree->model()->rowCount() > 0) mytree->model()->removeRows(0, mytree->model()->rowCount(), mytree->currentIndex().parent());
			}

			// output data
			// ... ...
		}
	}
	sssflag = mywidget->timeslider->value();
}

double res_volume(Mdl_Face FACE1, Mdl_Face FACE2, QList<Mdl_Tetrahedra> TetList)
{
	double det_kongxi = 0;

	double AREA = trilist_area(FACE1);

	det_kongxi = det_volume(TetList) / AREA;

	return det_kongxi;
}

double trilist_area(Mdl_Face FACE)
{
	double area = 0;

	for (int i = 0; i < FACE.triList.size(); ++i)
	{
		double a = DIST(FACE.triList[i].p1, FACE.triList[i].p2);
		double b = DIST(FACE.triList[i].p2, FACE.triList[i].p3);
		double c = DIST(FACE.triList[i].p1, FACE.triList[i].p3);
		double s = (a + b + c) / 2.0; 
		area += sqrt(s * (s - a) * (s - b) * (s - c));
	}

	return area;
}

double det_volume(QList<Mdl_Tetrahedra> TetList)
{
	double det_v = 0;

	for (int i = 0; i < TetList.size(); ++i)
	{
		det_v += ((TetList[i].Top_Tet_Porosity - TetList[i].Tet_Porosity)/(1 + TetList[i].Tet_Porosity)) * TetList[i].Get_Tet_Volume();
	}

	return det_v;
}


