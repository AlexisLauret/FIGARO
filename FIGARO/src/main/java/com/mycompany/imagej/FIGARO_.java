package com.mycompany.imagej;

import ij.*;
import ij.plugin.PlugIn;
import ij.gui.GenericDialog;
import ij.io.SaveDialog;
import ij.text.TextWindow;


import java.util.Arrays;
import Jama.*;
import java.lang.Math;

import com.itextpdf.layout.Canvas;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.AreaBreak;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Image;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.AreaBreakType;
import com.itextpdf.layout.property.TextAlignment;
import com.itextpdf.layout.property.UnitValue;
import com.itextpdf.kernel.pdf.PdfWriter;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import com.itextpdf.io.image.ImageDataFactory;
import com.itextpdf.kernel.font.PdfFont;
import com.itextpdf.kernel.font.PdfFontFactory;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.kernel.geom.Rectangle;
import com.itextpdf.kernel.pdf.PdfDocument;
import com.itextpdf.kernel.pdf.PdfPage;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * A plugin to approximate a 3D image with a gaussian thanks to the FIGARO algorithm
 * GRAY8, GRAY16, GRAY32 images.
 *
 * @author Alexis Lauret
 */


public class FIGARO_ implements PlugIn {
	protected ImagePlus image;

	// image property members
	private int width;
	private int height;
	private int depth;

	// plugin parameters
	public double dx;
	public double dy;
	public double dz;
	
	public double lambex;
	public double lambem;
	public double bead;
	public double refr;
	public double NA;

	
	
	//The main function. Automatically launched when one launch the plugin in FIGARO.
	@Override
	public void run(String arg) {
	    
		image=IJ.getImage();
		
		
		if (showDialog()) {

		    try {
				process(image);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		}
	}

	
	//The input dialog
	private boolean showDialog() {
		GenericDialog gd = new GenericDialog("Inputs");
		gd.setSize(500,800);
		
		gd.addMessage("Please enter :\n-The x, y and z resolution in µm\n-The excitation and emission wavelength in nm\n-The numerical aperture\n-The refractive index\n-The bead diameter in µm");
		
		gd.addNumericField("dx", 0.050, 3);
		gd.addNumericField("dy", 0.050, 3);
		gd.addNumericField("dz", 0.100, 3);
		gd.addNumericField("excitation wavelength", 800., 0);
		gd.addNumericField("emission wavelength", 515., 0);
		gd.addNumericField("NA", 1.00, 2);
		gd.addNumericField("n", 1.33, 2);
		gd.addNumericField("bead diameter", 0.2, 1);
		

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		
		dx = gd.getNextNumber();
		dy = gd.getNextNumber();
		dz = gd.getNextNumber();
		lambex=gd.getNextNumber();
		lambem=gd.getNextNumber();
		NA=gd.getNextNumber();
		refr=gd.getNextNumber();
		bead=gd.getNextNumber();
		
		return true;
	}
	
	//This is the processing of the image 
	public void process(ImagePlus image) throws IOException {
		
		//We change the image into data here
		

		ImageStack imagestack=image.getImageStack();

		width = image.getWidth();
		height = image.getHeight();
		depth = image.getNSlices();
		
		double [][][] Y=new double[width][height][depth];

		TextWindow txtwindow=new TextWindow("FIGARO algorithm","Organisms :",1000,700);
		txtwindow.append("CVN, CentraleSupelec, Paris Saclay University, France");
		txtwindow.append("XLIM, CNRS, Limoges, France");
		txtwindow.append("\nInvolved persons :");
		txtwindow.append("Emilie Chouzenoux, Tsz Kit Lau, Claire Lefort, Jean-Christophe Pesquet");
		txtwindow.append("Developer : Alexis Lauret");
		txtwindow.append("\nRelated article : \nT. K. Lau, E. Chouzenoux, C. Lefort and J.-C. Pesquet.  \nOptimal Multivariate Gaussian Fitting for PSF Modeling in Two-Photon Microscopy. \nIn Proceedings of the IEEE International Symposium on Biomedical Imaging (ISBI 2018), \nWashington DC, 4 - 7 April 2018");
		txtwindow.append("\nPlease quote this article when using this plugin.");
		txtwindow.append("\n\nsize data ="+width+"x"+height+"x"+depth);
		
		
		double Ymax=0;
		double Ymin=0;
		int[] imax={0,0,0};
		for (int i=0; i < width; i++) {
			for (int j=0; j < height; j++) {
				for (int k=0; k < depth; k++) {
					Y[i][j][k]=imagestack.getVoxel(i, j, k);
					if (Y[i][j][k]>Ymax) {
						Ymax=Y[i][j][k];
						imax[0]= i;
						imax[1]= j;
						imax[2]= k;
					}
					if (Y[i][j][k]<Ymin) {
						Ymin=Y[i][j][k];
					}
				}
			}
		}
		
		
		double [][][] Yfull=new double[width][height][depth];
		for (int i=0; i < width; i++) {
			for (int j=0; j < height; j++) {
				for (int k=0; k < depth; k++) {
					Yfull[i][j][k]=imagestack.getVoxel(i, j, k);
				}
			}
		}
		
		
	
		for (int i=0; i < width; i++) {
			for (int j=0; j < height; j++) {
				for (int k=0; k < depth; k++) {
					Y[i][j][k]=Y[i][j][k]/Ymax;
				}
			}
		}
		
		//These are the arrays of the coordinates
		double[] x1= new double[width];
		double[] x2= new double[height];
		double[] x3= new double[depth];
		
		for (int i=0; i < width; i++) {
			x1[i]=dx*i;
		}
		
		for (int j=0; j < height; j++) {
			x2[j]=dy*j;
		}
		
		for (int k=0; k < depth; k++) {
			x3[k]=dz*k;
		}
		
		//estimation={estimA, estimB, estimmu, estimC, estimP}
		
		int[] estimation= {1,1,1,1,1};
		double d=dx*dy*dz;
		double dd =  1./d;
		double A0=Ymin;
		double B0=1.;
		double[] Mu0=new double[] {x1[imax[0]],x2[imax[1]],x3[imax[2]]};
		double[][] C0=new double[][] {{dx,0.,0.},{0.,dy,0.},{0.,0.,dz}};
		double[][][] P0=Y;
		
		//warm start strategy for C
		
		//The text window acts like a console for ImageJ
		txtwindow.append("WARM START");
		Object[] O=Gauss3dFIT(Y,0.01/d,d,dd,x1,x2,x3,A0,B0,P0,Mu0,C0,0,estimation,10,txtwindow);
		C0=(double[][]) O[4];
		
		
		//search for optimal lambda
		double sigmaN=0.;
		for (int k=0; k<depth; k++) {
			double[][] Yw=new double[width][height];
			for (int i=0; i < width; i++) {
				for (int j=0; j < height; j++) {
					Yw[i][j]=Y[i][j][k];
				}
			}
			Yw=wiener2(Yw);
			double[][] Ymean= new double[width][height];
			for (int i=0; i < width; i++) {
				for (int j=0; j < height; j++) {
					Ymean[i][j]=Y[i][j][k]-Yw[i][j];
				}
			}
			sigmaN=sigmaN+std(Ymean);
		}
		sigmaN=sigmaN/depth;
		
		txtwindow.append("sigmaN="+sigmaN);
		txtwindow.append("------------------------");
		txtwindow.append("");
		
		double lambda = 1./d ;  
		boolean loopchi2  = true;
		
		//start chi2 loop
		double llambdad_min = -3. - Math.log10(d);
		double llambdad_max = 3.  - Math.log10(d);
		

		double tolambda = 1e-2;
		double tolchi2 = 1e-2;
		int loopct = 0;
		int nloop = 20;
		
		while(loopchi2) {
			loopct=loopct+1;
			double llambda = (llambdad_min + llambdad_max)/2;
		    lambda = Math.pow(10,llambda);
		    
		    txtwindow.append("TEST: lambda*d ="+lambda*d);
		    
		    O = Gauss3dFIT(Y,lambda,d,dd,x1,x2,x3,A0,B0,P0,Mu0,C0,0,estimation,10000,txtwindow);
		    double[] LS=(double[]) O[6];
		    double Chi2fit = 2.*LS[LS.length-1]/(Math.pow(sigmaN,2));
		    txtwindow.append("Chi2fit/N ="+Chi2fit/(width*height*depth));
		    		
		    if (Chi2fit/(width*height*depth) > 1) {
		        llambdad_max = llambda;
		    }
		    else {
		        llambdad_min = llambda;
		    }
		    
		    if ((llambdad_max-llambdad_min)/2 < tolambda) {
		    	
		    	txtwindow.append("break chi2 loop (tollambda)");
		    	
		        loopchi2 = false;
		    }
		    
		    if (loopct > nloop) {
		    	
		    	txtwindow.append("break chi2 loop (nloop)");
		    	
		        loopchi2 = false;
		    }
		    
		    if (Math.abs(Chi2fit/(width*height*depth) - 1) < tolchi2) {
		    	
		    	txtwindow.append("break chi2 loop (tolchi2)");
		    	
		        loopchi2 = false;
		    }
		    
		    txtwindow.append("------------------------");
		    txtwindow.append("");
		}
		
		//start optimization
		
		
		O=Gauss3dFIT(Y,lambda,d,dd,x1,x2,x3,A0,B0,P0,Mu0,C0,0,estimation,10000, txtwindow);
		double Ae=(double) O[0];
		double Be=(double) O[1];
		double[][][] Pe=(double[][][]) O[2];
		double[] LSlist=(double[]) O[6];
		double LS=LSlist[LSlist.length-1];
		double[] Mue=(double[]) O[3];
		double[][] Ce=(double[][]) O[4];
		double[][][] Ge=(double[][][]) GaussianDistrib3D(x1,x2,x3,Mue,Ce)[0];
		for (int i=0; i < width; i++) {
			for (int j=0; j < height; j++) {
				for (int k=0; k < depth; k++) {
					Ge[i][j][k]=Ae+Be*Ge[i][j][k];
				}
			}
		}
		
		
		//Position of the middle
		int[] MuInt=new int[] {0,0,0};
		double var=Math.abs(Mue[0]-x1[0]);
		for (int i=0; i < width; i++) {
			if (Math.abs(Mue[0]-x1[i])<var) {
				var=Math.abs(Mue[0]-x1[i]);
				MuInt[0]=i;
			}
		}
		var=Math.abs(Mue[1]-x2[0]);
		for (int j=0; j < height; j++) {
			if (Math.abs(Mue[1]-x2[j])<var) {
				var=Math.abs(Mue[1]-x2[j]);
				MuInt[1]=j;
			}
		}
		var=Math.abs(Mue[2]-x3[0]);
		for (int k=0; k < depth; k++) {
			if (Math.abs(Mue[2]-x3[k])<var) {
				var=Math.abs(Mue[2]-x3[k]);
				MuInt[2]=k;
			}
		}
		
		
		Matrix CeM=new Matrix(Ce);
		Matrix UeM=CeM.svd().getU();
		
		//Radius
		double[] Radius=new double[3];
		Radius=CeM.times(1./(2.*Math.log(2))).inverse().eig().getRealEigenvalues();
		for (int i=0; i<3; i++) {
			Radius[i]=Math.sqrt(Radius[i]);
		}
		int mainAxisPlace=0;
		if (Radius[1]>Radius[0]) {
			mainAxisPlace=1;
		}
		if (Radius[2]>Radius[mainAxisPlace]) {
			mainAxisPlace=2;
		}
		double[][] Ue=UeM.getArray();
		double[] Ue3=new double[] {Ue[0][2],Ue[1][2],Ue[2][2]};
		
		//Phi and Theta
		
		double phi = 180./Math.PI*atan2d(Ue3[1],Ue3[0]);
		double theta = 180./Math.PI*atan2d( Math.sqrt(Math.pow(Ue3[0],2) + Math.pow(Ue3[1],2)), Ue3[2]);


		txtwindow.append("");
		txtwindow.append("Radius="+Arrays.toString(Radius));
		txtwindow.append("Phi="+Double.toString(phi));
		txtwindow.append("Theta="+Double.toString(theta));
		txtwindow.append("Position of the middle: i="+MuInt[0]+" j="+MuInt[1]+" k="+MuInt[2]);
		txtwindow.append("Residual least square="+LS);
		

		
		
		ImagePlus ImgGe=IJ.createImage("Ge", width, height, depth, 32);
		ImagePlus ImgPe=IJ.createImage("Pe", width, height, depth, 32);
		
		ImageStack ImgStackGe=ImgGe.getImageStack();
		ImageStack ImgStackPe=ImgPe.getImageStack();
		
		for (int i=0; i < width; i++) {
			for (int j=0; j < height; j++) {
				for (int k=0; k < depth; k++) {
					ImgStackPe.setVoxel(i, j, k, Pe[i][j][k]);
					ImgStackGe.setVoxel(i, j, k, Ge[i][j][k]);
				}
			}
		}
		

		//We create the images of the profiles. The format is RGB to point out the middle in red.
		
		ImagePlus originalXY=IJ.createImage("originalXY", width, height, 1, 24);
		ImagePlus originalXZ=IJ.createImage("originalXZ", width, depth, 1, 24);
		ImagePlus originalYZ=IJ.createImage("originalYZ", height, depth, 1, 24);
		ImagePlus PeXY=IJ.createImage("PeXY with ellipse", width, height, 1, 24);
		ImagePlus PeXZ=IJ.createImage("PeXZ with ellipse", width, depth, 1, 24);
		ImagePlus PeYZ=IJ.createImage("PeYZ with ellipse", height, depth, 1, 24);
		ImagePlus PeXY2=IJ.createImage("PeXY", width, height, 1, 24);
		ImagePlus PeXZ2=IJ.createImage("PeXZ", width, depth, 1, 24);
		ImagePlus PeYZ2=IJ.createImage("PeYZ", height, depth, 1, 24);
		
		
		ImageStack ImgStackOriginalXY=originalXY.getImageStack();
		ImageStack ImgStackOriginalXZ=originalXZ.getImageStack();
		ImageStack ImgStackOriginalYZ=originalYZ.getImageStack();
		ImageStack ImgStackPeXY=PeXY.getImageStack();
		ImageStack ImgStackPeXZ=PeXZ.getImageStack();
		ImageStack ImgStackPeYZ=PeYZ.getImageStack();
		ImageStack ImgStackPeXY2=PeXY2.getImageStack();
		ImageStack ImgStackPeXZ2=PeXZ2.getImageStack();
		ImageStack ImgStackPeYZ2=PeYZ2.getImageStack();
		
		for (int i=0; i < width; i++) {
			for (int j=0; j < height; j++) {
				for (int k=0; k < depth; k++) {
					ImgStackOriginalXY.setVoxel(i, j, 0, Math.floor(Y[i][j][MuInt[2]]*255)*65793);
					if (Y[i][j][MuInt[2]]>=1) {
						ImgStackOriginalXY.setVoxel(i, j, 0, 16777215);
					}
					ImgStackOriginalXZ.setVoxel(i, k, 0, Math.floor(Y[i][MuInt[1]][k]*255)*65793);
					if (Y[i][MuInt[1]][k]>=1) {
						ImgStackOriginalXZ.setVoxel(i, k, 0, 16777215);
					}
					ImgStackOriginalYZ.setVoxel(j, k, 0, Math.floor(Y[MuInt[0]][j][k]*255)*65793);
					if (Y[MuInt[0]][j][k]>=1) {
						ImgStackOriginalYZ.setVoxel(j, k, 0, 16777215);
					}
					ImgStackPeXY.setVoxel(i, j, 0, Math.floor(Pe[i][j][MuInt[2]]*255)*65793);
					ImgStackPeXY2.setVoxel(i, j, 0, Math.floor(Pe[i][j][MuInt[2]]*255)*65793);
					
					if (Pe[i][j][MuInt[2]]>=1) {
						ImgStackPeXY.setVoxel(i, j, 0, 16777215);
						ImgStackPeXY2.setVoxel(i, j, 0, 16777215);
					}
					ImgStackPeXZ.setVoxel(i, k, 0, Math.floor(Pe[i][MuInt[1]][k]*255)*65793);
					ImgStackPeXZ2.setVoxel(i, k, 0, Math.floor(Pe[i][MuInt[1]][k]*255)*65793);
					if (Pe[i][MuInt[1]][k]>=1) {
						ImgStackPeXZ.setVoxel(i, k, 0, 16777215);
						ImgStackPeXZ2.setVoxel(i, k, 0, 16777215);
					}
					ImgStackPeYZ.setVoxel(j, k, 0, Math.floor(Pe[MuInt[0]][j][k]*255)*65793);
					ImgStackPeYZ2.setVoxel(j, k, 0, Math.floor(Pe[MuInt[0]][j][k]*255)*65793);
					if (Pe[MuInt[0]][j][k]>=1) {
						ImgStackPeYZ.setVoxel(j, k, 0, 16777215);
						ImgStackPeYZ2.setVoxel(j, k, 0, 16777215);
					}
				}
			}
		}
		
		

		ImgStackPeXY.setVoxel(MuInt[0], MuInt[1], 0, 16711680);
		ImgStackPeXZ.setVoxel(MuInt[0], MuInt[2], 0, 16711680);
		ImgStackPeYZ.setVoxel(MuInt[1], MuInt[2], 0, 16711680);
		
		
		//We change to BufferedImage to add an ellipse on some profiles and to prepare the change to the iText class of Image, to build the PDF
	    
	    BufferedImage awtImOriginalXY =toBufferedImage((java.awt.Image) awtScaledImage(originalXY)[0]);
	    BufferedImage awtImOriginalXZ = toBufferedImage((java.awt.Image) awtScaledImage(originalXZ)[0]);
	    BufferedImage awtImOriginalYZ = toBufferedImage((java.awt.Image) awtScaledImage(originalYZ)[0]);
	    BufferedImage awtImPeXY = toBufferedImage((java.awt.Image) awtScaledImage(PeXY)[0]);
	    BufferedImage awtImPeXZ = toBufferedImage((java.awt.Image) awtScaledImage(PeXZ)[0]);
	    BufferedImage awtImPeYZ = toBufferedImage((java.awt.Image) awtScaledImage(PeYZ)[0]);
	    
	    
	    float zoomXY=(float) awtScaledImage(PeXY)[1];
	    float zoomXZ=(float) awtScaledImage(PeXZ)[1];
	    float zoomYZ=(float) awtScaledImage(PeYZ)[1];
	    
	    
	    //Here we paint the ellipses
	    
	    Graphics2D gXY=(Graphics2D) awtImPeXY.getGraphics();
	    Graphics2D gXZ=(Graphics2D) awtImPeXZ.getGraphics();
	    Graphics2D gYZ=(Graphics2D) awtImPeYZ.getGraphics();
	    gXY.setColor(Color.BLUE);
	    gXZ.setColor(Color.BLUE);
	    gYZ.setColor(Color.BLUE);
	    int XYWidth=awtImPeXY.getWidth();
	    int XZWidth=awtImPeXZ.getWidth();
	    int YZWidth=awtImPeYZ.getWidth();
	    int XYHeight=awtImPeXY.getHeight();
	    int XZHeight=awtImPeXZ.getHeight();
	    int YZHeight=awtImPeYZ.getHeight();
	    
	    
	    double XYaEllipse=ellipseLengths(Ce[0][0],Ce[0][1],Ce[1][1],-2*Math.log(2))[0]*zoomXY;
	    double XZaEllipse=ellipseLengths(Ce[0][0],Ce[0][2],Ce[2][2],-2*Math.log(2))[0]*zoomXZ;
	    double YZaEllipse=ellipseLengths(Ce[1][1],Ce[1][2],Ce[2][2],-2*Math.log(2))[0]*zoomYZ;
	    double XYbEllipse=ellipseLengths(Ce[0][0],Ce[0][1],Ce[1][1],-2*Math.log(2))[1]*zoomXY;
	    double XZbEllipse=ellipseLengths(Ce[0][0],Ce[0][2],Ce[2][2],-2*Math.log(2))[1]*zoomXZ;
	    double YZbEllipse=ellipseLengths(Ce[1][1],Ce[1][2],Ce[2][2],-2*Math.log(2))[1]*zoomYZ;
	    double XYthetaEllipse=ellipseAngle(Ce[0][0],Ce[0][1],Ce[1][1]);
	    double XZthetaEllipse=ellipseAngle(Ce[0][0],Ce[0][2],Ce[2][2]);
	    double YZthetaEllipse=ellipseAngle(Ce[1][1],Ce[1][2],Ce[2][2]);
	    
	    
	    int x1XY=(int)  Math.round((((double) MuInt[0])/width*XYWidth + XYaEllipse*Math.cos(XYthetaEllipse)/dx))+1;
	    int x1XZ=(int)  Math.round((((double) MuInt[0])/width*XZWidth + XZaEllipse*Math.cos(XZthetaEllipse)/dx))+1;
	    int x1YZ=(int)  Math.round((((double) MuInt[1])/height*YZWidth + YZaEllipse*Math.cos(YZthetaEllipse)/dy))+1;
	    int y1XY=(int)  Math.round((((double) MuInt[1])/height*XYHeight + XYaEllipse*Math.sin(XYthetaEllipse)/dy))+1;
	    int y1XZ=(int)  Math.round((((double) MuInt[2])/depth*XZHeight + XZaEllipse*Math.sin(XZthetaEllipse)/dz))+1;
	    int y1YZ=(int)  Math.round((((double) MuInt[2])/depth*YZHeight + YZaEllipse*Math.sin(YZthetaEllipse)/dz))+1;
	    
	    for (int i=1;i<=360;i++) {
	    	int x2XY=(int)  Math.round((((double) MuInt[0])/width*XYWidth + (XYaEllipse*Math.cos(XYthetaEllipse)*Math.cos((double) i/360*2*Math.PI)-XYbEllipse*Math.sin(XYthetaEllipse)*Math.sin((double) i/360*2*Math.PI))/dx))+1;
		    int x2XZ=(int)  Math.round((((double) MuInt[0])/width*XZWidth + (XZaEllipse*Math.cos(XZthetaEllipse)*Math.cos((double) i/360*2*Math.PI)-XZbEllipse*Math.sin(XZthetaEllipse)*Math.sin((double) i/360*2*Math.PI))/dx))+1;
		    int x2YZ=(int)  Math.round((((double) MuInt[1])/height*YZWidth + (YZaEllipse*Math.cos(YZthetaEllipse)*Math.cos((double) i/360*2*Math.PI)-YZbEllipse*Math.sin(YZthetaEllipse)*Math.sin((double) i/360*2*Math.PI))/dy))+1;
		    int y2XY=(int)  Math.round((((double) MuInt[1])/height*XYHeight + (XYaEllipse*Math.sin(XYthetaEllipse)*Math.cos((double) i/360*2*Math.PI)+XYbEllipse*Math.cos(XYthetaEllipse)*Math.sin((double) i/360*2*Math.PI))/dy))+1;
		    int y2XZ=(int)  Math.round((((double) MuInt[2])/depth*XZHeight + (XZaEllipse*Math.sin(XZthetaEllipse)*Math.cos((double) i/360*2*Math.PI)+XZbEllipse*Math.cos(XZthetaEllipse)*Math.sin((double) i/360*2*Math.PI))/dz))+1;
		    int y2YZ=(int)  Math.round((((double) MuInt[2])/depth*YZHeight + (YZaEllipse*Math.sin(YZthetaEllipse)*Math.cos((double) i/360*2*Math.PI)+YZbEllipse*Math.cos(YZthetaEllipse)*Math.sin((double) i/360*2*Math.PI))/dz))+1;
	    
		    gXY.drawLine(x1XY,y1XY,x2XY,y2XY);
		    gXZ.drawLine(x1XZ,y1XZ,x2XZ,y2XZ);
		    gYZ.drawLine(x1YZ,y1YZ,x2YZ,y2YZ);
		    
		    x1XY=x2XY;
		    x1XZ=x2XZ;
		    x1YZ=x2YZ;
		    y1XY=y2XY;
		    y1XZ=y2XZ;
		    y1YZ=y2YZ;
		    
	    }
	    
	    
	    //Change to the iText Image class
	    
	    Image ImOriginalXY = new Image(ImageDataFactory.create(awtImOriginalXY,null));
	    Image ImOriginalXZ = new Image(ImageDataFactory.create(awtImOriginalXZ,null));
	    Image ImOriginalYZ = new Image(ImageDataFactory.create(awtImOriginalYZ,null));
	    Image ImPeXY = new Image(ImageDataFactory.create(awtImPeXY,null));
	    Image ImPeXZ = new Image(ImageDataFactory.create(awtImPeXZ,null));
	    Image ImPeYZ = new Image(ImageDataFactory.create(awtImPeYZ,null));
	    
	    SaveDialog sd = new SaveDialog("Save the report to...", "FIGARO report for " + image.getTitle(), ".pdf");
		String path = sd.getDirectory() + sd.getFileName();
		
		try {
			
			//Create the PDF doc
			PdfWriter writer=new PdfWriter(path);
			PdfDocument pdf = new PdfDocument(writer);
		    Document document = new Document(pdf);
		    
		    PdfFont normalFont = PdfFontFactory.createFont("C:\\Windows\\Fonts\\arial.ttf", "Identity-H", true);
		    document.setFont(normalFont);
		    
		    
		    //Three column to put the profiles into
		    float offSet = 36;
		    float columnWidth = (PageSize.A4.getWidth() - offSet * 2) / 3;
		    float columnHeight = PageSize.A4.getHeight()*2f/3; 
		    Rectangle[] columns = {
		        new Rectangle(offSet, 28, columnWidth, columnHeight),
		        new Rectangle(offSet + columnWidth , 28, columnWidth, columnHeight),
		        new Rectangle(offSet + 2*columnWidth, 28, columnWidth, columnHeight)};
		    
		    
		    document.add(bigTitle("FIGARO report for "+ image.getTitle()));
		    
		    PdfPage page = pdf.getFirstPage();
		    
		    PdfCanvas pdfCanvas = new PdfCanvas(page);
		    pdfCanvas.stroke();
		    Canvas canvas1 = new Canvas(pdfCanvas, pdf, columns[0]);
		    
		    Canvas canvas2 = new Canvas(pdfCanvas, pdf, columns[1]);
		    
		    Canvas canvas3 = new Canvas(pdfCanvas, pdf, columns[2]);
		    
		    
		    Paragraph p0=new Paragraph("Visual illustration of raw data and restored data");
		    p0.setFontSize(13f);
		    p0.setBold();
		    p0.setPaddingBottom(10f);
		    
		    Paragraph p1=new Paragraph("Profile XY");
		    
		    Paragraph p2=new Paragraph("");
		    p2.setPaddingTop(10f);
		    p2.setPaddingBottom(10f);
		    
		    Paragraph p3=new Paragraph("Profile XZ");
		    
		    Paragraph p4=new Paragraph("");
		    p4.setPaddingTop(10f);
		    p4.setPaddingBottom(10f);
		    
		    Paragraph p5=new Paragraph("Profile YZ");
		    
		    Paragraph p6=new Paragraph("");
		    p6.setPaddingTop(10f);
		    p6.setPaddingBottom(10f);
		    
		    
	        p2.add(ImOriginalXY);
	        p2.add("    ");
	        p2.add(ImPeXY);
	        p4.add(ImOriginalXZ);
	        p4.add("    ");
	        p4.add(ImPeXZ);
	        p6.add(ImOriginalYZ);
	        p6.add("    ");
	    	p6.add(ImPeYZ);
	    	
	    	
	    	document.add(p0);
	    	canvas1.add(p1);
	    	canvas1.add(p2);
	    	canvas2.add(p3);
	    	canvas2.add(p4);
	    	canvas3.add(p5);
	    	canvas3.add(p6);
	    	
	    	
	    	document.add(new AreaBreak(AreaBreakType.NEXT_PAGE));
	        
	        
	        Paragraph p7=new Paragraph("Quantitative values");
		    p7.setFontSize(13f);
		    p7.setBold();
		    p7.setPaddingBottom(10f);
		    
		    document.add(p7);
		    
		    Paragraph p8=new Paragraph("Approximation: y(u)=a+b*p(u)+w(u), y being the data, a, b and p estimated parameters, w the noise and u the position.");
		    p8.setFontSize(11f);
		    p8.setPaddingBottom(10f);
		    
		    document.add(p8);
		    
		    //Results table
		    Table table=new Table(new float[] {6f,1f,1f,1f,1f,1f,1f});
		    Table table2=new Table(new float[] {6f,1f,1f,1f,1f,1f,1f});
		    table.setWidth(UnitValue.createPercentValue(100));
		    table2.setWidth(UnitValue.createPercentValue(100));
		    
		    table.addCell(new Cell(1,1).add(new Paragraph("Excitation wavelength (nm)")));
		    table.addCell(new Cell(1,6).add(new Paragraph(Double.toString(lambex))));
		    table.addCell(new Cell(1,1).add(new Paragraph("Emission wavelength (nm)")));
		    table.addCell(new Cell(1,6).add(new Paragraph(Double.toString(lambem))));
		    table.addCell(new Cell(1,1).add(new Paragraph("Numerical Aperture")));
		    table.addCell(new Cell(1,6).add(new Paragraph(Double.toString(NA))));
		    table.addCell(new Cell(1,1).add(new Paragraph("Bead diameter (\u03BCm)")));
		    table.addCell(new Cell(1,6).add(new Paragraph(Double.toString(bead))));
		    table.addCell(new Cell(1,1).add(new Paragraph("Theoretical deep axis FWHM (\u03BCm) [Diaspro 2002]")));
		    
		    
    	    double diaspro1=0.7*lambem/(1000*NA);
    	    double diaspro2=2.3*refr*lambem/(1000*Math.pow(NA,2));
    		
		    table.addCell(new Cell(1,6).add(new Paragraph(Double.toString((double)((int)(diaspro2*1000))/1000))));
		    table.addCell(new Cell(1,1).add(new Paragraph("Theoretical lateral FWHM (\u03BCm) [Diaspro 2002]")));
		    
		    table.addCell(new Cell(1,6).add(new Paragraph(Double.toString((double)((int)(diaspro1*1000))/1000))));
		    table.addCell(new Cell(1,1).add(new Paragraph("Theorical lateral FWHM (\u03BCm) [Cox 2004]")));
		   
		    double FWHM=0.5*lambex/(1000*Math.sqrt(2)*NA);
		    table.addCell(new Cell(1,6).add(new Paragraph(Double.toString((double)((int)(FWHM*1000))/1000))));
		    
		    
		    
		    table2.addCell(new Cell(2,1).add(new Paragraph("Position of bead center")));
		    table2.addCell(new Cell(1,3).add(new Paragraph("Voxel position")));
		    table2.addCell(new Cell(1,3).add(new Paragraph("Distance position (\u03BCm)")));
		    table2.addCell(new Cell(1,1).add(new Paragraph(Integer.toString( MuInt[0]))));
		    table2.addCell(new Cell(1,1).add(new Paragraph(Integer.toString( MuInt[1]))));
		    table2.addCell(new Cell(1,1).add(new Paragraph(Integer.toString( MuInt[2]))));
		    table2.addCell(new Cell(1,1).add(new Paragraph(Double.toString( (double)((int)(MuInt[0]*dx*1000))/1000))));
		    table2.addCell(new Cell(1,1).add(new Paragraph(Double.toString( (double)((int)(MuInt[1]*dy*1000))/1000))));
		    table2.addCell(new Cell(1,1).add(new Paragraph(Double.toString( (double)((int)(MuInt[2]*dz*1000))/1000))));
		    table2.addCell(new Cell(1,1).add(new Paragraph("Background parameter a")));
		    table2.addCell(new Cell(1,6).add(new Paragraph(Double.toString((double)((int)(Ae*1000))/1000))));
		    table2.addCell(new Cell(2,1).add(new Paragraph("Extremal values")));
		    table2.addCell(new Cell(1,3).add(new Paragraph("Min")));
		    table2.addCell(new Cell(1,3).add(new Paragraph("Max")));
		    table2.addCell(new Cell(1,3).add(new Paragraph(Double.toString(Ymin))));
		    table2.addCell(new Cell(1,3).add(new Paragraph(Double.toString(Ymax))));
		    table2.addCell(new Cell(1,1).add(new Paragraph("Scale parameter b")));
		    table2.addCell(new Cell(1,6).add(new Paragraph(Double.toString((double)((int)(Be*1000))/1000))));
		    table2.addCell(new Cell(2,1).add(new Paragraph("Estimated FWHM (\u03BCm)")));
		    table2.addCell(new Cell(1,2).add(new Paragraph("Deep axis")));
		    table2.addCell(new Cell(1,2).add(new Paragraph("Horizontal axis")));
		    table2.addCell(new Cell(1,2).add(new Paragraph("Vertical axis")));
		    table2.addCell(new Cell(1,2).add(new Paragraph(Double.toString((double)((int)(2*Radius[mainAxisPlace]*1000))/1000))));
		    table2.addCell(new Cell(1,2).add(new Paragraph(Double.toString((double)((int)(2*Radius[(mainAxisPlace+1)%3]*1000))/1000))));
		    table2.addCell(new Cell(1,2).add(new Paragraph(Double.toString((double)((int)(2*Radius[(mainAxisPlace+2)%3]*1000))/1000))));
		    table2.addCell(new Cell(2,1).add(new Paragraph("Inclination angles in degrees")));
		    table2.addCell(new Cell(1,3).add(new Paragraph("Phi")));
		    table2.addCell(new Cell(1,3).add(new Paragraph("Theta")));
		    table2.addCell(new Cell(1,3).add(new Paragraph(Double.toString((double)((int)(phi*1000))/1000))));
		    table2.addCell(new Cell(1,3).add(new Paragraph(Double.toString((double)((int)(theta*1000))/1000))));
		    
		    document.add(table);
		    
		    Paragraph p9=new Paragraph("\nEstimated values resulting from FIGARO");
		    p9.setFontSize(11f);
		    p9.setPaddingBottom(10f);
		    
		    document.add(p9);
		    
		    document.add(table2);
		    
		    document.add(new AreaBreak(AreaBreakType.NEXT_PAGE));
		    
		    Paragraph p10=new Paragraph("Correlation matrix of the gaussian approximation");
		    p10.setFontSize(13f);
		    p10.setBold();
		    p10.setPaddingBottom(10f);
		    
		    document.add(p10);
		    
		    Paragraph p11=new Paragraph("\n( "+(double)((int)(Ce[0][0]*1000))/1000+"    "+(double)((int)(Ce[0][1]*1000))/1000+"    "+(double)((int)(Ce[0][2]*1000))/1000+"\n\n"+" "+(double)((int)(Ce[1][0]*1000))/1000+"    "+(double)((int)(Ce[1][1]*1000))/1000+"    "+(double)((int)(Ce[1][2]*1000))/1000+"\n\n"+" "+(double)((int)(Ce[2][0]*1000))/1000+"    "+(double)((int)(Ce[2][1]*1000))/1000+"    "+(double)((int)(Ce[2][2]*1000))/1000+" )");
		    document.add(p11);
		    
		    document.close();
		    
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			}
		
		ImagePlus NewPeXY=new ImagePlus("PeXY with ellipse",awtImPeXY);
		ImagePlus NewPeXZ=new ImagePlus("PeXZ with ellipse",awtImPeXZ);
		ImagePlus NewPeYZ=new ImagePlus("PeYZ with ellipse",awtImPeYZ);
	
		//show the images
		ImgGe.show();
		ImgPe.show();
		originalXY.show();
		originalXZ.show();
		originalYZ.show();
		PeXY2.show();
		PeXZ2.show();
		PeYZ2.show();
		NewPeXY.show();
		NewPeXZ.show();
		NewPeYZ.show();
		
		
		//Save the images
		SaveDialog sdpe = new SaveDialog("Save Pe", "Approximation Pe for " + image.getTitle(), ".tif");
		String pathPe = sdpe.getDirectory() + sdpe.getFileName();
		IJ.save(ImgPe,pathPe);
		
		SaveDialog sdge = new SaveDialog("Save Ge", "Gaussian Ge for " + image.getTitle(), ".tif");
		String pathGe = sdge.getDirectory() + sdge.getFileName();
		IJ.save(ImgGe,pathGe);

		SaveDialog sdOrXY = new SaveDialog("Save original XY profile", "XY profile for " + image.getTitle(), ".tif");
		String pathOrXY = sdOrXY.getDirectory() + sdOrXY.getFileName();
		IJ.save(originalXY,pathOrXY);
		
		SaveDialog sdOrXZ = new SaveDialog("Save original XZ profile", "XZ profile for " + image.getTitle(), ".tif");
		String pathOrXZ = sdOrXZ.getDirectory() + sdOrXZ.getFileName();
		IJ.save(originalXZ,pathOrXZ);
		
		SaveDialog sdOrYZ = new SaveDialog("Save original YZ profile", "YZ profile for " + image.getTitle(), ".tif");
		String pathOrYZ = sdOrYZ.getDirectory() + sdOrYZ.getFileName();
		IJ.save(originalYZ,pathOrYZ);
		
		SaveDialog sdPeXY2 = new SaveDialog("Save XY Pe profile", "XY Pe profile for " + image.getTitle(), ".tif");
		String pathPeXY2 = sdPeXY2.getDirectory() + sdPeXY2.getFileName();
		IJ.save(PeXY2,pathPeXY2);
		
		SaveDialog sdPeXZ2 = new SaveDialog("Save XZ Pe profile", "XZ Pe profile for " + image.getTitle(), ".tif");
		String pathPeXZ2 = sdPeXZ2.getDirectory() + sdPeXZ2.getFileName();
		IJ.save(PeXZ2,pathPeXZ2);
		
		SaveDialog sdPeYZ2 = new SaveDialog("Save YZ Pe profile", "YZ Pe profile for " + image.getTitle(), ".tif");
		String pathPeYZ2 = sdPeYZ2.getDirectory() + sdPeYZ2.getFileName();
		IJ.save(PeYZ2,pathPeYZ2);
		
		SaveDialog sdPeXY = new SaveDialog("Save XY Pe profile with ellipse", "XY Pe profile with ellipse for " + image.getTitle(), ".tif");
		String pathPeXY = sdPeXY.getDirectory() + sdPeXY.getFileName();
		IJ.save(NewPeXY,pathPeXY);
		
		SaveDialog sdPeXZ = new SaveDialog("Save XZ Pe profile with ellipse", "XZ Pe profile with ellipse for " + image.getTitle(), ".tif");
		String pathPeXZ = sdPeXZ.getDirectory() + sdPeXZ.getFileName();
		IJ.save(NewPeXZ,pathPeXZ);
		
		SaveDialog sdPeYZ = new SaveDialog("Save YZ Pe profile with ellipse", "YZ Pe profile with ellipse for " + image.getTitle(), ".tif");
		String pathPeYZ = sdPeYZ.getDirectory() + sdPeYZ.getFileName();
		IJ.save(NewPeYZ,pathPeYZ);
		
		
		
	}
	
	
	//The optimisation algorithm
	public Object[] Gauss3dFIT(double[][][] Y, double lamb, double d, double dd, double[] x1, double[] x2, double[] x3, double A0, double B0, double[][][] P0, double[] Mu0, double[][] C0, int display, int[] estimation, int NbIt, TextWindow txtwindow) {
		
		txtwindow.append("----- START GAUSSIAN FITTING ALGORITHM -----");
		
		double gamB = 1.;
		double gamA = 1.;
		double gamC = 1.;
		double gamM = 1.;
		double gamP = 1.;

		double A   = A0;
		double B   = B0;
		double[][] C   = C0;
		double[] Mu  = Mu0;
		double[][][] P   = P0;

		double Abar = A;
		double Bbar = B;
		double[]Mubar = Mu;
		double[][] Cbar = C;
		double[][][] Pbar = P;
		
		double[] F=new double[1];
		double[] LS=new double[1];
		double[] KL=new double[1];
		double[] Asave=new double[1];
		double[] Bsave=new double[1];
		double[][] Musave=new double[1][3];
		
		double[] L=GaussianFitCriterion(Y, lamb, x1, x2, x3, A, B, P, Mu, C, d);
		F[0]=L[0];
		LS[0]=L[1];
		KL[0]=L[2];
		Asave[0]=A0;
		Bsave[0]=B0;
		Musave[0]=Mu0;
		
		double stop = 1e-10;
		double stopp = 1e-5;
		double nu  = 0;
		double[][][] Gn =(double[][][]) GaussianDistrib3D(x1, x2, x3 , Mu, C)[0];
		for (int i=0; i < width; i++) {
			for (int j=0; j < height; j++) {
				for (int k=0; k < depth; k++) {
					Gn[i][j][k]=A+B*Gn[i][j][k];
				}
			}
		}
		for (int n=1; n<NbIt+1; n++) {
			
			//update P
			if (estimation[4]==1) {
				if (lamb>0) {
	                Object[] FP = ProxFP(Pbar, nu, gamP, lamb, d, dd, C, x1, x2, x3, Mu, A, B, Y);
	                P=(double[][][]) FP[0];
	                nu=(double) FP[1];
				}
	            else {
	                P=ProxFPlambda0(Pbar, gamP, dd, Y, A, B);
	                }
				Pbar=P;
	            }
			
			
			//update A
			if (estimation[0]==1){
		        A = ProxFA(Abar,gamA,B,P,Y);
		        Abar = A;
			}
			
			//update B
			if(estimation[1]==1) {
		        B = ProxFB(Bbar,gamB,A,P,Y);
		        Bbar = B;
			}
			
			//update Mu
			if (estimation[2]==1) {
		         
		        Mu = ProxFmu(Mubar, gamM, x1, x2, x3, C, lamb, P, d);
		        Mubar = Mu;  
			}
			
			//update C
			if(estimation[3]==1) {
		        C = ProxFC(Cbar, x1, x2, x3, gamC, P, Mu, lamb, d) ;
		        Cbar = C;
			}
			
			double[] A2save=new double[n+1];
			System.arraycopy(Asave, 0, A2save, 0, n);
			A2save[n] = A;
			Asave=A2save;
			double[] B2save=new double[n+1];
			System.arraycopy(Bsave, 0, B2save, 0, n);
			B2save[n] = B;
			Bsave=B2save;
			double[][] Mu2save=new double[n+1][3];
			System.arraycopy(Musave, 0, Mu2save, 0, n);
			Mu2save[n] = Mu;
			Musave=Mu2save;
		    L=GaussianFitCriterion(Y, lamb, x1, x2, x3, A, B, P, Mu, C, d);
		    double[] F2=new double[n+1];
			System.arraycopy(F, 0, F2, 0, n);
			F2[n] = L[0];
			F=F2;
			double[] LS2=new double[n+1];
			System.arraycopy(LS, 0, LS2, 0, n);
			LS2[n] = L[1];
			LS=LS2;
			double[] KL2=new double[n+1];
			System.arraycopy(KL, 0, KL2, 0, n);
			KL2[n] = L[2];
			KL=KL2;
		    if (Math.abs(F[n]-F[n-1])< stop*Math.abs(F[n-1])) {
		    	
		    	txtwindow.append("Nb iter ="+n+": STOP (F does not move anymore)");
		    	
		    	break;
		    }
		    double[][][] Ghold=new double[width][height][depth];
		    for (int i=0; i < width; i++) {
				for (int j=0; j < height; j++) {
					for (int k=0; k < depth; k++) {
						Ghold[i][j][k]=Gn[i][j][k];
					}
				}
			}
		    Gn=(double[][][]) GaussianDistrib3D(x1, x2, x3, Mu, C)[0];
		    for (int i=0; i < width; i++) {
				for (int j=0; j < height; j++) {
					for (int k=0; k < depth; k++) {
						Gn[i][j][k]=A+B*Gn[i][j][k];
					}
				}
			}
		    
		    double maxest=0.;
		    for (int i=0; i<4; i++) {
		    	if (estimation[i]>maxest) {
		    		maxest= estimation[i];
		    	}
		    }
		    if(maxest==1) {
		    	double norm1=0.;
		    	double norm2=0.;
		    	for (int i=0; i < width; i++) {
					for (int j=0; j < height; j++) {
						for (int k=0; k < depth; k++) {
							norm1=norm1+Math.pow(Gn[i][j][k]-Ghold[i][j][k], 2);
							norm2=norm2+Math.pow(Ghold[i][j][k], 2);
						}
					}
				}
		    	norm1=Math.sqrt(norm1);
		    	norm2=Math.sqrt(norm2);
		        double error = norm1/norm2;
		        if(error<stopp) {
		        	
		        	txtwindow.append("Nb iter = "+n+": STOP (G does not move anymore)");
		        	
		            break;
		    	}
			}
		}
		return new Object[] {A,B,P,Mu,C,F,LS,KL};
	}
	
	public double[] GaussianFitCriterion(double[][][] Y, double lamb, double[] x1, double[] x2, double[] x3, double A, double B, double[][][] P, double[] m, double[][] C, double d) {
		double a=0.;		
		for (int i=0; i < width; i++) {
			for (int j=0; j < height; j++) {
				for (int k=0; k < depth; k++) {
					a=a+Math.pow(Y[i][j][k]-A-B*P[i][j][k],2);
				}
			}
		}		
		double LS = 1./2.*a;	
		double[][][] G = (double[][][]) GaussianDistrib3D(x1, x2, x3,m,C)[0];
		double KL=d*(DKL(P,P)-DKL(G,P));
		double F = LS + lamb * KL;
		return new double[] {F,LS,KL};
	}
	
	public double DKL(double[][][] x, double[][][] y) {
		double epsilon=1e-30;
		double f=0.;
		for (int i=0; i<x.length; i++) {
			for (int j=0; j<x[0].length; j++) {
				for (int k=0; k<x[0][0].length; k++) {
					if (x[i][j][k]>epsilon) {
						f=f+y[i][j][k]*Math.log(x[i][j][k]);
					}
				}
			}
		}
		return f;
	}
	
	public Object[] GaussianDistrib3D(double[] x1, double[] x2, double[] x3, double[] mu, double[][] C) {
		Matrix CMatrix = new Jama.Matrix(C);
		double detS=1./CMatrix.det();
		double cte = 1./Math.sqrt(Math.pow((2*Math.PI),3)*detS) ;
		double[][][] G=new double[width][height][depth];
		for (int i=0; i < width; i++) {
			for (int j=0; j < height; j++) {
				for (int k=0; k < depth; k++) {
					Matrix a=new Jama.Matrix(new double[][] {mu});
					Matrix b=new Jama.Matrix(new double[][] {{x1[i],x2[j],x3[k]}});
					b.minusEquals(a);
					Matrix b2=b.transpose();
					G[i][j][k]=cte*Math.exp(-1./2.*b.times(CMatrix.times(b2)).get(0,0));
				}
			}
		}
		return new Object[] {G,cte};
	}

	public double ProxFA(double Abar, double gamA, double B, double[][][] P, double[][][] Y) {
		double a=0.;		
		for (int i=0; i < width; i++) {
			for (int j=0; j < height; j++) {
				for (int k=0; k < depth; k++) {
					a=a+Y[i][j][k]-B*P[i][j][k];
				}
			}
		}
		double A=(Abar+gamA*a)/(1.+gamA*width*height*depth);
		return A;
	}
	
	public double ProxFB(double Bbar, double gamB, double A, double[][][] P, double[][][] Y) {
		double a=0.;
		double b=0.;		
		for (int i=0; i < width; i++) {
			for (int j=0; j < height; j++) {
				for (int k=0; k < depth; k++) {
					a=a+Y[i][j][k]*P[i][j][k]-A*P[i][j][k];
					b=b+Math.pow(P[i][j][k], 2);
				}
			}
		}
		double B=(Bbar+gamB*a)/(1.+gamB*b);
		return Math.max(0., B);
	}
	
	public double[][] ProxFC(double[][] Cbar, double[] x1, double[] x2, double[] x3, double gamC, double[][][] P, double[] mu, double lambda, double d) {
		double sump=0.;
		for (int i=0; i < width; i++) {
			for (int j=0; j < height; j++) {
				for (int k=0; k < depth; k++) {
					sump=sump+P[i][j][k];
				}
			}
		}
		Matrix S=new Matrix(new double[][] {{0,0,0},{0,0,0},{0,0,0}});
		Matrix CMatrix=new Matrix(Cbar);
		Matrix muMatrix=new Matrix(new double[][] {{mu[0]},{mu[1]},{mu[2]}});
		for (int i=0; i < width; i++) {
			for (int j=0; j < height; j++) {
				for (int k=0; k < depth; k++) {
					Matrix X=new Matrix(new double[][] {{x1[i]},{x2[j]},{x3[k]}});
					X.minusEquals(muMatrix);
					S.plusEquals(X.times(X.transpose()).times(P[i][j][k]));
				}
			}
		}
		
		
		S=S.times(1./2.*gamC*lambda*d);
		Matrix M=CMatrix.minus(S);
		double m0=sump*gamC*lambda*d/2.;
		Matrix K=M.times(M).plus(new Matrix(new double[][] {{1,0,0},{0,1,0},{0,0,1}}).times(4.*m0));
		EigenvalueDecomposition E=K.eig();
		Matrix D=E.getD();
		D.set(0, 0, Math.sqrt(D.get(0,0)));
		D.set(1, 1, Math.sqrt(D.get(1,1)));
		D.set(2, 2, Math.sqrt(D.get(2,2)));
		Matrix V=E.getV();
		M.plusEquals(V.times(D.times(V.transpose())));
		M.timesEquals(1./2);
		double[][] A=M.getArray();
		return A;
	}
	
	public double[] ProxFmu(double[] Mubar, double gamM, double[] x1, double[] x2, double[] x3, double[][] C, double lambda, double[][][] P, double d) {
		double sump=0.;
		for (int i=0; i < width; i++) {
			for (int j=0; j < height; j++) {
				for (int k=0; k < depth; k++) {
					sump=sump+P[i][j][k];
				}
			}
		}
		Matrix CMatrix=new Matrix(C);
		Matrix A=new Matrix(new double[][] {{1,0,0},{0,1,0},{0,0,1}});
		A.plusEquals(CMatrix.times(gamM*lambda*d*sump));
		Matrix B=new Matrix(new double[][] {{0},{0},{0}});
		for (int i=0; i < width; i++) {
			for (int j=0; j < height; j++) {
				for (int k=0; k < depth; k++) {
					B.plusEquals(CMatrix.times(new Matrix(new double[][] {{x1[i]},{x2[j]},{x3[k]}})).times(P[i][j][k]));
				}
			}
		}
		B=B.times(gamM*lambda*d);
		B.plusEquals(new Matrix(new double[][] {{Mubar[0]},{Mubar[1]},{Mubar[2]}}));
		B=A.solve(B);
		double[] Mu=new double[] {B.get(0,0),B.get(1,0),B.get(2,0)};
		return Mu;
	}
	
	public double[][][] ProxFPlambda0(double[][][] Pbar, double gam, double dd, double[][][] Y, double A, double B) {
		double N=width*height*depth;
		double Lagrange=0.;
		for (int x=0; x < width; x++) {
			for (int y=0; y < height; y++) {
				for (int z=0; z < depth; z++) {
					Lagrange=Lagrange+Pbar[x][y][z]+(Y[x][y][z]-A)*B*gam;
				}
			}
		}
		Lagrange=1./N*(Lagrange-(1./dd)*(1.+Math.pow(B,2)*gam));
		
		double[][][] P=new double[width][height][depth];
		for (int x=0; x < width; x++) {
			for (int y=0; y < height; y++) {
				for (int z=0; z < depth; z++) {
					P[x][y][z]=(Pbar[x][y][z]-Lagrange+(Y[x][y][z]-A)*B*gam)/(1.+Math.pow(B,2)*gam);
				}
			}
		}
		return P;
	}
	
	public Object[] ProxFP(double[][][] Pbar, double nu0, double gamP, double lambda, double d, double dd, double[][] C, double[] x1, double[] x2, double[] x3, double[] Mu, double A, double B, double[][][] Y) {
		double[][][] CC=new double[width][height][depth];
		Matrix CMatrix=new Matrix(C);
		for (int i=0; i < width; i++) {
			for (int j=0; j < height; j++) {
				for (int k=0; k < depth; k++) {
					Matrix a=new Jama.Matrix(new double[][] {Mu});
					Matrix b=new Jama.Matrix(new double[][] {{x1[i],x2[j],x3[k]}});
					b.minusEquals(a);
					Matrix b2=b.transpose();
					CC[i][j][k]=1./2.*b.times(CMatrix.times(b2)).get(0,0)-1./2.*Math.log(CMatrix.det())+(3./2.)*Math.log(2.*Math.PI);
				}
			}
		}
		double thet = 1.;  
		double tol = 1e-9;
		int nmax=10;
		double rho3 = (thet * gamP * Math.pow(B, 2)+ 1.)/(thet * gamP * lambda * d) ;
		Object[] O=proxEntropyNu(Pbar, gamP, lambda, thet, d, CC, B, A, Y, nu0, dd);
		double phi=(double) O[0];
		double[][][] w=(double[][][]) O[1];
		double sumw1=0.;
		double cst = -1./(gamP * thet * lambda * d * rho3);
		double[] nuv=new double[nmax];
		double[] enuv=new double[nmax];
		for (int i=0; i < width; i++) {
			for (int j=0; j < height; j++) {
				for (int k=0; k < depth; k++) {
					sumw1=sumw1+(1.-1./(w[i][j][k]+1));
				}
			}
		}
		double dphi=cst*sumw1;
		nuv[0] = nu0 - phi/dphi;
		enuv[0]=Math.abs(nuv[0]-nu0);
		int index = 1;
		while (enuv[index-1] >= tol && index <nmax) {
			O = proxEntropyNu(Pbar, gamP, lambda, thet, d, CC, B, A, Y, nuv[index-1], dd);
			phi =(double) O[0];
	    	w=(double[][][]) O[1];
	    	sumw1=0.;
	    	for (int i=0; i < width; i++) {
	    		for (int j=0; j < height; j++) {
	    			for (int k=0; k < depth; k++) {
	    				sumw1=sumw1+(1.-1./(w[i][j][k]+1.));
	    			}
	    		}
	    	}
	    	dphi = cst * sumw1;
	    	nuv[index] = nuv[index-1] - phi/dphi;
	    	enuv[index] = Math.abs(nuv[index]-nuv[index-1]);
	    	index = index+1;
		}
		double nustar = nuv[index-1];
		double[][][] P=(double[][][]) proxEntropyNu(Pbar, gamP, lambda, thet, d, CC, B, A, Y, nustar, dd)[1];
		for (int i=0; i < width; i++) {
    		for (int j=0; j < height; j++) {
    			for (int k=0; k < depth; k++) {
    				P[i][j][k]=P[i][j][k]/rho3;
    			}
    		}
    	}
		return new Object[] {P,nustar};
	}
	
	public Object[] proxEntropyNu(double[][][] Pbar, double gam, double lamb, double thet, double d, double[][][] CC, double B, double A, double[][][] Y, double nu, double dd){
		double rho1 = (thet * gam * Math.pow(B,2)+1.)/(thet * gam * lamb * d);
		double phi=0.;
		double[][][] P=new double[width][height][depth];
		for (int i=0; i < width; i++) {
			for (int j=0; j < height; j++) {
				for (int k=0; k < depth; k++) {
					double tau=Math.log(rho1)-1.-CC[i][j][k]+(Pbar[i][j][k]+gam*B*(Y[i][j][k]-A)-nu)/(thet * gam * lamb * d);
					double limit=200.;
					if (tau<limit) {
						double s=LambertW(Math.exp(tau));
						P[i][j][k]=s;
						phi=phi+s;
					}
					else {
						double s=tau-Math.log(tau);
						P[i][j][k]=s;
						phi=phi+s;
					}
				}
			}
		}
		phi=phi/rho1-dd;
		return new Object[] {phi, P};
	}
	
	public double LambertW(double x) {
		double v=0.;
		double w=1.;
		while (Math.abs(w-v)/Math.abs(w)>1e-8) {
			v = w;
			double e = Math.exp(w);
			double f = w*e - x;  
			w = w - f/((e*(w+1.) - (w+2.)*f/(2.*w+2.)));
		}
		return w;
	}
	
	
	//this function is std(A(:))
	public double std(double[][] A) {
		int N1=A.length;
		int N2=A[0].length;
		double[] B=new double[N1*N2];
		for (int i=0; i<N1; i++) {
			for (int j=0; j<N2; j++) {
				B[j+N2*i]=A[i][j];
			}
		}
		double u=0.;
		for (int i=0; i<N1*N2; i++) {
			u=u+B[i];
		}
		u=u/(N1*N2);
		double V=0.;
		for (int i=0; i<N1*N2; i++) {
			V=V+Math.pow(B[i]-u,2);
		}
		V=Math.sqrt(V/(N1*N2-1));
		return V;
	}
	
	//Cross product
	public double[] cross(double[] A, double[] B){
		return new double[] {A[1]*B[2]-A[2]*B[1],A[2]*B[0]-A[0]*B[2],A[0]*B[1]-A[1]*B[0]};
	}
	
	
	public double dot(double[] A, double[] B){
		double S=0.;
		for (int i=0; i<A.length; i++) {
			S=S+A[i]*B[i];
		}
		return S;
	}
	
	
	public double atan2d(double x,double y) {
		double f = 0;
		if (x == 0) {
			if (y >= 0) {
				f=Math.PI/2;
			}
			    
			else {
				f=3./2*Math.PI;
			}
		}
		    
		else if (x > 0){
				if (y>=0) {
					f=Math.atan(y/x);
				}
				else {
					f=Math.atan(y/x)+2*Math.PI;
				}
		}
		else {
			if (y==0) {
				f=180;
			}
			else {
				f=Math.atan(y/x)+Math.PI;
			}
		}
		
		return f;
	}
	
	//norm of a vector
	public double norm(double[] A){
		double S=0.;
		for (int i=0; i<A.length; i++) {
			S=S+Math.pow(A[i],2);
		}
		S=Math.sqrt(S);
		return S;
	}

	//wiener2 filter
	public double[][] wiener2(double[][] I){
		
		int n=I.length;
		int m=I[0].length;
		double[][] FilteredI= new double[n][m];
		double[][] uMatrix=new double[n][m];
		double[][] o2Matrix=new double[n][m];
		double o2mean=0;
		
		for (int i=0; i < n; i++) {
			for (int j=0; j < m; j++) {
				if(i==0 && j==0) {
					double u=(I[0][0]+I[1][0]+I[0][1]+I[1][1])/9;
					double u2=Math.pow(u,2);
					double o2=(Math.pow(I[0][0],2)+Math.pow(I[1][0],2)+Math.pow(I[0][1],2)+Math.pow(I[1][1],2)-9*u2)/9;
					o2mean=o2mean+o2;
					uMatrix[i][j]=u;
					o2Matrix[i][j]=o2;
				}
				else if (i==0 && j==m-1) {
					double u=(I[0][j]+I[0][j-1]+I[1][j-1]+I[1][j])/9;
					double u2=Math.pow(u,2);
					double o2=(Math.pow(I[0][j],2)+Math.pow(I[1][j-1],2)+Math.pow(I[0][j-1],2)+Math.pow(I[1][j],2)-9*u2)/9;
					o2mean=o2mean+o2;
					uMatrix[i][j]=u;
					o2Matrix[i][j]=o2;
				}
				else if (i==n-1 && j==m-1) {
					double u=(I[i][j]+I[i][j-1]+I[i-1][j-1]+I[i-1][j])/9;
					double u2=Math.pow(u,2);
					double o2=(Math.pow(I[i][j],2)+Math.pow(I[i][j-1],2)+Math.pow(I[i-1][j-1],2)+Math.pow(I[i-1][j],2)-9*u2)/9;
					o2mean=o2mean+o2;
					uMatrix[i][j]=u;
					o2Matrix[i][j]=o2;
				}
				else if (i==n-1 && j==0) {
					double u=(I[i][0]+I[i-1][0]+I[i][1]+I[i-1][1])/9;
					double u2=Math.pow(u,2);
					double o2=(Math.pow(I[i][0],2)+Math.pow(I[i-1][0],2)+Math.pow(I[i][1],2)+Math.pow(I[i-1][1],2)-9*u2)/9;
					o2mean=o2mean+o2;
					uMatrix[i][j]=u;
					o2Matrix[i][j]=o2;
				}
				else if (i==0) {
					double u=(I[0][j]+I[0][j+1]+I[0][j-1]+I[1][j]+I[1][j+1]+I[1][j-1])/9;
					double u2=Math.pow(u,2);
					double o2=(Math.pow(I[0][j],2)+Math.pow(I[0][j-1],2)+Math.pow(I[0][j+1],2)+Math.pow(I[1][j],2)+Math.pow(I[1][j+1],2)+Math.pow(I[1][j-1],2)-9*u2)/9;
					o2mean=o2mean+o2;
					uMatrix[i][j]=u;
					o2Matrix[i][j]=o2;
				}
				else if (i==n-1) {
					double u=(I[i][j]+I[i][j+1]+I[i][j-1]+I[i-1][j]+I[i-1][j+1]+I[i-1][j-1])/9;
					double u2=Math.pow(u,2);
					double o2=(Math.pow(I[i][j],2)+Math.pow(I[i][j-1],2)+Math.pow(I[i][j+1],2)+Math.pow(I[i-1][j],2)+Math.pow(I[i-1][j+1],2)+Math.pow(I[i-1][j-1],2)-9*u2)/9;
					o2mean=o2mean+o2;
					uMatrix[i][j]=u;
					o2Matrix[i][j]=o2;
				}
				else if (j==0) {
					double u=(I[i][j]+I[i+1][j]+I[i-1][j]+I[i][j+1]+I[i+1][j+1]+I[i-1][j+1])/9;
					double u2=Math.pow(u,2);
					double o2=(Math.pow(I[i][j],2)+Math.pow(I[i+1][j],2)+Math.pow(I[i-1][j],2)+Math.pow(I[i][j+1],2)+Math.pow(I[i-1][j+1],2)+Math.pow(I[i+1][j+1],2)-9*u2)/9;
					o2mean=o2mean+o2;
					uMatrix[i][j]=u;
					o2Matrix[i][j]=o2;
				}
				else if (j==m-1) {
					double u=(I[i][j]+I[i+1][j]+I[i-1][j]+I[i][j-1]+I[i-1][j-1]+I[i+1][j-1])/9;
					double u2=Math.pow(u,2);
					double o2=(Math.pow(I[i][j],2)+Math.pow(I[i+1][j],2)+Math.pow(I[i-1][j],2)+Math.pow(I[i][j-1],2)+Math.pow(I[i-1][j-1],2)+Math.pow(I[i+1][j-1],2)-9*u2)/9;
					o2mean=o2mean+o2;
					uMatrix[i][j]=u;
					o2Matrix[i][j]=o2;
				}
				else {
					double u=(I[i][j]+I[i][j+1]+I[i][j-1]+I[i-1][j]+I[i-1][j+1]+I[i-1][j-1]+I[i+1][j]+I[i+1][j+1]+I[i+1][j-1])/9;
					double u2=Math.pow(u,2);
					double o2=(Math.pow(I[i][j],2)+Math.pow(I[i][j-1],2)+Math.pow(I[i][j+1],2)+Math.pow(I[i-1][j],2)+Math.pow(I[i-1][j+1],2)+Math.pow(I[i-1][j-1],2)+Math.pow(I[i+1][j],2)+Math.pow(I[i+1][j+1],2)+Math.pow(I[i+1][j-1],2)-9*u2)/9;
					o2mean=o2mean+o2;
					uMatrix[i][j]=u;
					o2Matrix[i][j]=o2;
				}
				
			}
		}
		
		o2mean=o2mean/(n*m);
		for (int i=0; i < n; i++) {
			for (int j=0; j < m; j++) {
				FilteredI[i][j]=uMatrix[i][j]+Math.max((o2Matrix[i][j]-o2mean), 0.)*(I[i][j]-uMatrix[i][j])/Math.max(o2Matrix[i][j], o2mean);
			}
		}
		
		return FilteredI;
	}
	
	

	public void showAbout() {
		IJ.showMessage("FIGARO",
			"Fitting Gaussians with Proximal Optimization"
		);
	}
	  
	//From ImagePlus to awt.Image scaled to fit in the report (unsure if always scale properly)
	public Object[] awtScaledImage(ImagePlus image)
	  {
	    
	    java.awt.Image img=image.getImage();
	    BufferedImage bimg=toBufferedImage(img);
	    
	    int imheight=bimg.getHeight();
	    int imwidth=bimg.getWidth();
	    float zoom=0;
	    if (750*imwidth>=75*imheight) {
	    	zoom=75f/imwidth;
	    }
	    else {
	    	zoom=750f/imheight;
	    }
	    
	    java.awt.Image img2=img.getScaledInstance((int) Math.floor(image.getWidth()*zoom),(int) Math.floor(image.getHeight()*zoom),java.awt.Image.SCALE_DEFAULT);
	    
	    
	    return new Object[] {img2,zoom};
	  }
	
	
	 public Paragraph bigTitle(String title)
	  {
	    Paragraph reportTitle = new Paragraph();
	    reportTitle.add(title);
	    reportTitle.setTextAlignment(TextAlignment.CENTER);
	    reportTitle.setBold();
	    reportTitle.setFontSize(16f);
	    reportTitle.setFixedLeading(20.0F);
	    reportTitle.setPaddingTop(15.0F);
	    reportTitle.setPaddingBottom(15.0F);
	    return reportTitle;
	  }
	 
	 //awt.Image to BufferedImage
	 public static BufferedImage toBufferedImage(java.awt.Image img)
	 {
	     if (img instanceof BufferedImage)
	     {
	         return (BufferedImage) img;
	     }

	     // Create a buffered image with transparency
	     BufferedImage bimage = new BufferedImage(img.getWidth(null), img.getHeight(null), BufferedImage.TYPE_INT_ARGB);

	     // Draw the image on to the buffered image
	     Graphics2D bGr = bimage.createGraphics();
	     bGr.drawImage(img, 0, 0, null);
	     bGr.dispose();

	     // Return the buffered image
	     return bimage;
	 }
	 
	 //From a,b,c in cartesian equation of an ellipse to angle for parametric equation
	 public static double ellipseAngle(double a, double b, double c) {
		 if (a==c) {
			 return Math.PI/4;
		 }
		 else if (b==0 && a<c) {
			 return 0.;
		 }
		 else if (b==0 && a>c) {
			 return Math.PI/2;
		 }
		 else if (b!=0 && a<c) {
			 return 1./2*Math.atan(2*b/(a-c));
		 }
		 else {
			 return Math.PI/2+1./2*Math.atan(2*b/(a-c));
		 }
		 
	 }
	 
	//From a,b,c in cartesian equation of an ellipse to axis lengths for parametric equation
	 public static double[] ellipseLengths(double a, double b, double c, double g) {
		 double A=Math.sqrt(2*(g*Math.pow(b,2)-a*c*g)/((Math.pow(b,2)-a*c)*(Math.sqrt(Math.pow(a-c,2)+4*Math.pow(b,2))-a-c)));
		 double B=Math.sqrt(2*(g*Math.pow(b,2)-a*c*g)/((Math.pow(b,2)-a*c)*(-Math.sqrt(Math.pow(a-c,2)+4*Math.pow(b,2))-a-c)));
	 
		 return new double[] {A,B};
	 }
	 
	 
	 
	 //this main function only serves a test purpose. It automatically launches ImageJ and the plugin on one of your image specified by the path
	 
	 public static void main(String[] args) {
			// set the plugins.dir property to make the plugin appear in the Plugins menu
			Class<?> clazz = FIGARO_.class;
			String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
			String pluginsDir = url.substring("file:".length(), url.length() - clazz.getName().length() - ".class".length());
			System.setProperty("plugins.dir", pluginsDir);

			// start ImageJ
			new ImageJ();
			
			//Use your own path to test 
			ImagePlus image = IJ.openImage("C:\\Users\\Flash\\eclipse-workspace\\cvn\\image\\0.2um_mycrop24.tif");
			
			image.show();

			// run the plugin
			IJ.runPlugIn(clazz.getName(), "");
		}


}